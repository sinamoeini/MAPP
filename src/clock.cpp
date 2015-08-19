#include "clock.h"
#include "thermo_dynamics.h"
#include "atom_types.h"
#include "ff.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock::Clock(MAPP* mapp):InitPtrs(mapp)
{
    ns_alloc=0;
    if(forcefield==NULL)
        error->abort("ff should be "
        "initiated before clock");
    
    if(mapp->mode!=DMD_mode)
        error->abort("clock works only "
        "for md mode");
    
    char** args;
    int narg=mapp->parse_line((char*)
    "Time FE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    
    fe_idx=1;
    stress_idx=2;
    time_idx=0;
    thermo=new ThermoDynamics(mapp,narg,args);
    for(int i=0;i<narg;i++)
        delete [] args[i];
    delete [] args;
    
    cdof_n=atoms->find_exist("cdof");
    
    int dim=atoms->dimension;
    if(dim)
    {
        CREATE1D(nrgy_strss,dim*(dim+1)/2+1);
        ns_alloc=1;
    }
    
    
    c_n=atoms->find("c");
    c_d_n=atoms->find("c_d");
    
    int c_dim=atoms->vectors[c_n]->dim;
    dof_lcl=atoms->natms*c_dim;
    
    int* c;
    atoms->vectors[c_n]->ret(c);
    
    int tmp_dof=dof_lcl;
    for(int idof=0;idof<dof_lcl;idof++)
        if(c[idof]<0)
            tmp_dof--;
    
    MPI_Allreduce(&tmp_dof,&dof_tot,1,MPI_INT,MPI_SUM,world);
    
    CREATE1D(y_0,dof_lcl);
    CREATE1D(a,dof_lcl);
    CREATE1D(g,dof_lcl);
    CREATE1D(g0,dof_lcl);
    CREATE1D(c0,dof_lcl);
    CREATE1D(h,dof_lcl);
    
    //defaults for solver
    min_gamma=0.0;
    gamma_red=0.5;
    slope=0.4;
    max_iter=50;
    m_tol=1.0e-9;
    a_tol=1.0e-6;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock::~Clock()
{
    if(dof_lcl)
    {
        delete [] y_0;
        delete [] a;
        delete [] g;
        delete [] g0;
        delete [] c0;
        delete [] h;
    }
    delete thermo;
    if(ns_alloc)
        delete [] nrgy_strss;
}

/*--------------------------------------------
 rectify
 --------------------------------------------*/
void Clock::rectify(type0* f)
{
    if(cdof_n==-1)
        return;
    atoms->vectors[cdof_n]->ret(cdof);
    
    int tot=(atoms->natms)*(atom_types->no_types);
    for(int i=0;i<tot;i++) if(cdof[i]==1) f[i]=0.0;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock::solve_n_err(type0& cost,type0& err)
{type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    
    type0 gamma,max_gamma=1.0;
    type0 inner,tmp0;
    type0 err_lcl;
    type0 ratio;
    type0 g0_g0,g_g,g_g0,g_h;
    type0 curr_cost,ideal_cost,init_cost;
    int chk;
    
    /* beginning of pre-conditioning */
    memcpy(c,y_0,dof_lcl*sizeof(type0));
    atoms->update_ph(c_n);
    /* end of pre-conditioning */
    
    /* find the steepest descent direction and cost */
    curr_cost=forcefield->g_calc_timer(1,beta,a,g,nrgy_strss);
    rectify(g);
    
    /* set the first trajectory */
    memcpy(h,g,dof_lcl*sizeof(type0));
    
    /* calculate g.h g_0.g_0 */
    inner=0.0;
    for(int i=0;i<dof_lcl;i++)
        if(c[i]>=0.0)
            inner+=g[i]*g[i];
    g0_g0=0.0;
    MPI_Allreduce(&inner,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
    g_h=g0_g0;
    int iter=0;
    
    
    while(curr_cost>m_tol*static_cast<type0>(dof_tot)
          && iter<max_iter && max_gamma>min_gamma)
    {
        
        memcpy(g0,g,dof_lcl*sizeof(type0));
        memcpy(c0,c,dof_lcl*sizeof(type0));
        
        
        /*
         beginning of finding maximum back tracking length
         */
        gamma=-1.0;
        
        for(int i=0;i<dof_lcl;i++)
        {
            if(c0[i]>=0.0)
            {
                if(gamma==-1.0)
                {
                    if(h[i]>0.0)
                    {
                        gamma=(1.0-c0[i])/h[i];
                    }
                    else if(h[i]<0.0)
                    {
                        gamma=(0.0-c0[i])/h[i];
                    }
                }
                else
                {
                    if(h[i]>0.0)
                    {
                        gamma=MIN((1.0-c0[i])/h[i],gamma);
                    }
                    else if(h[i]<0.0)
                    {
                        gamma=MIN((0.0-c0[i])/h[i],gamma);
                    }
                }
                
            }
        }
        
        if(gamma==-1.0)
            gamma=INFINITY;
        
        
        MPI_Allreduce(&gamma,&max_gamma,1,MPI_TYPE0,MPI_MIN,world);
        max_gamma*=0.999;
        /*
         end of finding maximum back tracking length
         */
        
        
        
        chk=1;
        init_cost=curr_cost;
        
        while(chk && max_gamma>min_gamma)
        {
            
            /*
             update c in h direction
             */
            for(int i=0;i<dof_lcl;i++)
                if(c0[i]>=0.0)
                    c[i]=c0[i]+max_gamma*h[i];
            atoms->update_ph(c_n);
            
            
            /*
             find the resultant cost
             */
            curr_cost=forcefield->g_calc_timer(0,beta,a,g,nrgy_strss);
            
            /* compare ideal cost and current cost */
            ideal_cost=init_cost-slope*max_gamma*g_h;
            if(curr_cost<ideal_cost)
                chk=0;
            max_gamma*=gamma_red;
            
            /* if the gamma dwarfs min_gamma back track all the way */
            if(max_gamma<=min_gamma)
            {
                
                memcpy(c,c0,dof_lcl*sizeof(type0));
                atoms->update_ph(c_n);
                
            }
        }
        
        curr_cost=forcefield->g_calc_timer(1,beta,a,g,nrgy_strss);
        rectify(g);
        /* if the back track was successful */
        if(chk==0)
        {
            /* calculate g_1.g_0 */
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
                if(c[i]>=0.0)
                    inner+=g[i]*g0[i];
            g_g0=0.0;
            MPI_Allreduce(&inner,&g_g0,1,MPI_TYPE0,MPI_SUM,world);
            
            /* calculate g_1.g_1 */
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
                if(c[i]>=0.0)
                    inner+=g[i]*g[i];
            g_g=0.0;
            MPI_Allreduce(&inner,&g_g,1,MPI_TYPE0,MPI_SUM,world);
            
            /* calculate (g_1.g_1-g_1.g_0)/g_0.g_0 */
            ratio=(g_g-g_g0)/g0_g0;
            
            /* update g_0.g_0 */
            g0_g0=g_g;
            
            /* calculate g_h */
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
            {
                if(c[i]>=0.0)
                {
                    h[i]*=ratio;
                    h[i]+=g[i];
                    inner+=h[i]*g[i];
                }
            }
            g_h=0.0;
            MPI_Allreduce(&inner,&g_h,1,MPI_TYPE0,MPI_SUM,world);
            
            /* if g_h is negative start from the begining */
            if(g_h<0.0)
            {
                memcpy(h,g,dof_lcl*sizeof(type0));
                g_h=g_g;
            }
        }
        
        iter++;
    }
    
    rectify(c_d);
    
    err_lcl=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        if(c[i]>=0.0)
        {
            tmp0=(y_0[i]-c[i]);
            err_lcl+=tmp0*tmp0;
        }
    }
    
    err=0.0;
    MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
    err=sqrt(err/static_cast<type0>(dof_tot))/a_tol;
    err*=err_prefac;
    

    cost=curr_cost/(m_tol*static_cast<type0>(dof_tot));
}

