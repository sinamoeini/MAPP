#include "clock_mbdf.h"
#include "neighbor.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_mbdf::Clock_mbdf(MAPP* mapp,int narg
,char** arg):Clock(mapp)
{
    min_gamma=0.0;
    gamma_red=0.8;
    slope=0.4;
    max_iter=50;
    max_order=6;
    m_tol=1.0e-10;
    a_tol=1.0e-6;
    e_tol=0.0;
    min_del_t=1.0e-12;
    max_del_t=1.0e4;
    initial_del_t=-1.0;
    
    if(narg<3)
        error->abort("clock mbdf should at least have 1 arguement");
    
    no_steps=atoi(arg[2]);
    if(no_steps<0)
        error->abort("no of steps for clock mbdf should be greater than 0");
    
    if(narg>3)
    {
        if((narg-3)%2!=0)
            error->abort("every keyword in clock mbdf should be followed by it's value");
        int iarg=3;
        while(iarg<narg)
        {
            if(strcmp(arg[iarg],"min_gamma")==0)
            {
                iarg++;
                min_gamma=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"red_gamma")==0)
            {
                iarg++;
                gamma_red=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"slope")==0)
            {
                iarg++;
                slope=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_iter")==0)
            {
                iarg++;
                max_iter=atoi(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_order")==0)
            {
                iarg++;
                max_order=atoi(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"m_tol")==0)
            {
                iarg++;
                m_tol=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"a_tol")==0)
            {
                iarg++;
                a_tol=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"e_tol")==0)
            {
                iarg++;
                e_tol=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"min_del_t")==0)
            {
                iarg++;
                min_del_t=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_del_t")==0)
            {
                iarg++;
                max_del_t=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"initial_del_t")==0)
            {
                iarg++;
                initial_del_t=atof(arg[iarg]);
                if(initial_del_t<=0.0)
                    error->abort("initial_del_t in clock mbdf should be greater than 0.0");
                iarg++;
            }
            else
                error->abort("unknown keyword in clock mbdf: %s",arg[iarg]);
        }
    }
    
    if(min_gamma<0.0)
        error->abort("min_gamma in clock mbdf should be greater than 0.0");
    if(gamma_red<=0.0 || gamma_red>=1.0)
        error->abort("gamma_red in clock mbdf should be between 0.0 & 1.0");
    if(slope<=0.0 || slope>=1.0)
        error->abort("slope in clock mbdf should be between 0.0 & 1.0");
    if(max_iter<=0)
        error->abort("max_iter in clock mbdf should be greater than 0");
    if(max_order<=0)
        error->abort("max_order in clock mbdf should be greater than 0");
    if(m_tol<=0.0)
        error->abort("m_tol in clock mbdf should be greater than 0.0");
    if(a_tol<=0.0)
        error->abort("a_tol in clock mbdf should be greater than 0.0");
    if(e_tol<0.0)
        error->abort("e_tol in clock mbdf should be equal or greater than 0.0");
    if(min_del_t<=0.0)
        error->abort("min_del_t in clock mbdf should be greater than 0.0");
    if(max_del_t<=0.0)
        error->abort("max_del_t in clock mbdf should be greater than 0.0");
    if(min_del_t>=max_del_t)
        error->abort("max_del_t in clock mbdf should be greater than min_del_t");

    
    
    c_n=atoms->find("c");
    c_d_n=atoms->find("c_d");
    
    int c_dim=atoms->vectors[c_n].dim;
    dof_lcl=atoms->natms*c_dim;
    MPI_Allreduce(&dof_lcl,&dof_tot,1,MPI_INT,MPI_SUM,world);
    
    CREATE1D(dy,dof_lcl);
    CREATE1D(y_0,dof_lcl);
    CREATE1D(a,dof_lcl);
    CREATE1D(g,dof_lcl);
    CREATE1D(g0,dof_lcl);
    CREATE1D(c0,dof_lcl);
    CREATE1D(h,dof_lcl);
    
    CREATE1D(t,max_order+1);
    CREATE1D(alpha_y,max_order+1);
    CREATE1D(dalpha_y,max_order+1);
    CREATE1D(y,max_order+1);
    for(int i=0;i<max_order+2;i++)
        CREATE1D(y[i],dof_lcl);
    
    CREATE1D(alph_err,max_order+2);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_mbdf::~Clock_mbdf()
{
    
    for(int i=0;i<max_order+1;i++)
        if(dof_lcl) delete [] y[i];
    
    delete [] alph_err;
    delete [] y;
    delete [] t;
    delete [] alpha_y;
    delete [] dalpha_y;

    
    if(dof_lcl)
    {
        delete [] dy;
        delete [] y_0;
        delete [] a;
        delete [] g0;
        delete [] c0;
        delete [] g;
        delete [] h;
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_mbdf::init()
{
    type0* c;
    type0* c_d;
    
    
    int f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<type0>(0,atoms->vectors[0].dim,"f");
    
    int id_n=atoms->find("id");
    if(cdof_n>-1)
    {
        int dof_n=atoms->find("dof");
        vecs_comm=new VecLst(mapp,6,0,c_n,c_d_n,cdof_n,dof_n,id_n);
        
    }
    else
    {
        vecs_comm=new VecLst(mapp,4,0,c_n,c_d_n,id_n);
    }
    
    vecs_comm->add_update(0);
    atoms->reset_comm(vecs_comm);
    forcefield->init();
    atoms->ph_setup(1,vecs_comm);
    
    neighbor->init();
    neighbor->create_list(0,1);
    atoms->store_0();
    forcefield->create_2nd_neigh_lst();
    
    
    
    forcefield->force_calc(1,nrgy_strss);
    
    
    thermo->start_force_time();
    forcefield->c_d_calc();
    thermo->stop_force_time();
    
    atoms->vectors[c_n].ret(c);
    atoms->vectors[c_d_n].ret(c_d);
    rectify(c_d);
    if(initial_del_t<0.0)
    {
        type0 sum,sum_lcl;
        sum_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
            sum_lcl+=c_d[i]*c_d[i];
        sum=0.0;
        MPI_Allreduce(&sum_lcl,&sum,1,MPI_TYPE0,MPI_SUM,world);
        sum=sqrt(sum/static_cast<type0>(dof_tot));
        initial_del_t=MAX(a_tol/sum,min_del_t);
        initial_del_t=MIN(initial_del_t,max_del_t);
    }
    
    t[0]=0.0;
    t[1]=-initial_del_t;
    memcpy(y[0],c,dof_lcl*sizeof(type0));
    memcpy(y[1],c,dof_lcl*sizeof(type0));
    memcpy(dy,c_d,dof_lcl*sizeof(type0));
    
    thermo->update(fe_idx,nrgy_strss[0]);
    thermo->update(stress_idx,6,&nrgy_strss[1]);
    thermo->update(time_idx,0.0);
    
    if(write!=NULL)
        write->init();
    thermo->init();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_mbdf::fin()
{
    
    if(write!=NULL)
        write->fin();
    
    forcefield->fin();
    neighbor->fin();
    
    delete vecs_comm;
    thermo->fin();
    
    
    atoms->x2s(atoms->natms);
}
/*--------------------------------------------
 init
 --------------------------------------------*/
int Clock_mbdf::interpolate(type0 del_t,int q)
{
    type0 tmp0,tmp1,tmp2;
    type0 curr_t=t[0]+del_t;
    type0 c0=0.0;
    type0 err_prefac0,err_prefac1;
    
    
    for(int i=0;i<q+1;i++)
    {
        c0+=1.0/(curr_t-t[i]);
        
        tmp0=tmp1=1.0;
        tmp2=0.0;
        for(int j=0;j<q+1;j++)
            if(i!=j)
                tmp0*=(curr_t-t[j])/(t[i]-t[j]);
        
        alpha_y[i]=tmp0;
        dalpha_y[i]=tmp0;

    }

    for(int i=0;i<q+1;i++)
    {
        tmp0=c0-1.0/(curr_t-t[i]);
        dalpha_y[i]*=tmp0;
    }
    
    tmp0=0.0;
    for(int i=0;i<q;i++)
        tmp0+=1.0/static_cast<type0>(i+1);
    beta=del_t/tmp0;

    int idof=0;
    int ret_val=1;
    int tot_ret_val;
    while(idof<dof_lcl && ret_val==1)
    {
        y_0[idof]=0.0;
        a[idof]=0.0;
        for(int j=0;j<q+1;j++)
        {
            a[idof]+=beta*dalpha_y[j]*y[j][idof];
            y_0[idof]+=alpha_y[j]*y[j][idof];
        }
        a[idof]-=y_0[idof];
        
        if(y_0[idof]<0.0 || y_0[idof]>1.0)
        {
            ret_val=-1;
        }
        
        idof++;
    }
    MPI_Allreduce(&ret_val,&tot_ret_val,1,MPI_INT,MPI_MIN,world);
    
    err_prefac0=(curr_t-t[0])/(curr_t-t[q]);
    err_prefac1=err_prefac0;
    for(int i=0;i<q;i++)
    {
        err_prefac1+=(curr_t-t[0])/(curr_t-t[i])
        -1.0/static_cast<type0>(i+1);
    }
    err_prefac=MAX(err_prefac0,fabs(err_prefac1));
    
    return tot_ret_val;

}
/*--------------------------------------------
 run
 --------------------------------------------*/
void Clock_mbdf::run()
{
    if(no_steps==0)
        return;
    
    type0* c;
    atoms->vectors[c_n].ret(c);
    type0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    type0* tmp_y;
    
    type0 del_t=initial_del_t,del_t_tmp;
    type0 err1,cost;
    type0 ratio;
    int ord=1;
    int chk;
    int initial_phase;
    int istep;
    
    initial_phase=1;
    const_stps=0;
    eq_ratio=1.0;
    istep=0;

    while (eq_ratio>=1.0 && istep <no_steps)
    {
        chk=interpolate(del_t,ord);
        while(chk==-1)
        {
            if(initial_phase) initial_phase=0;
            del_t*=0.8;
            del_t=MAX(del_t,min_del_t);
            chk=interpolate(del_t,ord);
        }

        cost=solve(del_t,ord);
        while(err>=1.0 || cost>=1.0)
        {
            if(initial_phase) initial_phase=0;
            err1=MAX(err,cost);
            ratio=pow(0.5/err1,1.0/static_cast<type0>(ord+1));
            
            if(ratio<0.5)
                ratio=0.5;
            else if(ratio>0.9)
                ratio=0.9;
            
            del_t_tmp=del_t*ratio;
            if(del_t_tmp<min_del_t)
                del_t=min_del_t;
            else
                del_t=del_t_tmp;
            
            
            chk=interpolate(del_t,ord);
            
            while(chk==-1)
            {
                del_t*=0.8;
                del_t=MAX(del_t,min_del_t);
                chk=interpolate(del_t,ord);
            }
            
            cost=solve(del_t,ord);
            
        }
        
        
        del_t_tmp=del_t;
        ord_dt(del_t,ord,initial_phase);
        
        tmp_y=y[max_order];
        for(int i=max_order;i>0;i--)
        {
            t[i]=t[i-1];
            y[i]=y[i-1];
        }
        
        y[0]=tmp_y;
        t[0]+=del_t_tmp;
        
        memcpy(y[0],c,dof_lcl*sizeof(type0));
        memcpy(dy,c_d,dof_lcl*sizeof(type0));
        
        
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step() || istep+1==no_steps)
        {
            thermo->start_force_time();
            forcefield->force_calc(1,nrgy_strss);
            thermo->stop_force_time();
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            thermo->update(time_idx,t[0]);
        }
        
        
        step_no++;
        istep++;
    }

    
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
type0 Clock_mbdf::step_size(type0 del_t,int ord)
{
    
    
    type0* c;
    atoms->vectors[c_n].ret(c);
    type0 tmp,ratio,tot_ratio;
    type0 max_h;
    
    
    ratio=INFINITY;
    
    
    for(int i=0;i<dof_lcl;i++)
    {
        if(dy[i]>0.0)
        {
            tmp=(1.0-c[i])/dy[i];
            ratio=MIN(ratio,tmp);
        }
        else if(dy[i]<0.0)
        {
            tmp=(0.0-c[i])/dy[i];
            ratio=MIN(ratio,tmp);
        }
    }
    MPI_Allreduce(&ratio,&tot_ratio,1,MPI_TYPE0,MPI_MIN,world);
    max_h=0.99*tot_ratio;
    
    type0 opt_h=del_t*pow(0.5/err,1.0/static_cast<type0>(ord+1));
    
    if(max_h<opt_h)
        opt_h=max_h;
    
    opt_h=MAX(min_del_t,opt_h);
    
    ratio=opt_h/del_t;
    
    if(ratio>2.0)
        ratio=1.0;
    else if (ratio<0.5)
        ratio=0.5;
    
    opt_h*=ratio;
    
    if(opt_h>=max_del_t)
    {
        opt_h=max_del_t;
        
    }
    else if(opt_h<=min_del_t)
    {
        opt_h=min_del_t;
    }
    ratio=opt_h/del_t;
    
    
    return ratio;
}
/*--------------------------------------------
 init
 --------------------------------------------*/

type0 Clock_mbdf::solve(type0 del_t,int ord)
{
    type0* c;
    atoms->vectors[c_n].ret(c);
    type0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    type0 gamma,max_gamma=1.0;
    type0 inner,tmp0,tmp1;
    type0 err_lcl;
    type0 ratio;
    type0 g0_g0,g_g,g_g0,g_h;
    type0 curr_cost,ideal_cost,cost;
    int chk;
    
    /*
    memcpy(c,y_0,dof_lcl*sizeof(type0));
    thermo->start_comm_time();
    atoms->update(c_n);
    thermo->stop_comm_time();
    */
    
    
    type0 tot_ratio;
    ratio=1.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp0=dy[i]*del_t+y[0][i];
        if(tmp0>1.0)
            ratio=MIN((1.0-y[0][i])/(tmp0-y[0][i]),ratio);
        else if(tmp0<0.0)
            ratio=MIN((y[0][i]-0.0)/(y[0][i]-tmp0),ratio);
    }
    MPI_Allreduce(&ratio,&tot_ratio,1,MPI_TYPE0,MPI_MIN,world);
    for(int i=0;i<dof_lcl;i++)
        c[i]=y[0][i]+dy[i]*del_t*tot_ratio;
    thermo->start_comm_time();
    atoms->update(c_n);
    thermo->stop_comm_time();
    
    
    
    thermo->start_force_time();
    curr_cost=forcefield->g_calc(0,beta,a,g);
    rectify(g);
    thermo->stop_force_time();
    /*
     type0 sum_lcl;
     type0 sum_tot;
     sum_lcl=0.0;
     sum_tot=0.0;
     for(int i=0;i<dof_lcl;i++)
     sum_lcl+=g[i];
     MPI_Allreduce(&sum_lcl,&sum_tot,1,MPI_TYPE0,MPI_SUM,world);
     for(int i=0;i<dof_lcl;i++)
     g[i]-=correc_fac*sum_tot;
     rectify(g);
     */
    memcpy(h,g,dof_lcl*sizeof(type0));
    
    inner=0.0;
    for(int i=0;i<dof_lcl;i++)
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
        
        gamma=0.99;
        
        for(int i=0;i<dof_lcl;i++)
        {
            tmp0=c0[i]+h[i];
            
            if(tmp0>1.0)
            {
                gamma=MIN((1.0-c0[i])/(tmp0-c0[i]),gamma);
            }
            else if(tmp0<0.0)
            {
                gamma=MIN((c0[i]-0.0)/(c0[i]-tmp0),gamma);
            }
        }
        MPI_Allreduce(&gamma,&max_gamma,1,MPI_TYPE0,MPI_MIN,world);
        
        max_gamma*=0.9;
        chk=1;
        
        cost=curr_cost;
        
        while(chk && max_gamma>min_gamma)
        {
            for(int i=0;i<dof_lcl;i++)
                c[i]=c0[i]+max_gamma*h[i];
            
            thermo->start_comm_time();
            atoms->update(c_n);
            thermo->stop_comm_time();
            
            thermo->start_force_time();
            curr_cost=forcefield->g_calc(0,beta,a,g);
            rectify(g);
            thermo->stop_force_time();
            /*
             sum_lcl=0.0;
             sum_tot=0.0;
             for(int i=0;i<dof_lcl;i++)
             sum_lcl+=g[i];
             MPI_Allreduce(&sum_lcl,&sum_tot,1,MPI_TYPE0,MPI_SUM,world);
             for(int i=0;i<dof_lcl;i++)
             g[i]-=correc_fac*sum_tot;
             rectify(g);
             */
            ideal_cost=cost-slope*max_gamma*g_h;
            if(curr_cost<ideal_cost)
                chk=0;
            max_gamma*=gamma_red;
            if(max_gamma<=min_gamma)
            {
                memcpy(c,c0,dof_lcl*sizeof(type0));
                
                thermo->start_comm_time();
                atoms->update(c_n);
                thermo->stop_comm_time();
                
                thermo->start_force_time();
                curr_cost=forcefield->g_calc(1,beta,a,g);
                thermo->stop_force_time();
            }
        }
        
        if(chk==0)
        {
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
                inner+=g[i]*g0[i];
            g_g0=0.0;
            MPI_Allreduce(&inner,&g_g0,1,MPI_TYPE0,MPI_SUM,world);
            
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
                inner+=g[i]*g[i];
            g_g=0.0;
            MPI_Allreduce(&inner,&g_g,1,MPI_TYPE0,MPI_SUM,world);
            
            ratio=(g_g-g_g0)/g0_g0;
            
            g0_g0=g_g;
            
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
            {
                h[i]*=ratio;
                h[i]+=g[i];
                inner+=h[i]*g[i];
            }
            g_h=0.0;
            MPI_Allreduce(&inner,&g_h,1,MPI_TYPE0,MPI_SUM,world);
            
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
    tmp1=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp0=(y_0[i]-c[i]);
        tmp1+=c_d[i]*c_d[i];
        err_lcl+=tmp0*tmp0;
    }
    
    err=0.0;
    MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
    err=sqrt(err/static_cast<type0>(dof_tot))/a_tol;
    err*=err_prefac;
    
    MPI_Allreduce(&tmp1,&eq_ratio,1,MPI_TYPE0,MPI_SUM,world);
    eq_ratio=sqrt(eq_ratio/static_cast<type0>(dof_tot))/e_tol;
    
    
    for(int i=0;i<dof_lcl;i++)
        if(c[i]<0.0 || c[i]>1.0)
            error->abort("c exceeded the domain");
    
    return curr_cost/(m_tol*static_cast<type0>(dof_tot));
}

/*--------------------------------------------
 factorial
 --------------------------------------------*/
type0 Clock_mbdf::fac(int no)
{
    int ans=1;
    for(int i=1;i<no+1;i++)
        ans*=i;
    return static_cast<type0>(ans);
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
void Clock_mbdf::ord_dt(type0& del_t,int& q
,int init_ph)
{
    if(init_ph)
    {
        del_t=MIN(del_t*2.0,max_del_t);
        if(q<max_order)
            q++;
        return;
    }
    
    type0 curr_t=t[0]+del_t;
    type0 tmp0,tmp1;
    type0 terkm2=0.0,terkm1=0.0,terk=0.0,terkp1=0.0;
    type0 tmp_err,ratio,del_t_tmp;
    
    int terkm2_flag,terkm1_flag,terk_flag,terkp1_flag
    ,tmp_q=q;
    
    terk_flag=1;
    terkp1_flag=0;
    terkm1_flag=0;
    terkm2_flag=0;
    
    if(q>2)
        terkm2_flag=1;
    if(q>1)
        terkm1_flag=1;
    
    
    if(const_stps==q+2 && q<max_order)
        terkp1_flag=1;
    
    if(terkm2_flag)
    {
        tmp1=fac(q-1)*pow(del_t,q-1);
        for(int i=0;i<q;i++)
        {
            if(i!=0)
            {
                tmp0=t[i-1]-curr_t;
                for(int j=0;j<q-1;j++)
                {
                    if(j!=i-1)
                        tmp0*=t[i-1]-t[j];
                }
            }
            else
            {
                tmp0=1.0;
                for(int j=0;j<q-1;j++)
                    tmp0*=curr_t-t[j];
            }
            
            alph_err[i]=tmp1/tmp0;
        }
        terkm2=err_est(q);
    }
    if(terkm1_flag)
    {
        tmp1=fac(q)*pow(del_t,q);
        for(int i=0;i<q+1;i++)
        {
            if(i!=0)
            {
                tmp0=t[i-1]-curr_t;
                for(int j=0;j<q;j++)
                {
                    if(j!=i-1)
                        tmp0*=t[i-1]-t[j];
                }
            }
            else
            {
                tmp0=1.0;
                for(int j=0;j<q;j++)
                    tmp0*=curr_t-t[j];
            }
            
            alph_err[i]=tmp1/tmp0;
        }
        terkm1=err_est(q+1);
    }
    if(terk_flag)
    {
        tmp1=fac(q+1)*pow(del_t,q+1);
        for(int i=0;i<q+2;i++)
        {
            if(i!=0)
            {
                tmp0=t[i-1]-curr_t;
                for(int j=0;j<q+1;j++)
                {
                    if(j!=i-1)
                        tmp0*=t[i-1]-t[j];
                }
            }
            else
            {
                tmp0=1.0;
                for(int j=0;j<q+1;j++)
                    tmp0*=curr_t-t[j];
            }
            
            alph_err[i]=tmp1/tmp0;
        }
        terk=err_est(q+2);
    }
    if(terkp1_flag)
    {
        tmp1=1.0;
        for(int i=0;i<q+3;i++)
        {
            if(i!=0)
            {
                tmp0=t[i-1]-curr_t;
                for(int j=0;j<q+2;j++)
                {
                    if(j!=i-1)
                        tmp0*=t[i-1]-t[j];
                }
            }
            else
            {
                tmp0=1.0;
                for(int j=0;j<q+2;j++)
                    tmp0*=curr_t-t[j];
            }
            
            alph_err[i]=tmp1/tmp0;
        }
        terkp1=err_est(q+3);
    }
    
    
    
    if(q>2)
    {
        if(MAX(terkm1,terkm2)<=terk)
        {
            tmp_q=q-1;
            tmp_err=terkm1;
        }
        else
        {
            if(terkp1<terk && terkp1_flag)
            {
                tmp_q=q+1;
                tmp_err=terkp1;
            }
            else
                tmp_err=terk;
        }
    }
    else if(q==2)
    {
        if(terkm1<=0.5*terk)
        {
            tmp_q=q-1;
            tmp_err=terkm1;
        }
        else if(terk<terkm1 && terkp1_flag)
        {
            if(terkp1<terk)
            {
                tmp_q=q+1;
                tmp_err=terkp1;
            }
            else
                tmp_err=terk;
        }
        else
            tmp_err=terk;
    }
    else
    {
        if(terkp1<0.5*terk && terkp1_flag)
        {
            tmp_q=q+1;
            tmp_err=terkp1;
        }
        else
            tmp_err=terk;
    }
    
    
    ratio=pow(0.5/tmp_err,1.0/static_cast<type0>(q+1));
    
    if(ratio>=2.0)
        const_stps=0;
    else if(ratio<=0.5)
    {
        ratio=0.5;
        const_stps=0;
    }
    else
    {
        ratio=1.0;
        if(tmp_q==q)
            const_stps++;
    }
    
    del_t_tmp=ratio*del_t;
    if(del_t_tmp>max_del_t)
        del_t=max_del_t;
    else if (del_t_tmp<min_del_t)
        del_t=min_del_t;
    else
        del_t=del_t_tmp;
    q=tmp_q;
    
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
type0 Clock_mbdf::err_est(int q)
{
    type0* c;
    atoms->vectors[c_n].ret(c);
    
    type0 tmp0,err_lcl,err_tot;
    err_lcl=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp0=alph_err[0]*c[i];
        for(int j=0;j<q-1;j++)
            tmp0+=alph_err[j+1]*y[j][i];
        err_lcl+=tmp0*tmp0;
    }
    MPI_Allreduce(&err_lcl,&err_tot,1,MPI_TYPE0,MPI_SUM,world);
    err_tot=sqrt(err_tot/static_cast<type0>(dof_tot))/a_tol;
    return err_tot;
}


