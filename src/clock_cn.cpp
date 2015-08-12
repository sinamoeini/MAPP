#include <stdlib.h>
#include "clock_cn.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_cn::Clock_cn(MAPP* mapp,int narg
                   ,char** arg):Clock(mapp)
{
    min_del_t=1.0e-12;
    max_del_t=1.0e4;
    initial_del_t=-1.0;
    max_t=1.0e7;
    
    
    if(narg>2)
    {
        if(narg%2!=0)
            error->abort("every keyword in clock bdf should be followed by it's value");
        int iarg=2;
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
            else if(strcmp(arg[iarg],"max_step")==0)
            {
                iarg++;
                max_step=atoi(arg[iarg]);
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
                    error->abort("initial_del_t in clock bdf should be greater than 0.0");
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_t")==0)
            {
                iarg++;
                max_t=atof(arg[iarg]);
                iarg++;
            }
            else
                error->abort("unknown keyword in clock bdf: %s",arg[iarg]);
        }
    }
    
    if(min_gamma<0.0)
        error->abort("min_gamma in clock bdf should be greater than 0.0");
    if(gamma_red<=0.0 || gamma_red>=1.0)
        error->abort("gamma_red in clock bdf should be between 0.0 & 1.0");
    if(slope<=0.0 || slope>=1.0)
        error->abort("slope in clock bdf should be between 0.0 & 1.0");
    if(max_step<=0)
        error->abort("max_step in clock bdf should be greater than 0");
    if(max_iter<=0)
        error->abort("max_iter in clock bdf should be greater than 0");
    if(max_order<=0)
        error->abort("max_order in clock bdf should be greater than 0");
    if(m_tol<=0.0)
        error->abort("m_tol in clock bdf should be greater than 0.0");
    if(a_tol<=0.0)
        error->abort("a_tol in clock bdf should be greater than 0.0");
    if(min_del_t<=0.0)
        error->abort("min_del_t in clock bdf should be greater than 0.0");
    if(max_del_t<=0.0)
        error->abort("max_del_t in clock bdf should be greater than 0.0");
    if(max_t<=0.0)
        error->abort("max_t in clock bdf should be greater than 0.0");
    
    
    if(initial_del_t!=-1.0 && initial_del_t>max_t)
        error->abort("max_t in clock bdf should be greater than initial_del_t");
    if(min_del_t>max_t)
        error->abort("max_t in clock bdf should be greater than min_del_t");
    if(max_del_t>max_t)
        error->abort("max_t in clock bdf should be greater than max_del_t");
    if(min_del_t>=max_del_t)
        error->abort("max_del_t in clock bdf should be greater than min_del_t");

    

    
    CREATE1D(t,2);
    CREATE1D(dy,2);
    for(int i=0;i<2;i++)
        CREATE1D(dy[i],dof_lcl);
    
    CREATE1D(y,dof_lcl);

    

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_cn::~Clock_cn()
{    
    if(dof_lcl)
    {
        for(int i=0;i<2;i++)
            delete [] dy[i];
        
        delete [] y;
    }
    
    delete [] t;
    delete [] dy;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_cn::init()
{
    old_skin=atoms->skin;
    old_comm_mode=atoms->comm_mode;
    atoms->chng_skin(0.0);
    atoms->comm_mode=COMM_MODE_5;
    
    type0* c;
    type0* c_d;
    
    int f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<type0>(0,atoms->vectors[0]->dim,"f");
    
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
    atoms->init(vecs_comm);
    
    
    forcefield->create_2nd_neigh_lst_timer();
    forcefield->force_calc_timer(1,nrgy_strss);
    forcefield->c_d_calc_timer(0,nrgy_strss);
    
    
    thermo->update(fe_idx,nrgy_strss[0]);
    thermo->update(stress_idx,6,&nrgy_strss[1]);
    thermo->update(time_idx,0.0);
    
    thermo->init();
    
    if(write!=NULL)
        write->init();
    
    
    t[0]=0.0;
    
    atoms->vectors[c_n]->ret(c);
    atoms->vectors[c_d_n]->ret(c_d);
    
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
    
    t[1]=-initial_del_t;
    memcpy(dy[0],c_d,dof_lcl*sizeof(type0));
    memcpy(dy[1],c_d,dof_lcl*sizeof(type0));
    memcpy(y,c,dof_lcl*sizeof(type0));

    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_cn::fin()
{
    
    if(write!=NULL)
        write->fin();
    thermo->fin();
    atoms->fin();
    delete vecs_comm;
    
    atoms->chng_skin(old_skin);
    atoms->comm_mode=old_comm_mode;
}

/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_cn::run()
{
    if(max_step==0)
        return;
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    type0* tmp_dy;
    type0 del_t=initial_del_t,del_t_tmp,err=0.0,max_err;
    type0 ratio=1.0,cost;


    del_t=initial_del_t;
    int err_chk;
    int istep;
    
    istep=0;
    while (istep<max_step && t[0]<max_t)
    {
        err_chk=1;
        
        while (err_chk)
        {
            
            interpolate(del_t);
            solve_n_err(cost,err);
            
            if(err<1.0 && cost<1.0)
                err_chk=0;

            if(err_chk)
            {
                max_err=MAX(err,cost);
                ratio=pow(0.5/max_err,1.0/3.0);
                
                del_t_tmp=del_t*ratio;
                if(del_t_tmp<min_del_t)
                    del_t=min_del_t;
                else
                    del_t=del_t_tmp;

            }
        }
        
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step() || istep==max_step-1 || t[0]+del_t==max_t)
        {
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            thermo->update(time_idx,t[0]+del_t);
        }
        del_t_tmp=del_t;
        ord_dt(del_t,err);
        
        t[1]=t[0];
        t[0]+=del_t_tmp;
        tmp_dy=dy[1];
        dy[1]=dy[0];
        dy[0]=tmp_dy;
        memcpy(dy[0],c_d,dof_lcl*sizeof(type0));
        memcpy(y,c,dof_lcl*sizeof(type0));
        

        step_no++;
        istep++;
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_cn::ord_dt(type0& del_t,type0 err)
{
    type0 ratio,del_t_tmp0,del_t_tmp1;
    ratio=pow(0.5/err,1.0/3.0);
    
    if(ratio<0.5)
        ratio=0.5;
    else if(ratio>2.0)
        ratio=2.0;
    
    del_t_tmp0=del_t*ratio;
    del_t_tmp1=del_t;
    
    if(del_t_tmp0<min_del_t)
        del_t=min_del_t;
    else
        del_t=del_t_tmp0;
    
    if(del_t+del_t_tmp1+t[0]>max_t)
        del_t=max_t-t[0]-del_t_tmp1;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_cn::interpolate(type0& del_t)
{
    type0 c0,c1;
    
    type0* c;
    atoms->vectors[c_n]->ret(c);
    
    int idof;
    int err_chk=1;
    int err_chk_lcl;
    
    while(err_chk)
    {
        c0=del_t+0.5*del_t*del_t/(t[0]-t[1]);
        c1=-0.5*del_t*del_t/(t[0]-t[1]);
        
        err_chk_lcl=0;
        idof=0;
        while(idof<dof_lcl && err_chk_lcl==0)
        {
            if(c[idof]>=0.0)
            {
                y_0[idof]=y[idof]+c0*dy[0][idof]+c1*dy[1][idof];
                a[idof]=-y[idof]-0.5*del_t*dy[0][idof];
                if(y_0[idof]<0.0 || y_0[idof]>1.0)
                {
                    err_chk_lcl=1;
                }
            }
            else
                y_0[idof]=c[idof];
            
            
            idof++;
        }
        beta=del_t*0.5;
        
        MPI_Allreduce(&err_chk_lcl,&err_chk,1,MPI_INT,MPI_MAX,world);
        
        if(err_chk)
        {
            del_t*=0.8;
            del_t=MAX(del_t,min_del_t);
        }
        err_prefac=(1.0/3.0)*(del_t/(del_t+t[0]-t[1]));
        
    }
    
}
