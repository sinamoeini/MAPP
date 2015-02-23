#include "clock_fe.h"
#include "neighbor.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
#define SMALL 1.0e-18
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_fe::Clock_fe(MAPP* mapp,int narg
                   ,char** arg):Clock(mapp)
{
    e_tol=0.0;
    a_tol=1.0e-6;
    min_del_t=1.0e-12;
    max_del_t=1.0e4;
    if(narg<3)
        error->abort("clock fe should at least have 1 arguement");
    
    no_steps=atoi(arg[2]);
    if(no_steps<0)
        error->abort("no of steps for clock fe should be greater than 0");
    
    if(narg>3)
    {
        if((narg-3)%2!=0)
            error->abort("every keyword in clock fe should be followed by it's value");
        int iarg=3;
        while(iarg<narg)
        {
            
            if(strcmp(arg[iarg],"a_tol")==0)
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
            else
                error->abort("unknown keyword in clock fe: %s",arg[iarg]);
        }
    }
    
    
    if(a_tol<=0.0)
        error->abort("a_tol in clock fe should be greater than 0.0");
    if(e_tol<0.0)
        error->abort("e_tol in clock fe should be equal or greater than 0.0");
    if(min_del_t<=0.0)
        error->abort("min_del_t in clock fe should be greater than 0.0");
    if(max_del_t<=0.0)
        error->abort("max_del_t in clock fe should be greater than 0.0");
    if(min_del_t>=max_del_t)
        error->abort("max_del_t in clock fe should be greater than min_del_t");
    
    
    
    
    c_n=atoms->find("c");
    c_d_n=atoms->find("c_d");
    
    int c_dim=atoms->vectors[c_n]->dim;
    dof_lcl=atoms->natms*c_dim;
    MPI_Allreduce(&dof_lcl,&dof_tot,1,MPI_INT,MPI_SUM,world);
    
    CREATE1D(y,dof_lcl);
    CREATE1D(dy,dof_lcl);
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_fe::~Clock_fe()
{
    if(dof_lcl)
    {
        delete [] y;
        delete [] dy;
    }
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_fe::init()
{
    int f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<type0>(0,atoms->vectors[0]->dim,"f");
    
    c_n=atoms->find_exist("c");
    c_d_n=atoms->find("c_d");
    
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

    forcefield->create_2nd_neigh_lst();
    
    
    
    forcefield->force_calc(1,nrgy_strss);
    
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    
    
    
    thermo->start_force_time();
    forcefield->c_d_calc();
    rectify(c_d);
    thermo->stop_force_time();
    
    memcpy(y,c,dof_lcl*sizeof(type0));
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
void Clock_fe::fin()
{
    if(write!=NULL)
        write->fin();
    
    forcefield->fin();
    neighbor->fin();
    
    delete vecs_comm;
    
    thermo->fin();
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void Clock_fe::run()
{
    if(no_steps==0)
        return;
    type0 curr_t=0.0;
    
    type0 del_t;
    
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    del_t=min_del_t;
    eq_ratio=1.0;
    int istep=0;
    
    while (eq_ratio>=1.0 && istep <no_steps)
    {
        solve(del_t);
        
        
        curr_t+=del_t;
        
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
            thermo->update(time_idx,curr_t);
        }
        
        
        ord_dt(del_t);

        
        thermo->start_force_time();
        forcefield->c_d_calc();
        rectify(c_d);
        thermo->stop_force_time();
        
        memcpy(y,c,dof_lcl*sizeof(type0));
        memcpy(dy,c_d,dof_lcl*sizeof(type0));
        
        istep++;
        step_no++;
    }
    
}

/*--------------------------------------------
 run
 --------------------------------------------*/
void Clock_fe::solve(type0& del_t)
{
    
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    
    type0 ratio,tot_ratio,tmp0,tmp1,err_lcl;
    ratio=1.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp0=y[i]+dy[i]*del_t;
        if(tmp0>1.0)
            ratio=MIN((1.0-y[i])/(tmp0-y[i]),ratio);
        else if(tmp0<0.0)
            ratio=MIN((y[i]-0.0)/(y[i]-tmp0),ratio);
    }
    
    MPI_Allreduce(&ratio,&tot_ratio,1,MPI_TYPE0,MPI_MIN,world);
    del_t*=ratio;
    
    err=1.0;
    while (err>=1.0)
    {
        for(int i=0;i<dof_lcl;i++)
            c[i]=y[i]+0.5*del_t*dy[i];
        
        thermo->start_comm_time();
        atoms->update_ph(c_n);
        thermo->stop_comm_time();
        
        thermo->start_force_time();
        forcefield->c_d_calc();
        rectify(c_d);
        thermo->stop_force_time();
        
        err_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
        {
            tmp0=c_d[i]-dy[i];
            err_lcl+=tmp0*tmp0;
        }
        err=0.0;
        MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
        err=0.5*del_t*sqrt(err/static_cast<type0>(dof_tot))/a_tol;
        
        
        if(del_t==min_del_t && err>=1.0)
            err=0.0;
        
        if(err>=1.0)
        {
            del_t*=0.9/err;
            if (ratio<=0.5)
                ratio=0.5;
            del_t*=ratio;
            del_t=MAX(del_t,min_del_t);
        }
    }

    tmp1=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp1+=c_d[i]*c_d[i];
        c[i]=y[i]+del_t*dy[i];
    }

    thermo->start_comm_time();
    atoms->update_ph(c_n);
    thermo->stop_comm_time();
    
    MPI_Allreduce(&tmp1,&eq_ratio,1,MPI_TYPE0,MPI_SUM,world);
    eq_ratio=sqrt(eq_ratio/static_cast<type0>(dof_tot))/e_tol;
 

}
/*--------------------------------------------
 run
 --------------------------------------------*/
void Clock_fe::ord_dt(type0& del_t)
{
    type0 ratio;
    ratio=0.9/err;
    
    if(ratio>=2.0)
        ratio=2.0;
    else if (ratio<=0.5)
        ratio=0.5;
    
    del_t*=ratio;
    del_t=MAX(del_t,min_del_t);
    del_t=MIN(del_t,max_del_t);
}

