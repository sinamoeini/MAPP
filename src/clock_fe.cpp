#include <stdlib.h>
#include "clock_fe.h"
#include "ff.h"
#include "error.h"
#include "memory.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_fe::Clock_fe(MAPP* mapp,int nargs
,char** args):ClockExplicit(mapp)
{
    if(nargs>2)
    {
        if(nargs%2!=0)
            error->abort("every keyword in clock adams should be followed by it's value");
        int iarg=2;
        while(iarg<nargs)
        {
            if(strcmp(args[iarg],"max_step")==0)
            {
                iarg++;
                max_step=atoi(args[iarg]);
                iarg++;
            }
            else if(strcmp(args[iarg],"max_t")==0)
            {
                iarg++;
                max_t=atof(args[iarg]);
                iarg++;
            }
            else if(strcmp(args[iarg],"a_tol")==0)
            {
                iarg++;
                a_tol=atof(args[iarg]);
                iarg++;
            }
            else if(strcmp(args[iarg],"min_del_t")==0)
            {
                iarg++;
                min_del_t=atof(args[iarg]);
                iarg++;
            }
            else if(strcmp(args[iarg],"initial_del_t")==0)
            {
                iarg++;
                initial_del_t=atof(args[iarg]);
                if(initial_del_t<=0.0)
                    error->abort("initial_del_t in clock fe should be greater than 0.0");
                iarg++;
            }
            else
                error->abort("unknown keyword in clock fe: %s",args[iarg]);
        }
    }
        

    if(max_step<0)
        error->abort("max_step in clock fe should be greater than 0");
    if(a_tol<=0.0)
        error->abort("a_tol in clock fe should be greater than 0.0");
    if(min_del_t<=0.0)
        error->abort("min_del_t in clock fe should be greater than 0.0");
    if(max_t<=0.0)
        error->abort("max_t in clock fe should be greater than 0.0");
    
    
    if(initial_del_t!=-1.0 && initial_del_t>max_t)
        error->abort("max_t in clock fe should be greater than initial_del_t");
    if(min_del_t>max_t)
        error->abort("max_t in clock fe should be greater than min_del_t");

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_fe::~Clock_fe()
{

    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void Clock_fe::allocate()
{
    vecs_1=new Vec<type0>*[2];
    vecs_1[0]=new Vec<type0>(atoms,c_dim);
    vecs_1[1]=new Vec<type0>(atoms,c_dim);
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Clock_fe::reset()
{
    y=vecs_1[0]->begin();
    dy=vecs_1[1]->begin();
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void Clock_fe::deallocate()
{
    delete vecs_1[1];
    delete vecs_1[0];
    delete [] vecs_1;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_fe::init()
{
    ClockExplicit::init();
    allocate();
    
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    rectify(c_d);

    
    if(initial_del_t<0.0)
    {
        type0 sum=forcefield_dmd->ddc_norm_timer()/sqrt(static_cast<type0>(dof_tot));
        initial_del_t=MIN(sqrt(4.0*a_tol/sum),1.0e-3*max_t);
    }
    
    init_stp_adj(initial_del_t);
    
    
    memcpy(y,c,dof_lcl*sizeof(type0));
    memcpy(dy,c_d,dof_lcl*sizeof(type0));
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_fe::fin()
{
    ClockExplicit::fin();
    deallocate();
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void Clock_fe::run()
{
    if(max_step==0)
        return;
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    
    type0 del_t=initial_del_t,del_t_tmp;
    type0 err=0.0;
    int err_chk;
    int istep;
    
    istep=0;
    while (istep<max_step && tot_t<max_t)
    {
        
        err_chk=1;
        while (err_chk)
        {
            interpolate_n_err(err,del_t);
            
            if(err<1.0)
                err_chk=0;
            fail_stp_adj(err,del_t);
        }
        
        max_succ_dt=MAX(max_succ_dt,del_t);
        tot_t+=del_t;

        
        
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step()|| istep==max_step-1 || tot_t>=max_t)
        {
            forcefield_dmd->force_calc_timer(true);
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            thermo->update(time_idx,tot_t);
        }
        
        del_t_tmp=del_t;
        ord_dt(del_t,err);
        
        memcpy(y,c,dof_lcl*sizeof(type0));
        memcpy(dy,c_d,dof_lcl*sizeof(type0));
        
        istep++;
        step_no++;
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void Clock_fe::interpolate_n_err(type0& err,type0& del_t)
{
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    
    int idof;
    int err_chk=1;
    type0 r,r_lcl,err_lcl;
    while (err_chk)
    {
        r=r_lcl=1.0;
        idof=0;
        while(idof<dof_lcl && r_lcl==1.0)
        {
            if(y[idof]>=0.0)
            {
                c[idof]=y[idof]+0.5*del_t*dy[idof];
                
                if(c[idof]<0.0 ||c[idof]>1.0)
                {
                    while (c[idof]<0.0)
                    {
                        r_lcl=-y[idof]/(0.5*(r_lcl*del_t));
                        c[idof]=y[idof]+0.5*(r_lcl*del_t)*dy[idof];
                    }
                    
                    while (c[idof]>1.0)
                    {
                        r_lcl=(1.0-y[idof])/(0.5*(r_lcl*del_t));
                        c[idof]=y[idof]+0.5*(r_lcl*del_t)*dy[idof];
                    }
                }
            }
            else
                c[idof]=y[idof];
            
            idof++;
        }
        MPI_Allreduce(&r_lcl,&r,1,MPI_TYPE0,MPI_MIN,world);
        if(r!=1.0)
        {
            if(max_t-tot_t<=2.0*min_del_t)
            {
                error->abort("reached minimum del_t (%e)",del_t);
            }
            else
            {
                if(del_t==min_del_t)
                {
                    error->abort("reached minimum del_t (%e)",del_t);
                }
                else
                {
                    if(r*del_t<min_del_t)
                        del_t=min_del_t;
                    else if(r*del_t>max_t-tot_t-min_del_t)
                        del_t=max_t-tot_t-min_del_t;
                    else
                        del_t*=r;
                }
            }
            intp_rej++;
            continue;
        }
        
        atoms->update(mapp->c);
        forcefield_dmd->dc_timer();
        
        err_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
        {
            if(c[i]>=0.0)
                err_lcl+=(c_d[i]-dy[i])*(c_d[i]-dy[i]);
        }
        
        err=0.0;
        MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
        err=0.5*del_t*sqrt(err/static_cast<type0>(dof_tot))/a_tol;
        
        if(err>=1.0)
        {
            intg_rej++;
            return;
        }
        r=r_lcl=1.0;
        idof=0;
        while(idof<dof_lcl && r_lcl==1.0)
        {
            if(y[idof]>=0.0)
            {
                c[idof]=y[idof]+0.5*del_t*(dy[idof]+c_d[idof]);
                
                if(c[idof]<0.0 ||c[idof]>1.0)
                {
                    while (c[idof]<0.0)
                    {
                        r_lcl=-y[idof]/(0.5*(del_t*r_lcl));
                        c[idof]=y[idof]+0.5*(r_lcl*del_t)*(dy[idof]+c_d[idof]);
                    }
                    
                    while (c[idof]>1.0)
                    {
                        r_lcl=(1.0-y[idof])/(0.5*(del_t*r_lcl));
                        c[idof]=y[idof]+0.5*(r_lcl*del_t)*(dy[idof]+c_d[idof]);
                    }
                }
            }
            else
                c[idof]=y[idof];
            
            idof++;
        }
        MPI_Allreduce(&r_lcl,&r,1,MPI_TYPE0,MPI_MIN,world);
        if(r!=1.0)
        {
            if(max_t-tot_t<=2.0*min_del_t)
            {
                error->abort("reached minimum del_t (%e)",del_t);
            }
            else
            {
                if(del_t==min_del_t)
                {
                    error->abort("reached minimum del_t (%e)",del_t);
                }
                else
                {
                    if(r*del_t<min_del_t)
                        del_t=min_del_t;
                    else if(r*del_t>max_t-tot_t-min_del_t)
                        del_t=max_t-tot_t-min_del_t;
                    else
                        del_t*=r;
                }
            }
            intp_rej++;
        }
        else
        {
            err_chk=0;
            atoms->update(mapp->c);
            forcefield_dmd->dc_timer();
            
        }
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void Clock_fe::ord_dt(type0& del_t,type0 err)
{
    type0 r;
    r=0.9/err;
    
    if(r>=2.0)
        r=2.0;
    else if(r<=0.5)
        r=0.5;
    else
        r=1.0;
    
    if(r*del_t>max_t-tot_t)
        del_t=max_t-tot_t;
    else
    {
        if(max_t-tot_t<=2.0*min_del_t)
        {
            del_t=2.0*min_del_t;
        }
        else
        {
            if(r*del_t<min_del_t)
                del_t=min_del_t;
            else if(r*del_t>=max_t-tot_t-min_del_t)
                del_t=max_t-tot_t-min_del_t;
            else
            {
                del_t*=r;
            }
        }
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
inline void Clock_fe::init_stp_adj(type0& del_t)
{
    if(del_t>max_t-tot_t)
        del_t=max_t-tot_t;
    else
    {
        if(max_t-tot_t<=2.0*min_del_t)
        {
            del_t=max_t-tot_t;
        }
        else
        {
            if(del_t<min_del_t)
                del_t=min_del_t;
            else if(del_t>=max_t-tot_t-min_del_t)
                del_t=max_t-tot_t-min_del_t;
        }
    }
}
/*--------------------------------------------
 step addjustment after failure
 --------------------------------------------*/
inline void Clock_fe::fail_stp_adj(type0 err,type0& del_t)
{
    if(max_t-tot_t<=2.0*min_del_t)
    {
        error->abort("reached minimum order & del_t (%e)",del_t);
    }
    else
    {
        if(del_t==min_del_t)
        {
            
            error->abort("reached minimum order & del_t (%e)",del_t);
        }
        else
        {
            type0 r=0.9/err;
            r=MAX(r,0.5);            

            if(r*del_t<min_del_t)
                del_t=min_del_t;
            else if(r*del_t>max_t-tot_t-min_del_t)
                del_t=max_t-tot_t-min_del_t;
            else
                del_t*=r;
        }
    }
}

