#include <stdlib.h>
#include "dmd_fe.h"
#include "ff.h"
#include "error.h"
#include "memory.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMD_fe::DMD_fe(MAPP* mapp,int nargs
,char** args):DMDExplicit(mapp)
{
    if(nargs>2)
    {
        if(nargs%2!=0)
            error->abort("every keyword in dmd adams should be followed by it's value");
        int iarg=2;
        while(iarg<nargs)
        {
            if(strcmp(args[iarg],"max_step")==0)
            {
                iarg++;
                max_step=atoi(args[iarg]);
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
                    error->abort("initial_del_t in dmd fe should be greater than 0.0");
                iarg++;
            }
            else
                error->abort("unknown keyword in dmd fe: %s",args[iarg]);
        }
    }
        

    if(max_step<0)
        error->abort("max_step in dmd fe should be greater than 0");
    if(a_tol<=0.0)
        error->abort("a_tol in dmd fe should be greater than 0.0");
    if(min_del_t<=0.0)
        error->abort("min_del_t in dmd fe should be greater than 0.0");

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
DMD_fe::~DMD_fe()
{

    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD_fe::allocate()
{
    vecs_1=new Vec<type0>*[2];
    vecs_1[0]=new Vec<type0>(atoms,c_dim);
    vecs_1[1]=new Vec<type0>(atoms,c_dim);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD_fe::reset()
{
    DMDExplicit::reset();
    y=vecs_1[0]->begin();
    dy=vecs_1[1]->begin();
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD_fe::deallocate()
{
    delete vecs_1[1];
    delete vecs_1[0];
    delete [] vecs_1;
}
/*--------------------------------------------

 --------------------------------------------*/
type0 DMD_fe::est_dt()
{
    type0 sum=forcefield_dmd->ddc_norm_timer()/sqrt(nc_dofs);
    type0 del_t=MIN(sqrt(4.0*a_tol/sum),1.0e-3*(max_t-tot_t));
    init_stp_adj(del_t);
    return del_t;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_fe::init()
{
    DMDExplicit::init();
    allocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_fe::fin()
{
    DMDExplicit::fin();
    deallocate();
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void DMD_fe::run()
{
    if(max_step==0)
        return;
    
    type0 del_t,del_t_tmp;
    type0 err=0.0;
    int istep;
    bool min_run=false;
    istep=0;
    while (istep<max_step && tot_t<max_t)
    {
        reset();
        type0* c=mapp->c->begin();
        type0* c_d=mapp->c_d->begin();
        del_t=est_dt();
        
        memcpy(y,c,ncs*sizeof(type0));
        memcpy(dy,c_d,ncs*sizeof(type0));
        
        while (istep<max_step && tot_t<max_t && !min_run)
        {
            err=1.0;
            while(err>=1.0)
            {
                interpolate_n_err(err,del_t);
                if(err<1.0)
                    continue;
                
                fail_stp_adj(err,del_t);
                memcpy(c,y,ncs*sizeof(type0));
                memcpy(c_d,dy,ncs*sizeof(type0));
            }
            
            min_run=decide_min(istep,del_t);
            if(min_run) continue;
            
            del_t_tmp=del_t;
            ord_dt(del_t,err);
            
            memcpy(y,c,ncs*sizeof(type0));
            memcpy(dy,c_d,ncs*sizeof(type0));
        }
        
        if(!min_run) continue;
        do_min();
        min_run=false;
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void DMD_fe::interpolate_n_err(type0& err,type0& del_t)
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
        while(idof<ncs && r_lcl==1.0)
        {
            if(y[idof]>=0.0)
            {
                c[idof]=y[idof]+0.5*del_t*dy[idof];
                
                if(c[idof]<0.0 || c[idof]>1.0)
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
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
                err_lcl+=(c_d[i]-dy[i])*(c_d[i]-dy[i]);
        }
        
        err=0.0;
        MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
        err=0.5*del_t*sqrt(err/nc_dofs)/a_tol;
        
        if(err>=1.0)
        {
            intg_rej++;
            return;
        }
        r=r_lcl=1.0;
        idof=0;
        while(idof<ncs && r_lcl==1.0)
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
void DMD_fe::ord_dt(type0& del_t,type0 err)
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
inline void DMD_fe::init_stp_adj(type0& del_t)
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
inline void DMD_fe::fail_stp_adj(type0 err,type0& del_t)
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
