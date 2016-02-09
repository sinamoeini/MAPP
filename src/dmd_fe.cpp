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
            else if(strcmp(args[iarg],"dt_min")==0)
            {
                iarg++;
                dt_min=atof(args[iarg]);
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
    if(dt_min<=0.0)
        error->abort("dt_min in dmd fe should be greater than 0.0");

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
 destructor
 --------------------------------------------*/
void DMD_fe::deallocate()
{
    delete vecs_1[1];
    delete vecs_1[0];
    delete [] vecs_1;
}
/*--------------------------------------------
 restart a simulation
 --------------------------------------------*/
void DMD_fe::restart(type0& del_t,int& q)
{
    reset();
    y=vecs_1[0]->begin();
    dy=vecs_1[1]->begin();
    
    type0 sum=forcefield_dmd->ddc_norm_timer()/sqrt(nc_dofs);
    del_t=MIN(sqrt(4.0*a_tol/sum),1.0e-3*(t_fin-t_cur));
    
    if(del_t>t_fin-t_cur)
        del_t=t_fin-t_cur;
    else
    {
        if(t_fin-t_cur<=2.0*dt_min)
        {
            del_t=t_fin-t_cur;
        }
        else
        {
            if(del_t<dt_min)
                del_t=dt_min;
            else if(del_t>=t_fin-t_cur-dt_min)
                del_t=t_fin-t_cur-dt_min;
        }
    }
    
    memcpy(y,mapp->c->begin(),ncs*sizeof(type0));
    memcpy(dy,mapp->c_d->begin(),ncs*sizeof(type0));
}
/*--------------------------------------------
 store the vectors
 --------------------------------------------*/
void DMD_fe::store_vecs(type0 del_t)
{
    memcpy(y,mapp->c->begin(),ncs*sizeof(type0));
    memcpy(dy,mapp->c_d->begin(),ncs*sizeof(type0));
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
            if(t_fin-t_cur<=2.0*dt_min)
            {
                error->abort("reached minimum del_t (%e)",del_t);
            }
            else
            {
                if(del_t==dt_min)
                {
                    error->abort("reached minimum del_t (%e)",del_t);
                }
                else
                {
                    if(r*del_t<dt_min)
                        del_t=dt_min;
                    else if(r*del_t>t_fin-t_cur-dt_min)
                        del_t=t_fin-t_cur-dt_min;
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
            if(t_fin-t_cur<=2.0*dt_min)
            {
                error->abort("reached minimum del_t (%e)",del_t);
            }
            else
            {
                if(del_t==dt_min)
                {
                    error->abort("reached minimum del_t (%e)",del_t);
                }
                else
                {
                    if(r*del_t<dt_min)
                        del_t=dt_min;
                    else if(r*del_t>t_fin-t_cur-dt_min)
                        del_t=t_fin-t_cur-dt_min;
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
    
    if(r*del_t>t_fin-t_cur)
        del_t=t_fin-t_cur;
    else
    {
        if(t_fin-t_cur<=2.0*dt_min)
        {
            del_t=2.0*dt_min;
        }
        else
        {
            if(r*del_t<dt_min)
                del_t=dt_min;
            else if(r*del_t>=t_fin-t_cur-dt_min)
                del_t=t_fin-t_cur-dt_min;
            else
            {
                del_t*=r;
            }
        }
    }
}
