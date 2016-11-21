#include <stdlib.h>
#include <limits>
#include "atoms.h"
#include "dmd_bdf_y.h"
#include "ff.h"
#include "error.h"
#include "memory.h"
#include "write.h"
#include "thermo_dynamics.h"
#include "cmd.h"
#include "MAPP.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMD_bdf_y::DMD_bdf_y(int nargs
,char** args):DMDImplicit()
{
    q_max=5;
    
    char* dmd_style=NULL;
    
    Pattern cmd;
    
    /*----------------------------*/
    cmd.cmd("dmd");
    cmd.add_var(dmd_style,"style");
    cmd.add_vdesc(0,"defines the style of dmd");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=VLogics("eq","bdf");
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("max_step");
    cmd.add_var(max_step,"nstep");
    cmd.add_vdesc(0,"defines maximum number of steps");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=VLogics("gt",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("max_iter");
    cmd.add_var(max_iter,"niter");
    cmd.add_vdesc(0,"defines maximum number of iterations");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=VLogics("gt",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("q_max");
    cmd.add_var(q_max,"q");
    cmd.add_vdesc(0,"defines maximum order of bdf");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=VLogics("gt",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("a_tol");
    cmd.add_var(a_tol,"tol");
    cmd.add_vdesc(0,"defines absolute tolerance in local truncation error for performing integration");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=VLogics("gt",0.0)*VLogics("le",5);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("dt_min");
    cmd.add_var(dt_min,"dt");
    cmd.add_vdesc(0,"defines minimum time step");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=VLogics("gt",0.0);
    /*------------------------------------------------------------------------------------*/

    cmd.scan(args,nargs);
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
DMD_bdf_y::~DMD_bdf_y()
{
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_bdf_y::allocate()
{
    vecs_1=new Vec<type0>*[q_max+2];
    for(int ivec=0;ivec<q_max+2;ivec++)
        vecs_1[ivec]=new Vec<type0>(atoms,c_dim);
    
    CREATE1D(alpha_y,q_max+1);
    CREATE1D(dalpha_y,q_max+1);
    CREATE1D(alph_err,q_max+2);
    CREATE1D(t,q_max+1);
    CREATE1D(y,q_max+1);
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_bdf_y::deallocate()
{
    
    delete [] y;
    delete [] t;
    delete [] alph_err;
    delete [] dalpha_y;
    delete [] alpha_y;

    for(int ivec=0;ivec<q_max+2;ivec++)
        delete vecs_1[ivec];
    delete [] vecs_1;
}
/*--------------------------------------------
 restart a simulation
 --------------------------------------------*/
inline void DMD_bdf_y::restart()
{
    reset();
    for(int i=0;i<q_max+1;i++)
        y[i]=vecs_1[i]->begin();
    dy=vecs_1[q_max+1]->begin();
}
/*--------------------------------------------
 store the vectors
 --------------------------------------------*/
inline void DMD_bdf_y::start()
{
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    type0 sum,sum_lcl;
    sum_lcl=0.0;
    for(int i=0;i<ncs;i++)
        if(c[i]>=0.0)
            sum_lcl+=c_d[i]*c_d[i];
    sum=0.0;
    MPI_Allreduce(&sum_lcl,&sum,1,MPI_TYPE0,MPI_SUM,world);
    sum=sqrt(sum/nc_dofs);
    dt=MIN(2.0*a_tol/sum,1.0e-3*(t_fin-t_cur));
    
    
    if(dt>t_fin-t_cur)
        dt=t_fin-t_cur;
    else
    {
        if(t_fin-t_cur<=2.0*dt_min)
        {
            dt=t_fin-t_cur;
        }
        else
        {
            if(dt<dt_min)
                dt=dt_min;
            else if(dt>=t_fin-t_cur-dt_min)
                dt=t_fin-t_cur-dt_min;
        }
    }
    
    memcpy(y[0],mapp->c->begin(),ncs*sizeof(type0));
    memcpy(y[1],mapp->c->begin(),ncs*sizeof(type0));
    memcpy(dy,mapp->c_d->begin(),ncs*sizeof(type0));
    for(int i=0;i<q_max+1;i++) t[i]=0.0;
    t[0]=0.0;
    t[1]=-dt;
    q=1;
}
/*--------------------------------------------
 store the vectors
 --------------------------------------------*/
inline void DMD_bdf_y::update_for_next()
{
    
    type0* tmp_y=y[q_max];
    for(int i=q_max;i>0;i--)
    {
        t[i]=t[i-1]-dt;
        y[i]=y[i-1];
    }
    y[0]=tmp_y;
    memcpy(y[0],mapp->c->begin(),ncs*sizeof(type0));
    memcpy(dy,mapp->c_d->begin(),ncs*sizeof(type0));
    
    dt=dt_new;
    q+=dq;
}
/*--------------------------------------------
 interpolation fail
 --------------------------------------------*/
void DMD_bdf_y::interpolate_fail()
{
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    type0 sum,sum_lcl;
    sum_lcl=0.0;
    for(int i=0;i<ncs;i++)
        if(c[i]>=0.0)
            sum_lcl+=c_d[i]*c_d[i];
    sum=0.0;
    MPI_Allreduce(&sum_lcl,&sum,1,MPI_TYPE0,MPI_SUM,world);
    sum=sqrt(sum/nc_dofs);
    dt=MIN(2.0*a_tol/sum,1.0e-3*(t_fin-t_cur));
    
    
    if(dt>t_fin-t_cur)
        dt=t_fin-t_cur;
    else
    {
        if(t_fin-t_cur<=2.0*dt_min)
        {
            dt=t_fin-t_cur;
        }
        else
        {
            if(dt<dt_min)
                dt=dt_min;
            else if(dt>=t_fin-t_cur-dt_min)
                dt=t_fin-t_cur-dt_min;
        }
    }
    
    memcpy(y[0],mapp->c->begin(),ncs*sizeof(type0));
    memcpy(y[1],mapp->c->begin(),ncs*sizeof(type0));
    memcpy(y_0,mapp->c->begin(),ncs*sizeof(type0));
    memcpy(dy,mapp->c_d->begin(),ncs*sizeof(type0));
    for(int i=0;i<q_max+1;i++) t[i]=0.0;
    t[0]=0.0;
    t[1]=-dt;
    q=1;
    
    beta=dt;
    beta_inv=1.0/dt;
    
    for(int i=0;i<ncs;i++)
        if(c[i]>=0.0)
            a[i]-=y_0[i];

}
/*--------------------------------------------
 init
 --------------------------------------------*/
inline bool DMD_bdf_y::interpolate()
{
    type0 k0=0.0,tmp0;
    
    
    for(int i=0;i<q+1;i++)
    {
        k0+=1.0/(1.0-t[i]/dt);
        tmp0=1.0;
        for(int j=0;j<q+1;j++)
            if(i!=j)
                tmp0*=(dt-t[j])/(t[i]-t[j]);
        
        alpha_y[i]=tmp0;
    }
    
    for(int i=0;i<q+1;i++)
    {
        tmp0=k0-1.0/(1.0-t[i]/dt);
        dalpha_y[i]=alpha_y[i]*tmp0/dt;
    }
    
    tmp0=0.0;
    for(int i=0;i<q;i++)
        tmp0+=1.0/static_cast<type0>(i+1);
    beta=dt/tmp0;
    beta_inv=tmp0/dt;
    
    bool xcd_dom=false;
    type0* c=mapp->c->begin();
    
    for(int i=0;i<ncs && !xcd_dom;i++)
    {
        if(c[i]>=0.0)
        {
            y_0[i]=0.0;
            a[i]=0.0;
            for(int j=0;j<q+1;j++)
            {
                a[i]+=dalpha_y[j]*y[j][i];
                y_0[i]+=alpha_y[j]*y[j][i];
            }
            a[i]*=beta;
            a[i]-=y_0[i];
            
            if(y_0[i]<0.0 || y_0[i]>1.0)
            {
                xcd_dom=true;
            }
        }
        else
            y_0[i]=c[i];
    }
    int err_dom_lcl=0,err_dom;
    if(xcd_dom)
        err_dom_lcl=1;
    MPI_Allreduce(&err_dom_lcl,&err_dom,1,MPI_INT,MPI_MAX,world);
    
    if(err_dom)
        return false;
    
    return true;

}
/*--------------------------------------------
 step addjustment after success
 --------------------------------------------*/
inline void DMD_bdf_y::err_fac_calc()
{
    type0 err_prefac0=dt/(dt-t[q]);
    type0 err_prefac1=err_prefac0;
    for(int i=0;i<q;i++)
    {
        err_prefac1+=dt/(dt-t[i])
        -1.0/static_cast<type0>(i+1);
    }
    err_fac=MAX(err_prefac0,fabs(err_prefac1));
}
/*--------------------------------------------
 step addjustment after success
 --------------------------------------------*/
inline void DMD_bdf_y::ord_dt(type0& r)
{
    dq=0;
    if(intg_rej==0 && intp_rej==0 && solve_rej==0)
    {
        r=2.0;
        if(q<q_max) dq=1;
        return;
    }
    
    
    type0 terkm2=0.0,terkm1=0.0,terk=0.0,terkp1=0.0;
    type0 est;
    
    bool terkm2_flag,terkm1_flag,terkp1_flag;
    
    terkp1_flag=false;
    terkm1_flag=false;
    terkm2_flag=false;
    
    if(q>2)
        terkm2_flag=true;
    if(q>1)
        terkm1_flag=true;
    
    if(const_q>q+1 && const_dt>q+1 && q<q_max)
        terkp1_flag=true;
    
    
    if(terkm2_flag)
        terkm2=err_est(q-1);
    
    if(terkm1_flag)
        terkm1=err_est(q);
    
    terk=err_est(q+1);
    
    if(terkp1_flag)
        terkp1=err_est(q+2);
    
    if(q>2)
    {
        if(MAX(terkm1,terkm2)<=terk)
        {
            est=terkm1/static_cast<type0>(q);
            dq=-1;
        }
        else
        {
            if(terkp1<terk && terkp1_flag)
            {
                est=terkp1/static_cast<type0>(q+2);
                dq=1;
            }
            else
                est=err;
        }
    }
    else if(q==2)
    {
        if(terkm1<=0.5*terk)
        {
            est=terkm1/static_cast<type0>(q);
            dq=-1;
        }
        else if(terkp1_flag && terk<terkm1)
        {
            if(terkp1<terk)
            {
                est=terkp1/static_cast<type0>(q+2);
                dq=1;
            }
            else
                est=terk/static_cast<type0>(q+1);
        }
        else
            est=terk/static_cast<type0>(q+1);
    }
    else
    {
        if(terkp1_flag && terkp1<0.5*terk)
        {
            est=terkp1/static_cast<type0>(q+2);
            dq=1;
        }
        else
            est=terk/static_cast<type0>(q+1);
    }
    
    r=pow(0.5/err,1.0/static_cast<type0>(q+1));
    if(dq!=0)
    {
        r*=pow(est,1.0/static_cast<type0>(q+dq+1))/pow(terk/static_cast<type0>(q+1),1.0/static_cast<type0>(q+1));
    }
    
    r=pow(0.5/est,1.0/static_cast<type0>(q+dq+1));
    
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
inline type0 DMD_bdf_y::err_est(int q_)
{
    type0 tmp0,err_lcl,err,k0;
    type0* c=mapp->c->begin();
    
    k0=1.0;
    for(int i=0;i<q_;i++)
        k0*=static_cast<type0>(i+1)/(1.0-t[i]/dt);
    
    for(int i=0;i<q_;i++)
    {
        tmp0=1.0;
        for(int j=0;j<q_;j++)
            if(j!=i)
                tmp0*=(dt-t[j])/(t[i]-t[j]);
        
        alph_err[i]=-tmp0;
    }
    
    err_lcl=0.0;
    for(int i=0;i<ncs;i++)
    {
        if(c[i]>=0.0)
        {
            tmp0=c[i];
            for(int j=0;j<q_;j++)
                tmp0+=alph_err[j]*y[j][i];
            err_lcl+=tmp0*tmp0;
        }
    }
    MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
    err=sqrt(err/nc_dofs)/a_tol;
    err*=fabs(k0);
    return err;
}
/*--------------------------------------------
 error calculation
 --------------------------------------------*/
void DMD_bdf_y::err_calc()
{
    type0 tmp0,err_lcl=0.0,max_dy_lcl=0.0,max_dy;
    type0 c_d_norm_lcl=0.0;
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    for(int i=0;i<ncs;i++)
    {
        if(c[i]>=0.0)
        {
            c_d_norm_lcl+=c_d[i]*c_d[i];
            tmp0=(y_0[i]-c[i]);
            max_dy_lcl=MAX(max_dy_lcl,fabs(c_d[i]));
            err_lcl+=tmp0*tmp0;
        }
    }
    MPI_Allreduce(&max_dy_lcl,&max_dy,1,MPI_TYPE0,MPI_MAX,world);
    MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
    MPI_Allreduce(&c_d_norm_lcl,&c_d_norm,1,MPI_TYPE0,MPI_SUM,world);
    c_d_norm=sqrt(c_d_norm/nc_dofs)/a_tol;
    err=sqrt(err/nc_dofs)/a_tol;
    err*=err_fac;
}

