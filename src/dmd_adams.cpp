#include <stdlib.h>
#include <limits>
#include "dmd_adams.h"
#include "ff.h"
#include "error.h"
#include "memory.h"
#include "write.h"
#include "thermo_dynamics.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMD_adams::DMD_adams(MAPP* mapp,int nargs
,char** args):DMDImplicit(mapp)
{
    q_max=12;
    
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
            else if(strcmp(args[iarg],"max_iter")==0)
            {
                iarg++;
                max_iter=atoi(args[iarg]);
                iarg++;
            }
            else if(strcmp(args[iarg],"q_max")==0)
            {
                iarg++;
                q_max=atoi(args[iarg]);
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
                error->abort("unknown keyword in dmd adams: %s",args[iarg]);
        }
    }

    if(max_step<0)
        error->abort("max_step in dmd adams should be greater than 0");
    if(max_iter<=0)
        error->abort("max_iter in dmd adams should be greater than 0");
    if(q_max<=0)
        error->abort("q_max in dmd adams should be greater than 0");
    if(a_tol<=0.0)
        error->abort("a_tol in dmd adams should be greater than 0.0");
    if(dt_min<=0.0)
        error->abort("dt_min in dmd adams should be greater than 0.0");
    
    int n=1+q_max/2;
    CREATE1D(xi,n);
    CREATE1D(wi,n);
    for(int i=0;i<n;i++)
    {
        CREATE1D(xi[i],i+1);
        CREATE1D(wi[i],i+1);
    }
    
    XMath* xmath= new XMath();
    for(int i=0;i<n;i++)
    {
        xmath->quadrature_lg(i+1,xi[i],wi[i]);
        for(int j=0;j<i+1;j++)
        {
            xi[i][j]*=0.5;
            xi[i][j]+=0.5;
            wi[i][j]*=0.5;
        }
    }
    delete xmath;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
DMD_adams::~DMD_adams()
{

    int n=1+q_max/2;
    for(int i=0;i<n;i++)
    {
        delete [] xi[i];
        delete [] wi[i];
    }
    if(n)
    {
        delete [] xi;
        delete [] wi;
    }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD_adams::allocate()
{
    vecs_1=new Vec<type0>*[q_max+2];
    for(int ivec=0;ivec<q_max+2;ivec++)
        vecs_1[ivec]=new Vec<type0>(atoms,c_dim);
    
    CREATE1D(t,q_max);
    CREATE1D(dy,q_max);
    CREATE1D(alpha_dy,q_max);
    CREATE1D(dalpha_dy,q_max);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD_adams::deallocate()
{
    if(q_max)
    {
        delete [] dy;
        delete [] t;
        delete [] alpha_dy;
        delete [] dalpha_dy;
    }
    
    for(int ivec=0;ivec<q_max+2;ivec++)
        delete vecs_1[ivec];
    delete [] vecs_1;
}
/*--------------------------------------------
 restart a simulation
 --------------------------------------------*/
inline void DMD_adams::restart()
{
    reset();
    for(int i=0;i<q_max;i++)
        dy[i]=vecs_1[i]->begin();
    y=vecs_1[q_max]->begin();
    e_n=vecs_1[q_max+1]->begin();
}
/*--------------------------------------------
 start the vectors
 --------------------------------------------*/
inline void DMD_adams::start()
{
    type0 sum=forcefield_dmd->ddc_norm_timer()/sqrt(nc_dofs);
    dt=MIN(sqrt(2.0*a_tol/sum),1.0e-3*(t_fin-t_cur));
    type0 dt_max_lcl=numeric_limits<type0>::infinity(),dt_max;
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    for(int i=0;i<ncs;i++)
    {
        if(c[i]>=0.0)
        {
            if(c_d[i]>0.0)
                dt_max_lcl=MIN((1.0-c[i])/c_d[i],dt_max_lcl);
            else if(c_d[i]<0.0)
                dt_max_lcl=MIN(-c[i]/c_d[i],dt_max_lcl);
        }
    }
    MPI_Allreduce(&dt_max_lcl,&dt_max,1,MPI_TYPE0,MPI_MIN,world);
    dt=MIN(dt,dt_max*0.9999999999999);
    
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
    
    memcpy(y,mapp->c->begin(),ncs*sizeof(type0));
    memcpy(dy[0],mapp->c_d->begin(),ncs*sizeof(type0));
    for(int i=0;i<q_max;i++) t[i]=0.0;
    q=1;
}
/*--------------------------------------------
 store the vectors
 --------------------------------------------*/
inline void DMD_adams::update_for_next()
{
    type0* tmp_dy=dy[q_max-1];
    for(int i=q_max-1;i>0;i--)
    {
        t[i]=t[i-1]-dt;
        dy[i]=dy[i-1];
    }
    dy[0]=tmp_dy;
    memcpy(y,mapp->c->begin(),ncs*sizeof(type0));
    memcpy(dy[0],mapp->c_d->begin(),ncs*sizeof(type0));
    for(int i=0;i<ncs;i++)
        if(y[i]>=0.0)
            e_n[i]=y[i]-y_0[i];
    
    dt=dt_new;
    q+=dq;
}
/*--------------------------------------------
 interpolation fail
 --------------------------------------------*/
void DMD_adams::interpolate_fail()
{
    
    type0 dt_max_lcl=inf,dt_max;
    type0 dt_max_lcl_;
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    
    for(int i=0;i<ncs;i++)
    {
        
        if(c[i]>=0.0)
        {
            if(c_d[i]>0.0)
            {
                dt_max_lcl_=(1.0-c[i])/c_d[i];
                if(dt_max_lcl_<dt_max_lcl)
                {
                    while(c[i]+c_d[i]*dt_max_lcl_>1.0)
                        dt_max_lcl_=nextafter(dt_max_lcl_,0.0);
                    
                    dt_max_lcl=dt_max_lcl_;
                }
                
            }
            else if(c_d[i]<0.0)
            {
                dt_max_lcl_=-c[i]/c_d[i];
                if(dt_max_lcl_<dt_max_lcl)
                {
                    while(c[i]+c_d[i]*dt_max_lcl_<0.0)
                        dt_max_lcl_=nextafter(dt_max_lcl_,0.0);
                    
                    dt_max_lcl=dt_max_lcl_;
                }
            }
        }
    }
    
    
    MPI_Allreduce(&dt_max_lcl,&dt_max,1,MPI_TYPE0,MPI_MIN,world);
    
    dt=dt_max+eps;
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
    
    q=1;
    type0 k0,k1=0.0,k2=0.0,y0,dy0;
    int n;
    type0 tmp0,tmp1;
    
    for(int j=0;j<q;j++)
        alpha_dy[j]=0.0;
    
    n=1+q/2;
    for(int i=0;i<n;i++)
    {
        tmp0=wi[n-1][i];
        for(int j=0;j<q;j++)
            tmp0*=xi[n-1][i]-t[j]/dt;
        
        for(int j=0;j<q;j++)
            alpha_dy[j]+=tmp0/(xi[n-1][i]-t[j]/dt);
        
    }
    
    for(int i=0;i<q;i++)
    {
        tmp0=1.0;
        tmp1=1.0;
        for(int j=0;j<q;j++)
            if(i!=j)
            {
                tmp0*=t[i]/dt-t[j]/dt;
                tmp1*=1.0-t[j]/dt;
            }
        alpha_dy[i]*=dt/tmp0;
        dalpha_dy[i]=tmp1/tmp0;
    }
    dalpha_y=0.0;
    alpha_y=1.0;
    
    k0=1.0;
    for(int i=0;i<q-1;i++)
        k0*=1.0-t[i]/dt;
    
    k1=0.0;
    k2=0.0;
    n=1+q/2;
    for(int i=0;i<n;i++)
    {
        tmp0=wi[n-1][i];
        for(int j=0;j<q-1;j++)
            tmp0*=xi[n-1][i]-t[j]/dt;
        k1+=tmp0;
        k2+=tmp0*(xi[n-1][i]-1.0);
    }
    
    beta=dt*k1/k0;
    beta_inv=k0/(dt*k1);
    
    bool xcd_dom=false;
    for(int i=0;i<ncs && !xcd_dom;i++)
    {
        if(y[i]>=0.0)
        {
            dy0=dalpha_y*y[i];
            y0=alpha_y*y[i];
            for(int j=0;j<q;j++)
            {
                dy0+=dalpha_dy[j]*dy[j][i];
                y0+=alpha_dy[j]*dy[j][i];
            }
            
            
            if(y0<=0.0)
            {
                y0=0.0;
                if(dy0>0.0)
                    dy0=0.0;
            }
            if(y0>1.0)
            {
                y0=1.0;
                if(dy0<0.0)
                    dy0=0.0;
                
            }
            
            y_0[i]=y0;
            a[i]=-y0/beta+dy0;
        }
        else
            y_0[i]=y[i];
    }

    
}
/*--------------------------------------------
 calculate the coefficients
 --------------------------------------------*/
inline bool DMD_adams::interpolate()
{
    type0 k0,k1=0.0,k2=0.0;
    int n;
    type0 tmp0,tmp1;
    
    for(int j=0;j<q;j++)
        alpha_dy[j]=0.0;
    
    n=1+q/2;
    for(int i=0;i<n;i++)
    {
        tmp0=wi[n-1][i];
        for(int j=0;j<q;j++)
            tmp0*=xi[n-1][i]-t[j]/dt;
        
        for(int j=0;j<q;j++)
            alpha_dy[j]+=tmp0/(xi[n-1][i]-t[j]/dt);
        
    }
    
    for(int i=0;i<q;i++)
    {
        tmp0=1.0;
        tmp1=1.0;
        for(int j=0;j<q;j++)
            if(i!=j)
            {
                tmp0*=t[i]/dt-t[j]/dt;
                tmp1*=1.0-t[j]/dt;
            }
        alpha_dy[i]*=dt/tmp0;
        dalpha_dy[i]=tmp1/tmp0;
    }
    dalpha_y=0.0;
    alpha_y=1.0;
    
    
    k0=1.0;
    for(int i=0;i<q-1;i++)
        k0*=1.0-t[i]/dt;
    
    k1=0.0;
    k2=0.0;
    n=1+q/2;
    for(int i=0;i<n;i++)
    {
        tmp0=wi[n-1][i];
        for(int j=0;j<q-1;j++)
            tmp0*=xi[n-1][i]-t[j]/dt;
        k1+=tmp0;
        k2+=tmp0*(xi[n-1][i]-1.0);
    }
    
    beta=dt*k1/k0;
    beta_inv=k0/(dt*k1);
    
    bool xcd_dom=false;
    for(int i=0;i<ncs && !xcd_dom;i++)
    {
        if(y[i]>=0.0)
        {
            a[i]=dalpha_y*y[i];
            y_0[i]=alpha_y*y[i];
            for(int j=0;j<q;j++)
            {
                a[i]+=dalpha_dy[j]*dy[j][i];
                y_0[i]+=alpha_dy[j]*dy[j][i];
            }
            
            a[i]-=y_0[i]/beta;
            if(y_0[i]<0.0 || y_0[i]>1.0)
                xcd_dom=true;
        }
        else
            y_0[i]=y[i];
    }
    int err_dom,err_dom_lcl=0;
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
inline void DMD_adams::err_fac_calc()
{
    type0 k1=0.0,k2=0.0,tmp0;
    int n=1+q/2;
    for(int i=0;i<n;i++)
    {
        tmp0=wi[n-1][i];
        for(int j=0;j<q-1;j++)
            tmp0*=xi[n-1][i]-t[j]/dt;
        k1+=tmp0;
        k2+=tmp0*(xi[n-1][i]-1.0);
    }
    err_fac=fabs((k2/k1)/(1.0-t[q-1]/dt));
}
/*--------------------------------------------
 step addjustment after success
 --------------------------------------------*/
inline void DMD_adams::ord_dt(type0& r)
{
    type0 lo_r,hi_r;
    
    r=pow(0.5/err,1.0/static_cast<type0>(q+1));
    ratio_calc(lo_r,hi_r);
    
    dq=0;
    if(hi_r>lo_r && hi_r>r)
    {
        dq=1;
        r=hi_r;
    }
    else if(lo_r>hi_r && lo_r>r)
    {
        dq=-1;
        r=lo_r;
    }
}
/*--------------------------------------------
 find the new del_t and order
 --------------------------------------------*/
inline void DMD_adams::ratio_calc(type0& lo_ratio,type0& hi_ratio)
{
    type0* c=mapp->c->begin();
    type0 lo_err,hi_err,lo_err_lcl,hi_err_lcl;
    type0 tmp0;
    type0 k1;
    int n;
    
    lo_ratio=hi_ratio=0.0;
    
    k1=0.0;
    n=1+q/2;
    for(int i=0;i<n;i++)
    {
        tmp0=wi[n-1][i];
        for(int j=0;j<q-1;j++)
            tmp0*=xi[n-1][i]-t[j]/dt;
        k1+=tmp0;
    }
    
    
    if(q>1)
    {
        type0 a0,a1;
        type0 k0,k4;
        
        k0=1.0;
        for(int i=0;i<q-1;i++)
            k0*=1.0-t[i]/dt;
        
        k4=0.0;
        n=1+(q-1)/2;
        for(int i=0;i<n;i++)
        {
            tmp0=wi[n-1][i]*(xi[n-1][i]-1.0);
            for(int j=0;j<q-2;j++)
                tmp0*=xi[n-1][i]-t[j]/dt;
            k4+=tmp0;
        }
        
        a0=fabs(k4/k1);
        a1=k1*dt/(k0*(1.0-t[q-1]/dt));
        
        lo_err_lcl=0.0;
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
            {
                tmp0=c[i]-y_0[i];
                for(int j=0;j<q;j++)
                    tmp0+=a1*(1.0-t[j]/dt)*dalpha_dy[j]*dy[j][i];
                
                lo_err_lcl+=tmp0*tmp0;
            }
        }
        
        lo_err=0.0;
        MPI_Allreduce(&lo_err_lcl,&lo_err,1,MPI_TYPE0,MPI_SUM,world);
        lo_err=sqrt(lo_err/nc_dofs)/a_tol;
        lo_err*=a0;
        lo_ratio=pow(0.5/lo_err,1.0/static_cast<type0>(q));
    }
    
    if(q<q_max && t[1]!=0.0)
    {
        type0 a0,a1;
        type0 k2,k3,tmp1;
        
        k2=0.0;
        k3=0.0;
        n=1+q/2;
        for(int i=0;i<n;i++)
        {
            tmp0=wi[n-1][i]*(xi[n-1][i]-1.0);
            tmp1=wi[n-1][i];
            for(int j=0;j<q-1;j++)
            {
                tmp0*=xi[n-1][i]-t[j]/dt;
                tmp1*=xi[n-1][i]-1.0+t[j+1]/t[1];
            }
            k2+=tmp0;
            k3+=tmp1;
        }
        
        a0=fabs((k2/k1)/((1.0-t[q-1]/dt)*static_cast<type0>(q+1)));
        a1=-(k1/k3)*(1.0-t[q-1]/dt)/(t[q]/t[1]);
        for(int i=0;i<q+1;i++)
            a1*=-dt/t[1];
        
        hi_err_lcl=0.0;
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
            {
                tmp0=c[i]-y_0[i]+a1*e_n[i];
                hi_err_lcl+=tmp0*tmp0;
            }
        }
        hi_err=0.0;
        MPI_Allreduce(&hi_err_lcl,&hi_err,1,MPI_TYPE0,MPI_SUM,world);
        hi_err=sqrt(hi_err/nc_dofs)/a_tol;
        hi_err*=a0;
        hi_ratio=pow(0.5/hi_err,1.0/static_cast<type0>(q+2));
        
    }
    
}

/*--------------------------------------------
 error calculation
 --------------------------------------------*/
void DMD_adams::err_calc()
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

