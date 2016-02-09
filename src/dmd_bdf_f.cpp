#include <stdlib.h>
#include <limits>
#include "dmd_bdf_f.h"
#include "ff.h"
#include "error.h"
#include "memory.h"
#include "write.h"
#include "thermo_dynamics.h"
#define BDF_VC 0
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMD_bdf_f::DMD_bdf_f(MAPP* mapp,int nargs
,char** args):DMDImplicit(mapp)
{
    q_max=6;
    
    if(nargs>2)
    {
        if(nargs%2!=0)
            error->abort("every keyword in dmd bdf_f should be followed by it's value");
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
                error->abort("unknown keyword in dmd bdf_f: %s",args[iarg]);
        }
    }
    
    if(max_step<0)
        error->abort("max_step in dmd bdf_f should be greater than 0");
    if(max_iter<=0)
        error->abort("max_iter in dmd bdf_f should be greater than 0");
    if(q_max<=0)
        error->abort("max_order in dmd bdf_f should be greater than 0");
    if(a_tol<=0.0)
        error->abort("a_tol in dmd bdf_f should be greater than 0.0");
    if(dt_min<=0.0)
        error->abort("dt_min in dmd bdf_f should be greater than 0.0");

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
DMD_bdf_f::~DMD_bdf_f()
{
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_bdf_f::allocate()
{
    vecs_1=new Vec<type0>*[2];
    vecs_1[0]=new Vec<type0>(atoms,c_dim*(q_max+1));
    vecs_1[1]=new Vec<type0>(atoms,c_dim);
    
    CREATE1D(t,q_max+1);
    CREATE1D(l,q_max+1);
    
    A=new type0*[q_max+1];
    *A=new type0[(q_max+2)*(q_max+1)/2];
    for(int i=1;i<q_max+1;i++)
        A[i]=A[i-1]+q_max+1-i;
    
    A_bar=new type0*[q_max+1];
    *A_bar=new type0[(q_max+2)*(q_max+1)/2];
    for(int i=1;i<q_max+1;i++)
        A_bar[i]=A_bar[i-1]+q_max-i;
    
    for(int i=0;i<q_max+1;i++)
        A[0][i]=A[i][i]=1.0;
    
    for(int i=2;i<q_max+1;i++)
        for(int j=1;j<i;j++)
            A[j][i]=A[j][i-1]+A[j-1][i-1];

}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_bdf_f::deallocate()
{
    delete [] *A;
    delete [] A;
    delete [] *A_bar;
    delete [] A_bar;
    
    delete [] t;
    delete [] l;
    delete vecs_1[0];
    delete vecs_1[1];
    delete [] vecs_1;
}
/*--------------------------------------------
 restart a simulation
 --------------------------------------------*/
inline void DMD_bdf_f::restart()
{
    reset();
    z=vecs_1[0]->begin();
    dy=vecs_1[1]->begin();
}
/*--------------------------------------------
 start the vectors
 --------------------------------------------*/
inline void DMD_bdf_f::start()
{
    type0 sum=forcefield_dmd->ddc_norm_timer()/sqrt(nc_dofs);
    dt=MIN(sqrt(2.0*a_tol/sum),1.0e-3*(t_fin-t_cur));
    
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
    
    type0* z_=z;
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    for(int i=0;i<ncs;i++)
    {
        if(c[i]>=0.0)
        {
            z_[0]=c[i];
            z_[1]=c_d[i];
        }
        else
        {
            z_[0]=0.0;
            z_[1]=0.0;
        }
        
        z_+=q_max+1;
    }
    
    for(int i=0;i<q_max+1;i++)
        t[i]=0.0;
    q=1;
    lo_avail=hi_avail=false;
}
/*--------------------------------------------
 store the vectors
 --------------------------------------------*/
inline void DMD_bdf_f::update_for_next()
{
    type0* c=mapp->c->begin();
    for(int i=0;i<ncs;i++)
    {
        if(c[i]>=0.0)
        {
            dy[i]=c[i]-y_0[i];
        }
        else
            dy[i]=0.0;
    }
    
    if(BDF_VC)
        update_z_vc();
    else
        update_z_flc();
    
    for(int i=q_max;i>0;i--)
        t[i]=t[i-1]-dt;
    
    dt=dt_new;
    q+=dq;
}
/*--------------------------------------------
 update A_bar
 --------------------------------------------*/
inline void DMD_bdf_f::update_A_bar()
{
    type0 dt_;
    for(int i=0;i<q+1;i++)
    {
        dt_=1.0;
        for(int j=i;j<q+1;j++)
        {
            A_bar[i][j]=A[i][j]*dt_;
            dt_*=dt;
        }
    }
}
/*--------------------------------------------

 --------------------------------------------*/
void DMD_bdf_f::beta_calc_flc()
{
    type0 iq=0.0;
    beta_inv=0.0;
    for(int i=0;i<q;i++,iq++)
        beta_inv+=1.0/(iq+1.0);
    beta_inv/=dt;
    beta=1.0/beta_inv;
}
/*--------------------------------------------

 --------------------------------------------*/
void DMD_bdf_f::beta_calc_vc()
{
    beta_inv=0.0;
    for(int i=0;i<q;i++)
        beta_inv+=1.0/(dt-t[i]);
    beta=1.0/beta_inv;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
inline bool DMD_bdf_f::interpolate()
{
    update_A_bar();
    if(BDF_VC)
        beta_calc_vc();
    else
        beta_calc_flc();
    
    type0* z_=z;
    type0 y0;
    bool xcd_dom=false;
    type0* c=mapp->c->begin();
    
    for(int i=0;i<ncs && !xcd_dom;i++)
    {
        if(c[i]>=0.0)
        {
            y0=0.0;
            for(int j=0;j<q+1;j++)
                y0+=A_bar[0][j]*z_[j];
            
            if(y0<0.0 || y0>1.0)
            {
                xcd_dom=true;
                continue;
            }
            y_0[i]=y0;
            a[i]=-y0*beta_inv;
            for(int j=1;j<q+1;j++)
                a[i]+=A_bar[1][j]*z_[j];
        }
        else
            y_0[i]=c[i];
        
        
        z_+=q_max+1;
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
inline void DMD_bdf_f::err_fac_calc()
{
    if(BDF_VC)
        calc_err_coef_vc();
    else
        calc_err_coef_flc();
}
/*--------------------------------------------
 step addjustment after success
 --------------------------------------------*/
inline void DMD_bdf_f::ord_dt(type0& r)
{
    eta_calc();
    if(eta[0]>eta[1] && eta[0]>eta[2])
    {
        dq=-1;
        r=eta[0];
    }
    else if(eta[2]>eta[0] && eta[2]>eta[1])
    {
        dq=1;
        r=eta[2];
    }
    else
    {
        dq=0;
        r=eta[1];
    }
}
/*--------------------------------------------
 error calculation
 --------------------------------------------*/
void DMD_bdf_f::err_calc()
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
/*--------------------------------------------
 update l array
 --------------------------------------------*/
void DMD_bdf_f::update_l(int q_)
{
    l[0]=1.0;
    for(int i=1;i<q_+1;i++)
        l[i]=0.0;

    type0 k0;

    for(int i=0;i<q_;i++)
    {
        k0=1.0/(dt-t[i]);
        for(int j=i+1;j>0;j--)
            l[j]+=k0*l[j-1];
        
    }
}
/*--------------------------------------------
 update Nordsieck
 --------------------------------------------*/
void DMD_bdf_f::update_z_vc()
{
    // VC-f
    
    
    if(dq==1)
    {
        update_l(q);
        
        type0* z_=z;
        type0 dy_;
        for(int i=0;i<ncs;i++)
        {
            dy_=dy[i];
            for(int j=0;j<q+1;j++)
            {
                for(int k=j+1;k<q+1;k++)
                    z_[j]+=A_bar[j][k]*z_[k];
                z_[j]+=l[j]*dy_;
            }
            z_[q+1]=0.0;
            
            z_+=q_max+1;
        }
        
    }
    else if(dq==-1)
    {
        update_l(q-2);
        type0 p_q_2=1.0;
        for(int i=0;i<q-2;i++)
            p_q_2*=dt-t[i];
        for(int i=2;i<q;i++)
            A_bar[i][q]-=p_q_2*l[i-2];
        
        type0 k0=1.0/(dt-t[q-1])+1.0/(dt-t[q-2]);
        
        l[q-1]=0.0;
        for(int i=q-1;i>0;i--)
            l[i]+=k0*l[i-1];
        
        type0* z_=z;
        type0 dy_;
        for(int i=0;i<ncs;i++)
        {
            dy_=dy[i];
            for(int j=0;j<q;j++)
            {
                for(int k=j+1;k<q+1;k++)
                    z_[j]+=A_bar[j][k]*z_[k];
                z_[j]+=l[j]*dy_;
            }
            
            z_+=q_max+1;
        }
        
    }
    else
    {
        update_l(q);
        
        type0* z_=z;
        type0 dy_;
        for(int i=0;i<ncs;i++)
        {
            dy_=dy[i];
            for(int j=0;j<q+1;j++)
            {
                for(int k=j+1;k<q+1;k++)
                    z_[j]+=A_bar[j][k]*z_[k];
                z_[j]+=l[j]*dy_;
            }
            
            z_+=q_max+1;
        }
        
    }

    
}
/*--------------------------------------------
 update Nordsieck
 --------------------------------------------*/
void DMD_bdf_f::update_z_flc()
{
    // FLC-f
    if(dq==1)
    {

        type0 k2=0.0,iq=0.0;
        for(int i=0;i<q;i++,iq++)
            k2+=1.0/(dt*(iq+1.0))-1.0/(dt-t[i]);
            
        
        update_l(q);
        
        
        
        type0* z_=z;
        type0 dy_;
        for(int i=0;i<ncs;i++)
        {
            dy_=dy[i];
            for(int j=0;j<q+1;j++)
            {
                for(int k=j+1;k<q+1;k++)
                    z_[j]+=A_bar[j][k]*z_[k];
                z_[j]+=l[j]*dy_;
            }
            z_[q+1]=0.0;
            
            z_+=q_max+1;
        }
        
    }
    else if(dq==-1)
    {
        
        update_l(q-2);
        type0 p_q_2=1.0;
        for(int i=0;i<q-2;i++)
            p_q_2*=dt-t[i];
        for(int i=2;i<q;i++)
            A_bar[i][q]-=p_q_2*l[i-2];

        
        type0 k1=0.0,iq=0.0;
        for(int i=0;i<q-2;i++,iq++)
            k1+=1.0/(dt*(iq+1.0))-1.0/(dt-t[i]);
        k1+=1.0/(dt*(iq+1.0))+1.0/(dt*(iq+2.0));
        
        
        l[q-1]=0.0;
        for(int i=q-1;i>0;i--)
            l[i]+=k1*l[i-1];
        
        type0* z_=z;
        type0 dy_;
        for(int i=0;i<ncs;i++)
        {
            dy_=dy[i];
            for(int j=0;j<q;j++)
            {
                for(int k=j+1;k<q+1;k++)
                    z_[j]+=A_bar[j][k]*z_[k];
                z_[j]+=l[j]*dy_;
            }
            
            z_+=q_max+1;
        }
        
    }
    else
    {
        update_l(q-1);

        type0 k0=0.0,iq=0.0;
        for(int i=0;i<q-1;i++,iq++)
            k0+=1.0/(dt*(iq+1.0))-1.0/(dt-t[i]);
        k0+=1.0/(dt*(iq+1.0));
        
        l[q]=0.0;
        for(int i=q;i>0;i--)
            l[i]+=k0*l[i-1];
        
        
        type0* z_=z;
        type0 dy_;
        for(int i=0;i<ncs;i++)
        {
            dy_=dy[i];
            for(int j=0;j<q+1;j++)
            {
                for(int k=j+1;k<q+1;k++)
                    z_[j]+=A_bar[j][k]*z_[k];
                z_[j]+=l[j]*dy_;
            }
            
            z_+=q_max+1;
        }
        
    }
    
    
}
/*--------------------------------------------
 calculate the error coefficent
 --------------------------------------------*/
void DMD_bdf_f::calc_err_coef_flc()
{
    type0 u_q=0.0,s_bar_q=0.0;
    type0 iq=0.0;
    
    for(int i=0;i<q;i++,iq++)
    {
        s_bar_q+=1.0/(iq+1.0);
        u_q+=1.0/(dt-t[i]);
    }
    
    
    u_q*=dt;
    
    u_q++;
    u_q-=s_bar_q;
    
    err_fac=fabs(u_q/(s_bar_q*(1.0+iq*u_q)));
    
    
    
    if(const_q>1 && q<q_max)
    {
        type0 r=1.0;
        for(int i=1;i<q;i++)
            r*=(t[i]-dt)/t[i];
        type0 u_p_q=0.0;
        for(int i=0;i<q;i++)
            u_p_q+=1.0/t[i+1];
        
        u_p_q*=t[1];
        u_p_q++;
        u_p_q-=s_bar_q;
        
        type0 tmp=(1.0-t[q]/dt)*(u_q-1.0/(iq+1.0)+dt/(dt-t[q]));
        tmp/=(1.0+iq*u_q)*(s_bar_q+1.0/(iq+1.0))*(iq+2.0);
        hi_err_fac[0]=fabs(tmp);
        hi_err_fac[1]=-dt*dt*r*(1.0+iq*u_q)/(t[1]*t[q]*(1.0+iq*u_p_q));
    }
    
    if(q>1)
    {
        type0 p_q=1.0;
        for(int i=0;i<q;i++,iq++)
            p_q*=dt-t[i];
        
        type0 tmp=u_q-dt/(dt-t[q-1])+1.0/iq;
        tmp/=s_bar_q-1.0/iq;
        tmp/=1.0-t[q-1]/dt;
        
        lo_err_fac[0]=fabs(tmp);
        lo_err_fac[1]=p_q;
    }

}
/*--------------------------------------------
 calculate the error coefficent
 --------------------------------------------*/
void DMD_bdf_f::calc_err_coef_vc()
{
    
    type0 s_q=0.0,r=1.0,iq=0.0;
    for(int i=0;i<q;i++)
        s_q+=1.0/(dt-t[i]);
    for(int i=1;i<q;i++,iq++)
        r*=(t[i]-dt)/t[i];
    
    
    err_fac=fabs(1.0/(s_q*(1.0+r)*dt));

    if(const_q>1 && q<q_max)
    {
        type0 r_p=1.0;
        for(int i=1;i<q;i++)
            r_p*=t[i+1]/(t[i+1]-t[1]);
        
        type0 tmp=(dt-t[q])/(dt*dt*(iq+2.0)*r);
        tmp/=s_q+1.0/(dt-t[q]);
        
        hi_err_fac[0]=fabs(tmp);
        hi_err_fac[1]=-dt*dt*r*(1.0+r)/(t[1]*t[q]*(1.0+r_p));
    }
    
    if(q>1)
    {
        type0 p_q=1.0;
        for(int i=0;i<q;i++,iq++)
            p_q*=dt-t[i];
        
        lo_err_fac[0]=fabs(1.0/((dt-t[q])*s_q-1.0));
        lo_err_fac[1]=p_q;
    }
    
}
/*--------------------------------------------
 eta calculation
 --------------------------------------------*/
void DMD_bdf_f::eta_calc()
{
    type0* c=mapp->c->begin();
    type0 norm_lcl=0.0;
    type0 tmp,norm,iq=static_cast<type0>(q);
    eta[1]=pow(0.5/err,1.0/(iq+1.0));
    
    eta[0]=eta[2]=0.0;
    
    if(const_q>1 && q<q_max)
    {
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
            {
                tmp=c[i]-y_0[i]+hi_err_fac[1]*dy[i];
                norm_lcl+=tmp*tmp;
            }
        }
        
        MPI_Allreduce(&norm_lcl,&tmp,1,MPI_TYPE0,MPI_SUM,world);
        norm=hi_err_fac[0]*sqrt(tmp/nc_dofs)/a_tol;
        eta[2]=pow(0.5/norm,1.0/(iq+2.0));
    }
    
    if(q>1)
    {
        type0* z_=z;
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
            {
                tmp=c[i]-y_0[i]+lo_err_fac[1]*z_[q];
                norm_lcl+=tmp*tmp;
            }
            
            z_+=q_max+1;
        }
        MPI_Allreduce(&norm_lcl,&tmp,1,MPI_TYPE0,MPI_SUM,world);
        norm=lo_err_fac[0]*sqrt(tmp/nc_dofs)/a_tol;
        eta[0]=pow(0.5/norm,1.0/(iq));
    }
}



