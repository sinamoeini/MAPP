#include <stdlib.h>
#include <limits>
#include "dmd_bdf.h"
#include "ff.h"
#include "error.h"
#include "memory.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMD_bdf::DMD_bdf(MAPP* mapp,int nargs
,char** args):DMDImplicit(mapp)
{
    max_order=6;
    
    if(nargs>2)
    {
        if(nargs%2!=0)
            error->abort("every keyword in dmd bdf should be followed by it's value");
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
            else if(strcmp(args[iarg],"max_order")==0)
            {
                iarg++;
                max_order=atoi(args[iarg]);
                iarg++;
            }
            else if(strcmp(args[iarg],"m_tol")==0)
            {
                iarg++;
                m_tol=atof(args[iarg]);
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
                    error->abort("initial_del_t in dmd bdf should be greater than 0.0");
                iarg++;
            }
            else
                error->abort("unknown keyword in dmd bdf: %s",args[iarg]);
        }
    }

    if(max_step<0)
        error->abort("max_step in dmd bdf should be greater than 0");
    if(max_iter<=0)
        error->abort("max_iter in dmd bdf should be greater than 0");
    if(max_order<=0)
        error->abort("max_order in dmd bdf should be greater than 0");
    if(m_tol<=0.0)
        error->abort("m_tol in dmd bdf should be greater than 0.0");
    if(a_tol<=0.0)
        error->abort("a_tol in dmd bdf should be greater than 0.0");
    if(min_del_t<=0.0)
        error->abort("min_del_t in dmd bdf should be greater than 0.0");
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
DMD_bdf::~DMD_bdf()
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD_bdf::allocate()
{
    DMDImplicit::allocate();
    vecs_1=new Vec<type0>*[max_order+2];
    for(int ivec=0;ivec<max_order+2;ivec++)
        vecs_1[ivec]=new Vec<type0>(atoms,c_dim);
    
    
    CREATE1D(t,max_order);
    CREATE1D(y,max_order);
    CREATE1D(alpha_y,max_order);
    CREATE1D(dalpha_y,max_order);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD_bdf::reset()
{
    DMDImplicit::reset();    
    for(int i=0;i<max_order;i++)
        y[i]=vecs_1[i]->begin();
    dy=vecs_1[max_order]->begin();
    e_n=vecs_1[max_order+1]->begin();
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD_bdf::deallocate()
{
    if(max_order)
    {
        delete [] y;
        delete [] t;
        delete [] alpha_y;
        delete [] dalpha_y;
    }
    
    for(int ivec=0;ivec<max_order+2;ivec++)
        delete vecs_1[ivec];
    delete [] vecs_1;
    
    DMDImplicit::deallocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
type0 DMD_bdf::est_dt()
{
    type0 sum=forcefield_dmd->ddc_norm_timer()/sqrt(nc_dofs);
    type0 del_t=MIN(sqrt(2.0*a_tol/sum),1.0e-3*(max_t-tot_t));
    init_stp_adj(del_t);
    return del_t;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_bdf::init()
{
    DMDImplicit::init();
    allocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_bdf::fin()
{
    DMDImplicit::fin();
    deallocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_bdf::run()
{
    if(max_step==0)
        return;
    
    type0 del_t,del_t_tmp;
    type0 cost,err=0.0;
    int q;
    int istep;
    bool min_run=false;
    
    istep=0;
    while (istep<max_step && tot_t<max_t)
    {
        reset();
        type0* c=mapp->c->begin();
        type0* c_d=mapp->c_d->begin();
        type0* tmp_y;
        del_t=est_dt();
        
        memcpy(y[0],c,ncs*sizeof(type0));
        memcpy(dy,c_d,ncs*sizeof(type0));
        for(int i=0;i<max_order;i++) t[i]=0.0;
        q=1;
        while (istep<max_step && tot_t<max_t && !min_run)
        {
            err=1.0;
            while(err>=1.0)
            {
                interpolate(del_t,q);
                solve_n_err(cost,err);
                if(err<1.0)
                    continue;
                
                fail_stp_adj(err,del_t,q);
                memcpy(c,y[0],ncs*sizeof(type0));
                memcpy(c_d,dy,ncs*sizeof(type0));
            }
            
            max_succ_q=MAX(max_succ_q,q);


            min_run=decide_min(istep,del_t);
            if(min_run) continue;
            
            del_t_tmp=del_t;
            ord_dt(del_t,q,err);
            
            tmp_y=y[max_order-1];
            for(int i=max_order-1;i>0;i--)
            {
                t[i]=t[i-1]-del_t_tmp;
                y[i]=y[i-1];
            }
            y[0]=tmp_y;
            memcpy(y[0],c,ncs*sizeof(type0));
            memcpy(dy,c_d,ncs*sizeof(type0));
        }
        
        if(!min_run) continue;
        do_min();
        min_run=false;
    }
}
/*--------------------------------------------
 calculate the coefficients
 --------------------------------------------*/
void DMD_bdf::interpolate(type0& del_t,int& q)
{
    type0 tmp0,tmp1,k0,k1,k2;
    
    int idof;
    int err_chk=1;
    int err_chk_lcl;
    
    while(err_chk)
    {
        k0=1.0;
        k1=k2=0.0;
        for(int i=1;i<q;i++)
        {
            tmp0=del_t/t[i];
            tmp1=t[i]/del_t;
            k0*=1.0-tmp0;
            k1+=tmp0;
            k2+=1.0/(1.0-tmp1);
        }
        alpha_y[0]=k0*(1.0+k1);
        for(int i=1;i<q;i++)
        {
            tmp1=del_t/t[i];
            tmp0=tmp1*tmp1;
            for(int j=1;j<q;j++)
                if(i!=j)
                    tmp0*=(del_t-t[j])/(t[i]-t[j]);
            alpha_y[i]=tmp0;
        }
        dalpha_y[0]=k0*(k1+k2+k1*k2)/del_t;
        for(int i=1;i<q;i++)
            dalpha_y[i]=alpha_y[i]*(2.0+k2-1.0/(1.0-t[i]/del_t))/del_t;
        
        
        alpha_dy=k0*del_t;
        dalpha_dy=k0*(1.0+k2);
        
        beta=del_t/(1.0+k2);
        
        err_chk_lcl=0;
        idof=0;
        while(idof<ncs && err_chk_lcl==0)
        {
            if(y[0][idof]>=0.0)
            {
                a[idof]=beta*dalpha_dy*dy[idof];
                y_0[idof]=alpha_dy*dy[idof];
                for(int j=0;j<q;j++)
                {
                    a[idof]+=beta*dalpha_y[j]*y[j][idof];
                    y_0[idof]+=alpha_y[j]*y[j][idof];
                }
                
                a[idof]-=y_0[idof];
                if(y_0[idof]<0.0 || y_0[idof]>1.0)
                {
                    err_chk_lcl=1;
                }
            }
            else
                y_0[idof]=y[0][idof];
            
            
            idof++;
        }
        
        MPI_Allreduce(&err_chk_lcl,&err_chk,1,MPI_INT,MPI_MAX,world);
        if(err_chk)
        {
            if(q>1)
            {
                q--;
            }
            else
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
                        type0 r_lcl=1.0,r;
                        for(int i=0;i<ncs;i++)
                        {
                            if(y[0][i]>=0.0)
                            {
                                tmp0=y[0][i]+r_lcl*del_t*dy[i];
                                if(tmp0>1.0)
                                    r_lcl=0.999*((1.0-y[0][i])/(del_t*dy[i]));
                                if(tmp0<0.0)
                                    r_lcl=-0.999*((y[0][i])/(del_t*dy[i]));
                            }
                        }
                        
                        MPI_Allreduce(&r_lcl,&r,1,MPI_TYPE0,MPI_MIN,world);
                        
                        if(r*del_t<min_del_t)
                            del_t=min_del_t;
                        else if(r*del_t>max_t-tot_t-min_del_t)
                            del_t=max_t-tot_t-min_del_t;
                        else
                            del_t*=r;
                    }
                }
            }
        }
        err_prefac=fabs(-1.0/((1.0+k0)*(1.0+k2)));
        
        if(err_chk)
            intp_rej++;
    }

}
/*--------------------------------------------
 find the new del_t and order
 --------------------------------------------*/
void DMD_bdf::ord_dt(type0& del_t,int& q
,type0 err)
{
    type0* c=mapp->c->begin();
    type0 lo_r,hi_r,r;
    
    r=pow(0.5/err,1.0/static_cast<type0>(q+1));
    ratio_calc(q,del_t,lo_r,hi_r);
    
    if(hi_r>lo_r && hi_r>r)
    {
        q++;
        r=hi_r;
    }
    else if(lo_r>hi_r && lo_r>r)
    {
        q--;
        r=lo_r;
    }
    
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
    
    for(int i=0;i<ncs;i++)
        if(c[i]>=0.0)
            e_n[i]=c[i]-y_0[i];

}
/*--------------------------------------------
 error calculation
 --------------------------------------------*/
void DMD_bdf::ratio_calc(int q,type0 del_t,
type0& lo_ratio,type0& hi_ratio)
{
    type0* c=mapp->c->begin();
    type0 lo_err,hi_err,lo_err_lcl,hi_err_lcl;
    type0 tmp0,tmp1;
    type0 k0,k2;
    
    lo_ratio=hi_ratio=0.0;
    
    k0=1.0;
    k2=0.0;
    for(int i=1;i<q;i++)
    {
        tmp0=del_t/t[i];
        tmp1=t[i]/del_t;
        k0*=1.0-tmp0;
        k2+=1.0/(1.0-tmp1);
    }
    
    
    if(q>1)
    {
        type0 a0;
        a0=fabs(-1.0/(static_cast<type0>(q-1)*(1.0+k2)*(1.0-t[q-1]/del_t)));
        
        lo_err_lcl=0.0;
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
            {
                tmp0=c[i];
                tmp0-=k0*y[0][i];
                for(int j=1;j<q;j++)
                    tmp0-=(t[j]/del_t)*alpha_y[j]*y[j][i];
                
                lo_err_lcl+=tmp0*tmp0;
            }
        }
        
        lo_err=0.0;
        MPI_Allreduce(&lo_err_lcl,&lo_err,1,MPI_TYPE0,MPI_SUM,world);
        lo_err=sqrt(lo_err/nc_dofs)/a_tol;
        lo_err*=a0;
        lo_ratio=pow(0.5/lo_err,1.0/static_cast<type0>(q));
    }
    
    if(q<max_order && t[1]!=0.0)
    {
        type0 a0,a1;
        type0 k3;
        
        k3=1.0;
        for(int i=0;i<q-1;i++)
            k3*=1.0/(1.0-t[1]/t[i+2]);
        
        a0=fabs(-(1.0-t[q]/del_t)/(static_cast<type0>((q+1)*(q+2))*(1.0+k0)*(1.0+k2)));
        a1=k0*(1.0+k0)*del_t*del_t/((1.0+k3)*t[1]*t[q]);
        
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
 step addjustment after failure
 --------------------------------------------*/
inline void DMD_bdf::fail_stp_adj(type0 err,type0& del_t,int& q)
{
    if(max_t-tot_t<=2.0*min_del_t)
    {
        if(q>1)
            q--;
        else
            error->abort("reached minimum order & del_t (%e)",del_t);
    }
    else
    {
        if(del_t==min_del_t)
        {
            if(q>1)
                q--;
            else
                error->abort("reached minimum order & del_t (%e)",del_t);
        }
        else
        {
            type0 r=pow(0.5/err,1.0/static_cast<type0>(q+1));
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
/*--------------------------------------------
 init
 --------------------------------------------*/
inline void DMD_bdf::init_stp_adj(type0& del_t)
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

