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
    max_order=12;
    
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
                    error->abort("initial_del_t in dmd adams should be greater than 0.0");
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
    if(max_order<=0)
        error->abort("max_order in dmd adams should be greater than 0");
    if(m_tol<=0.0)
        error->abort("m_tol in dmd adams should be greater than 0.0");
    if(a_tol<=0.0)
        error->abort("a_tol in dmd adams should be greater than 0.0");
    if(min_del_t<=0.0)
        error->abort("min_del_t in dmd adams should be greater than 0.0");
    
    int n=1+max_order/2;
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

    int n=1+max_order/2;
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
 
 --------------------------------------------*/
void DMD_adams::reset()
{
    DMDImplicit::reset();
    for(int i=0;i<max_order;i++)
        dy[i]=vecs_1[i]->begin();
    y=vecs_1[max_order]->begin();
    e_n=vecs_1[max_order+1]->begin();
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD_adams::allocate()
{
    DMDImplicit::allocate();
    
    vecs_1=new Vec<type0>*[max_order+2];
    for(int ivec=0;ivec<max_order+2;ivec++)
        vecs_1[ivec]=new Vec<type0>(atoms,c_dim);
    
    CREATE1D(t,max_order);
    CREATE1D(dy,max_order);
    CREATE1D(alpha_dy,max_order);
    CREATE1D(dalpha_dy,max_order);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD_adams::deallocate()
{
    if(max_order)
    {
        delete [] dy;
        delete [] t;
        delete [] alpha_dy;
        delete [] dalpha_dy;
    }
    
    for(int ivec=0;ivec<max_order+2;ivec++)
        delete vecs_1[ivec];
    delete [] vecs_1;
    
    DMDImplicit::deallocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_adams::init()
{
    DMDImplicit::init();
    allocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_adams::fin()
{
    deallocate();    
    DMDImplicit::fin();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_adams::run()
{
    if(max_step==0)
        return;

    type0 del_t=initial_del_t,del_t_tmp;
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
        type0* tmp_dy;
        int iter=-1,iter_0=-1;
        iter_dcr_cntr=-1;        
        del_t=est_dt();
        
        memcpy(y,c,ncs*sizeof(type0));
        memcpy(dy[0],c_d,ncs*sizeof(type0));
        for(int i=0;i<max_order;i++) t[i]=0.0;
        q=1;
        while (istep<max_step && tot_t<max_t && !min_run)
        {
            err=1.0;
            cost=1.0;
            while (MAX(cost,err)>=1.0)
            {
                interpolate(del_t,q);
                if(atoms->my_p==0) printf("%d %e: ",q,del_t);
                iter=solve_n_err(cost,err);
                if(MAX(cost,err)<1.0)
                    continue;
                if(cost>=1.0)
                    iter_0=-1;
             
                atoms->update(mapp->c);
            }
            
            if(iter_0!=-1 && (iter<iter_0 || iter==0))
                iter_dcr_cntr++;
            else
                iter_dcr_cntr=0;
            
            iter_0=iter;
            
            max_succ_q=MAX(max_succ_q,q);

            min_run=decide_min(istep,del_t);
            if(min_run) continue;
            
            del_t_tmp=del_t;
            ord_dt(err,del_t,q);
            
            tmp_dy=dy[max_order-1];
            for(int i=max_order-1;i>0;i--)
            {
                t[i]=t[i-1]-del_t_tmp;
                dy[i]=dy[i-1];
            }
            dy[0]=tmp_dy;
            memcpy(y,c,ncs*sizeof(type0));
            memcpy(dy[0],c_d,ncs*sizeof(type0));
        }

        if(!min_run) continue;
        do_min();
        min_run=false;
    }
}
/*--------------------------------------------

 --------------------------------------------*/
inline type0 DMD_adams::est_dt()
{
    type0 sum=forcefield_dmd->ddc_norm_timer()/sqrt(nc_dofs);
    type0 del_t=MIN(sqrt(2.0*a_tol/sum),1.0e-3*(max_t-tot_t));
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
    return del_t;
}
/*--------------------------------------------
 calculate the coefficients
 --------------------------------------------*/
inline void DMD_adams::interpolate(type0& del_t,int& q)
{
    int n;
    type0 tmp0,tmp1;
    type0 k0,k1,k2;
    
    int idof;
    int err_chk=1;
    int err_chk_lcl;
    
    while(err_chk)
    {
        for(int j=0;j<q;j++)
            alpha_dy[j]=0.0;
        
        n=1+q/2;
        for(int i=0;i<n;i++)
        {
            tmp0=wi[n-1][i];
            for(int j=0;j<q;j++)
                tmp0*=xi[n-1][i]-t[j]/del_t;
            
            for(int j=0;j<q;j++)
                alpha_dy[j]+=tmp0/(xi[n-1][i]-t[j]/del_t);
            
        }
        
        for(int i=0;i<q;i++)
        {
            tmp0=1.0;
            tmp1=1.0;
            for(int j=0;j<q;j++)
                if(i!=j)
                {
                    tmp0*=t[i]/del_t-t[j]/del_t;
                    tmp1*=1.0-t[j]/del_t;
                }
            alpha_dy[i]*=del_t/tmp0;
            dalpha_dy[i]=tmp1/tmp0;
        }
        dalpha_y=0.0;
        alpha_y=1.0;
        
        
        k0=1.0;
        for(int i=0;i<q-1;i++)
            k0*=1.0-t[i]/del_t;
        
        k1=0.0;
        k2=0.0;
        n=1+q/2;
        for(int i=0;i<n;i++)
        {
            tmp0=wi[n-1][i];
            for(int j=0;j<q-1;j++)
                tmp0*=xi[n-1][i]-t[j]/del_t;
            k1+=tmp0;
            k2+=tmp0*(xi[n-1][i]-1.0);
        }
        
        beta=del_t*k1/k0;
        
        err_chk_lcl=0;
        idof=0;
        while(idof<ncs && err_chk_lcl==0)
        {
            if(y[idof]>=0.0)
            {
                a[idof]=beta*dalpha_y*y[idof];
                y_0[idof]=alpha_y*y[idof];
                for(int j=0;j<q;j++)
                {
                    a[idof]+=beta*dalpha_dy[j]*dy[j][idof];
                    y_0[idof]+=alpha_dy[j]*dy[j][idof];
                }
                
                a[idof]-=y_0[idof];
                if(y_0[idof]<0.0 || y_0[idof]>1.0)
                {
                    err_chk_lcl=1;
                }
            }
            else
                y_0[idof]=y[idof];
            
            
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
                            if(y[i]>=0.0)
                            {
                                tmp0=y[i]+r_lcl*del_t*dy[0][i];
                                if(tmp0>1.0)
                                    r_lcl=0.999*((1.0-y[i])/(del_t*dy[0][i]));
                                if(tmp0<0.0)
                                    r_lcl=-0.999*((y[i])/(del_t*dy[0][i]));
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
        err_prefac=fabs((k2/k1)/(1.0-t[q-1]/del_t));
        if(err_chk)
            intp_rej++;
    }
}
/*--------------------------------------------
 find the new del_t and order
 --------------------------------------------*/
inline void DMD_adams::ord_dt(type0 err,
type0& del_t,int& q)
{
    type0* c=mapp->c->begin();
    type0 lo_r,hi_r,r;
    int tmp_q=q;

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
        
    if(iter_dcr_cntr<iter_dcr_thrsh)
        r=MIN(r,1.0);
    
    if(r>=2.0)
        r=2.0;
    else if(r<=0.5)
        r=0.5;
    else
        r=1.0;
    
    if(r>1.0)
        iter_dcr_cntr=0;
    
    bool const_stp_chk=false;
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
                if(r==1.0)
                    const_stp_chk=true;
            }
        }
    }
    
    if(const_stp_chk && tmp_q==q)
        const_stps++;
    else
        const_stps=0;
    
    for(int i=0;i<ncs;i++)
        if(c[i]>=0.0)
            e_n[i]=c[i]-y_0[i];
}
/*--------------------------------------------
 find the new del_t and order
 --------------------------------------------*/
inline void DMD_adams::ratio_calc(int q,type0 del_t,
type0& lo_ratio,type0& hi_ratio)
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
            tmp0*=xi[n-1][i]-t[j]/del_t;
        k1+=tmp0;
    }
    
    
    if(q>1)
    {
        type0 a0,a1;
        type0 k0,k4;
        
        k0=1.0;
        for(int i=0;i<q-1;i++)
            k0*=1.0-t[i]/del_t;
        
        k4=0.0;
        n=1+(q-1)/2;
        for(int i=0;i<n;i++)
        {
            tmp0=wi[n-1][i]*(xi[n-1][i]-1.0);
            for(int j=0;j<q-2;j++)
                tmp0*=xi[n-1][i]-t[j]/del_t;
            k4+=tmp0;
        }
        
        a0=fabs(k4/k1);
        a1=k1*del_t/(k0*(1.0-t[q-1]/del_t));
        
        lo_err_lcl=0.0;
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
            {
                tmp0=c[i]-y_0[i];
                for(int j=0;j<q;j++)
                    tmp0+=a1*(1.0-t[j]/del_t)*dalpha_dy[j]*dy[j][i];
                
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
                tmp0*=xi[n-1][i]-t[j]/del_t;
                tmp1*=xi[n-1][i]-1.0+t[j+1]/t[1];
            }
            k2+=tmp0;
            k3+=tmp1;
        }
        
        a0=fabs((k2/k1)/((1.0-t[q-1]/del_t)*static_cast<type0>(q+1)));
        a1=-(k1/k3)*(1.0-t[q-1]/del_t)/(t[q]/t[1]);
        for(int i=0;i<q+1;i++)
            a1*=-del_t/t[1];
        
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
inline void DMD_adams::fail_stp_adj(type0 err,type0 m_err,type0& del_t,int& q)
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
            type0 r=pow(0.5/MAX(err,m_err),1.0/static_cast<type0>(q+1));
            
            r=MAX(r,0.1);
            r=MIN(r,0.9);
            
            if(r*del_t<min_del_t)
                del_t=min_del_t;
            else if(r*del_t>max_t-tot_t-min_del_t)
                del_t=max_t-tot_t-min_del_t;
            else
                del_t*=r;
        }
    }
}

