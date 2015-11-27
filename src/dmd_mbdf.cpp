#include <stdlib.h>
#include <limits>
#include "dmd_mbdf.h"
#include "ff.h"
#include "error.h"
#include "memory.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMD_mbdf::DMD_mbdf(MAPP* mapp,int nargs
,char** args):DMDImplicit(mapp)
{
    pre_cond=1;
    max_order=5;
    
    if(nargs>2)
    {
        if(nargs%2!=0)
            error->abort("every keyword in dmd mbdf should be followed by it's value");
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
                    error->abort("initial_del_t in dmd mbdf should be greater than 0.0");
                iarg++;
            }
            else
                error->abort("unknown keyword in dmd mbdf: %s",args[iarg]);
        }
    }
    
    if(max_step<0)
        error->abort("max_step in dmd mbdf should be greater than 0");
    if(max_iter<=0)
        error->abort("max_iter in dmd mbdf should be greater than 0");
    if(max_order<=0)
        error->abort("max_order in dmd mbdf should be greater than 0");
    if(m_tol<=0.0)
        error->abort("m_tol in dmd mbdf should be greater than 0.0");
    if(a_tol<=0.0)
        error->abort("a_tol in dmd mbdf should be greater than 0.0");
    if(min_del_t<=0.0)
        error->abort("min_del_t in dmd mbdf should be greater than 0.0");
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
DMD_mbdf::~DMD_mbdf()
{
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_mbdf::allocate()
{
    DMDImplicit::allocate();
    
    vecs_1=new Vec<type0>*[max_order+2];
    for(int ivec=0;ivec<max_order+2;ivec++)
        vecs_1[ivec]=new Vec<type0>(atoms,c_dim);
    
    CREATE1D(alpha_y,max_order+1);
    CREATE1D(dalpha_y,max_order+1);
    CREATE1D(alph_err,max_order+2);
    CREATE1D(t,max_order+1);
    CREATE1D(y,max_order+1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD_mbdf::reset()
{
    DMDImplicit::reset();
    for(int i=0;i<max_order+1;i++)
        y[i]=vecs_1[i]->begin();
    dy=vecs_1[max_order+1]->begin();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_mbdf::deallocate()
{
    
    delete [] y;
    delete [] t;
    delete [] alph_err;
    delete [] dalpha_y;
    delete [] alpha_y;

    for(int ivec=0;ivec<max_order+2;ivec++)
        delete vecs_1[ivec];
    delete [] vecs_1;
    
    DMDImplicit::deallocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
type0 DMD_mbdf::est_dt()
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
    type0 del_t=MIN(2.0*a_tol/sum,1.0e-3*(max_t-tot_t));
    init_stp_adj(del_t);
    return del_t;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_mbdf::init()
{
    DMDImplicit::init();
    allocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_mbdf::fin()
{
    DMDImplicit::fin();
    deallocate();
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void DMD_mbdf::run()
{
    if(max_step==0)
        return;

    type0 del_t=initial_del_t,del_t_tmp;
    type0 cost,err;
    int q;
    int initial_phase,const_stps;
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
        memcpy(y[1],c,ncs*sizeof(type0));
        memcpy(dy,c_d,ncs*sizeof(type0));
        for(int i=0;i<max_order+1;i++) t[i]=0.0;
        t[0]=0.0;
        t[1]=-del_t;
        q=1;
        initial_phase=1;
        const_stps=0;
        
        while (istep<max_step && tot_t<max_t && !min_run)
        {
            err=1.0;
            while (err>=1.0)
            {
                interpolate(initial_phase,const_stps,del_t,q);
                solve_n_err(cost,err);
                err=MAX(err,cost);
                if(err<1.0)
                    continue;
                
                fail_stp_adj(initial_phase,const_stps,err,del_t,q);
                memcpy(c,y[0],ncs*sizeof(type0));
                memcpy(c_d,dy,ncs*sizeof(type0));
            }
            
            max_succ_q=MAX(max_succ_q,q);
            
            min_run=decide_min(istep,del_t);
            if(min_run) continue;

            del_t_tmp=del_t;
            ord_dt(initial_phase,const_stps,del_t,q);
            
            tmp_y=y[max_order];
            for(int i=max_order;i>0;i--)
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
 init
 --------------------------------------------*/
void DMD_mbdf::interpolate(int& init_phase,int& const_stps,type0& del_t,int& q)
{
    type0 tmp0;
    
    type0 k0;
    type0 err_prefac0,err_prefac1;
    
    type0* c=mapp->c->begin();
    
    int idof;
    int err_chk_lcl;
    int err_chk=1;
    
    
    while(err_chk)
    {
        k0=0.0;
        
        
        for(int i=0;i<q+1;i++)
        {
            k0+=1.0/(1.0-t[i]/del_t);
            //k0+=1.0/(del_t-t[i]);
            
            tmp0=1.0;
            for(int j=0;j<q+1;j++)
                if(i!=j)
                    tmp0*=(del_t-t[j])/(t[i]-t[j]);
            
            alpha_y[i]=tmp0;
        }
        
        for(int i=0;i<q+1;i++)
        {
            tmp0=k0-1.0/(1.0-t[i]/del_t);
            dalpha_y[i]=alpha_y[i]*tmp0/del_t;
        }
        
        tmp0=0.0;
        for(int i=0;i<q;i++)
            tmp0+=1.0/static_cast<type0>(i+1);
        beta=del_t/tmp0;
        
        
        
        
        err_chk_lcl=0;
        idof=0;
        while(idof<ncs && err_chk_lcl==0)
        {
            if(c[idof]>=0.0)
            {
                y_0[idof]=0.0;
                a[idof]=0.0;
                for(int j=0;j<q+1;j++)
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
                y_0[idof]=c[idof];
            
            
            idof++;
        }
        
        MPI_Allreduce(&err_chk_lcl,&err_chk,1,MPI_INT,MPI_MAX,world);
        
        if(err_chk)
        {
            const_stps=0;
            init_phase=0;
            
            if(q>1)
            {
                q--;
            }
            else
            {
                memcpy(y[1],y[0],ncs*sizeof(type0));
                memcpy(y_0,y[0],ncs*sizeof(type0));
                for(int i=0;i<ncs;i++)
                    if(c[i]>=0.0)
                        a[i]=-y_0[i];
                err_chk=0;
                
            }
        }
        
        err_prefac0=(del_t)/(del_t-t[q]);
        err_prefac1=err_prefac0;
        for(int i=0;i<q;i++)
        {
            err_prefac1+=(del_t)/(del_t-t[i])
            -1.0/static_cast<type0>(i+1);
        }
        err_prefac=MAX(err_prefac0,fabs(err_prefac1));
        
        if(err_chk)
            intp_rej++;
    }
    
    for(int i=0;i<ncs;i++)
    {
        if(y[0][i]>=0.0)
        {
            y_1[i]=y[0][i]+del_t*dy[i];
            
            if(y_1[i]<0.0 || y_1[i]>1.0)
                error->abort("exceeded the domain %e dt",y_1[i],del_t);
        }
        else
            y_1[i]=y[0][i];
    }
}
/*--------------------------------------------
 step addjustment after success
 --------------------------------------------*/
void DMD_mbdf::ord_dt(int& init_phase,int& const_stps,type0& del_t,int& q)
{
    if(init_phase)
    {
        type0 r=2.0;
        if(q<max_order)
            q++;
        
        if(r*del_t>max_t-tot_t)
        {
            init_phase=0;
            del_t=max_t-tot_t;
        }
        else
        {
            
            if(max_t-tot_t<=2.0*min_del_t)
            {
                init_phase=0;
                del_t=2.0*min_del_t;
            }
            else
            {
                if(r*del_t<min_del_t)
                {
                    init_phase=0;
                    del_t=min_del_t;
                }
                else if(r*del_t>=max_t-tot_t-min_del_t)
                {
                    init_phase=0;
                    del_t=max_t-tot_t-min_del_t;
                }
                else
                    del_t*=r;
            }
        }
        
        
        type0 rr_lcl=1.0,rr,tmp;
        type0* c=mapp->c->begin();
        type0* c_d=mapp->c_d->begin();
        
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
            {
                tmp=c[i]+rr_lcl*del_t*c_d[i];
                if(tmp>1.0)
                    rr_lcl=0.999*((1.0-c[i])/(del_t*c_d[i]));
                if(tmp<0.0)
                    rr_lcl=-0.999*((c[i])/(del_t*c_d[i]));
            }
        }
        
        MPI_Allreduce(&rr_lcl,&rr,1,MPI_TYPE0,MPI_MIN,world);
        if(rr!=1.0)
        {
            init_phase=0;
            if(max_t-tot_t<=2.0*min_del_t)
                error->abort("reached minimum order & del_t (%e)",del_t);
            else
            {
                if(rr*del_t<min_del_t)
                    error->abort("reached minimum order & del_t (%e)",del_t);
                else if(rr*del_t>=max_t-tot_t-min_del_t)
                    del_t=max_t-tot_t-min_del_t;
                else
                    del_t*=rr;
            }
        }
        return;
    }
    
    
    type0 terkm2=0.0,terkm1=0.0,terk=0.0,terkp1=0.0;
    type0 est;
    
    int terkm2_flag,terkm1_flag,terk_flag,terkp1_flag
    ,tmp_q=q;
    
    terk_flag=1;
    terkp1_flag=0;
    terkm1_flag=0;
    terkm2_flag=0;
    
    if(q>2)
        terkm2_flag=1;
    if(q>1)
        terkm1_flag=1;
    
    
    if(const_stps==q+2 && q<max_order)
        terkp1_flag=1;
    
    if(terkm2_flag)
        terkm2=err_est(q-1,del_t);
    
    if(terkm1_flag)
        terkm1=err_est(q,del_t);
    
    if(terk_flag)
        terk=err_est(q+1,del_t);
    
    if(terkp1_flag)
        terkp1=err_est(q+2,del_t);
    
    tmp_q=q;
    
    if(q>2)
    {
        if(MAX(terkm1,terkm2)<=terk)
        {
            est=terkm1/static_cast<type0>(q);
            q--;
        }
        else
        {
            if(terkp1<terk && terkp1_flag)
            {
                est=terkp1/static_cast<type0>(q+2);
                q++;
            }
            else
            {
                est=terk/static_cast<type0>(q+1);
            }
        }
    }
    else if(q==2)
    {
        if(terkm1<=0.5*terk)
        {
            est=terkm1/static_cast<type0>(q);
            q--;
        }
        else if(terkp1_flag && terk<terkm1)
        {
            if(terkp1<terk)
            {
                est=terkp1/static_cast<type0>(q+2);
                q++;
            }
            else
            {
                est=terk/static_cast<type0>(q+1);
            }
        }
        else
        {
            est=terk/static_cast<type0>(q+1);
        }
    }
    else
    {
        if(terkp1_flag && terkp1<0.5*terk)
        {
            est=terkp1/static_cast<type0>(q+2);
            q++;
        }
        else
        {
            est=terk/static_cast<type0>(q+1);
        }
    }
    
    type0 r=pow(0.5/est,1.0/static_cast<type0>(q+1));
    int const_stp_chk=0;
    
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
                if(r==1.0)
                    const_stp_chk=1;
            }
        }
    }
    
    type0 rr=precond_rat_adj(del_t);
    
    if(rr!=1.0)
    {
        const_stp_chk=0;
        if(max_t-tot_t<=2.0*min_del_t)
            error->abort("reached minimum order & del_t (%e)",del_t);
        else
        {
            if(rr*del_t<min_del_t)
                error->abort("reached minimum order & del_t (%e)",del_t);
            else if(rr*del_t>=max_t-tot_t-min_del_t)
                del_t=max_t-tot_t-min_del_t;
            else
                del_t*=rr;
        }
    }
    
    if(const_stp_chk && tmp_q==q)
        const_stps++;
    else
        const_stps=0;

    
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
type0 DMD_mbdf::err_est(int q,type0 del_t)
{
    type0 tmp0,err_lcl,err,k0;
    type0* c=mapp->c->begin();
    
    k0=1.0;
    for(int i=0;i<q;i++)
        k0*=static_cast<type0>(i+1)/(1.0-t[i]/del_t);
    
    for(int i=0;i<q;i++)
    {
        tmp0=1.0;
        for(int j=0;j<q;j++)
            if(j!=i)
                tmp0*=(del_t-t[j])/(t[i]-t[j]);
        
        alph_err[i]=-tmp0;
    }
    
    err_lcl=0.0;
    for(int i=0;i<ncs;i++)
    {
        if(c[i]>=0.0)
        {
            tmp0=c[i];
            for(int j=0;j<q;j++)
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
 step addjustment after failure
 --------------------------------------------*/
inline void DMD_mbdf::fail_stp_adj(int& init_phase,int& const_stps,type0 err,type0& del_t,int& q)
{
    const_stps=0;
    init_phase=0;
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
/*--------------------------------------------
 init
 --------------------------------------------*/
inline void DMD_mbdf::init_stp_adj(type0& del_t)
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
    
    type0 rr=precond_rat_adj(del_t);
    
    if(rr!=1.0)
    {
        if(max_t-tot_t<=2.0*min_del_t)
            error->abort("reached minimum order & del_t (%e)",del_t);
        else
        {
            if(rr*del_t<min_del_t)
                error->abort("reached minimum order & del_t (%e)",del_t);
            else if(rr*del_t>=max_t-tot_t-min_del_t)
                del_t=max_t-tot_t-min_del_t;
            else
                del_t*=rr;
        }
    }
}
/*--------------------------------------------
 step addjustment after success
 --------------------------------------------*/
inline type0 DMD_mbdf::precond_rat_adj(type0 del_t)
{
    type0 r_lcl=1.0,r,tmp;
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    
    for(int i=0;i<ncs;i++)
    {
        if(c[i]>=0.0)
        {
            tmp=c[i]+r_lcl*del_t*c_d[i];
            if(tmp>1.0)
                r_lcl=0.999*((1.0-c[i])/(del_t*c_d[i]));
            if(tmp<0.0)
                r_lcl=-0.999*((c[i])/(del_t*c_d[i]));
        }
    }
    
    MPI_Allreduce(&r_lcl,&r,1,MPI_TYPE0,MPI_MIN,world);
    return r;
}
