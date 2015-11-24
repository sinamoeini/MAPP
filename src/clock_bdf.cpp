#include <stdlib.h>
#include <limits>
#include "clock_bdf.h"
#include "ff.h"
#include "error.h"
#include "memory.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_bdf::Clock_bdf(MAPP* mapp,int nargs
,char** args):ClockImplicit(mapp)
{
    max_order=6;
    
    if(nargs>2)
    {
        if(nargs%2!=0)
            error->abort("every keyword in clock bdf should be followed by it's value");
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
                    error->abort("initial_del_t in clock bdf should be greater than 0.0");
                iarg++;
            }
            else
                error->abort("unknown keyword in clock bdf: %s",args[iarg]);
        }
    }
    
    if(min_gamma<0.0)
        error->abort("min_gamma in clock bdf should be greater than 0.0");
    if(gamma_red<=0.0 || gamma_red>=1.0)
        error->abort("gamma_red in clock bdf should be between 0.0 & 1.0");
    if(slope<=0.0 || slope>=1.0)
        error->abort("slope in clock bdf should be between 0.0 & 1.0");
    if(max_step<0)
        error->abort("max_step in clock bdf should be greater than 0");
    if(max_iter<=0)
        error->abort("max_iter in clock bdf should be greater than 0");
    if(max_order<=0)
        error->abort("max_order in clock bdf should be greater than 0");
    if(m_tol<=0.0)
        error->abort("m_tol in clock bdf should be greater than 0.0");
    if(a_tol<=0.0)
        error->abort("a_tol in clock bdf should be greater than 0.0");
    if(min_del_t<=0.0)
        error->abort("min_del_t in clock bdf should be greater than 0.0");
    if(max_t<=0.0)
        error->abort("max_t in clock bdf should be greater than 0.0");
    
    
    if(initial_del_t!=-1.0 && initial_del_t>max_t)
        error->abort("max_t in clock bdf should be greater than initial_del_t");
    if(min_del_t>max_t)
        error->abort("max_t in clock bdf should be greater than min_del_t");
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_bdf::~Clock_bdf()
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void Clock_bdf::allocate()
{
    ClockImplicit::allocate();
    vecs_1=new Vec<type0>*[max_order+2];
    for(int ivec=0;ivec<max_order+2;ivec++)
        vecs_1[ivec]=new Vec<type0>(atoms,c_dim);
    
    
    CREATE1D(t,max_order);
    CREATE1D(y,max_order);
    CREATE1D(alpha_y,max_order);
    CREATE1D(dalpha_y,max_order);
    
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Clock_bdf::reset()
{
    for(int i=0;i<max_order;i++)
        y[i]=vecs_1[i]->begin();
    dy=vecs_1[max_order]->begin();
    e_n=vecs_1[max_order+1]->begin();
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void Clock_bdf::deallocate()
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
    
    ClockImplicit::deallocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_bdf::init()
{
    ClockImplicit::init();
    allocate();
    
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    rectify(c_d);
    
    if(initial_del_t<0.0)
    {
        type0 sum=forcefield_dmd->ddc_norm_timer()/sqrt(static_cast<type0>(dof_tot));
        initial_del_t=MIN(sqrt(2.0*a_tol/sum),1.0e-3*max_t);
    }
    
    init_stp_adj(initial_del_t);
    
    for(int i=0;i<max_order;i++)
        t[i]=0.0;
    memcpy(y[0],c,dof_lcl*sizeof(type0));
    memcpy(dy,c_d,dof_lcl*sizeof(type0));
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_bdf::fin()
{
    ClockImplicit::fin();
    deallocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_bdf::run()
{
    if(max_step==0)
        return;
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    
    type0* tmp_y;
    
    type0 del_t=initial_del_t,del_t_tmp;
    type0 cost,err=0.0;
    int q=1;
    int err_chk;
    int istep;
    bool thermo_flag;
    
    istep=0;
    while (istep<max_step && tot_t<max_t)
    {
        err_chk=1;
        while (err_chk)
        {
            interpolate(del_t,q);
            solve_n_err(cost,err);
            
            if(err<1.0 && cost<1.0)
                err_chk=0;
            
            if(err_chk)
                fail_stp_adj(MAX(cost,err),del_t,q);
        }
        
        max_succ_dt=MAX(max_succ_dt,del_t);
        max_succ_q=MAX(max_succ_q,q);
        tot_t+=del_t;
        
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        thermo_flag=(thermo->test_prev_step()|| istep==max_step-1 || tot_t>=max_t);
        

        del_t_tmp=del_t;
        if(min==NULL)
        {
            if(thermo_flag)
                forcefield_dmd->force_calc_timer(true);

            ord_dt(del_t,q,err);
        }
        else
        {
            type0 f_norm=min->calc_ave_f_norm();
            if(f_norm-init_f_norm>=f_tol)
            {
                forcefield_dmd->dynamic_flag=true;
                min->run();
                nmin++;
                forcefield_dmd->dynamic_flag=false;
                init_f_norm=min->calc_ave_f_norm();
                reset();
                c=mapp->c->begin();
                c_d=mapp->c_d->begin();
                neighbor_dmd->create_2nd_list();
                forcefield_dmd->dc_timer();
                type0 sum=forcefield_dmd->ddc_norm_timer()/sqrt(static_cast<type0>(dof_tot));
                type0 del_t_=MIN(sqrt(2.0*a_tol/sum),1.0e-3*(max_t-tot_t));
                init_stp_adj(del_t_);
                for(int i=0;i<max_order;i++)
                    t[i]=0.0;
                q=1;
            }
            else
                ord_dt(del_t,q,err);

        }


        
        tmp_y=y[max_order-1];
        for(int i=max_order-1;i>0;i--)
        {
            t[i]=t[i-1]-del_t_tmp;
            y[i]=y[i-1];
        }
        
        y[0]=tmp_y;
        memcpy(y[0],c,dof_lcl*sizeof(type0));
        memcpy(dy,c_d,dof_lcl*sizeof(type0));
        
        if(thermo_flag)
        {
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            thermo->update(time_idx,tot_t);
        }
        step_no++;
        istep++;
    }
}
/*--------------------------------------------
 calculate the coefficients
 --------------------------------------------*/
void Clock_bdf::interpolate(type0& del_t,int& q)
{
    type0 tmp0,tmp1,c0,c1,c2;

    int idof;
    int err_chk=1;
    int err_chk_lcl;
    
    while(err_chk)
    {
        c0=1.0;
        c1=c2=0.0;
        for(int i=1;i<q;i++)
        {
            tmp0=del_t/t[i];
            tmp1=t[i]/del_t;
            c0*=1.0-tmp0;
            c1+=tmp0;
            c2+=1.0/(1.0-tmp1);
        }
        alpha_y[0]=c0*(1.0+c1);
        for(int i=1;i<q;i++)
        {
            tmp1=del_t/t[i];
            tmp0=tmp1*tmp1;
            for(int j=1;j<q;j++)
                if(i!=j)
                    tmp0*=(del_t-t[j])/(t[i]-t[j]);
            alpha_y[i]=tmp0;
        }
        dalpha_y[0]=c0*(c1+c2+c1*c2)/del_t;
        for(int i=1;i<q;i++)
            dalpha_y[i]=alpha_y[i]*(2.0+c2-1.0/(1.0-t[i]/del_t))/del_t;
        
        
        alpha_dy=c0*del_t;
        dalpha_dy=c0*(1.0+c2);
        
        beta=del_t/(1.0+c2);
        
        err_chk_lcl=0;
        idof=0;
        while(idof<dof_lcl && err_chk_lcl==0)
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
                        for(int i=0;i<dof_lcl;i++)
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
        err_prefac=fabs(-1.0/((1.0+c0)*(1.0+c2)));
        
        if(err_chk)
            intp_rej++;
    }
}
/*--------------------------------------------
 find the new del_t and order
 --------------------------------------------*/
void Clock_bdf::ord_dt(type0& del_t,int& q
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
    
    for(int i=0;i<dof_lcl;i++)
        if(c[i]>=0.0)
            e_n[i]=c[i]-y_0[i];

}
/*--------------------------------------------
 error calculation
 --------------------------------------------*/
void Clock_bdf::ratio_calc(int q,type0 del_t,
type0& lo_ratio,type0& hi_ratio)
{
    type0* c=mapp->c->begin();
    type0 lo_err,hi_err,lo_err_lcl,hi_err_lcl;
    type0 tmp0,tmp1;
    type0 c0,c2;
    
    lo_ratio=hi_ratio=0.0;
    
    c0=1.0;
    c2=0.0;
    for(int i=1;i<q;i++)
    {
        tmp0=del_t/t[i];
        tmp1=t[i]/del_t;
        c0*=1.0-tmp0;
        c2+=1.0/(1.0-tmp1);
    }
    
    
    if(q>1)
    {
        type0 a0;
        a0=fabs(-1.0/(static_cast<type0>(q-1)*(1.0+c2)*(1.0-t[q-1]/del_t)));
        
        lo_err_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
        {
            if(c[i]>=0.0)
            {
                tmp0=c[i];
                tmp0-=c0*y[0][i];
                for(int j=1;j<q;j++)
                    tmp0-=(t[j]/del_t)*alpha_y[j]*y[j][i];
                
                lo_err_lcl+=tmp0*tmp0;
            }
        }
        
        lo_err=0.0;
        MPI_Allreduce(&lo_err_lcl,&lo_err,1,MPI_TYPE0,MPI_SUM,world);
        lo_err=sqrt(lo_err/static_cast<type0>(dof_tot))/a_tol;
        lo_err*=a0;
        lo_ratio=pow(0.5/lo_err,1.0/static_cast<type0>(q));
    }
    
    if(q<max_order && t[1]!=0.0)
    {
        type0 a0,a1;
        type0 c3;
        
        c3=1.0;
        for(int i=0;i<q-1;i++)
            c3*=1.0/(1.0-t[1]/t[i+2]);
        
        a0=fabs(-(1.0-t[q]/del_t)/(static_cast<type0>((q+1)*(q+2))*(1.0+c0)*(1.0+c2)));
        a1=c0*(1.0+c0)*del_t*del_t/((1.0+c3)*t[1]*t[q]);
        
        hi_err_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
        {
            if(c[i]>=0.0)
            {
                tmp0=c[i]-y_0[i]+a1*e_n[i];
                hi_err_lcl+=tmp0*tmp0;
            }
        }
        hi_err=0.0;
        MPI_Allreduce(&hi_err_lcl,&hi_err,1,MPI_TYPE0,MPI_SUM,world);
        hi_err=sqrt(hi_err/static_cast<type0>(dof_tot))/a_tol;
        hi_err*=a0;
        hi_ratio=pow(0.5/hi_err,1.0/static_cast<type0>(q+2));
    }
    
}
/*--------------------------------------------
 step addjustment after failure
 --------------------------------------------*/
inline void Clock_bdf::fail_stp_adj(type0 err,type0& del_t,int& q)
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
inline void Clock_bdf::init_stp_adj(type0& del_t)
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

