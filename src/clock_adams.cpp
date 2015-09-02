#include <stdlib.h>
#include "clock_adams.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_adams::Clock_adams(MAPP* mapp,int narg
,char** arg):Clock(mapp)
{
    
    max_order=12;
    min_del_t=1.0e-12;
    max_del_t=1.0e4;
    initial_del_t=-1.0;
    max_t=1.0e7;
    
    
    if(narg>2)
    {
        if(narg%2!=0)
            error->abort("every keyword in clock bdf should be followed by it's value");
        int iarg=2;
        while(iarg<narg)
        {
            if(strcmp(arg[iarg],"min_gamma")==0)
            {
                iarg++;
                min_gamma=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"red_gamma")==0)
            {
                iarg++;
                gamma_red=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"slope")==0)
            {
                iarg++;
                slope=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_step")==0)
            {
                iarg++;
                max_step=atoi(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_iter")==0)
            {
                iarg++;
                max_iter=atoi(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_order")==0)
            {
                iarg++;
                max_order=atoi(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"m_tol")==0)
            {
                iarg++;
                m_tol=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"a_tol")==0)
            {
                iarg++;
                a_tol=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"min_del_t")==0)
            {
                iarg++;
                min_del_t=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_del_t")==0)
            {
                iarg++;
                max_del_t=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"initial_del_t")==0)
            {
                iarg++;
                initial_del_t=atof(arg[iarg]);
                if(initial_del_t<=0.0)
                    error->abort("initial_del_t in clock bdf should be greater than 0.0");
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_t")==0)
            {
                iarg++;
                max_t=atof(arg[iarg]);
                iarg++;
            }
            else
                error->abort("unknown keyword in clock bdf: %s",arg[iarg]);
        }
    }
    
    if(min_gamma<0.0)
        error->abort("min_gamma in clock bdf should be greater than 0.0");
    if(gamma_red<=0.0 || gamma_red>=1.0)
        error->abort("gamma_red in clock bdf should be between 0.0 & 1.0");
    if(slope<=0.0 || slope>=1.0)
        error->abort("slope in clock bdf should be between 0.0 & 1.0");
    if(max_step<=0)
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
    if(max_del_t<=0.0)
        error->abort("max_del_t in clock bdf should be greater than 0.0");
    if(max_t<=0.0)
        error->abort("max_t in clock bdf should be greater than 0.0");
    
    
    if(initial_del_t!=-1.0 && initial_del_t>max_t)
        error->abort("max_t in clock bdf should be greater than initial_del_t");
    if(min_del_t>max_t)
        error->abort("max_t in clock bdf should be greater than min_del_t");
    if(max_del_t>max_t)
        error->abort("max_t in clock bdf should be greater than max_del_t");
    if(min_del_t>=max_del_t)
        error->abort("max_del_t in clock bdf should be greater than min_del_t");
    
    
    CREATE1D(y_0,dof_lcl);
    CREATE1D(e_n,dof_lcl);
    CREATE1D(y,dof_lcl);
    
    CREATE1D(t,max_order);
    CREATE1D(dy,max_order);
    for(int i=0;i<max_order;i++)
        CREATE1D(dy[i],dof_lcl);
    
    CREATE1D(alpha_dy,max_order);
    CREATE1D(dalpha_dy,max_order);
    
    int n=1+max_order/2;
    CREATE1D(xi,n);
    CREATE1D(wi,n);
    for(int i=0;i<n;i++)
    {
        CREATE1D(xi[i],i+1);
        CREATE1D(wi[i],i+1);
    }
    
    XMath* xmath= new XMath(mapp);
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
Clock_adams::~Clock_adams()
{
    for(int i=0;i<max_order;i++)
        if(dof_lcl)
            delete [] dy[i];
    
    if(max_order)
    {
        delete [] dy;
        delete [] t;
        delete [] alpha_dy;
        delete [] dalpha_dy;
    }
    
    if(dof_lcl)
    {
        delete [] e_n;
        delete [] y;
    }
    
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
 init
 --------------------------------------------*/
void Clock_adams::interpolate(type0& del_t,int& q)
{
    int n;
    type0 tmp0,tmp1;
    type0 c0,c1,c2;
    
    type0* c;
    atoms->vectors[c_n]->ret(c);
    
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
        
        
        c0=1.0;
        for(int i=0;i<q-1;i++)
            c0*=1.0-t[i]/del_t;
        
        c1=0.0;
        c2=0.0;
        n=1+q/2;
        for(int i=0;i<n;i++)
        {
            tmp0=wi[n-1][i];
            for(int j=0;j<q-1;j++)
                tmp0*=xi[n-1][i]-t[j]/del_t;
            c1+=tmp0;
            c2+=tmp0*(xi[n-1][i]-1.0);
        }
        
        beta=del_t*c1/c0;
        
        err_chk_lcl=0;
        idof=0;
        while(idof<dof_lcl && err_chk_lcl==0)
        {
            if(c[idof]>=0.0)
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
                y_0[idof]=c[idof];
            
            
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
                del_t*=0.8;
                del_t=MAX(del_t,min_del_t);
            }
        }
        err_prefac=fabs((c2/c1)/(1.0-t[q-1]/del_t));
    }
}
/*--------------------------------------------
 find the new del_t and order
 --------------------------------------------*/
void Clock_adams::ord_dt(type0& del_t,int& q
,type0 err)
{
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0 lo_ratio,hi_ratio,ratio;

    ratio=pow(0.5/err,1.0/static_cast<type0>(q+1));
    ratio_calc(q,del_t,lo_ratio,hi_ratio);
    
    if(hi_ratio>lo_ratio && hi_ratio>ratio)
    {
        q++;
        ratio=hi_ratio;
    }
    else if(lo_ratio>hi_ratio && lo_ratio>ratio)
    {
        q--;
        ratio=lo_ratio;
    }
    
    if(ratio>=2.0)
        ratio=2.0;
    else if(ratio<=0.5)
        ratio=0.5;
    else
        ratio=1.0;
    
    del_t=MAX(min_del_t,MIN(ratio*del_t,MIN(max_t-tot_t,max_del_t)));
}
/*--------------------------------------------
 find the new del_t and order
 --------------------------------------------*/
void Clock_adams::ratio_calc(int q,type0 del_t,
type0& lo_ratio,type0& hi_ratio)
{
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0 lo_err,hi_err,lo_err_lcl,hi_err_lcl;
    type0 tmp0;
    type0 c1;
    int n;
    
    lo_ratio=hi_ratio=0.0;

    c1=0.0;
    n=1+q/2;
    for(int i=0;i<n;i++)
    {
        tmp0=wi[n-1][i];
        for(int j=0;j<q-1;j++)
            tmp0*=xi[n-1][i]-t[j]/del_t;
        c1+=tmp0;
    }
    
    
    if(q>1)
    {
        type0 a0,a1;
        type0 c0,c4;
        
        c0=1.0;
        for(int i=0;i<q-1;i++)
            c0*=1.0-t[i]/del_t;
        
        c4=0.0;
        n=1+(q-1)/2;
        for(int i=0;i<n;i++)
        {
            tmp0=wi[n-1][i]*(xi[n-1][i]-1.0);
            for(int j=0;j<q-2;j++)
                tmp0*=xi[n-1][i]-t[j]/del_t;
            c4+=tmp0;
        }
        
        a0=fabs(c4/c1);
        a1=c1*del_t/(c0*(1.0-t[q-1]/del_t));
        
        lo_err_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
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
        lo_err=sqrt(lo_err/static_cast<type0>(dof_tot))/a_tol;
        lo_err*=a0;
        lo_ratio=pow(0.5/lo_err,1.0/static_cast<type0>(q));
    }
    
    if(q<max_order && t[1]!=0.0)
    {
        type0 a0,a1;
        type0 c2,c3,tmp1;
        
        c2=0.0;
        c3=0.0;
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
            c2+=tmp0;
            c3+=tmp1;
        }
        
        a0=fabs((c2/c1)/((1.0-t[q-1]/del_t)*static_cast<type0>(q+1)));
        a1*=-(c1/c3)*(1.0-t[q-1]/del_t)/(t[q]/t[1]);
        for(int i=0;i<q+1;i++)
            a1*=-del_t/t[1];
        
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
 init
 --------------------------------------------*/
void Clock_adams::init()
{
    old_skin=atoms->skin;
    old_comm_mode=atoms->comm_mode;
    atoms->chng_skin(0.0);
    atoms->comm_mode=COMM_MODE_5;
    
    type0* c;
    type0* c_d;
    
    int f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<type0>(0,atoms->vectors[0]->dim,"f");
    
    int id_n=atoms->find("id");
    if(cdof_n>-1)
    {
        int dof_n=atoms->find("dof");
        vecs_comm=new VecLst(mapp,6,0,c_n,c_d_n,cdof_n,dof_n,id_n);
        
    }
    else
    {
        vecs_comm=new VecLst(mapp,4,0,c_n,c_d_n,id_n);
    }
    
    vecs_comm->add_update(0);
    atoms->init(vecs_comm);
    
    
    forcefield->create_2nd_neigh_lst_timer();
    forcefield->force_calc_timer(1,nrgy_strss);
    forcefield->c_d_calc_timer(0,nrgy_strss);
    
    thermo->update(fe_idx,nrgy_strss[0]);
    thermo->update(stress_idx,6,&nrgy_strss[1]);
    thermo->update(time_idx,0.0);
    
    thermo->init();
    
    if(write!=NULL)
        write->init();
    
    atoms->vectors[c_n]->ret(c);
    atoms->vectors[c_d_n]->ret(c_d);
    
    rectify(c_d);
    if(initial_del_t<0.0)
    {
        type0 sum,sum_lcl;
        sum_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
            if(c[i]!=-1.0)
                sum_lcl+=c_d[i]*c_d[i];
        sum=0.0;
        MPI_Allreduce(&sum_lcl,&sum,1,MPI_TYPE0,MPI_SUM,world);
        sum=sqrt(sum/static_cast<type0>(dof_tot));
        initial_del_t=MAX(a_tol/sum,min_del_t);
        initial_del_t=MIN(initial_del_t,max_del_t);
        if(initial_del_t>max_t)
            initial_del_t=max_t;
    }
    
    for(int i=0;i<max_order;i++)
        t[i]=0.0;
    memcpy(y,c,dof_lcl*sizeof(type0));
    memcpy(dy[0],c_d,dof_lcl*sizeof(type0));
    
    tot_t=0.0;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_adams::fin()
{
    
    if(write!=NULL)
        write->fin();
    thermo->fin();
    atoms->fin();
    delete vecs_comm;
    
    atoms->chng_skin(old_skin);
    atoms->comm_mode=old_comm_mode;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_adams::run()
{
    if(max_step==0)
        return;
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    
    type0* tmp_dy;
    
    type0 del_t=initial_del_t,del_t_tmp,max_err;
    type0 ratio=1.0,cost,err=0.0;
    int q=1;
    int err_chk;
    int istep;
    
    
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
            {
                if(del_t==min_del_t)
                {
                    if (q>1)
                        q--;
                    else
                        error->abort("reached minimum order & del_t (%e)",min_del_t);
                }
                else
                {
                    max_err=MAX(err,cost);
                    ratio=pow(0.5/max_err,1.0/static_cast<type0>(q+1));
                    
                    if(ratio<=0.5)
                        ratio=0.5;
                    
                    del_t_tmp=del_t*ratio;
                    if(del_t_tmp<min_del_t)
                        del_t=min_del_t;
                    else
                        del_t=del_t_tmp;
                }
            }
        }
        
        tot_t+=del_t;
        
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step()|| istep==max_step-1 || tot_t>=max_t)
        {
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            thermo->update(time_idx,tot_t);
        }
        
        del_t_tmp=del_t;
        
        ord_dt(del_t,q,err);
        
        tmp_dy=dy[max_order-1];
        for(int i=max_order-1;i>0;i--)
        {
            t[i]=t[i-1]-del_t_tmp;
            dy[i]=dy[i-1];
        }
        
        dy[0]=tmp_dy;
        memcpy(dy[0],c_d,dof_lcl*sizeof(type0));
        memcpy(y,c,dof_lcl*sizeof(type0));
        
        step_no++;
        istep++;
    }
    
}
