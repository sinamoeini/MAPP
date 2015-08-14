#include <stdlib.h>
#include "clock_bdf.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
//test
#include "timer.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_bdf::Clock_bdf(MAPP* mapp,int narg
,char** arg):Clock(mapp)
{

    max_order=6;
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
    CREATE1D(dy,dof_lcl);
    
    CREATE1D(t,max_order);
    CREATE1D(y,max_order);
    for(int i=0;i<max_order;i++)
        CREATE1D(y[i],dof_lcl);
    
    CREATE1D(alpha_y,max_order);
    CREATE1D(dalpha_y,max_order);
    CREATE1D(err_coef,max_order);
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_bdf::~Clock_bdf()
{
    
    for(int i=0;i<max_order;i++)
        if(dof_lcl)
            delete [] y[i];
    
    if(max_order)
    {
        delete [] y;
        delete [] t;
        delete [] alpha_y;
        delete [] dalpha_y;
        delete [] err_coef;
    }
    
    if(dof_lcl)
    {
        delete [] e_n;
        delete [] dy;
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_bdf::init()
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
    
    t[0]=0.0;
    memcpy(y[0],c,dof_lcl*sizeof(type0));
    memcpy(dy,c_d,dof_lcl*sizeof(type0));

    tot_t=0.0;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_bdf::fin()
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
 calculate the coefficients
 --------------------------------------------*/
void Clock_bdf::interpolate(type0& del_t,int& q)
{
    type0 tmp0,tmp1,k0,k1,k2,k3,k4;

    type0* c;
    atoms->vectors[c_n]->ret(c);
    
    int idof;
    int err_chk=1;
    int err_chk_lcl;
    
    while(err_chk)
    {

        
        k0=k3=1.0;
        k1=k2=0.0;
        k4=2.0/del_t;
        for(int i=1;i<q;i++)
        {
            k3*=-1.0/t[i];
            tmp0=del_t/t[i];
            tmp1=tmp0*tmp0;
            for(int j=1;j<q;j++)
            {
                if(i!=j)
                    tmp1*=(del_t-t[j])/(t[i]-t[j]);
            }
            alpha_y[i]=tmp1;
            
            k0*=-(del_t-t[i])/t[i];
            k1-=1.0/t[i];
            k2+=1.0/(del_t-t[i]);
        }
        
        alpha_y[0]=k0*(1.0-del_t*k1);
        alpha_dy=k0*del_t;
        
        dalpha_dy=k0*(1.0+del_t*k2);
        dalpha_y[0]=k0*(k2-k1-del_t*k1*k2);
        
        for(int i=1;i<q;i++)
            dalpha_y[i]=alpha_y[i]*(k4+k2-1.0/(del_t-t[i]));
        

        tmp0=0.0;
        for(int i=0;i<q;i++)
            tmp0+=1.0/(del_t-t[i]);
        beta=1.0/tmp0;
        
        err_chk_lcl=0;
        idof=0;
        while(idof<dof_lcl && err_chk_lcl==0)
        {
            if(c[idof]>=0.0)
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
        err_prefac=fabs(-1.0/((del_t*k2+1.0)*(1.0+k0)));
    }

    
}
/*--------------------------------------------
 find the new del_t and order
 --------------------------------------------*/
void Clock_bdf::ord_dt(type0& del_t,int& q
,int istep,type0 err)
{
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0 lo_err_lcl,hi_err_lcl;
    type0 tmp0,tmp1;   
    type0 lo_ratio,hi_ratio,ratio,del_t_tmp0,del_t_tmp1;
    type0 hi_err=0.0,lo_err=0.0;
    int lo_ord_avail=0,hi_ord_avail=0;
    
    
    if(istep!=0 && q<max_order)
        hi_ord_avail=1;
    if(q>1)
        lo_ord_avail=1;
    
    
    if(hi_ord_avail)
    {
        type0 c0,c1,c2,c3;
        type0 a0,a1;
        
        c0=1.0/del_t;
        c1=-del_t/t[1];
        c2=c3=1.0;
        
        for(int i=1;i<q;i++)
        {
            c0+=1.0/(del_t-t[i]);
            c1*=-(del_t-t[i])/t[i+1];
            c2*=-(del_t-t[i])/t[i];
            c3*=-t[i+1]/(t[1]-t[i+1]);
        }
        
        a0=(del_t-t[q])/(del_t*c0*static_cast<type0>((q+1)*(q+2)));
        a1=a0;
        a0*=-1.0/(del_t*(1.0+c2));
        a1*=-c1/(t[1]*(1.0+c3));
        
        hi_err_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
        {
            if(c[i]>=0.0)
            {
                tmp0=a0*(c[i]-y_0[i])+a1*e_n[i];
                hi_err_lcl+=tmp0*tmp0;
            }
        }
        hi_err=0.0;
        MPI_Allreduce(&hi_err_lcl,&hi_err,1
        ,MPI_TYPE0,MPI_SUM,world);
        hi_err=sqrt(hi_err/static_cast<type0>(dof_tot))/a_tol;
        
    }
    
    if(lo_ord_avail)
    {
        type0 c0;
        type0 err_dy,err_e_n;

        c0=1.0;
        tmp0=0.0;
        for(int i=0;i<q;i++)
        {
            c0*=del_t-t[i];
            tmp0+=1.0/(del_t-t[i]);
        }
        tmp0*=-static_cast<type0>(q-1)*(del_t-t[q-1]);
        err_e_n=1.0/tmp0;
        c0*=err_e_n;
        
        tmp0=1.0;
        tmp1=0.0;
        for(int i=1;i<q;i++)
        {
            tmp0*=-t[i];
            tmp1-=1.0/t[i];
        }
        err_dy=c0/tmp0;
        err_coef[0]=-c0*tmp1/tmp0;
        
        
        for(int i=1;i<q;i++)
        {
            tmp0=t[i]*t[i];
            for(int j=1;j<q;j++)
                if(i!=j)
                    tmp0*=(t[i]-t[j]);
            err_coef[i]=c0/tmp0;
        }
        
        
        lo_err_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
        {
            if(c[i]>=0.0)
            {
                tmp0=err_e_n*(c[i]-y_0[i]);
                tmp0+=err_dy*dy[i];
                for(int j=0;j<q;j++)
                    tmp0+=err_coef[j]*y[j][i];
                
                lo_err_lcl+=tmp0*tmp0;
            }
        }
        
        lo_err=0.0;
        MPI_Allreduce(&lo_err_lcl,&lo_err,1,MPI_TYPE0,MPI_SUM,world);
        lo_err=sqrt(lo_err/static_cast<type0>(dof_tot))/a_tol;
    }
    
    for(int i=0;i<dof_lcl;i++)
        if(c[i]>=0.0)
            e_n[i]=c[i]-y_0[i];

    if(hi_ord_avail && lo_ord_avail)
    {
        lo_ratio=pow(0.5/lo_err,1.0/static_cast<type0>(q));
        hi_ratio=pow(0.5/hi_err,1.0/static_cast<type0>(q+2));
        ratio=pow(0.5/err,1.0/static_cast<type0>(q+1));
        
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
        
    }
    else if(hi_ord_avail==1&& lo_ord_avail==0)
    {
        hi_ratio=pow(0.5/hi_err,1.0/static_cast<type0>(q+2));
        ratio=pow(0.5/err,1.0/static_cast<type0>(q+1));
        
        if(hi_ratio>ratio)
        {
            q++;
            ratio=hi_ratio;
        }
    }
    else if(hi_ord_avail==0 && lo_ord_avail==1)
    {
        lo_ratio=pow(0.5/lo_err,1.0/static_cast<type0>(q));
        ratio=pow(0.5/err,1.0/static_cast<type0>(q+1));
        if(lo_ratio>ratio)
        {
            q--;
            ratio=lo_ratio;
        }
    }
    else
        ratio=pow(0.5/err,1.0/static_cast<type0>(q+1));
    
    
    if(ratio>=2.0)
        ratio=2.0;
    else if(ratio<=0.5)
        ratio=0.5;
    else
        ratio=1.0;
    

    del_t_tmp0=ratio*del_t;
    del_t_tmp1=del_t;
    if(del_t_tmp0>max_del_t)
        del_t=max_del_t;
    else if (del_t_tmp0<min_del_t)
        del_t=min_del_t;
    else
        del_t=del_t_tmp0;
    
    if(del_t+del_t_tmp1+tot_t>max_t)
        del_t=max_t-tot_t-del_t_tmp1;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_bdf::run()
{
    if(max_step==0)
        return;
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
     
    type0* tmp_y;
    
    type0 del_t=initial_del_t,del_t_tmp,max_err;
    type0 ratio=1.0,cost,err=0.0;
    int ord=1;
    int err_chk;
    int istep;
    

    istep=0;
    while (istep<max_step && tot_t<max_t)
    {
        err_chk=1;
        while (err_chk)
        {
            interpolate(del_t,ord);
            solve_n_err(cost,err);
            
            
            if(err<1.0 && cost<1.0)
                err_chk=0;
            
            if(err_chk)
            {
                if(del_t==min_del_t)
                {
                    if (ord>1)
                        ord--;
                    else
                        error->abort("Fucked up");
                }
                else
                {
                    max_err=MAX(err,cost);
                    ratio=pow(0.5/max_err,1.0/static_cast<type0>(ord+1));
                    
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

        if(thermo->test_prev_step()|| istep==max_step-1 || tot_t==max_t)
        {
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            thermo->update(time_idx,tot_t);
        }
        
        del_t_tmp=del_t;
        ord_dt(del_t,ord,istep,err);
        
        tmp_y=y[max_order-1];
        for(int i=max_order-1;i>0;i--)
        {
            t[i]=t[i-1]-del_t_tmp;
            y[i]=y[i-1];
        }
        
        y[0]=tmp_y;

        
        memcpy(y[0],c,dof_lcl*sizeof(type0));
        memcpy(dy,c_d,dof_lcl*sizeof(type0));
        
        step_no++;
        istep++;
    }
    
    
}
