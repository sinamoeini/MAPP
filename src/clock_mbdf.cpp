#include <stdlib.h>
#include "clock_mbdf.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_mbdf::Clock_mbdf(MAPP* mapp,int narg
                       ,char** arg):Clock(mapp)
{
    max_order=5;
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
    
    
    
    
    CREATE1D(dy,dof_lcl);
    
    CREATE1D(t,max_order+1);
    CREATE1D(alpha_y,max_order+1);
    CREATE1D(dalpha_y,max_order+1);
    CREATE1D(y,max_order+1);
    for(int i=0;i<max_order+1;i++)
        CREATE1D(y[i],dof_lcl);
    
    CREATE1D(alph_err,max_order+2);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_mbdf::~Clock_mbdf()
{
    
    for(int i=0;i<max_order+1;i++)
        if(dof_lcl)
            delete [] y[i];
    
    delete [] alph_err;
    delete [] y;
    delete [] t;
    delete [] alpha_y;
    delete [] dalpha_y;
    
    
    if(dof_lcl)
    {
        delete [] dy;
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_mbdf::init()
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
            if(c[i]>=0.0)
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
    t[1]=-initial_del_t;
    memcpy(y[0],c,dof_lcl*sizeof(type0));
    memcpy(y[1],c,dof_lcl*sizeof(type0));
    memcpy(dy,c_d,dof_lcl*sizeof(type0));
    
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_mbdf::fin()
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
void Clock_mbdf::interpolate(type0& del_t,int& q,int initial_phase)
{
    type0 tmp0,tmp1,tmp2;
    type0 curr_t;
    type0 c0;
    type0 err_prefac0,err_prefac1;
    
    type0* c;
    atoms->vectors[c_n]->ret(c);
    
    int idof;
    int err_chk_lcl;
    int err_chk=1;
    
    
    while(err_chk)
    {
        c0=0.0;
        curr_t=t[0]+del_t;
        
        for(int i=0;i<q+1;i++)
        {
            c0+=1.0/(curr_t-t[i]);
            
            tmp0=tmp1=1.0;
            tmp2=0.0;
            for(int j=0;j<q+1;j++)
                if(i!=j)
                    tmp0*=(curr_t-t[j])/(t[i]-t[j]);
            
            alpha_y[i]=tmp0;
            dalpha_y[i]=tmp0;
            
        }
        
        for(int i=0;i<q+1;i++)
        {
            tmp0=c0-1.0/(curr_t-t[i]);
            dalpha_y[i]*=tmp0;
        }
        
        tmp0=0.0;
        for(int i=0;i<q;i++)
            tmp0+=1.0/static_cast<type0>(i+1);
        beta=del_t/tmp0;
        
        
        
        
        err_chk_lcl=0;
        idof=0;
        while(idof<dof_lcl && err_chk_lcl==0)
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
            if(initial_phase)
            {
                initial_phase=0;
            }
            const_stps=0;
            if(q>1)
            {
                q--;
            }
            else
            {
                memcpy(y[1],y[0],dof_lcl*sizeof(type0));
                
            }
        }
        
        err_prefac0=(curr_t-t[0])/(curr_t-t[q]);
        err_prefac1=err_prefac0;
        for(int i=0;i<q;i++)
        {
            err_prefac1+=(curr_t-t[0])/(curr_t-t[i])
            -1.0/static_cast<type0>(i+1);
        }
        err_prefac=MAX(err_prefac0,fabs(err_prefac1));
    }
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
void Clock_mbdf::ord_dt(type0& del_t,int& q
                        ,int init_ph)
{
    type0 del_t_tmp0,del_t_tmp1;
    del_t_tmp1=del_t;
    if(init_ph)
    {
        del_t=MIN(del_t*2.0,max_del_t);
        if(q<max_order)
            q++;
        if(del_t+del_t_tmp1+t[0]>max_t)
            del_t=max_t-t[0]-del_t_tmp1;
        return;
    }
    
    type0 terkm2=0.0,terkm1=0.0,terk=0.0,terkp1=0.0;
    type0 tmp_err,ratio;
    
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
    
    
    
    if(q>2)
    {
        if(MAX(terkm1,terkm2)<=terk)
        {
            tmp_q=q-1;
            tmp_err=terkm1;
        }
        else
        {
            if(terkp1<terk && terkp1_flag)
            {
                tmp_q=q+1;
                tmp_err=terkp1;
            }
            else
                tmp_err=terk;
        }
    }
    else if(q==2)
    {
        if(terkm1<=0.5*terk)
        {
            tmp_q=q-1;
            tmp_err=terkm1;
        }
        else if(terk<terkm1 && terkp1_flag)
        {
            if(terkp1<terk)
            {
                tmp_q=q+1;
                tmp_err=terkp1;
            }
            else
                tmp_err=terk;
        }
        else
            tmp_err=terk;
    }
    else
    {
        if(terkp1<0.5*terk && terkp1_flag)
        {
            tmp_q=q+1;
            tmp_err=terkp1;
        }
        else
            tmp_err=terk;
    }
    
    
    ratio=pow(0.5/tmp_err,1.0/static_cast<type0>(q+1));
    
    if(ratio>=2.0)
    {
        const_stps=0;
        ratio=2.0;
    }
    else if(ratio<=0.5)
    {
        ratio=0.5;
        const_stps=0;
    }
    else
    {
        ratio=1.0;
        if(tmp_q==q)
            const_stps++;
    }
    
    del_t_tmp0=ratio*del_t;
    if(del_t_tmp0>max_del_t)
        del_t=max_del_t;
    else if (del_t_tmp0<min_del_t)
        del_t=min_del_t;
    else
        del_t=del_t_tmp0;
    
    q=tmp_q;
    if(del_t+del_t_tmp1+t[0]>max_t)
        del_t=max_t-t[0]-del_t_tmp1;
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void Clock_mbdf::run()
{
    if(max_step==0)
        return;
    
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    
    type0* tmp_y;
    
    type0 del_t=initial_del_t,del_t_tmp;
    type0 err,max_err,cost;
    type0 ratio;
    int ord=1;
    int err_chk;
    int initial_phase;
    int istep;
    
    initial_phase=1;
    const_stps=0;
    istep=0;
    
    while (istep<max_step && t[0]<max_t)
    {
        err_chk=1;
        
        while (err_chk)
        {
            
            interpolate(del_t,ord,initial_phase);
            solve_n_err(cost,err);
            
            if(err<1.0 && cost<1.0)
                err_chk=0;
            
            if(err_chk)
            {
                if(initial_phase) initial_phase=0;
                const_stps=0;
                max_err=MAX(err,cost);
                ratio=pow(0.5/max_err,1.0/static_cast<type0>(ord+1));
                
                if(ratio<0.5)
                    ratio=0.5;
                else if(ratio>0.9)
                    ratio=0.9;
                
                del_t_tmp=del_t*ratio;
                if(del_t_tmp<min_del_t)
                    del_t=min_del_t;
                else
                    del_t=del_t_tmp;
            }
        }
        
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step()|| istep==max_step-1 || t[0]+del_t==max_t)
        {
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            thermo->update(time_idx,t[0]+del_t);
        }
        del_t_tmp=del_t;
        ord_dt(del_t,ord,initial_phase);
        
        tmp_y=y[max_order];
        for(int i=max_order;i>0;i--)
        {
            t[i]=t[i-1];
            y[i]=y[i-1];
        }
        
        y[0]=tmp_y;
        t[0]+=del_t_tmp;
        
        memcpy(y[0],c,dof_lcl*sizeof(type0));
        memcpy(dy,c_d,dof_lcl*sizeof(type0));
        
        step_no++;
        istep++;
    }
    
    
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
type0 Clock_mbdf::err_est(int q,type0 dt)
{
    type0 tmp0,tmp1,err_lcl,err;
    type0* c;
    atoms->vectors[c_n]->ret(c);
    
    
    tmp0=1.0;
    for(int i=0;i<q;i++)
        tmp0*=static_cast<type0>(i+1)*dt;
    
    for(int i=0;i<q;i++)
    {
        tmp1=t[i]-(t[0]+dt);
        for(int j=0;j<q;j++)
            if(i!=j)
                tmp1*=t[i]-t[j];
        
        alph_err[i+1]=tmp0/tmp1;
    }
    
    tmp1=1.0;
    for(int i=0;i<q;i++)
        tmp1*=t[0]+dt-t[i];
    alph_err[0]=tmp0/tmp1;
    
    
    err_lcl=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        if(c[i]>=0.0)
        {
            tmp0=alph_err[0]*c[i];
            for(int j=0;j<q;j++)
                tmp0+=alph_err[j+1]*y[j][i];
            err_lcl+=tmp0*tmp0;
        }
    }
    MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
    err=sqrt(err/static_cast<type0>(dof_tot))/a_tol;
    return err;
    
}

