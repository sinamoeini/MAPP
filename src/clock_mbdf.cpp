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
    
    for(int i=0;i<max_order+1;i++)
        t[i]=0.0;
    
    t[0]=0.0;
    t[1]=-initial_del_t;
    memcpy(y[0],c,dof_lcl*sizeof(type0));
    memcpy(y[1],c,dof_lcl*sizeof(type0));
    memcpy(dy,c_d,dof_lcl*sizeof(type0));
    tot_t=0.0;
    
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
void Clock_mbdf::interpolate(type0& del_t,int& q)
{
    type0 tmp0;

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

        
        for(int i=0;i<q+1;i++)
        {
            c0+=1.0/(1.0-t[i]/del_t);
            //c0+=1.0/(del_t-t[i]);
            
            tmp0=1.0;
            for(int j=0;j<q+1;j++)
                if(i!=j)
                    tmp0*=(del_t-t[j])/(t[i]-t[j]);
            
            alpha_y[i]=tmp0;
        }
        
        for(int i=0;i<q+1;i++)
        {
            tmp0=c0-1.0/(1.0-t[i]/del_t);
            dalpha_y[i]=alpha_y[i]*tmp0/del_t;
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
            const_stps=0;
            if(q>1)
            {
                q--;
            }
            else
            {
                memcpy(y[1],y[0],dof_lcl*sizeof(type0));
                memcpy(y_0,y[0],dof_lcl*sizeof(type0));
                for(int i=0;i<dof_lcl;i++)
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
    }
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
void Clock_mbdf::ord_dt(type0& del_t,int& q)
{
    type0 terkm2=0.0,terkm1=0.0,terk=0.0,terkp1=0.0;
    type0 est,ratio;
    
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
    
    
    ratio=pow(0.5/est,1.0/static_cast<type0>(q+1));
    
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
    
    del_t=MAX(min_del_t,MIN(ratio*del_t,MIN(max_t-tot_t,max_del_t)));
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
    int q=1;
    int err_chk;
    int initial_phase;
    int istep;
    
    initial_phase=1;
    const_stps=0;
    istep=0;
    err=0.0;
    
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
                    if(initial_phase) initial_phase=0;
                    const_stps=0;
                    max_err=MAX(err,cost);
                    ratio=pow(0.5/max_err,1.0/static_cast<type0>(q+1));
                    
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
        
        if(initial_phase)
        {
            del_t=MIN(MIN(2.0*del_t,max_del_t),max_t-tot_t);
            if(q<max_order)
                q++;
        }
        else
            ord_dt(del_t,q);
        
        
        tmp_y=y[max_order];
        for(int i=max_order;i>0;i--)
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
/*--------------------------------------------
 max step size
 --------------------------------------------*/
type0 Clock_mbdf::err_est(int q,type0 del_t)
{
    type0 tmp0,err_lcl,err,c0;
    type0* c;
    atoms->vectors[c_n]->ret(c);
    
    
    c0=1.0;
    for(int i=0;i<q;i++)
        c0*=static_cast<type0>(i+1)/(1.0-t[i]/del_t);
    
    for(int i=0;i<q;i++)
    {
        tmp0=1.0;
        for(int j=0;j<q;j++)
            if(j!=i)
                tmp0*=(del_t-t[j])/(t[i]-t[j]);
        
        alph_err[i]=-tmp0;
    }
    
    err_lcl=0.0;
    for(int i=0;i<dof_lcl;i++)
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
    err=sqrt(err/static_cast<type0>(dof_tot))/a_tol;
    err*=fabs(c0);
    return err;
}

