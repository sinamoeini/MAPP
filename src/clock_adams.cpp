#include <stdlib.h>
#include "clock_adams.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
#define SMALL 1.0e-18
#define MAX_ITER 50
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
    CREATE1D(err_coef,max_order);
    
    int n=1+max_order/2;
    
    CREATE1D(xi,n);
    CREATE1D(wi,n);
    for(int i=0;i<n;i++)
    {
        CREATE1D(xi[i],i+1);
        CREATE1D(wi[i],i+1);
    }
    for(int i=0;i<n;i++)
        quad(i+1,xi[i],wi[i]);
    
    /*
    t[0]=12.0;
    t[1]=10.0;
    t[2]=9.0;
    t[3]=3.0;
    type0 del_t=1.0;
    int q=4;
    interpolate(del_t,q);
    for(int i=0;i<q;i++)
        printf("alpha_dy[%d] %lf \n",i,alpha_dy[i]);
    printf("\n\n");
    
    for(int i=0;i<q;i++)
        printf("dalpha_dy[%d] %lf \n",i,dalpha_dy[i]);

    
    error->abort("END");
     */
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
        delete [] err_coef;
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
    
    t[0]=0.0;
    memcpy(y,c,dof_lcl*sizeof(type0));
    memcpy(dy[0],c_d,dof_lcl*sizeof(type0));
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
    int ord=1;
    int err_chk;
    int istep;
    
    
    istep=0;
    while (istep<max_step && t[0]<max_t)
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
                max_err=MAX(err,cost);
                ratio=pow(0.5/max_err,1.0/static_cast<type0>(ord+1));
                
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
        ord_dt(del_t,ord,istep,err);
        
        tmp_dy=dy[max_order-1];
        for(int i=max_order-1;i>0;i--)
        {
            t[i]=t[i-1];
            dy[i]=dy[i-1];
        }
        
        dy[0]=tmp_dy;
        t[0]+=del_t_tmp;
        
        memcpy(y,c,dof_lcl*sizeof(type0));
        memcpy(dy[0],c_d,dof_lcl*sizeof(type0));
        
        step_no++;
        istep++;
    }

}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_adams::interpolate(type0& del_t,int& q)
{
    int n;
    type0 tmp0,tmp1,k0;
    
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
                tmp0*=t[0]+xi[n-1][i]*del_t-t[j];
            
            for(int j=0;j<q;j++)
                alpha_dy[j]+=tmp0/(t[0]+xi[n-1][i]*del_t-t[j]);
        }
        for(int i=0;i<q;i++)
        {
            tmp0=1.0;
            tmp1=1.0;
            for(int j=0;j<q;j++)
                if(i!=j)
                {
                    tmp0*=t[i]-t[j];
                    tmp1*=t[0]+del_t-t[j];
                }
            alpha_dy[i]*=del_t/tmp0;
            dalpha_dy[i]=tmp1/tmp0;
        }
        dalpha_y=0.0;
        alpha_y=1.0;
        
        n=1+q/2;
        beta=0.0;
        k0=0.0;
        for(int i=0;i<n;i++)
        {
            tmp0=wi[n-1][i];
            for(int j=0;j<q-1;j++)
                tmp0*=t[0]+xi[n-1][i]*del_t-t[j];
            beta+=tmp0;
            k0+=tmp0*(xi[n-1][i]-1.0);
        }
        k0/=beta;
        
        tmp0=1.0;
        for(int i=0;i<q-1;i++)
            tmp0*=t[0]+del_t-t[i];
        
        beta*=del_t/tmp0;
        
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
        err_prefac=fabs(del_t*k0/(t[0]+del_t-t[q-1]));

    }
}
/*--------------------------------------------
 find the new del_t and order
 --------------------------------------------*/
void Clock_adams::ord_dt(type0& del_t,int& q
,int istep,type0 err)
{
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0 lo_err_lcl,hi_err_lcl;
    type0 tmp0;
    type0 lo_ratio,hi_ratio,ratio,del_t_tmp0,del_t_tmp1;
    type0 hi_err=0.0,lo_err=0.0;
    int lo_ord_avail=0,hi_ord_avail=0,n;
    
    
    if(istep!=0 && q<max_order)
        hi_ord_avail=1;
    if(q>1)
        lo_ord_avail=1;
    
    
    if(hi_ord_avail)
    {
        type0 a0,a1;
        type0 c0,c1,c2,tmp1;
        
        n=1+(q-1)/2;
        c0=c1=c2=0.0;
        for(int i=0;i<n;i++)
        {
            tmp0=wi[n-1][i];
            tmp1=wi[n-1][i];
            for(int j=0;j<q-1;j++)
            {
                tmp0*=t[0]+xi[n-1][i]*del_t-t[j];
                tmp1*=t[1]+xi[n-1][i]*(t[0]-t[1])-t[j+1];
            }
            c0+=tmp0;
            c1+=tmp1;
            c2+=tmp0*(xi[n-1][i]-1.0)*(t[0]+xi[n-1][i]*del_t-t[q-1]);
        }

        a0=a1=c2/static_cast<type0>(q+1);
        a0*=1.0/(c0*(t[0]+del_t-t[q-1]));
        a1*=-del_t/(c1*(t[0]-t[1])*(t[0]-t[q]));
        
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
        MPI_Allreduce(&hi_err_lcl,&hi_err,1,MPI_TYPE0,MPI_SUM,world);
        hi_err=sqrt(hi_err/static_cast<type0>(dof_tot))/a_tol;
    }
    
    if(lo_ord_avail)
    {
        type0 err_e_n;
        type0 c0,c1;
        
        n=1+(q-1)/2;
        c0=c1=0.0;
        for(int i=0;i<n;i++)
        {
            tmp0=wi[n-1][i];
            for(int j=0;j<q-2;j++)
                tmp0*=t[0]+xi[n-1][i]*del_t-t[j];
            c0+=tmp0*(xi[n-1][i]-1.0);
            c1+=tmp0*(t[0]+xi[n-1][i]*del_t-t[q-2]);
        }
        err_e_n=del_t*c0/c1;
        
        for(int i=0;i<q;i++)
        {
            tmp0=1.0;
            for(int j=0;j<q;j++)
                if(i!=j)
                    tmp0*=t[i]-t[j];
            
            err_coef[i]=c0*del_t*del_t/tmp0;
        }

        
        lo_err_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
        {
            if(c[i]>=0.0)
            {
                tmp0=err_e_n*(c[i]-y_0[i]);
                for(int j=0;j<q;j++)
                    tmp0+=err_coef[j]*dy[j][i];
                
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
    
    if(del_t+del_t_tmp1+t[0]>max_t)
        del_t=max_t-t[0]-del_t_tmp1;
}

/*--------------------------------------------
 construct a legendre polynomial of degree n
 --------------------------------------------*/
void Clock_adams::quad(int n,type0* x,type0* w)
{

    int m=n/2+1;
    int err_chk,ord,icurs,jcurs;
    type0 a,u0,inv_u0,f,up,df,tmp0,tol;
    type0 ii,jj;
    type0* p_coef;
    type0* dp_coef;
    
    
    CREATE1D(p_coef,m);
    CREATE1D(dp_coef,m);

    
    
    up=0.0;
    for(int i=0;i<m;i++)
        p_coef[i]=dp_coef[i]=0.0;

    p_coef[0]=1.0;
    
    for(int i=1;i<n+1;i++)
    {
        ii=static_cast<type0>(i);
        if(i%2==0)
        {
            m=i/2+1;
            p_coef[m-1]=(2.0-1.0/ii)*p_coef[m-2];
            for(int j=m-2;j>0;j--)
            {
                jj=static_cast<type0>(j);
                p_coef[j]*=-(2.0*jj+1)/ii;
                p_coef[j]+=(1.0+(2*jj-1.0)/ii)*p_coef[j-1];
            }
            p_coef[0]*=-1.0/ii;
        }
        else
        {
            m=(i+1)/2;
            for(int j=0;j<m-1;j++)
            {
                jj=static_cast<type0>(j);
                p_coef[j]*=1.0+2.0*jj/ii;
                p_coef[j]-=2.0*(jj+1.0)*p_coef[j+1]/ii;
            }
            //2.0-1.0/ii=(1+2*(m-1)/i)
            p_coef[m-1]*=2.0-1.0/ii;
        }
    }
    
    
    
    
    m=n/2+1;
    
    for(int i=1;i<m;i++)
    {
        dp_coef[i-1]=p_coef[i]*static_cast<type0>(i);
    }
    
    
    tol=1.0e-14;
    ord=m;
    icurs=n-1;
    a=p_coef[m-1];
    
    for(int i=0;i<m;i++)
        p_coef[i]/=a;
    for(int i=1;i<m;i++)
        dp_coef[i-1]/=a;
    
    while (ord>1)
    {
        u0=1.0;
        err_chk=1;
        while (err_chk)
        {
            df=f=0.0;
            tmp0=1.0;
            for(int j=0;j<ord-1;j++)
            {
                f+=p_coef[j]*tmp0;
                df+=dp_coef[j]*tmp0;
                tmp0*=u0;
            }
            f+=p_coef[ord-1]*tmp0;
            if(abs(f)<tol)
                err_chk=0;
            if(err_chk)
                u0-=f/df;
        }

        x[icurs]=sqrt(u0);
        
        inv_u0=1.0/u0;
        p_coef[0]*=-inv_u0;
        for(int i=1;i<ord-1;i++)
        {
            p_coef[i]*=-inv_u0;
            p_coef[i]+=inv_u0*p_coef[i-1];
            dp_coef[i-1]=p_coef[i]*static_cast<type0>(i);
        }
        
        ord--;
        icurs--;
    }
    delete [] p_coef;
    delete [] dp_coef;
    
    
    if(n%2==0)
    {
        icurs++;
        jcurs=icurs-1;
        for(int i=icurs;i<n;i++)
        {
            tmp0=a;
            for(int j=icurs;j<n;j++)
            {
                if(i!=j)
                    tmp0*=x[i]*x[i]-x[j]*x[j];
                else
                    tmp0*=x[i]+x[j];
            }
            w[i]=2.0/(tmp0*tmp0*(1.0-x[i]*x[i]));
            w[jcurs]=w[i];
            x[jcurs]=-x[i];
            jcurs--;
        }
    }
    else
    {
        x[icurs]=0.0;
        icurs++;
        tmp0=a;
        for(int i=icurs;i<n;i++)
            tmp0*=-x[i]*x[i];
        w[icurs-1]=2.0/(tmp0*tmp0);
        
        jcurs=icurs-2;
        for(int i=icurs;i<n;i++)
        {
            tmp0=a*x[i];
            for(int j=icurs;j<n;j++)
            {
                if(i!=j)
                    tmp0*=x[i]*x[i]-x[j]*x[j];
                else
                    tmp0*=x[i]+x[j];
            }
            w[i]=2.0/(tmp0*tmp0*(1.0-x[i]*x[i]));
            w[jcurs]=w[i];
            x[jcurs]=-x[i];
            jcurs--;
        }
    }
    
    for(int i=0;i<n;i++)
    {
        x[i]*=0.5;
        x[i]+=0.5;
        w[i]*=0.5;
    }
    

}
