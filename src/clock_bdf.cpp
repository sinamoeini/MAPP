#include "clock_bdf.h"
#include "neighbor.h"
#include "neighbor.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_bdf::Clock_bdf(MAPP* mapp,int narg
                     ,char** arg):Clock(mapp)
{
    min_gamma=0.0;
    gamma_red=0.8;
    slope=0.4;
    max_iter=50;
    max_order=6;
    m_tol=1.0e-10;
    a_tol=1.0e-6;
    e_tol=0.0;
    min_del_t=1.0e-12;
    max_del_t=1.0e4;
    initial_del_t=-1.0;
    
    if(narg<3)
        error->abort("clock bdf should at least have 1 arguement");
    
    no_steps=atoi(arg[2]);
    if(no_steps<=0)
        error->abort("no of steps for clock bdf should be greater than 0");
    
    if(narg>3)
    {
        if((narg-3)%2!=0)
            error->abort("every keyword in clock bdf should be followed by it's value");
        int iarg=3;
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
            else if(strcmp(arg[iarg],"e_tol")==0)
            {
                iarg++;
                e_tol=atof(arg[iarg]);
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
            else
                error->abort("unknown keyword in clock bdf: %s",arg[iarg]);
        }
    }
    
    if(min_gamma<=0.0)
        error->abort("min_gamma in clock bdf should be greater than 0.0");
    if(gamma_red<=0.0 || gamma_red>=1.0)
        error->abort("gamma_red in clock bdf should be between 0.0 & 1.0");
    if(slope<=0.0 || slope>=1.0)
        error->abort("slope in clock bdf should be between 0.0 & 1.0");
    if(max_iter<=0)
        error->abort("max_iter in clock bdf should be greater than 0");
    if(max_order<=0)
        error->abort("max_order in clock bdf should be greater than 0");
    if(m_tol<=0.0)
        error->abort("m_tol in clock bdf should be greater than 0.0");
    if(a_tol<=0.0)
        error->abort("a_tol in clock bdf should be greater than 0.0");
    if(e_tol<0.0)
        error->abort("e_tol in clock bdf should be equal or greater than 0.0");
    if(min_del_t<=0.0)
        error->abort("min_del_t in clock bdf should be greater than 0.0");
    if(max_del_t<=0.0)
        error->abort("max_del_t in clock bdf should be greater than 0.0");
    if(min_del_t>=max_del_t)
        error->abort("max_del_t in clock bdf should be greater than min_del_t");
    
    
    c_n=atoms->find("c");
    c_d_n=atoms->find("c_d");
    
    int c_dim=atoms->vectors[c_n].dim;
    dof_lcl=atoms->natms*c_dim;
    MPI_Allreduce(&dof_lcl,&dof_tot,1,MPI_INT,MPI_SUM,world);
    
    CREATE1D(y_0,dof_lcl);
    CREATE1D(e_n,dof_lcl);
    CREATE1D(a,dof_lcl);
    CREATE1D(g,dof_lcl);
    CREATE1D(g0,dof_lcl);
    CREATE1D(c0,dof_lcl);
    CREATE1D(h,dof_lcl);
    CREATE1D(dy,dof_lcl);
    
    CREATE1D(t,max_order);
    CREATE1D(y,max_order);
    for(int i=0;i<max_order;i++)
        CREATE1D(y[i],dof_lcl);
    
    CREATE1D(alpha_y,max_order);
    CREATE1D(dalpha_y,max_order);
    CREATE1D(lwr_alpha,max_order);
    
    
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
        delete [] lwr_alpha;
    }
    
    if(dof_lcl)
    {
        delete [] y_0;
        delete [] e_n;
        delete [] a;
        delete [] g;
        delete [] g0;
        delete [] c0;
        delete [] h;
        delete [] dy;
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_bdf::init()
{
    type0* c;
    type0* c_d;
    
    int f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<type0>(0,atoms->vectors[0].dim,"f");
    
    if(cdof_n>-1)
    {
        int dof_n=atoms->find("dof");
        vecs_comm=new VecLst(mapp,5,0,c_n,c_d_n,cdof_n,dof_n);
        
    }
    else
    {
        vecs_comm=new VecLst(mapp,3,0,c_n,c_d_n);
    }

    
    vecs_comm->add_update(0);
    atoms->reset_comm(vecs_comm);
    forcefield->init();
    atoms->ph_setup(1,vecs_comm);
    
    neighbor->init();
    neighbor->create_list(0,1);
    atoms->store_0();
    forcefield->create_2nd_neigh_lst();
    

    
    forcefield->force_calc(1,nrgy_strss);
    
    thermo->update(fe_idx,nrgy_strss[0]);
    thermo->update(stress_idx,6,&nrgy_strss[1]);
    thermo->update(time_idx,0.0);
    
    if(write!=NULL)
        write->init();
    thermo->init();
    
    
    t[0]=0.0;
    thermo->start_force_time();
    forcefield->c_d_calc();
    thermo->stop_force_time();
    atoms->vectors[c_n].ret(c);
    atoms->vectors[c_d_n].ret(c_d);
    
    rectify(c_d);
    if(initial_del_t<0.0)
    {
        type0 sum,sum_lcl;
        sum_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
            sum_lcl+=c_d[i]*c_d[i];
        sum=0.0;
        MPI_Allreduce(&sum_lcl,&sum,1,MPI_TYPE0,MPI_SUM,world);
        sum=sqrt(sum/static_cast<type0>(dof_tot));
        initial_del_t=MAX(a_tol/sum,min_del_t);
        initial_del_t=MIN(initial_del_t,max_del_t);
    }
    
    memcpy(y[0],c,dof_lcl*sizeof(type0));
    memcpy(dy,c_d,dof_lcl*sizeof(type0));
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_bdf::fin()
{
    
    if(write!=NULL)
        write->fin();
    
    forcefield->fin();
    neighbor->fin();
    
    delete vecs_comm;
    thermo->fin();
    
    
    atoms->x2s(atoms->natms);
}
/*--------------------------------------------
 calculate the coefficients
 --------------------------------------------*/
int Clock_bdf::interpolate(type0 del_t,int q)
{
    type0 tmp0,tmp1,k0,k1,k2,k3,k4,curr_t;
    curr_t=t[0]+del_t;
    
    k0=k3=1.0;
    k1=k2=0.0;
    k4=2.0/(curr_t-t[0]);
    for(int i=1;i<q;i++)
    {
        k3*=1.0/(t[0]-t[i]);
        tmp0=(curr_t-t[0])/(t[i]-t[0]);
        tmp1=tmp0*tmp0;
        for(int j=1;j<q;j++)
        {
            if(i!=j)
                tmp1*=(curr_t-t[j])/(t[i]-t[j]);
        }
        alpha_y[i]=tmp1;
        
        k0*=(curr_t-t[i])/(t[0]-t[i]);
        k1+=1.0/(t[0]-t[i]);
        k2+=1.0/(curr_t-t[i]);
    }
    
    alpha_y[0]=k0*(1.0-(curr_t-t[0])*k1);
    alpha_dy_0=k0*(curr_t-t[0]);
    
    
    //fore derivatives
    dalpha_dy_0=k0*(1.0+(curr_t-t[0])*k2);
    dalpha_y[0]=k0*(k2-k1-(curr_t-t[0])*k1*k2);
    
    for(int i=1;i<q;i++)
        dalpha_y[i]=alpha_y[i]*(k4+k2-1.0/(curr_t-t[i]));
    
    
    
    tmp0=0.0;
    for(int i=0;i<q;i++)
        tmp0+=1.0/static_cast<type0>(i+1);
    beta=del_t/tmp0;
    
    
    int ret_val=1;
    int tot_ret_val;
    int idof=0;
    while(idof<dof_lcl && ret_val==1)
    {
        
        a[idof]=beta*dalpha_dy_0*dy[idof];
        y_0[idof]=alpha_dy_0*dy[idof];
        for(int j=0;j<q;j++)
        {
            a[idof]+=beta*dalpha_y[j]*y[j][idof];
            y_0[idof]+=alpha_y[j]*y[j][idof];
        }
        
        a[idof]-=y_0[idof];
        if(y_0[idof]<0.0 || y_0[idof]>1.0)
        {
            ret_val=-1;
        }
        
        idof++;
    }
    
    MPI_Allreduce(&ret_val,&tot_ret_val,1,MPI_INT,MPI_MIN,world);
    err_prefac=fabs(-1.0/((curr_t-t[0])*(k2+0.5*k4)*(1.0+k0)));
    return tot_ret_val;
}

/*--------------------------------------------
 find the new del_t and order
 --------------------------------------------*/
void Clock_bdf::ord_dt(type0& del_t,int& q
,int istep)
{
    
    type0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    type0* c;
    atoms->vectors[c_n].ret(c);
    type0 curr_t=t[0]+del_t;
    type0 lo_err_lcl,hi_err_lcl;
    type0 tmp0;
    type0 Q_n=1.0;
    type0 hi_prefac=0.0,lo_prefac=0.0;
    type0 c5,c6,c7,c8,c9,c10;
    
    type0 lo_ratio,hi_ratio,ratio,del_t_tmp;
    type0 hi_err=0.0,lo_err=0.0;
    int lo_ord_avail=0,hi_ord_avail=0;
    
    
    if(istep!=0 && q<max_order)
        hi_ord_avail=1;
    if(q>1)
        lo_ord_avail=1;
    
    
    if(hi_ord_avail)
    {
        c5=c6=1.0;
        c7=(curr_t-t[0])/(t[0]-t[1]);
        c8=1.0/(curr_t-t[0])+1.0/(curr_t-t[q]);
        
        for(int i=1;i<q;i++)
        {
            c5*=(curr_t-t[i])/(curr_t-t[0]);
            c6*=(t[0]-t[i+1])/(t[0]-t[1]);
            c7*=(curr_t-t[i])/(t[0]-t[i+1]);
            c8+=1.0/(curr_t-t[i]);
        }
        c5++;
        c6++;
        Q_n=(curr_t-t[0])*c7*c5/(c6*(t[0]-t[1]));
        hi_prefac=-(curr_t-t[q])/((curr_t-t[0])
        *(curr_t-t[0])*c8*c5*static_cast<type0>(q+2));
        
        
        hi_err_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
        {
            tmp0=((c[i]-y_0[i])-Q_n*e_n[i]);
            hi_err_lcl+=tmp0*tmp0;
        }
        hi_err=0.0;
        MPI_Allreduce(&hi_err_lcl,&hi_err,1
        ,MPI_TYPE0,MPI_SUM,world);
        hi_err=fabs(hi_prefac)*
        sqrt(hi_err/static_cast<type0>(dof_tot))/a_tol;
        
    }
    
    if(lo_ord_avail)
    {
        c9=0.0;
        c10=1.0;
        for(int i=0;i<q-1;i++)
        {
            c9+=1.0/(curr_t-t[i]);
            c10=curr_t-t[i];
        }
        lo_prefac=c10/c9;
        
        for(int i=0;i<q-1;i++)
        {
            tmp0=1.0/(curr_t-t[i]);
            lwr_alpha[i]=tmp0*tmp0*lo_prefac;
            
            for(int j=0;j<q-1;j++)
                if(i!=j)
                    lwr_alpha[i]*=1.0/(t[i]-t[j]);
        }
        lwr_alpha_dy=1.0/c9;
        lwr_alpha_y=-1.0;
        
        lo_err_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
        {
        
            tmp0=lwr_alpha_y*c[i];
            tmp0+=lwr_alpha_dy*c_d[i];
            for(int j=0;j<q-1;j++)
                tmp0+=lwr_alpha[j]*y[j][i];
            
            lo_err_lcl+=tmp0*tmp0;
        }
        
        lo_err=0.0;
        MPI_Allreduce(&lo_err_lcl,&lo_err,1,MPI_TYPE0,MPI_SUM,world);
        lo_err=sqrt(lo_err/static_cast<type0>(dof_tot))/a_tol;
        
    }
    
    for(int i=0;i<dof_lcl;i++)
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
    
    del_t_tmp=ratio*del_t;
    if(del_t_tmp>max_del_t)
        del_t=max_del_t;
    else if (del_t_tmp<min_del_t)
        del_t=min_del_t;
    else
        del_t=del_t_tmp;
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_bdf::run()
{
    type0* c;
    atoms->vectors[c_n].ret(c);
    type0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    type0* tmp_y;
    
    type0 del_t=initial_del_t,del_t_tmp,err1;
    type0 ratio=1.0,cost;
    int ord=1;
    int chk;
    int istep;
    
    eq_ratio=1.0;
    istep=0;
    while (eq_ratio>=1.0 && istep <no_steps)
    {
        chk=interpolate(del_t,ord);
        while(chk==-1)
        {

            del_t*=0.8;
            del_t=MAX(del_t,min_del_t);
            chk=interpolate(del_t,ord);
        }
        cost=solve(del_t,ord);
        while(err>=1.0 || cost>=1.0)
        {
            err1=MAX(err,cost);
            ratio=pow(0.5/err1,1.0/static_cast<type0>(ord+1));

            del_t_tmp=del_t*ratio;
            if(del_t_tmp<min_del_t)
                del_t=min_del_t;
            else
                del_t=del_t_tmp;
            
            
            chk=interpolate(del_t,ord);
            
            while(chk==-1)
            {
                del_t*=0.8;
                chk=interpolate(del_t,ord);
            }
            
            cost=solve(del_t,ord);
        }

        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step()|| istep==no_steps-1)
        {
            thermo->start_force_time();
            forcefield->force_calc(1,nrgy_strss);
            thermo->stop_force_time();
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            thermo->update(time_idx,t[0]+del_t);
        }
        
        del_t_tmp=del_t;

        ord_dt(del_t,ord,istep);
        
        
        tmp_y=y[max_order-1];
        for(int i=max_order-1;i>0;i--)
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
 init
 --------------------------------------------*/
type0 Clock_bdf::solve(type0 del_t,int q)
{
    
    type0* c;
    atoms->vectors[c_n].ret(c);
    type0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    type0 gamma,max_gamma=1.0;
    type0 inner,tmp0,tmp1;
    type0 err_lcl;
    type0 ratio;
    type0 g0_g0,g_g,g_g0,g_h;
    type0 curr_cost,ideal_cost,cost;
    int chk;
    
    
    memcpy(c,y_0,dof_lcl*sizeof(type0));
    thermo->start_comm_time();
    atoms->update(c_n);
    thermo->stop_comm_time();
    
    
    /*
    type0 tot_ratio;
    ratio=1.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp=dy[i]*del_t+y[0][i];
        if(tmp>1.0)
            ratio=MIN((1.0-y[0][i])/(tmp-y[0][i]),ratio);
        else if(tmp<0.0)
            ratio=MIN((y[0][i]-0.0)/(y[0][i]-tmp),ratio);
    }
    MPI_Allreduce(&ratio,&tot_ratio,1,MPI_TYPE0,MPI_MIN,world);
    for(int i=0;i<dof_lcl;i++)
        c[i]=y[0][i]+dy[i]*del_t*tot_ratio;
    thermo->start_comm_time();
    atoms->update(c_n);
    thermo->stop_comm_time();
    */
    
    
    thermo->start_force_time();
    curr_cost=forcefield->g_calc(0,beta,a,g);
    rectify(g);
    thermo->stop_force_time();
    /*
    type0 sum_lcl;
    type0 sum_tot;
    sum_lcl=0.0;
    sum_tot=0.0;
    for(int i=0;i<dof_lcl;i++)
        sum_lcl+=g[i];
    MPI_Allreduce(&sum_lcl,&sum_tot,1,MPI_TYPE0,MPI_SUM,world);
    for(int i=0;i<dof_lcl;i++)
        g[i]-=correc_fac*sum_tot;
    rectify(g);
     */
    memcpy(h,g,dof_lcl*sizeof(type0));
    
    inner=0.0;
    for(int i=0;i<dof_lcl;i++)
        inner+=g[i]*g[i];
    g0_g0=0.0;
    MPI_Allreduce(&inner,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
    g_h=g0_g0;
    
    int iter=0;
    
    while(curr_cost>m_tol*static_cast<type0>(dof_tot)
          && iter<max_iter && max_gamma>min_gamma)
    {
        memcpy(g0,g,dof_lcl*sizeof(type0));
        memcpy(c0,c,dof_lcl*sizeof(type0));
        
        gamma=0.99;
        
        for(int i=0;i<dof_lcl;i++)
        {
            tmp0=c0[i]+h[i];
            
            if(tmp0>1.0)
            {
                gamma=MIN((1.0-c0[i])/(tmp0-c0[i]),gamma);
            }
            else if(tmp0<0.0)
            {
                gamma=MIN((c0[i]-0.0)/(c0[i]-tmp0),gamma);
            }
        }
        MPI_Allreduce(&gamma,&max_gamma,1,MPI_TYPE0,MPI_MIN,world);
        
        max_gamma*=0.9;
        chk=1;
        
        cost=curr_cost;
            
        while(chk && max_gamma>min_gamma)
        {
            for(int i=0;i<dof_lcl;i++)
                c[i]=c0[i]+max_gamma*h[i];
            
            thermo->start_comm_time();
            atoms->update(c_n);
            thermo->stop_comm_time();
            
            thermo->start_force_time();
            curr_cost=forcefield->g_calc(0,beta,a,g);
            rectify(g);
            thermo->stop_force_time();
            /*
            sum_lcl=0.0;
            sum_tot=0.0;
            for(int i=0;i<dof_lcl;i++)
                sum_lcl+=g[i];
            MPI_Allreduce(&sum_lcl,&sum_tot,1,MPI_TYPE0,MPI_SUM,world);
            for(int i=0;i<dof_lcl;i++)
                g[i]-=correc_fac*sum_tot;
            rectify(g);
            */
            ideal_cost=cost-slope*max_gamma*g_h;
            if(curr_cost<ideal_cost)
                chk=0;
            max_gamma*=gamma_red;
            if(max_gamma<=min_gamma)
            {
                memcpy(c,c0,dof_lcl*sizeof(type0));
                
                thermo->start_comm_time();
                atoms->update(c_n);
                thermo->stop_comm_time();
                
                thermo->start_force_time();
                curr_cost=forcefield->g_calc(1,beta,a,g);
                thermo->stop_force_time();
            }
        }
        
        if(chk==0)
        {
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
                inner+=g[i]*g0[i];
            g_g0=0.0;
            MPI_Allreduce(&inner,&g_g0,1,MPI_TYPE0,MPI_SUM,world);
            
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
                inner+=g[i]*g[i];
            g_g=0.0;
            MPI_Allreduce(&inner,&g_g,1,MPI_TYPE0,MPI_SUM,world);
            
            ratio=(g_g-g_g0)/g0_g0;
            
            g0_g0=g_g;
            
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
            {
                h[i]*=ratio;
                h[i]+=g[i];
                inner+=h[i]*g[i];
            }
            g_h=0.0;
            MPI_Allreduce(&inner,&g_h,1,MPI_TYPE0,MPI_SUM,world);
            
            if(g_h<0.0)
            {
                memcpy(h,g,dof_lcl*sizeof(type0));
                g_h=g_g;
            }
        }
        
        iter++;
    }
    
    rectify(c_d);
    

    err_lcl=0.0;
    tmp1=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp0=(y_0[i]-c[i]);
        tmp1+=c_d[i]*c_d[i];
        err_lcl+=tmp0*tmp0;
    }
    
    err=0.0;
    MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
    err=sqrt(err/static_cast<type0>(dof_tot))/a_tol;
    err*=err_prefac;
    
    MPI_Allreduce(&tmp1,&eq_ratio,1,MPI_TYPE0,MPI_SUM,world);
    eq_ratio=sqrt(eq_ratio/static_cast<type0>(dof_tot))/e_tol;
    
    
    for(int i=0;i<dof_lcl;i++)
        if(c[i]<0.0 || c[i]>1.0)
            error->abort("c exceeded the domain");
    
    return curr_cost/(m_tol*static_cast<type0>(dof_tot));
    
}
