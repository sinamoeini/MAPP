#include "clock.h"
#include "thermo_dynamics.h"
#include "atom_types.h"
#include "write.h"
#include "error.h"
#include "memory.h"
#include "timer.h"
#include <limits>
#define INITIAL_STEP_MODE 1
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock::Clock(MAPP* mapp):InitPtrs(mapp)
{
    ns_alloc=0;
    if(forcefield==NULL)
        error->abort("ff should be "
        "initiated before clock");
    
    if(mapp->mode!=DMD_mode)
        error->abort("clock works only "
        "for md mode");
    
    char** args;
    int nargs=mapp->parse_line(
    "Time FE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    
    fe_idx=1;
    stress_idx=2;
    time_idx=0;
    thermo=new ThermoDynamics(mapp,nargs,args);
    for(int i=0;i<nargs;i++)
        delete [] args[i];
    if(nargs)
        delete [] args;
    
    
    int dim=atoms->dimension;
    if(dim)
    {
        CREATE1D(nrgy_strss,dim*(dim+1)/2+1);
        ns_alloc=1;
    }
    
    
    c_dim=mapp->c->dim;
    forcefield_dmd=dynamic_cast<ForceFieldDMD*>(forcefield);
    neighbor_dmd=dynamic_cast<Neighbor_dmd*>(neighbor);
    
    a_tol=sqrt(numeric_limits<type0>::epsilon());
    min_del_t=std::numeric_limits<type0>::epsilon();
    initial_del_t=-1.0;
    max_t=1.0e7;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock::~Clock()
{
    if(ns_alloc)
        delete [] nrgy_strss;
    delete thermo;
}
/*--------------------------------------------
 rectify
 --------------------------------------------*/
void Clock::rectify(type0* f)
{
    if(mapp->cdof==NULL)
        return;
    bool* cdof=mapp->cdof->begin();
    
    int tot=(atoms->natms)*(atom_types->no_types);
    for(int i=0;i<tot;i++) if(cdof[i]==1) f[i]=0.0;
}
/*--------------------------------------------
 default init
 --------------------------------------------*/
void Clock::init()
{
    old_skin=atoms->get_skin();
    atoms->set_skin(0.0);
    if(mapp->c_d==NULL)
        mapp->c_d=new Vec<type0>(atoms,c_dim,"c_d");
    if(mapp->f==NULL)
        mapp->f=new Vec<type0>(atoms,mapp->x->dim);
    
    vecs_comm=new VecLst(atoms);
    vecs_comm->add_updt(mapp->c);
    vecs_comm->add_updt(mapp->ctype);
    atoms->init(vecs_comm,false);
    
    
    neighbor_dmd->create_2nd_list();
    forcefield_dmd->force_calc_timer(1,nrgy_strss);
    forcefield_dmd->dc_timer();
    
    thermo->update(fe_idx,nrgy_strss[0]);
    thermo->update(stress_idx,6,&nrgy_strss[1]);
    thermo->update(time_idx,0.0);
    thermo->init();
    
    if(write!=NULL)
        write->init();
    
    type0* c=mapp->c->begin();
    dof_lcl=atoms->natms*c_dim;
    int tmp_dof=dof_lcl;
    for(int idof=0;idof<dof_lcl;idof++)
        if(c[idof]<0.0)
            tmp_dof--;
    MPI_Allreduce(&tmp_dof,&dof_tot,1,MPI_INT,MPI_SUM,world);
    
    tot_t=0.0;
}
/*--------------------------------------------
 default fin
 --------------------------------------------*/
void Clock::fin()
{
    if(write!=NULL)
        write->fin();
    thermo->fin();
    atoms->fin();
    print_stats();
    timer->print_stats();
    neighbor->print_stats();
    delete vecs_comm;
    atoms->set_skin(old_skin);
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ClockImplicit::ClockImplicit(MAPP* mapp):Clock(mapp)
{
    epsilon=std::numeric_limits<type0>::epsilon();
    golden=0.5+0.5*sqrt(5.0);
    //defaults for solver
    min_gamma=0.0;
    gamma_red=0.5;
    slope=0.4;
    max_iter=50;
    m_tol=sqrt(numeric_limits<type0>::epsilon());
    ls_mode=LS_BT;
    pre_cond=1;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ClockImplicit::~ClockImplicit()
{
}
/*--------------------------------------------
 solve the implicit equation
 --------------------------------------------*/
void ClockImplicit::solve_n_err(type0& cost,type0& err)
{
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    
    type0 gamma;
    type0 inner0,inner1,tmp0;
    type0 err_lcl;
    type0 ratio;
    type0 g0_g0,g_g,g_g0,g_h;
    type0 curr_cost;
    int iter,line_search_succ;
    prev_val=-1.0;
    
    if(pre_cond==1)
    {
        /* beginning of pre-conditioning */
        memcpy(c,y_0,dof_lcl*sizeof(type0));
        atoms->update(mapp->c);
        /* end of pre-conditioning */
    }
    else if(pre_cond==2)
    {
        /* beginning of pre-conditioning */
        memcpy(c,y_1,dof_lcl*sizeof(type0));
        atoms->update(mapp->c);
        /* end of pre-conditioning */
    }

    
    /* find the steepest descent direction and cost */
    curr_cost=forcefield_dmd->imp_cost_grad_timer(1,beta,a,g);
    rectify(g);
    

    
    
    /* set the first trajectory */
    memcpy(h,g,dof_lcl*sizeof(type0));
    
    /* calculate g.h g_0.g_0 */
    inner0=0.0;
    for(int i=0;i<dof_lcl;i++)
        if(c[i]>=0.0)
            inner0+=g[i]*g[i];
    g0_g0=0.0;
    MPI_Allreduce(&inner0,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
    g_h=g0_g0;

    iter=0;
    line_search_succ=1;
    
    while(curr_cost>m_tol*static_cast<type0>(dof_tot)
    && iter<max_iter && line_search_succ==1)
    {
        memcpy(g0,g,dof_lcl*sizeof(type0));
        memcpy(c0,c,dof_lcl*sizeof(type0));
        
        
        /* do the line search here */
        
        gamma=0.0;
        //test(gamma,curr_cost,-g_h);
        if(ls_mode==LS_GS)
            line_search_succ=line_search_gs(gamma,curr_cost,g_h);
        else if(ls_mode==LS_BT)
            line_search_succ=line_search_bt(gamma,curr_cost,g_h);
        else if(ls_mode==LS_BRENT)
            line_search_succ=line_search_brent(gamma,curr_cost,g_h);

        /* after a successful line search */
        if(curr_cost<=m_tol*static_cast<type0>(dof_tot)
        || line_search_succ==0)
            continue;
       
        curr_cost=forcefield_dmd->imp_cost_grad_timer(1,beta,a,g);

        rectify(g);

        /* calculate g_1.g_0 & g_1.g_1 */
        inner0=inner1=0.0;
        for(int i=0;i<dof_lcl;i++)
            if(c[i]>=0.0)
            {
                inner0+=g[i]*g0[i];
                inner1+=g[i]*g[i];
            }

        MPI_Allreduce(&inner0,&g_g0,1,MPI_TYPE0,MPI_SUM,world);
        MPI_Allreduce(&inner1,&g_g,1,MPI_TYPE0,MPI_SUM,world);

        
        /* calculate (g_1.g_1-g_1.g_0)/g_0.g_0 */
        ratio=(g_g-g_g0)/g0_g0;
        
        /* update g_0.g_0 */
        g0_g0=g_g;
        
        /* calculate g_h */
        inner0=0.0;
        for(int i=0;i<dof_lcl;i++)
        {
            if(c[i]>=0.0)
            {
                h[i]*=ratio;
                h[i]+=g[i];
                inner0+=h[i]*g[i];
            }
        }
        MPI_Allreduce(&inner0,&g_h,1,MPI_TYPE0,MPI_SUM,world);
        
        /* if g_h is negative start from the begining */
        if(g_h<0.0)
        {
            memcpy(h,g,dof_lcl*sizeof(type0));
            g_h=g_g;
        }
        
        
        iter++;
    }
    
    rectify(c_d);
    
    err_lcl=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        if(c[i]>=0.0)
        {
            tmp0=(y_0[i]-c[i]);
            err_lcl+=tmp0*tmp0;
        }
    }
    
    err=0.0;
    MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
    err=sqrt(err/static_cast<type0>(dof_tot))/a_tol;
    err*=err_prefac;

    cost=curr_cost/(m_tol*static_cast<type0>(dof_tot));
    if(iter)
    {
        if(cost<1.0)
            solve_acc++;
        else
            solve_rej++;
    }
    if(err>=1.0)
    {
        intg_rej++;
    }

}
/*--------------------------------------------
 find the the cost function given gamma
 --------------------------------------------*/
inline type0 ClockImplicit::cost_func(type0 gamma)
{
    type0* c=mapp->c->begin();
    
    for(int i=0;i<dof_lcl;i++)
        if(c0[i]>=0.0)
            c[i]=c0[i]+gamma*h[i];
    
    atoms->update(mapp->c);
    
    return forcefield_dmd->imp_cost_grad_timer(0,beta,a,g);
}
/*--------------------------------------------
 given the direction h do the lin seach
 --------------------------------------------*/
int ClockImplicit::line_search_gs(type0& a0,type0& fa0,type0 dfa0)
{
    type0 h_norm_lcl,h_norm;
    type0 max_a_lcl,max_a;
    type0 a1,fa1,a2,fa2;
    type0 r,q,ulim,u,fu;
    type0 tol,x0,x1,x2,x3,f1,f2;
    int calc;
    int max_iter_ls;

    max_iter_ls=20;
    calc=0;
    a2=0.0;

    
    /* beginning of finding the maximum gamma */
    max_a_lcl=numeric_limits<type0>::infinity();
    h_norm_lcl=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        if(c0[i]>=0.0)
        {
            h_norm_lcl+=h[i]*h[i];
            if(h[i]>0.0)
                max_a_lcl=MIN((1.0-c0[i])/h[i],max_a_lcl);
            else if(h[i]<0.0)
                max_a_lcl=MIN((0.0-c0[i])/h[i],max_a_lcl);
        }
    }
    
    MPI_Allreduce(&max_a_lcl,&max_a,1,MPI_TYPE0,MPI_MIN,world);
    MPI_Allreduce(&h_norm_lcl,&h_norm,1,MPI_TYPE0,MPI_SUM,world);
    max_a*=0.999;

    if(max_a<epsilon)
    {
        type0* c=mapp->c->begin();
        memcpy(c,c0,dof_lcl*sizeof(type0));
        atoms->update(mapp->c);
        return 0;
    }
    /* end of finding the maximum gamma */
    
    /* beginning of bracketing minimum */
    a0=0.0;
    a1=epsilon/sqrt(h_norm);
    if(a1>max_a)
        a1=max_a*epsilon;
    
    fa1=fa0;
    u=r=a1;
    while(fa1==fa0 && u<max_a)
    {
        a1=u;
        fa1=cost_func(a1);
        u+=r;
    }
    
    if(fa1>fa0)
    {
        type0* c=mapp->c->begin();
        memcpy(c,c0,dof_lcl*sizeof(type0));
        atoms->update(mapp->c);
        return 0;
    }
    
    fa2=fa1;
    while (fa2<=fa1)
    {
        a2=a1+golden*(a1-a0);
        if(a2>max_a)
        {
            type0* c=mapp->c->begin();
            memcpy(c,c0,dof_lcl*sizeof(type0));
            atoms->update(mapp->c);
            return 0;
        }
        fa2=cost_func(a2);
        
        if(fa2>fa1)
            continue;
        ulim=MIN(a1+(golden+2.0)*(a1-a0),max_a);
        
        r=(a1-a0)*(fa1-fa2);
        q=(a1-a2)*(fa1-fa0);
        
        u=0.5*a1+(a2*q-a0*r)/(2.0*(q-r));
        
        if(a1<u && u<a2)
        {
            fu=cost_func(u);
            if(fu<fa2)
            {
                a0=a1;
                a1=u;
                fa0=fa1;
                fa1=fu;
            }
            else if(fu>fa1)
            {
                a2=u;
                fa2=fu;
            }
            else
            {
                a0=a1;
                a1=a2;
                
                fa0=fa1;
                fa1=fa2;
            }
            
        }
        else if (a2<u)
        {
            u=MIN(u,ulim);
            fu=cost_func(u);
            
            a0=a1;
            a1=a2;
            a2=u;
            
            fa0=fa1;
            fa1=fa2;
            fa2=fu;
            
        }
        else
        {
            a0=a1;
            a1=a2;
            
            fa0=fa1;
            fa1=fa2;
        }
    }
    /* end of bracketing minimum */
    
    
    /* beginning of golden section */
    r=golden-1.0;
    q=2.0-golden;
    
    x0=a0;
    x3=a2;
    
    if(a2-a1>a1-a0)
    {
        x1=a1;
        x2=a1+q*(a2-a1);
        
        f1=fa1;
        f2=cost_func(x2);
        calc=2;
        
    }
    else
    {
        x1=a1+q*(a0-a1);
        x2=a1;
        
        f1=cost_func(x1);
        f2=fa1;
        calc=1;
    }
    
    tol=sqrt(2.0*epsilon);

    while(fabs(x3-x0)>tol*(fabs(x1)+fabs(x2)) && x3-x0>epsilon && max_iter_ls)
    {
        if(f2<f1)
        {
            x0=x1;
            x1=x2;
            x2=r*x2+q*x3;
            
            f1=f2;
            f2=cost_func(x2);
            calc=2;
        }
        else
        {
            x3=x2;
            x2=x1;
            x1=r*x1+q*x0;
            f2=f1;
            f1=cost_func(x1);
            calc=1;
        }
        
        max_iter_ls--;
    }
    
    /* end of golden section */
    
    if(f1<f2)
    {
        if(calc==2)
            f1=cost_func(x1);
        fa0=f1;
        a0=x1;
    }
    else
    {
        if(calc==1)
            f2=cost_func(x2);
        fa0=f2;
        a0=x2;
    }

    return 1;
}
/*--------------------------------------------
 given the direction h do the lin seach
 --------------------------------------------*/
int ClockImplicit::line_search_brent(type0& a0,type0& fa0,type0 dfa0)
{
    type0 h_norm_lcl,h_norm;
    type0 max_a_lcl,max_a;
    type0 a1,fa1,a2,fa2;
    type0 r,q,ulim,u,fu;
    type0 tol;
    type0 ax,bx,cx,p;
    type0 x,w,v,fx,fw,fv,xm,tol1,tol2,zeps,e,etemp,d,cgold;
    int max_iter_ls;
    int calc;
    
    zeps=tol=epsilon;
    cgold=0.38196601;
    max_iter_ls=20;
    e=0.0;
    d=0.0;
    calc=0;
    a2=0.0;
    
    /* beginning of finding the maximum gamma */
    max_a_lcl=numeric_limits<type0>::infinity();
    h_norm_lcl=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        if(c0[i]>=0.0)
        {
            h_norm_lcl+=h[i]*h[i];
            if(h[i]>0.0)
                max_a_lcl=MIN((1.0-c0[i])/h[i],max_a_lcl);
            else if(h[i]<0.0)
                max_a_lcl=MIN((0.0-c0[i])/h[i],max_a_lcl);
        }
    }
    
    MPI_Allreduce(&max_a_lcl,&max_a,1,MPI_TYPE0,MPI_MIN,world);
    MPI_Allreduce(&h_norm_lcl,&h_norm,1,MPI_TYPE0,MPI_SUM,world);
    max_a*=0.999;
    
    if(max_a<epsilon)
    {
        type0* c=mapp->c->begin();
        memcpy(c,c0,dof_lcl*sizeof(type0));
        atoms->update(mapp->c);
        return 0;
    }
    /* end of finding the maximum gamma */
    
    /* beginning of bracketing minimum */
    a0=0.0;
    a1=epsilon/sqrt(h_norm);
    if(a1>max_a)
        a1=max_a*epsilon;
    
    fa1=fa0;
    u=r=a1;
    while(fa1==fa0 && u<max_a)
    {
        a1=u;
        fa1=cost_func(a1);
        u+=r;
    }
    
    if(fa1>fa0)
    {
        type0* c=mapp->c->begin();
        memcpy(c,c0,dof_lcl*sizeof(type0));
        atoms->update(mapp->c);
        return 0;
    }
    
    fa2=fa1;
    while (fa2<=fa1)
    {
        a2=a1+golden*(a1-a0);
        if(a2>max_a)
        {
            type0* c=mapp->c->begin();
            memcpy(c,c0,dof_lcl*sizeof(type0));
            atoms->update(mapp->c);
            return 0;
        }
        fa2=cost_func(a2);
        
        if(fa2>fa1)
            continue;
        ulim=MIN(a1+(golden+2.0)*(a1-a0),max_a);
        
        r=(a1-a0)*(fa1-fa2);
        q=(a1-a2)*(fa1-fa0);
        
        u=0.5*a1+(a2*q-a0*r)/(2.0*(q-r));
        
        if(a1<u && u<a2)
        {
            fu=cost_func(u);
            if(fu<fa2)
            {
                a0=a1;
                a1=u;
                fa0=fa1;
                fa1=fu;
            }
            else if(fu>fa1)
            {
                a2=u;
                fa2=fu;
            }
            else
            {
                a0=a1;
                a1=a2;
                
                fa0=fa1;
                fa1=fa2;
            }
            
        }
        else if (a2<u)
        {
            u=MIN(u,ulim);
            fu=cost_func(u);
            
            a0=a1;
            a1=a2;
            a2=u;
            
            fa0=fa1;
            fa1=fa2;
            fa2=fu;
            
        }
        else
        {
            a0=a1;
            a1=a2;
            
            fa0=fa1;
            fa1=fa2;
        }
    }
    /* end of bracketing minimum */
    
    /* beginning of brent */
    ax=a0;
    bx=a1;
    cx=a2;
    
    a0=ax;
    a1=cx;
    
    x=w=v=bx;
    fx=fw=fv=fa1;
    
    for(int iter=0;iter<max_iter_ls;iter++)
    {
        xm=0.5*(a0+a1);
        tol1=tol*fabs(x)+zeps;
        tol2=2.0*tol1;
        
        if(fabs(x-xm)<=(tol2-0.5*(a1-a0)))
        {
            if(u!=x)
                cost_func(x);
            fa0=fx;
            a0=x;
            return 1;
        }
        
        if(fabs(e)>tol1)
        {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if(q>0.0) p=-p;
            q=fabs(q);
            etemp=e;
            e=d;
            if(fabs(p)>=fabs(0.5*q*etemp)
               || p<=q*(a0-x)
               || p>=q*(a1-x))
            {
                if(x>=xm)
                    e=(a0-x);
                else
                    e=a1-x;
                d=cgold*e;
            }
            else
            {
                d=p/q;
                u=x+d;
                if (u-a0<tol2 || a1-u<tol2)
                {
                    if(xm-x>=0.0)
                        d=fabs(tol1);
                    else
                        d=-fabs(tol1);
                }
            }
        }
        else
        {
            if(x>=xm)
                e=(a0-x);
            else
                e=a1-x;
            
            d=cgold*e;
        }
        
        if(fabs(d)>=tol1)
        {
            u=x+d;
        }
        else
        {
            if(d>=0.0)
            {
                u=x+fabs(tol1);
            }
            else
            {
                u=x-fabs(tol1);
            }
        }
        
        fu=cost_func(u);
        
        if(fu<=fx)
        {
            if(u>=x)
                a0=x;
            else
                a1=x;
            
            v=w;
            w=x;
            x=u;
            
            fv=fw;
            fw=fx;
            fx=fu;
        }
        else
        {
            if(u<x)
                a0=u;
            else
                a1=u;
            
            if(fu<=fw
               || w==x)
            {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            }
            else if(fu<=fv
                    || v==x
                    || v==w)
            {
                v=u;
                fv=fu;
            }
        }
    }
    
    
    if(u!=x)
        cost_func(x);
    fa0=fx;
    a0=x;
    
    /* beginning of brent */
    return 1;
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
int ClockImplicit::line_search_bt(type0& a0,type0& fa0,type0 dfa0)
{
    type0 u,a,max_a;

    
    /* beginning of finding the maximum gamma */
    a=numeric_limits<type0>::infinity();
    
    for(int i=0;i<dof_lcl;i++)
    {
        if(c0[i]>=0.0)
        {
            if(h[i]>0.0)
            {
                a=MIN((1.0-c0[i])/h[i],a);
            }
            else if(h[i]<0.0)
            {
                a=MIN((0.0-c0[i])/h[i],a);
            }
            
            
        }
    }
    
    MPI_Allreduce(&a,&max_a,1,MPI_TYPE0,MPI_MIN,world);
    max_a*=0.999;
    
    if(max_a<numeric_limits<type0>::epsilon())
    {
        type0* c=mapp->c->begin();
        memcpy(c,c0,dof_lcl*sizeof(type0));
        atoms->update(mapp->c);
        return 0;
    }
    /* end of finding the maximum gamma */
    
    
    
    type0 curr_cost,ideal_cost;
    
    curr_cost=ideal_cost=0.0;
    if(prev_val!=-1.0)
    {
        if(INITIAL_STEP_MODE==1)
        {
            u=MIN(max_a,-prev_val/dfa0);
        }
        else if (INITIAL_STEP_MODE==2)
        {
            u=MIN(max_a,-2.0*(fa0-prev_val)/dfa0);
            prev_val=fa0;
        }
        else
            u=max_a;
    }
    else
    {
        u=max_a;
    }
    
    while (curr_cost>=ideal_cost && u>min_gamma)
    {
        a=u;
        curr_cost=cost_func(a);
        ideal_cost=fa0-slope*a*dfa0;
        u=gamma_red*a;
    }
    
    if(curr_cost<ideal_cost)
    {
        a0=a;
        fa0=curr_cost;
        if(INITIAL_STEP_MODE==1)
            prev_val=-a0*dfa0;
        return 1;
    }

    
    type0* c=mapp->c->begin();
    memcpy(c,c0,dof_lcl*sizeof(type0));
    atoms->update(mapp->c);
    return 0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ClockImplicit::allocate()
{
    CREATE1D(y_0,dof_lcl);
    CREATE1D(y_1,dof_lcl);
    CREATE1D(a,dof_lcl);
    CREATE1D(g,dof_lcl);
    CREATE1D(g0,dof_lcl);
    CREATE1D(c0,dof_lcl);
    CREATE1D(h,dof_lcl);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ClockImplicit::deallocate()
{
    if(dof_lcl)
    {
        delete [] y_0;
        delete [] y_1;
        delete [] a;
        delete [] g;
        delete [] g0;
        delete [] c0;
        delete [] h;
    }
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void ClockImplicit::print_stats()
{
    if(atoms->my_p==0)
    {
        fprintf(output,"clock stats:\n");
        fprintf(output,"efficiancy fac: %e\n",tot_t/timer->tot_time);
        fprintf(output,"max timestep: %e\n",max_succ_dt);
        fprintf(output,"max    order: %d\n",max_succ_q);
        fprintf(output,"rejected integration   attempts: %d\n",intg_rej);
        fprintf(output,"rejected interpolation attempts: %d\n",intp_rej);
        fprintf(output,"total implicit equations: %d = accepted(%d) + rejected(%d)\n",solve_rej+solve_acc,solve_acc,solve_rej);
    }
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void ClockImplicit::init()
{
    Clock::init();
    max_succ_q=1;
    max_succ_dt=0.0;
    solve_rej=0;
    solve_acc=0;
    intg_rej=0;
    intp_rej=0;
    
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void ClockImplicit::fin()
{
    Clock::fin();
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
int ClockImplicit::test(type0 a0,type0 fa0,type0 dfa0)
{
    type0 u,a,max_a;
    int no=100;
    type0 frac;
    
    a=numeric_limits<type0>::infinity();
    
    for(int i=0;i<dof_lcl;i++)
    {
        if(c0[i]>=0.0)
        {
            if(h[i]>0.0)
            {
                a=MIN((1.0-c0[i])/h[i],a);
            }
            else if(h[i]<0.0)
            {
                a=MIN((0.0-c0[i])/h[i],a);
            }
            
            
        }
    }
    
    MPI_Allreduce(&a,&max_a,1,MPI_TYPE0,MPI_MIN,world);
    max_a=MIN(max_a,-1.1*fa0/dfa0);
    max_a*=0.999;
    
    if(max_a>1.0e-9)
    {

        frac=max_a/static_cast<type0>(no);
        u=0.0;
        printf("fa0 %24.22lf dfa0 %24.22lf max_a %e\n",fa0,dfa0,max_a);
        
        for(int i=0;i<no;i++)
        {
            printf("%e %e %e\n",u,cost_func(u)-fa0,u*dfa0);
            u+=frac;
        }
        
        type0* c=mapp->c->begin();
        memcpy(c,c0,dof_lcl*sizeof(type0));
        atoms->update(mapp->c);
    }
    return 0;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ClockExplicit::ClockExplicit(MAPP* mapp):
Clock(mapp)
{
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ClockExplicit::~ClockExplicit()
{
    
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void ClockExplicit::init()
{
    Clock::init();
    max_succ_dt=0.0;
    intg_rej=0;
    intp_rej=0;
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void ClockExplicit::fin()
{
    Clock::fin();
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void ClockExplicit::print_stats()
{
    if(atoms->my_p==0)
    {
        fprintf(output,"clock stats:\n");
        fprintf(output,"efficiancy fac: %e\n",tot_t/timer->tot_time);
        fprintf(output,"max timestep: %e\n",max_succ_dt);
        fprintf(output,"rejected integration   attempts: %d\n",intg_rej);
        fprintf(output,"rejected interpolation attempts: %d\n",intp_rej);
    }
}
