#include <stdlib.h>
#include <limits>
#include "dmd_inf.h"
#include "ff.h"
#include "write.h"
#include "error.h"
#include "memory.h"
#include "thermo_dynamics.h"
#include "atom_types.h"
#define INITIAL_STEP_MODE 1
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMD_inf::DMD_inf(MAPP* mapp,int nargs
,char** args):DMDImplicit(mapp)
{
    error->abort("this dmd style is not ready yet");
    epsilon=std::numeric_limits<type0>::epsilon();
    epsilon_sqrt=sqrt(epsilon);
    golden=0.5+0.5*sqrt(5.0);
    
    mod=MIN_MOD_EN_DC_PROJ;

    if(nargs>2)
    {
        int iarg=2;
        while(iarg<nargs)
        {

            if(strcmp(args[iarg],"max_step")==0)
            {
                iarg++;
                max_step=atoi(args[iarg]);
                iarg++;
            }
            else if(strcmp(args[iarg],"a_tol")==0)
            {
                iarg++;
                a_tol=atof(args[iarg]);
                iarg++;
            }
            else if(strcmp(args[iarg],"dc_proj")==0)
            {
                mod=MIN_MOD_EN_DC_PROJ;
                iarg++;
            }
            else if(strcmp(args[iarg],"en_proj")==0)
            {
                mod=MIN_MOD_FORM_EN;
                iarg++;
            }
            else if(strcmp(args[iarg],"dc_norm_proj")==0)
            {
                mod=MIN_MOD_DC_NORM;
                iarg++;
            }
            else
                error->abort("unknown keyword in dmd int: %s",args[iarg]);
        }
    }
    
    if(max_step<=0)
        error->abort("max_step in dmd bdf should be greater than 0");
    if(a_tol<=0.0)
        error->abort("a_tol in dmd bdf should be greater than 0.0");


}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
DMD_inf::~DMD_inf()
{
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD_inf::allocate()
{
    DMDImplicit::allocate();
    CREATE1D(c1,ncs);
    CREATE1D(g_orig,ncs);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD_inf::deallocate()
{

    if(ncs)
    {
        delete [] c1;
        delete [] g_orig;
    }
    DMDImplicit::deallocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_inf::init()
{
    
    DMDImplicit::init();
    allocate();

    ls_mode=LS_GS;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_inf::fin()
{
    DMDImplicit::fin();
    deallocate();
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void DMD_inf::run()
{
    type0* c=mapp->c->begin();
    
    type0 gamma;
    type0 inner0,inner1;
    type0 ratio;
    type0 g0_g0,g_g,g_g0,g_h;
    type0 curr_cost,cost,iter=0;
    int line_search_succ,err;
    prev_val=-1.0;
    
    /* find the steepest descent direction and cost */
    
    if(mod==MIN_MOD_EN_DC_PROJ)
    {
        curr_cost=forcefield_dmd->dc_en_proj_timer(1,h,g_h);
        rectify(h);
        inner0=0.0;
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
            {
                inner0+=h[i]*h[i];
            }
        }
        g0_g0=0.0;
        MPI_Allreduce(&inner0,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
    }
    else
    {
        if(mod==MIN_MOD_FORM_EN)
            curr_cost=forcefield_dmd->en_grad_timer(1,g_orig,g);
        else
            curr_cost=forcefield_dmd->dc_norm_grad_timer(1,g_orig,g);
        
        rectify(g);
        rectify(g_orig);
        
        /* set the first trajectory */
        memcpy(h,g,ncs*sizeof(type0));
        
        /* calculate g.h g_0.g_0 */
        inner0=0.0;
        for(int i=0;i<ncs;i++)
            if(c[i]>=0.0)
                inner0+=g_orig[i]*g[i];
        g0_g0=0.0;
        MPI_Allreduce(&inner0,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
        g_h=g0_g0;
    }
    
    iter=0;
    line_search_succ=1;
    err=0;
    while(err==0)
    {

        memcpy(c0,c,ncs*sizeof(type0));
        memcpy(g0,g,ncs*sizeof(type0));
        
        cost=curr_cost;
        
        if(write!=NULL)
            write->write();
        
        thermo->thermo_print();
        
        ls_prep();
        
        
        
        gamma=0.0;
        if(ls_mode==LS_GS)
            line_search_succ=line_search_gs(gamma,curr_cost,g_h);
        else if(ls_mode==LS_BT)
            line_search_succ=line_search_bt(gamma,curr_cost,g_h);
        else if(ls_mode==LS_BRENT)
            line_search_succ=line_search_brent(gamma,curr_cost,g_h);
       
        if(line_search_succ==0)
        {
            //test1(0.0,cost,g_h);
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            err=1;
            continue;
        }
        
        prev_val=gamma;
        
        
        if(mod==MIN_MOD_EN_DC_PROJ)
        {
            curr_cost=forcefield_dmd->dc_en_proj_timer(1,h,g_h);
            rectify(h);
            inner0=0.0;
            for(int i=0;i<ncs;i++)
            {
                if(c[i]>=0.0)
                {
                    inner0+=h[i]*h[i];
                }
            }
            g0_g0=0.0;
            MPI_Allreduce(&inner0,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
        }
        else
        {
            if(mod==MIN_MOD_FORM_EN)
                curr_cost=forcefield_dmd->en_grad_timer(1,g_orig,g);
            else
                curr_cost=forcefield_dmd->dc_norm_grad_timer(1,g_orig,g);
            rectify(g);
            rectify(g_orig);
            
            /* set the first trajectory */
            memcpy(h,g,ncs*sizeof(type0));
            
            /* calculate g.h g_0.g_0 */
            inner0=inner1=0.0;
            for(int i=0;i<ncs;i++)
                if(c[i]>=0.0)
                {
                    inner0+=g_orig[i]*g0[i];
                    inner1+=g_orig[i]*g[i];
                }
            
            MPI_Allreduce(&inner0,&g_g0,1,MPI_TYPE0,MPI_SUM,world);
            MPI_Allreduce(&inner1,&g_g,1,MPI_TYPE0,MPI_SUM,world);
            
            
            ratio=(g_g-g_g0)/g0_g0;
            
            /* update g_0.g_0 */
            g0_g0=g_g;
            
            /* calculate g_h */
            inner0=0.0;
            for(int i=0;i<ncs;i++)
            {
                if(c[i]>=0.0)
                {
                    h[i]*=ratio;
                    h[i]+=g[i];
                    inner0+=g_orig[i]*h[i];
                }
            }
            MPI_Allreduce(&inner0,&g_h,1,MPI_TYPE0,MPI_SUM,world);
            
            /* if g_h is negative start from the begining */
            if(g_h<0.0)
            {
                memcpy(h,g,ncs*sizeof(type0));
                g_h=g_g;
            }
        }
        
        if(iter+1==max_step)
            err=1;
        if(cost-curr_cost<a_tol)
            err=1;
        if(g0_g0==0.0)
            err=1;
        
        if(thermo->test_prev_step() || err)
        {
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
        }
        
        iter++;
        step_no++;

    }

}
/*--------------------------------------------
 find the the cost function given gamma
 --------------------------------------------*/
inline type0 DMD_inf::cost_func(type0 gamma)
{
    type0 tmp;
    type0* c=mapp->c->begin();
    
    if(gamma<max_a)
    {
        for(int i=0;i<ncs;i++)
            if(c0[i]>=0.0)
                c[i]=c0[i]+gamma*h[i];
    }
    else
    {
        for(int i=0;i<ncs;i++)
            if(c0[i]>=0.0)
                c[i]=c1[i];
    }
    
    atoms->update(mapp->c);
    
    if(mod==MIN_MOD_EN_DC_PROJ)
        return forcefield_dmd->dc_en_proj_timer(0,g,tmp);
    else if(mod==MIN_MOD_FORM_EN)
        return forcefield_dmd->en_grad_timer(0,g_orig,g);
    else
        return forcefield_dmd->dc_norm_grad_timer(0,g_orig,g);

}
/*--------------------------------------------
 find the the cost function given gamma
 --------------------------------------------*/
inline type0 DMD_inf::dcost_func(type0 gamma,type0& df)
{
    type0* c=mapp->c->begin();
    
    if(gamma<max_a)
    {
        for(int i=0;i<ncs;i++)
            if(c0[i]>=0.0)
                c[i]=c0[i]+gamma*h[i];
    }
    else
    {
        for(int i=0;i<ncs;i++)
            if(c0[i]>=0.0)
                c[i]=c1[i];
    }
    
    atoms->update(mapp->c);
    
    type0 cost=forcefield_dmd->en_grad_timer(1,g_orig,g);
    type0 inner0=0.0;
    for(int i=0;i<ncs;i++)
    {
        if(c[i]>=0.0)
        {
            inner0+=g_orig[i]*h[i];
        }
    }
    MPI_Allreduce(&inner0,&df,1,MPI_TYPE0,MPI_SUM,world);
    
    return cost;
}
/*--------------------------------------------
 find the the cost function given gamma
 --------------------------------------------*/
inline void DMD_inf::ls_prep()
{
    type0 h_norm_lcl;
    type0 max_a_lcl;
    max_a_lcl=numeric_limits<type0>::infinity();
    h_norm_lcl=0.0;
    for(int i=0;i<ncs;i++)
    {
        if(c0[i]>=0.0)
        {
            //printf("c0[%d] h[%d]: %32.20lf  %32.20lf\n",i,i,c0[i],h[i]);
            h_norm_lcl+=h[i]*h[i];
            if(h[i]>0.0)
            {
                c1[i]=(1.0-c0[i])/h[i];
                max_a_lcl=MIN(c1[i],max_a_lcl);
            }
            else if(h[i]<0.0)
            {
                c1[i]=(0.0-c0[i])/h[i];
                max_a_lcl=MIN(c1[i],max_a_lcl);
            }
        }
    }
    MPI_Allreduce(&max_a_lcl,&max_a,1,MPI_TYPE0,MPI_MIN,world);
    MPI_Allreduce(&h_norm_lcl,&h_norm,1,MPI_TYPE0,MPI_SUM,world);
    h_norm=sqrt(h_norm);
    
    if(h_norm==0.0)
        error->abort("h_norm is zero");
    
    
    for(int i=0;i<ncs;i++)
    {
        if(c0[i]>=0.0)
        {
            h_norm_lcl+=h[i]*h[i];
            if(h[i]>0.0)
            {
                if(c1[i]>max_a)
                    c1[i]=c0[i]+max_a*h[i];
                else
                    c1[i]=1.0;
            }
            else if(h[i]<0.0)
            {
                if(c1[i]>max_a)
                    c1[i]=c0[i]+max_a*h[i];
                else
                    c1[i]=0.0;
            }
            else
                c1[i]=c0[i];
        }
    }
}
/*--------------------------------------------
 given the direction h do the lin seach
 --------------------------------------------*/
int DMD_inf::line_search_gs(type0& a0,type0& fa0,type0 dfa0)
{
    type0 a1,fa1,a2,fa2;
    type0 r,q,ulim,u,fu;
    type0 tol,x0,x1,x2,x3,f1,f2;
    int calc;
    int max_iter_ls;
    
    max_iter_ls=20;
    calc=0;
    a2=0.0;
    
    
    /* beginning of finding the maximum gamma */
    if(max_a<epsilon)
    {
        fa2=cost_func(max_a);
        if(fa2<=fa0)
        {
            a0=max_a;
            fa0=fa2;
            return 1;
        }
        type0* c=mapp->c->begin();
        memcpy(c,c0,ncs*sizeof(type0));
        atoms->update(mapp->c);
        return 0;
    }
    /* end of finding the maximum gamma */
    
    /* beginning of bracketing minimum */
    if(dfa0!=numeric_limits<type0>::infinity())
    {
        if(prev_val<0.0)
            a1=MIN(MAX(epsilon_sqrt/fabs(dfa0),epsilon_sqrt/h_norm),max_a);
        else
        {
            a1=MIN(epsilon_sqrt/fabs(dfa0),0.1*prev_val);
        }
    }
    else
    {
        a1=MIN(epsilon_sqrt*sqrt(epsilon_sqrt)/fabs(dfa0),max_a);
    }
    

    
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
        fa2=cost_func(max_a);
        if(fa2<=fa0)
        {
            a0=max_a;
            fa0=fa2;
            return 1;
        }
        
        if(dfa0*a1<epsilon_sqrt)
        {
            a1=MIN(epsilon_sqrt*sqrt(epsilon_sqrt)/fabs(dfa0),max_a);
            fa1=cost_func(a1);
        }
        
        if(fa1>=fa0)
        {
            type0* c=mapp->c->begin();
            memcpy(c,c0,ncs*sizeof(type0));
            atoms->update(mapp->c);
            return 0;
        }
    }
    
    fa2=fa1;
    int end=0;
    while (fa2<=fa1 && end==0)
    {
        a2=a1+golden*(a1-a0);
        if(a2>max_a)
        {
            a2=max_a;
            fa2=cost_func(a2);
            if(fa2<=fa0)
            {
                a0=max_a;
                fa0=fa2;
                return 1;
            }
            end=1;
            continue;
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
int DMD_inf::line_search_brent(type0& a0,type0& fa0,type0 dfa0)
{
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
    if(max_a<epsilon)
    {
        fa2=cost_func(max_a);
        if(fa2<=fa0)
        {
            a0=max_a;
            fa0=fa2;
            return 1;
        }
        type0* c=mapp->c->begin();
        memcpy(c,c0,ncs*sizeof(type0));
        atoms->update(mapp->c);
        return 0;

    }
    /* end of finding the maximum gamma */
    
    /* beginning of bracketing minimum */
    if(dfa0!=numeric_limits<type0>::infinity())
    {
        if(prev_val<0.0)
            a1=MIN(epsilon_sqrt/fabs(dfa0),1.0e-10*max_a);
        else
            a1=MIN(1.0e-12/fabs(dfa0),0.25*prev_val);
    }
    else
    {
        a1=1.0e-12*max_a;
    }
    
    
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
        fa2=cost_func(max_a);
        if(fa2<=fa0)
        {
            a0=max_a;
            fa0=fa2;
            return 1;
        }
        type0* c=mapp->c->begin();
        memcpy(c,c0,ncs*sizeof(type0));
        atoms->update(mapp->c);
        return 0;
    }
    
    fa2=fa1;
    int end=0;
    while (fa2<=fa1 && end==0)
    {
        a2=a1+golden*(a1-a0);
        if(a2>max_a)
        {
            a2=max_a;
            fa2=cost_func(a2);
            if(fa2<=fa0)
            {
                a0=max_a;
                fa0=fa2;
                return 1;
            }
            end=1;
            continue;
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
int DMD_inf::line_search_bt(type0& a0,type0& fa0,type0 dfa0)
{
    type0 u,a=0.0,fa2;
    
    
    /* beginning of finding the maximum gamma */
    if(max_a<epsilon)
    {
        fa2=cost_func(max_a);
        if(fa2<=fa0)
        {
            a0=max_a;
            fa0=fa2;
            return 1;
        }
        type0* c=mapp->c->begin();
        memcpy(c,c0,ncs*sizeof(type0));
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
    memcpy(c,c0,ncs*sizeof(type0));
    atoms->update(mapp->c);
    return 0;
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
int DMD_inf::test1(type0 a0,type0 fa0,type0 dfa0)
{

    dfa0*=-1.0;
    type0* c=mapp->c->begin();
   
    /*
    type0* mu;
    atoms->vectors[atoms->find("mu")]->ret(mu);
    
    for(int i=0;i<ncs;i++)
    {
        
        if(c0[i]>=0.0)
        {
            printf("%d c %e h %e mu %24.22lf \n",i,c[i],h[i],mu[i]);
        }
    }
     */
    
    


    int no=100;
    type0 u,frac,tmp0;
    
     frac=2.0*max_a*1.0e-5/static_cast<type0>(no);
     u=0.0;
     printf("fa0 %24.22lf dfa0 %e max_a %e dfa0*max_a %e hnorm %e\n",fa0,dfa0,max_a,dfa0*max_a,h_norm);
     
     for(int i=0;i<no;i++)
     {
         tmp0=cost_func(u);
         printf("%e %e %e \n",u,cost_func(u)-fa0,u*dfa0);
         u+=frac;
     }

     memcpy(c,c0,ncs*sizeof(type0));
     atoms->update(mapp->c);
    

    return 0;

}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
int DMD_inf::test0()
{


    
    return 0;
}

