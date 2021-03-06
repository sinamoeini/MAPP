#ifdef LS_Style
    LSStyle(LineSearch_brent,brent)
#else
#ifndef __MAPP__ls_brent__
#define __MAPP__ls_brent__

#include "ls.h"
namespace MAPP_NS
{
    template<class Func>
    class LineSearch_brent : public LineSearch<Func>
    {
    private:
        Func* func;
        type0 epsilon;
        type0 tol,zeps;
        int max_iter;
        int ls_brack_timer_mode;
        int ls_timer_mode;
        bool brack;
    protected:
    public:
        LineSearch_brent(int,char**);
        ~LineSearch_brent();
        void init(Func*);
        int line_min(type0&,type0&,int);
    };
    
}
#include <stdlib.h>
#include "cmd.h"
/*--------------------------------------------
 constructor
 --------------------------------------------*/
template<class Func>
LineSearch_brent<Func>::
LineSearch_brent(int nargs,char** args)
:LineSearch<Func>()
{
    epsilon=LineSearch<Func>::epsilon;
    
    max_iter=5;
    tol=sqrt(2.0*epsilon);
    zeps=epsilon*1.0e-3;
    brack=true;

    Pattern cmd;
    
    char* ls_style=NULL;
    char* bracket=NULL;
    
    cmd.cmd("ls");
    cmd.add_var(ls_style);
    cmd.add_vlog(0)=VLogics("eq","brent");

    cmd.cmd("tol");
    cmd.add_var(tol);
    cmd.add_vlog(0)=VLogics("gt",0.0);
    
    cmd.cmd("zeps");
    cmd.add_var(zeps);
    cmd.add_vlog(0)=VLogics("ge",0.0);
    
    cmd.cmd("max_iter");
    cmd.add_var(max_iter);
    cmd.add_vlog(0)=VLogics("gt",0);
    
    cmd.cmd("bracket");
    cmd.add_var(bracket);
    cmd.add_vlog(0)=VLogics("eq","yes")+VLogics("eq","no");
    
    cmd.scan(args,nargs);
    
    if(bracket&&!strcmp(bracket,"no"))
        brack=false;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<class Func>
LineSearch_brent<Func>::
~LineSearch_brent()
{
    
}
/*--------------------------------------------

 --------------------------------------------*/
template<class Func>
void LineSearch_brent<Func>::
init(Func* func_)
{
    LineSearch<Func>::init(func_);
    func=func_;
}
/*--------------------------------------------
 minimize line
 --------------------------------------------*/
template<class Func>
int LineSearch_brent<Func>::line_min(type0& nrgy
,type0& alpha,int init_flag)
{
    type0& prev_val=LineSearch<Func>::prev_val;
    type0 a,b;
    type0 fa;
    type0 max_a,h_norm;
    type0 cgold,gold;
    type0 dfa;

    type0 d=0.0,etemp,fu,fv,fw,fx;
    type0 r,q,p,tol1,tol2,u,v,w,x,xm;
    type0 e=0.0;
    u=fu=0.0;
    
    gold=0.61803399;
    cgold=0.38196601;
    
    func->ls_prep(dfa,h_norm,max_a);
    if(max_a==0.0 && dfa==0.0)
        return LS_F_GRAD0;
    
    if(dfa>=0.0)
        return LS_F_DOWNHILL;
    
    a=0.0;
    fa=nrgy;
    
    
    type0 fb;
    
    if(brack)
    {
        int chk_bracket=LineSearch<Func>::bracket(dfa,max_a,a,x,b,fa,fx,fb);
        
        if(chk_bracket!=B_S)
        {
            func->F_reset();
            return chk_bracket;
        }
    }
    else
    {
        a=0.0;
        b=max_a;
        x=cgold*(b-a);
        fx=func->F(x);
    }
    
    w=v=x;
    fw=fv=fx;

    for(int iter=0;iter<max_iter;iter++)
    {
        xm=0.5*(a+b);
        tol1=tol*fabs(x)+zeps;
        tol2=2.0*tol1;
        
        if(fabs(x-xm)<=(tol2-0.5*(b-a)))
        {
            if(u!=x)
                func->F(x);
            nrgy=fx;
            alpha=x;
            prev_val=alpha;
            return LS_S;
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
               || p<=q*(a-x)
               || p>=q*(b-x))
            {
                if(x>=xm)
                    e=(a-x);
                else
                    e=b-x;
                d=cgold*e;
            }
            else
            {
                d=p/q;
                u=x+d;
                if(u-a<tol2 || b-u<tol2)
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
                e=(a-x);
            else
                e=b-x;
            
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
        
        fu=func->F(u);
        
        if(fu<=fx)
        {
            if(u>=x)
                a=x;
            else
                b=x;
            
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
                a=u;
            else
                b=u;
            
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

    //make sure that the result is less than initial value
    if(fa<fx)
    {
        func->F_reset();
        alpha=0.0;
        return LS_MIN_ALPHA;
    }
    
    if(u!=x)
        func->F(x);
    
    nrgy=fx;
    alpha=x;
    
    prev_val=alpha;
    return LS_S;
}


#endif
#endif
