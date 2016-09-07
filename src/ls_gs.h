#ifdef LS_Style
    LSStyle(LineSearch_goldensection,golden)
#else
#ifndef __MAPP__ls_gs__
#define __MAPP__ls_gs__

#include "ls.h"
namespace MAPP_NS
{
    template<class Func>
    class LineSearch_goldensection : public LineSearch<Func>
    {
    private:
        Func* func;
        type0 tol;
        type0 max_iter;
        bool brack;
    protected:
    public:
        LineSearch_goldensection(int,char**);
        virtual void init(Func*);
        ~LineSearch_goldensection();
        int line_min(type0&,type0&,int);
    };
}
#include <stdlib.h>
#include "cmd.h"
/*--------------------------------------------
 constructor
 --------------------------------------------*/
template<class Func>
LineSearch_goldensection<Func>::
LineSearch_goldensection(int nargs,char** args):
LineSearch<Func>()
{
    tol=sqrt(LineSearch<Func>::epsilon);
    max_iter=5;
    brack=true;
    
    
    Pattern cmd;
    char* ls_style=NULL;
    char* bracket=NULL;
    
    cmd.cmd("ls");
    cmd.add_var(ls_style);
    cmd.add_vlog(0)=vlogic("eq","golden");
    
    cmd.cmd("tol");
    cmd.add_var(tol);
    cmd.add_vlog(0)=vlogic("gt",0.0);
    
    cmd.cmd("max_iter");
    cmd.add_var(max_iter);
    cmd.add_vlog(0)=vlogic("gt",0);
    
    cmd.cmd("bracket");
    cmd.add_var(bracket);
    cmd.add_vlog(0)=vlogic("eq","yes")+vlogic("eq","no");
    
    cmd.scan(args,nargs);
    
    if(bracket!=NULL && !strcmp(bracket,"no"))
        brack=false;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<class Func>
LineSearch_goldensection<Func>::
~LineSearch_goldensection()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class Func>
void LineSearch_goldensection<Func>::
init(Func* func_)
{
    LineSearch<Func>::init(func_);
    func=func_;
}
/*--------------------------------------------
 minimize line
 --------------------------------------------*/
template<class Func>
int LineSearch_goldensection<Func>::line_min(type0& nrgy
,type0& alpha,int init_flag)
{
    int calc=0;
    type0 a,b,c;
    type0 fa,fb,fc;
    type0 max_a,h_norm;
    type0 x0,x1,x2,x3,f1,f2;
    type0 dfa;
    type0 gold,cgold;
    type0& prev_val=LineSearch<Func>::prev_val;
    type0 epsilon=LineSearch<Func>::epsilon;
    
    
    gold=0.61803399;
    cgold=0.38196601;
    
    func->ls_prep(dfa,h_norm,max_a);
    
    LineSearch<Func>::h_norm=h_norm;
    
    if(max_a==0.0 && dfa==0.0)
        return LS_F_GRAD0;
    
    if(dfa>=0.0)
        return LS_F_DOWNHILL;
    
    a=0.0;
    fa=nrgy;
    
    
    if(brack)
    {
        int chk_bracket=LineSearch<Func>::bracket(dfa,max_a,a,b,c,fa,fb,fc);
        if(chk_bracket!=B_S)
        {
            func->F_reset();
            return chk_bracket;
        }
        x0=a;
        x3=c;
        
        
        if(c-b>b-a)
        {
            x1=b;
            x2=b+cgold*(c-b);
            
            f1=fb;
            f2=func->F(x2);
            calc=2;
        }
        else
        {
            x1=b+cgold*(a-b);
            x2=b;
            
            f1=func->F(x1);
            f2=fb;
            calc=1;
        }
        
        
        int iter=max_iter;
        while(x3-x0>tol*(x1+x2) && x3-x0>epsilon && iter)
        {
            if(f2<f1)
            {
                x0=x1;
                x1=x2;
                x2=gold*x2+cgold*x3;
                
                f1=f2;
                f2=func->F(x2);
                calc=2;
            }
            else
            {
                x3=x2;
                x2=x1;
                x1=gold*x1+cgold*x0;
                f2=f1;
                f1=func->F(x1);
                calc=1;
            }
            iter--;
        }
    }
    else
    {
        type0 delta,f0,f3;
        f1=f2=0.0;
        x0=0.0;
        f0=fa;
        x3=max_a;
        f3=func->F(x3);
        
        int iter=max_iter;
        
        bool set_left=false;
        bool set_right=false;
        while(x3-x0>tol*(x1+x2) && x3-x0>epsilon && iter)
        {
            delta=cgold*(x3-x0);
            if(!set_left)
            {
                x1=x0+delta;
                f1=func->F(x1);
                calc=1;
                set_left=1;
            }
            
            if(!set_right)
            {
                x2=x3-delta;
                f2=func->F(x2);
                calc=2;
                set_right=1;
            }
            
            
            if((f0<f2 && f0<f3) || (f1<f2 && f1<f3))
            {
                set_left=0;
                x3=x2;
                x2=x1;
                f3=f2;
                f2=f1;
            }
            else
            {
                set_right=0;
                x0=x1;
                x1=x2;
                f0=f1;
                f1=f2;
            }
            
            iter--;
        }
        
        if(!set_left)
        {
            if(f0<f1)
            {
                if(x0==0.0)
                {
                    func->F_reset();
                    alpha=0.0;
                    return LS_MIN_ALPHA;
                }
                
                x1=x0;
                f1=func->F(x1);
                calc=1;
            }
        }
        if(!set_right)
        {
            if(f3<f2)
            {
                
                x2=x3;
                f2=func->F(x2);
                calc=2;
            }
        }
    }
    
    if(f1<f2)
    {
        if(calc==2)
            f1=func->F(x1);
        nrgy=f1;
        alpha=x1;
    }
    else
    {
        if(calc==1)
            f2=func->F(x2);
        nrgy=f2;
        alpha=x2;
    }

    prev_val=-dfa*alpha;
    return LS_S;
}

#endif
#endif
