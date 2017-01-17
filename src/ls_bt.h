#ifdef LS_Style
    LSStyle(LineSearch_backtrack,btrack)
#else
#ifndef __MAPP__ls_bt__
#define __MAPP__ls_bt__

#include "ls.h"
namespace MAPP_NS
{
    template<class Func>
    class LineSearch_backtrack : public LineSearch<Func>
    {
    private:
        Func* func;
    protected:
        type0 c,rho,min_alpha;
    public:
        LineSearch_backtrack(int,char**);
        virtual void init(Func*);
        ~LineSearch_backtrack();
        int line_min(type0&,type0&,int);
    };
}
#include <stdlib.h>
#include "cmd.h"
/*--------------------------------------------
 constructor
 --------------------------------------------*/
template<class Func>
LineSearch_backtrack<Func>::
LineSearch_backtrack(int nargs,char** args):
LineSearch<Func>()
{
    min_alpha=0.0;
    c=0.4;
    rho=0.5;
    
    Pattern cmd;

    char* ls_style=NULL;
    cmd.cmd("ls");
    cmd.add_var(ls_style);
    cmd.add_vlog(0)=VLogics("eq","bt");
    
    cmd.cmd("min_alpha");
    cmd.add_var(min_alpha);
    cmd.add_vlog(0)+=VLogics("ge",0.0);
    
    cmd.cmd("c");
    cmd.add_var(c);
    cmd.add_vlog(0)+=VLogics("gt",0.0)*VLogics("lt",1.0);
    
    cmd.cmd("rho");
    cmd.add_var(rho);
    cmd.add_vlog(0)+=VLogics("gt",0.0)*VLogics("lt",1.0);
    
    cmd.scan(args,nargs);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<class Func>
LineSearch_backtrack<Func>::~LineSearch_backtrack()
{
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<class Func>
void LineSearch_backtrack<Func>::
init(Func* func_)
{
    LineSearch<Func>::init(func_);
    func=func_;
}
/*--------------------------------------------
 minimize line
 --------------------------------------------*/
template<class Func>
int LineSearch_backtrack<Func>::line_min(type0& nrgy
,type0& alpha,int init_flag)
{
    type0 max_a,h_norm;
    type0 dfa,init_energy;
    type0 current_energy,ideal_energy;
    type0& prev_val=LineSearch<Func>::prev_val;
    
    init_energy=nrgy;
    
    func->ls_prep(dfa,h_norm,max_a);

    if(max_a==0.0 && dfa==0.0)
        return LS_F_GRAD0;
    
    if(dfa>=0.0)
        return LS_F_DOWNHILL;
    
    if(init_flag==0)
    {
        max_a=MIN(max_a,1.0);
    }
    else if(init_flag==1)
    {
        if(prev_val>0.0 && dfa<0.0)
            max_a=MIN(max_a,-prev_val/dfa);
        prev_val=-dfa;
    }
    else if(init_flag==2)
    {
        if(prev_val>0.0 && dfa<0.0)
            max_a=MIN(max_a,MIN(1.0,2.02*(prev_val+nrgy)/dfa));
        prev_val=-nrgy;
    }
    
    if(max_a<=min_alpha)
        return LS_MIN_ALPHA;
    
    
    alpha=max_a;
    while (1)
    {
        ideal_energy=nrgy+alpha*c*dfa;
        current_energy=func->F(alpha);
        if(current_energy<=ideal_energy)
        {
            if(init_flag==1)
                prev_val*=alpha;
            nrgy=current_energy;
            return LS_S;
        }
        alpha*=rho;
        
        if(alpha<=min_alpha)
        {
            nrgy=init_energy;
            func->F_reset();
            return LS_MIN_ALPHA;
        }
    }
}
#endif
#endif
