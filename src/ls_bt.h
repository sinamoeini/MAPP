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
        Error*& error;
    protected:
        type0 c,rho,min_alpha;
    public:
        LineSearch_backtrack(MAPP*,int,char**);
        virtual void init(Func*);
        ~LineSearch_backtrack();
        int line_min(type0&,type0&,int);
    };
}

#include <stdlib.h>
/*--------------------------------------------
 constructor
 --------------------------------------------*/
template<class Func>
LineSearch_backtrack<Func>::
LineSearch_backtrack(MAPP* mapp,int nargs,char** args):
LineSearch<Func>(mapp),
error(mapp->error)
{
    min_alpha=0.0;
    c=0.4;
    rho=0.5;
    
    if(nargs>2)
    {
        if(nargs%2!=0)
            error->abort("every keyword in ls backtrack should be followed by it's value");
        int iarg=2;
        while(iarg<nargs)
        {
            if(strcmp(args[iarg],"min_alpha")==0)
            {
                iarg++;
                min_alpha=atof(args[iarg]);
                iarg++;
            }
            else if(strcmp(args[iarg],"c")==0)
            {
                iarg++;
                c=atof(args[iarg]);
                iarg++;
            }
            else if(strcmp(args[iarg],"rho")==0)
            {
                iarg++;
                rho=atof(args[iarg]);
                iarg++;
            }
            else
                error->abort("unknown keyword in ls backtrack: %s",args[iarg]);
        }
    }
    
    if(min_alpha<0.0)
        error->abort("min_alpha in ls backtrack should be greater than 0.0");
    if(c<=0.0 || c>=1.0)
        error->abort("c in ls backtrack should be between 0.0 & 1.0");
    if(rho<=0.0 || rho>=1.0)
        error->abort("rho in ls backtrack should be between 0.0 & 1.0");
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