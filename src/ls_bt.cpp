#include "ls_bt.h"
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearch_backtrack::
LineSearch_backtrack(MAPP* mapp,int narg,char** arg)
:LineSearch(mapp)
{
    min_alpha=0.0;
    c=0.4;
    rho=0.5;
    
    if(narg>2)
    {
        if(narg%2!=0)
            error->abort("every keyword in ls backtrack should be followed by it's value");
        int iarg=2;
        while(iarg<narg)
        {
            if(strcmp(arg[iarg],"max_dx")==0)
            {
                iarg++;
                max_dx=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"min_alpha")==0)
            {
                iarg++;
                min_alpha=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"c")==0)
            {
                iarg++;
                c=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"rho")==0)
            {
                iarg++;
                rho=atof(arg[iarg]);
                iarg++;
            }
            else
                error->abort("unknown keyword in ls backtrack: %s",arg[iarg]);
        }
    }
    
    if(max_dx<=0.0)
        error->abort("max_dx in ls backtrack should be greater than 0.0");
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
LineSearch_backtrack::~LineSearch_backtrack()
{
    
}
/*--------------------------------------------
 minimize line
 --------------------------------------------*/
int LineSearch_backtrack::line_min(type0& nrgy
,type0& alpha,int init_flag)
{
    
    type0 max_a,h_norm;
    type0 f_h;
    
    type0 current_energy,ideal_energy;
    
    
    
    init_manip(f_h,h_norm,max_a);
    if(max_a==0.0 && f_h==0.0)
        return LS_F_GRAD0;
    
    if(f_h<=0.0)
        return LS_F_DOWNHILL;
    
    if(init_flag==0)
    {
        max_a=MIN(max_a,1.0);
    }
    else if(init_flag==1)
    {
        if(prev_val!=0.0 && f_h>0.0)
            max_a=MIN(max_a,prev_val/f_h);
        prev_val=f_h;
    }
    else if(init_flag==2)
    {
        if(prev_val!=0.0 && f_h>0.0)
            max_a=MIN(max_a,MIN(1.0,2.02*(prev_val-nrgy)/f_h));
        prev_val=nrgy;
    }
    
    if(max_a<=min_alpha)
        return LS_MIN_ALPHA;
    
    
    alpha=max_a;
    while (1)
    {
        ideal_energy=nrgy-alpha*c*f_h;
        current_energy=energy(alpha);
        
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
            nrgy=energy(0.0);
            return LS_MIN_ALPHA;
        }
    }
    
    
}
