#include "ls_gs.h"
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearch_goldensection::
LineSearch_goldensection(MAPP* mapp,int narg,char** arg)
:LineSearch(mapp)
{
    tol=sqrt(2.0*epsilon);
    
    if(narg>2)
    {
        if(narg%2!=0)
            error->abort("every keyword in ls goldensection should be followed by it's value");
        int iarg=2;
        while(iarg<narg)
        {
            if(strcmp(arg[iarg],"max_dx")==0)
            {
                iarg++;
                max_dx=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"tol")==0)
            {
                iarg++;
                tol=atof(arg[iarg]);
                iarg++;
            }
            else
                error->abort("unknown keyword in ls goldensection: %s",arg[iarg]);
        }
    }
    
    if(max_dx<=0.0)
        error->abort("max_dx in ls goldensection should be greater than 0.0");
    if(tol<=0.0)
        error->abort("tol in ls goldensection should be greater than 0.0");
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
LineSearch_goldensection::
~LineSearch_goldensection()
{
    
}
/*--------------------------------------------
 minimize line
 --------------------------------------------*/
int LineSearch_goldensection::line_min(type0& nrgy
,type0& alpha,int init_flag)
{
    int calc=0;
    type0 a,b,c;
    type0 fa,fb,fc;
    type0 max_a,h_norm;
    type0 x0,x1,x2,x3,f1,f2;
    type0 r,q;
    type0 f_h;
    r=golden-1.0;
    q=2.0-golden;
    
    init_manip(f_h,h_norm,max_a);
    
    a=0.0;
    fa=nrgy;
    
    int chk_bracket=bracket(f_h,max_a,a,b,c,fa,fb,fc);
    if(chk_bracket!=B_S)
    {
        nrgy=energy(0.0);
        return chk_bracket;
    }
    
    x0=a;
    x3=c;
    
    if(fabs(c-b)>fabs(b-a))
    {
        x1=b;
        x2=b+q*(c-b);
        
        f1=fb;
        f2=energy(x2);
        calc=2;
    }
    else
    {
        x1=b+q*(a-b);
        x2=b;
        
        f1=energy(x1);
        f2=fb;
        calc=1;
    }
    
    
    while(fabs(x3-x0)>tol*(fabs(x1)+fabs(x2)))
    {
        if(f2<f1)
        {
            x0=x1;
            x1=x2;
            x2=r*x2+q*x3;
            
            f1=f2;
            f2=energy(x2);
            calc=2;
        }
        else
        {
            x3=x2;
            x2=x1;
            x1=r*x1+q*x0;
            f2=f1;
            f1=energy(x1);
            calc=1;
        }
    }
    
    if(f1<f2)
    {
        if(calc==2)
            f1=energy(x1);
        nrgy=f1;
        alpha=x1;
    }
    else
    {
        if(calc==1)
            f2=energy(x2);
        nrgy=f2;
        alpha=x2;
    }
    return LS_S;
}
