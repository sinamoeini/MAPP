#include "ls_brent.h"
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearch_brent::
LineSearch_brent(MAPP* mapp,int narg,char** arg)
:LineSearch(mapp)
{
    max_iter=50;
    tol=sqrt(2.0*epsilon);
    zeps=epsilon*1.0e-3;
    zeps=tol;
    
    if(narg>2)
    {
        if(narg%2!=0)
            error->abort("every keyword in ls brent should be followed by it's value");
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
            else if(strcmp(arg[iarg],"zeps")==0)
            {
                iarg++;
                zeps=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_iter")==0)
            {
                iarg++;
                max_iter=atoi(arg[iarg]);
                iarg++;
            }
            else
                error->abort("unknown keyword in ls brent: %s",arg[iarg]);
        }
    }
    
    if(max_dx<=0.0)
        error->abort("max_dx in ls brent should be greater than 0.0");
    if(tol<=0.0)
        error->abort("tol in ls brent should be greater than 0.0");
    if(zeps<0.0)
        error->abort("zeps in ls brent should be greater than 0.0");
    if(max_iter<=0)
        error->abort("max_iter in ls brent should be greater than 0");
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
LineSearch_brent::
~LineSearch_brent()
{
    
}
/*--------------------------------------------
 minimize line
 --------------------------------------------*/
int LineSearch_brent::line_min(type0& nrgy
,type0& alpha,int init_flag)
{
    type0 a,b,c;
    type0 fa,fb,fc;
    type0 max_a,h_norm;
    type0 cgold,gold;
    type0 dfa;
    type0 ax,bx,cx;
    type0 d=0.0,etemp,fu,fv,fw,fx;
    type0 r,q,p,tol1,tol2,u,v,w,x,xm;
    type0 e=0.0;
    u=fu=0.0;
    
    gold=0.61803399;
    cgold=0.38196601;
    
    init_manip(dfa,h_norm,max_a);
    if(max_a==0.0 && dfa==0.0)
        return LS_F_GRAD0;
    
    if(dfa>=0.0)
        return LS_F_DOWNHILL;
    
    a=0.0;
    fa=nrgy;
    
    int chk_bracket=bracket(dfa,max_a,a,b,c,fa,fb,fc);
    if(chk_bracket!=B_S)
    {
        reset();
        return chk_bracket;
    }


    
    ax=a;
    bx=b;
    cx=c;
    
    a=ax;
    b=cx;
    
    x=w=v=bx;
    fx=fw=fv=fb;
    
    for(int iter=0;iter<max_iter;iter++)
    {
        xm=0.5*(a+b);
        tol1=tol*fabs(x)+zeps;
        tol2=2.0*tol1;
        
        if(fabs(x-xm)<=(tol2-0.5*(b-a)))
        {
            energy(x);
            nrgy=fx;
            alpha=x;
            return  LS_S;
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
                if (u-a < tol2 || b-u < tol2)
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
        
        fu=energy(u);
        
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
    
    alpha=u;
    nrgy=fu;
    prev_val=-dfa*alpha;
    return LS_S;
}

