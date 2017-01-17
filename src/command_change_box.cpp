#include "command_change_box.h"
#include "atoms.h"
#include "xmath.h"
#include "error.h"
#include "memory.h"
#include <stdlib.h>
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_change_box::Command_change_box(int nargs,char** args)
{
    
    type0 A[__dim__][__dim__];
    
    int iarg=1;
    
    if(strcmp(args[iarg],"strain")==0)
    {
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                A[i][j]=0.0;
    }
    else if(strcmp(args[iarg],"dilation")==0)
    {
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                A[i][j]=1.0;
    }
    else if(strcmp(args[iarg],"equal")==0)
    {
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                A[i][j]=atoms->H[i][j];
    }
    else
        Error::abort("unknown keyword for change_box: %s",args[iarg]);
    
    iarg++;
    
    
    int icmp,jcmp;
    while(iarg<nargs)
    {
        if(sscanf(args[iarg],"H[%d][%d]",&icmp,&jcmp)==2)
        {
            if(icmp<0 || icmp>=__dim__ || jcmp<0 || jcmp>=__dim__)
                Error::abort("wrong component in change_box for H[%d][%d]",icmp,jcmp);
            iarg++;
            
            if(nargs-iarg<1)
                Error::abort("H[%d][%d] in change_box followed by 1 arguement",icmp,jcmp);

            
            if(icmp<=jcmp)
                A[jcmp][icmp]=atof(args[iarg]);
            else
                A[icmp][jcmp]=atof(args[iarg]);
            iarg++;
            
        }
        else
            Error::abort("unknown keyword for change_box: %s",args[iarg]);
    }    
    
    atoms->x2s_lcl();
    
    if(strcmp(args[1],"strain")==0)
    {
        for(int i=0;i<__dim__;i++)
            A[i][i]+=1.0;
        type0 C[__dim__][__dim__];

        
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                C[i][j]=0.0;
        
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                for(int k=0;k<__dim__;k++)
                    C[i][j]+=atoms->H[i][k]*A[k][j];
        
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                atoms->H[i][j]=C[i][j];
        
    }
    else if(strcmp(args[1],"dilation")==0)
    {
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                atoms->H[i][j]*=A[i][j];
    }
    else if(strcmp(args[1],"equal")==0)
    {
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
                atoms->H[i][j]=A[i][j];
    }
    
    if(__dim__==3)
        M3INV_LT(atoms->H,atoms->B);
    else
        XMath::invert_lower_triangle(atoms->H,atoms->B);
    
    atoms->s2x_lcl();
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_change_box::~Command_change_box()
{
}
