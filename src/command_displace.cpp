#include "command_displace.h"
#include "memory.h"
#include "error.h"
#include "atoms.h"
#include "group.h"
#include <stdlib.h>
using namespace MAPP_NS;
enum {disp_x,disp_s};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_displace::Command_displace(int nargs,char** args)
{
    int mode=-1;
    int dim=__dim__;
    if(nargs<3+dim)
        error->abort("incorrect displace command");
    if(strcmp(args[1],"x")==0)
        mode=disp_x;
    else if(strcmp(args[1],"s")==0)
        mode=disp_s;
    else
        error->abort("unknown keyword %s",args[1]);
        
    type0* disp;
    CREATE1D(disp,dim);
    for(int iarg=2;iarg<dim+2;iarg++)
        disp[iarg-2]=atof(args[iarg]);
    
    int iarg=dim+2;
    
    if(strcmp(args[iarg],"group")==0)
    {
    }
    else
        error->abort("unknown keyword %s",args[iarg]);
    iarg++;
    
    int ngrps=nargs-iarg;
    int** grp_idx=new int*[ngrps];
    int* grp_sz=new int[ngrps];
    
    for(int i=0;i<ngrps;i++)
    {
        Group* grp=groups->find_grp(args[iarg]);
        grp->get_idx(grp_sz[i],grp_idx[i]);
        iarg++;
    }
    
    type0* x=atoms->x->begin();
    type0* xi;
    int x_dim=atoms->x->dim;
    
    if(mode==disp_s)
        V3M_LT(disp,atoms->H,disp);
    
    for(int igrp=0;igrp<ngrps;igrp++)
    {
        for(int i=0;i<grp_sz[igrp];i++)
        {
            xi=x+x_dim*grp_idx[igrp][i];
            for(int idim=0;idim<dim;idim++)
                xi[idim]+=disp[idim];
        }
    }

    if(ngrps)
    {
        delete [] grp_idx;
        delete [] grp_sz;
    }
    delete [] disp;
    
    atoms->reset();
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_displace::~Command_displace()
{
    
}
