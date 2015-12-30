#include "command_displace.h"
#include "memory.h"
#include "error.h"
#include "atoms.h"
#include "group.h"
#include <stdlib.h>
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_displace::Command_displace(MAPP* mapp
,int nargs,char** args):InitPtrs(mapp)
{
    int dim=atoms->dimension;
    if(nargs<2+dim)
        error->abort("incorrect displace command");
    
    type0* disp;
    CREATE1D(disp,dim);
    for(int iarg=1;iarg<dim+1;iarg++)
        disp[iarg-1]=atof(args[iarg]);
    
    int iarg=dim+1;
    
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
    
    type0* x=mapp->x->begin();
    type0* xi;
    int x_dim=mapp->x->dim;
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
