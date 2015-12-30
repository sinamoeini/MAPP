#include "command_fix.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "atom_types.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_fix::Command_fix(MAPP* mapp
,int nargs,char** args):InitPtrs(mapp)
{
    if(strcmp(args[1],"release")==0)
    {
        if(nargs!=2)
            error->abort("incorrect fix command");
        
        if(mapp->x_dof!=NULL)
        {
            delete mapp->x_dof;
            mapp->x_dof=NULL;
        }
        
        if(mapp->c_dof!=NULL)
        {
            delete mapp->c_dof;
            mapp->c_dof=NULL;
        }
        
        return;
    }
    
    int iarg;
    int dim=atoms->dimension;
    int no_types=atom_types->no_types;
    byte* x_dof;
    byte* alpha_dof;
    byte* c_dof;
    
    CREATE1D(x_dof,dim);
    CREATE1D(alpha_dof,no_types);
    CREATE1D(c_dof,no_types);
    for(int i=0;i<dim;i++)
        x_dof[i]=true;
    for(int i=0;i<no_types;i++)
        alpha_dof[i]=true;
    for(int i=0;i<no_types;i++)
        c_dof[i]=true;
    
    
    
    int c_comp;
    int alpha_comp;
    bool grps_started=false;
    iarg=1;
    while(iarg<nargs && !grps_started)
    {
        if(strcmp(args[iarg],"x")==0)
            x_dof[0]=false;
        else if(strcmp(args[iarg],"y")==0)
            x_dof[1]=false;
        else if(strcmp(args[iarg],"z")==0)
            x_dof[2]=false;
        else if(sscanf(args[iarg],"c[%d]",&c_comp)==1)
        {
            if(mapp->mode!=DMD_mode)
                error->abort("c[%d] can only be used in dmd mode",c_comp);
            if(c_comp>=no_types || c_comp<0)
                error->abort("invalid c component %d",c_comp);
            c_dof[c_comp]=false;
        }
        else if(sscanf(args[iarg],"alpha[%d]",&alpha_comp)==1)
        {
            if(mapp->mode!=DMD_mode)
                error->abort("alpha[%d] can only be used in dmd mode",alpha_comp);
            if(alpha_comp>=no_types || alpha_comp<0)
                error->abort("invalid alpha component %d",alpha_comp);
            alpha_dof[alpha_comp]=false;
        }
        else
        {
            grps_started=true;
            continue;
        }
        iarg++;
    }
    
    if(iarg>=nargs)
        error->abort("at least 2 file have to be provided for fix command");
    
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
    


    bool x_xst=false;
    bool alpha_xst=false;
    bool c_xst=false;
    
    for(int i=0;i<dim && x_xst==0;i++)
        if(!x_dof[i])
            x_xst=true;
    
    for(int i=0;i<no_types && alpha_xst==0;i++)
        if(!alpha_dof[i])
            alpha_xst=true;
    
    for(int i=0;i<no_types && c_xst==0;i++)
        if(!c_dof[i])
            c_xst=true;
    
    int x_dim=0;
    int c_dim=0;
    int natms=atoms->natms;
    if(x_xst || alpha_xst)
    {
        x_dim=mapp->x->dim;
        if(mapp->x_dof==NULL)
        {
            mapp->x_dof=new Vec<bool>(atoms,x_dim);
            bool* dof=mapp->x_dof->begin();
            for(int i=0;i<natms*x_dim;i++)
                dof[i]=true;
        }
        
        if(alpha_xst)
            c_dim=mapp->c->dim;
    }
    if(c_xst)
    {
        c_dim=mapp->c->dim;
        if(mapp->c_dof==NULL)
        {
            mapp->c_dof=new Vec<bool>(atoms,c_dim);
            bool* cdof=mapp->c_dof->begin();
            for(int i=0;i<natms*c_dim;i++)
                cdof[i]=true;
        }
    }
    
    if(x_xst)
    {
        bool* dof=mapp->x_dof->begin();
        int icurs;
        
        for(int igrp=0;igrp<ngrps;igrp++)
        {
            for(int i=0;i<grp_sz[igrp];i++)
            {
                icurs=grp_idx[igrp][i]*x_dim;
                for(int k=0;k<dim;k++)
                {
                    if(!x_dof[k])
                        dof[icurs+k]=false;
                }
            }
        }
    }
    
    if(alpha_xst)
    {
        dmd_type* type=mapp->ctype->begin();
        bool* dof=mapp->x_dof->begin();
        int icurs;
        int jcurs;
        
        for(int igrp=0;igrp<ngrps;igrp++)
        {
            for(int i=0;i<grp_sz[igrp];i++)
            {
                icurs=grp_idx[igrp][i]*x_dim+dim;
                jcurs=grp_idx[igrp][i]*c_dim;
                for(int k=0;k<c_dim;k++)
                {
                    if(!alpha_dof[type[jcurs+k]])
                        dof[icurs+k]=false;
                }
            }
        }
    }
    
    if(c_xst)
    {
        dmd_type* type=mapp->ctype->begin();
        bool* cdof=mapp->c_dof->begin();
        int icurs;
        
        for(int igrp=0;igrp<ngrps;igrp++)
        {
            for(int i=0;i<grp_sz[igrp];i++)
            {
                icurs=grp_idx[igrp][i]*c_dim;
                for(int k=0;k<c_dim;k++)
                {
                    if(!c_dof[type[icurs+k]])
                        cdof[icurs+k]=false;
                }
            }
        }
    }
    
    if(dim)
        delete [] x_dof;
    if(no_types)
    {
        delete [] alpha_dof;
        delete [] c_dof;
    }
    
    if(ngrps)
    {
        delete [] grp_idx;
        delete [] grp_sz;
    }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_fix::~Command_fix()
{
}



