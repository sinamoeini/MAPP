#include "command_fix.h"
#include "memory.h"
#include "error.h"
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
        
        if(mapp->dof!=NULL)
        {
            delete mapp->dof;
            mapp->dof=NULL;
        }
        
        if(mapp->cdof!=NULL)
        {
            delete mapp->cdof;
            mapp->cdof=NULL;
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
        x_dof[i]=false;
    for(int i=0;i<no_types;i++)
        alpha_dof[i]=false;
    for(int i=0;i<no_types;i++)
        c_dof[i]=false;
    
    
    
    int c_comp;
    int alpha_comp;
    int files_started=0;
    iarg=1;
    while(iarg<nargs && files_started==0)
    {
        if(strcmp(args[iarg],"x")==0)
            x_dof[0]=true;
        else if(strcmp(args[iarg],"y")==0)
            x_dof[1]=true;
        else if(strcmp(args[iarg],"z")==0)
            x_dof[2]=true;
        else if(sscanf(args[iarg],"c[%d]",&c_comp)==1)
        {
            if(mapp->mode!=DMD_mode)
                error->abort("c[%d] can only be used in dmd mode",c_comp);
            if(c_comp>=no_types || c_comp<0)
                error->abort("invalid c component %d",c_comp);
            c_dof[c_comp]=true;
        }
        else if(sscanf(args[iarg],"alpha[%d]",&alpha_comp)==1)
        {
            if(mapp->mode!=DMD_mode)
                error->abort("alpha[%d] can only be used in dmd mode",alpha_comp);
            if(alpha_comp>=no_types || alpha_comp<0)
                error->abort("invalid alpha component %d",alpha_comp);
            alpha_dof[alpha_comp]=true;
        }
        else
        {
            files_started=1;
            continue;
        }
        iarg++;
    }
    
    if(iarg>=nargs)
        error->abort("at least 1 file have to be provided for fix command");
    
    FILE* fp=NULL;
    char* line;
    
    int natms=atoms->natms;
    int tot_natms=atoms->tot_natms;
    int* id=mapp->id->begin();
    int* list=NULL;
    int list_size=0;
    int list_cpcty=0;
    int list_grow=10;
    int list_size_tot=0;
    int iatm;
    
    CREATE1D(line,MAXCHAR);
    while (iarg<nargs)
    {
        if(atoms->my_p==0)
        {
            fp=fopen(args[iarg],"r");
            if(fp==NULL)
                error->abort("file %s not found",args[iarg]);
            fgets(line,MAXCHAR,fp);
        }
        
        while(mapp->read_line(fp,line)!=-1)
        {
            nargs=mapp->hash_remover(line);
            if(nargs==0)
                continue;
            
            if(sscanf(line,"%d",&iatm)!=1)
                error->abort("every line of file %s can contain only 1argument",args[iarg]);
            if(iatm<0 || iatm>=tot_natms)
                error->abort("the id of atoms in file %s should be between 0 & %d",args[iarg],tot_natms);
            int i=0;
            while(i<natms && id[i]!=iatm) i++;
            if(i!=natms)
            {
                if(list_size+1>list_cpcty)
                {
                    GROW(list,list_size,list_size+1+list_grow);
                    list_cpcty=list_size+1+list_grow;
                }
                
                list[list_size++]=i;
            }
            list_size_tot++;
        }
        
        if(atoms->my_p==0)
            fclose(fp);
        
        iarg++;
    }
    delete [] line;
    

    int x_xst=0;
    int alpha_xst=0;
    int c_xst=0;
    
    for(int i=0;i<dim && x_xst==0;i++)
        if(x_dof[i] && list_size_tot)
            x_xst=1;
    
    for(int i=0;i<no_types && alpha_xst==0;i++)
        if(alpha_dof[i]  && list_size_tot)
            alpha_xst=1;
    
    for(int i=0;i<no_types && c_xst==0;i++)
        if(c_dof[i]  && list_size_tot)
            c_xst=1;
    
    int x_dim=0;
    int c_dim=0;
    
    if(x_xst || alpha_xst)
    {
        x_dim=mapp->x->dim;
        if(mapp->dof==NULL)
        {
            mapp->dof=new Vec<bool>(atoms,x_dim);
            bool* dof=mapp->dof->begin();
            for(int i=0;i<natms*x_dim;i++)
                dof[i]=false;
        }
        
        if(alpha_xst)
            c_dim=mapp->c->dim;
    }
    if(c_xst)
    {
        c_dim=mapp->c->dim;
        if(mapp->cdof==NULL)
        {
            mapp->cdof=new Vec<bool>(atoms,c_dim);
            bool* cdof=mapp->cdof->begin();
            for(int i=0;i<natms*c_dim;i++)
                cdof[i]=false;
        }
    }
    
    if(x_xst)
    {
        bool* dof=mapp->dof->begin();
        int icurs;
        for(int i=0;i<list_size;i++)
        {
            icurs=list[i]*x_dim;
            for(int k=0;k<dim;k++)
            {
                if(x_dof[k])
                    dof[icurs+k]=true;
            }
        }
    }
    
    if(alpha_xst)
    {
        dmd_type* type=mapp->ctype->begin();
        bool* dof=mapp->dof->begin();
        int icurs;
        int jcurs;
        for(int i=0;i<list_size;i++)
        {
            icurs=list[i]*x_dim+dim;
            jcurs=list[i]*c_dim;
            for(int k=0;k<c_dim;k++)
            {
                if(alpha_dof[type[jcurs+k]])
                    dof[icurs+k]=true;
            }
        }
    }
    
    if(c_xst)
    {
        dmd_type* type=mapp->ctype->begin();
        bool* cdof=mapp->cdof->begin();
        int icurs;
        for(int i=0;i<list_size;i++)
        {
            icurs=list[i]*c_dim;
            for(int k=0;k<c_dim;k++)
            {
                if(c_dof[type[icurs+k]])
                    cdof[icurs+k]=true;
            }
        }
    }
    
    if(list_cpcty)
        delete [] list;
    if(dim)
        delete [] x_dof;
    if(no_types)
    {
        delete [] alpha_dof;
        delete [] c_dof;
    }
    

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_fix::~Command_fix()
{
}



