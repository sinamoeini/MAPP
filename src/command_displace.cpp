#include "command_displace.h"
#include "memory.h"
#include "error.h"
#include "atoms.h"
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
    
    type0* x=mapp->x->begin();
    int x_dim=mapp->x->dim;

    for(int i=0;i<list_size;i++)
    {
        iatm=list[i];
        for(int idim=0;idim<dim;idim++)
            x[i*x_dim+idim]+=disp[idim];
    }
    
    delete [] line;
    delete [] disp;
    
    atoms->reset();
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_displace::~Command_displace()
{
    
}
