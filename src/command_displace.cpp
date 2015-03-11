#include "command_displace.h"
#include <stdlib.h>
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_displace::Command_displace(MAPP* mapp
,int narg,char** args):InitPtrs(mapp)
{
    if((narg-1)%4!=0)
        error->abort("wrong displace command");
    int dim=atoms->dimension;
    
    int no_traj=(narg-1)/(1+dim);
    int iarg=1;
    int no_atoms=0;
    char* line;
    CREATE1D(line,MAXCHAR);
    
    FILE* fp=NULL;
    int* list;
    type0* disp;
    CREATE1D(disp,dim);
    
    id_n=atoms->find("id");
    
    
    
    
    for(int i=0;i<no_traj;i++)
    {
        if(atoms->my_p_no==0)
        {
            fp=fopen(args[iarg],"r");
            if(fp==NULL)
                error->abort("file %s not found",args[iarg]);
            
            fgets(line,MAXCHAR,fp);
            sscanf(line,"%d",&no_atoms);
            
        }
        
        MPI_Bcast(&no_atoms,1,MPI_INT,0,world);
        CREATE1D(list,no_atoms);
        if(atoms->my_p_no==0)
        {
            for(int iatm=0;iatm<no_atoms;iatm++)
            {
                fgets(line,MAXCHAR,fp);
                sscanf(line,"%d",&list[iatm]);
            }
            
            
            fclose(fp);
        }
        MPI_Bcast(list,no_atoms,MPI_INT,0,world);
        
        iarg++;
        
        for(int j=0;j<dim;j++)
            disp[j]=atof(args[iarg++]);
        

        
        move(list,no_atoms,disp);

        if(no_atoms)
            delete [] list;
    }

    
    delete [] line;
    delete [] disp;
    
    int* tmp_list;
    int no_vecs=atoms->no_vecs;
    CREATE1D(tmp_list,no_vecs);
    for(int i=0;i<no_vecs;i++)
        tmp_list[i]=i;
    
    VecLst* vecs_comm=new VecLst(mapp,tmp_list,no_vecs);
    delete [] tmp_list;
    
    atoms->x2s(atoms->natms);
    atoms->xchng_cmplt(vecs_comm);
    atoms->s2x(atoms->natms);
    
    
    delete vecs_comm;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_displace::~Command_displace()
{
    
}
/*--------------------------------------------
 move
 --------------------------------------------*/
void Command_displace::move(int* list,
int no_atoms,type0* disp)
{
    for(int i=0;i<no_atoms;i++)
        if(list[i]<0 || list[i]>atoms->tot_natms-1)
            error->abort("invalid atom id %d",list[i]);
    int* id;
    atoms->vectors[id_n]->ret(id);
    type0* x;
    atoms->vectors[0]->ret(x);
    int x_dim=atoms->vectors[0]->dim;
    int dim=atoms->dimension;
    
    int natms=atoms->natms;
    for(int i=0;i<natms;i++)
        for(int j=0;j<no_atoms;j++)
            if(id[i]==list[j])
                for(int k=0;k<dim;k++)
                    x[i*x_dim+k]+=disp[k];
    
    
}
