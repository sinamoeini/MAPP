#include "group.h"
#include "memory.h"
#include "atoms.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Group::Group(MAPP* mapp,char* name
,char* file_name):InitPtrs(mapp)
{
    FILE* fp=NULL;
    char* line;
    int lngth;
    group_natms=0;
 
    
    lngth=static_cast<int>(strlen(name))+1;
    CREATE1D(group_name,lngth);
    memcpy(group_name,name,lngth*sizeof(char));
    
    
    CREATE1D(line,MAXCHAR);
    
    if(atoms->my_p_no==0)
    {
        fp=fopen(file_name,"r");
        if(fp==NULL)
            error->abort("file %s not found",file_name);
    }
    
    group_natms=line_read(fp,line,file_name);
    
    if(group_natms<1 || group_natms>atoms->tot_natms)
        error->abort("number of atoms in id "
        "file %s should be between 1 & %d"
        ,file_name,atoms->tot_natms);

    

    CREATE1D(group_id_list,group_natms);
    
    for(int iatm=0;iatm<group_natms;iatm++)
    {
        group_id_list[iatm]=line_read(fp,line,file_name);
        if(group_id_list[iatm]<0 || group_id_list[iatm]>atoms->tot_natms-1)
            error->abort("atom id should be between 0 & %d",atoms->tot_natms-1);
    }
    
    
    for(int iatm=0;iatm<group_natms;iatm++)
        for(int jatm=iatm+1;jatm<group_natms;jatm++)
            if(group_id_list[iatm]==group_id_list[jatm])
                error->abort("duplicate atom id in file %s",file_name);
        
    
    if(MAXCHAR)
        delete [] line;
    
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Group::Group(MAPP* mapp,char* name
,Region* region,int inc):InitPtrs(mapp)
{
    int lngth;
    int group_natms_lcl;
    type0* x;
    int* id;
    int x_dim;
    int icmp;
    int pos;
    int* id_buff;
    int id_buff_size;
    int id_buff_grow=100;
    int group_natms_lcl_;
    int pos_;
        
    group_natms=0;
    
    lngth=static_cast<int>(strlen(name))+1;
    CREATE1D(group_name,lngth);
    memcpy(group_name,name,lngth*sizeof(char));
    

    atoms->vectors[0]->ret(x);
    atoms->vectors[atoms->find("id")]->ret(id);
    x_dim=atoms->vectors[0]->dim;
    icmp=0;
    group_natms_lcl=0;
    
    id_buff=NULL;
    id_buff_size=0;
    
    for(int iatm=0;iatm<atoms->natms;iatm++)
    {
        if(region->belong(&x[icmp])==inc)
        {
            if(group_natms_lcl==id_buff_size)
            {
                GROW(id_buff,id_buff_size,id_buff_size+id_buff_grow);
                id_buff_size+=id_buff_grow;
            }
            id_buff[group_natms_lcl]=id[iatm];
            group_natms_lcl++;
        }
        icmp+=x_dim;
    }
    
    MPI_Allreduce(&group_natms_lcl,&group_natms,1,MPI_INT,MPI_SUM,world);
    
    MPI_Scan(&group_natms_lcl,&pos,1,MPI_INT,MPI_SUM,world);
    
    pos-=group_natms_lcl;
    
    CREATE1D(group_id_list,group_natms);
    
    memcpy(&group_id_list[pos],id_buff,group_natms_lcl*sizeof(int));

    if(id_buff_size)
        delete [] id_buff;
    
    
    group_natms_lcl_=0;
    pos_=0;
    
    for(int iproc=0;iproc<atoms->tot_p;iproc++)
    {
        if(iproc==atoms->my_p_no)
        {
            pos_=pos;
            group_natms_lcl_=group_natms_lcl;
        }
        
        MPI_Bcast(&pos_,1,MPI_INT,iproc,world);
        MPI_Bcast(&group_natms_lcl_,1,MPI_INT,iproc,world);
        MPI_Bcast(&group_id_list[pos_],group_natms_lcl_,MPI_INT,iproc,world);
    }
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Group::~Group()
{
    delete [] group_name;
    if(group_natms)
        delete [] group_id_list;
}
/*--------------------------------------------
 read line and broadcast
 --------------------------------------------*/
int inline Group::line_read(FILE* file
,char*& line,char* file_name)
{
    int eof=0;
    int lngth;
    int icurs;
    int narg;

    while (narg==0)
    {
        if(atoms->my_p_no==0)
            if(feof(file))
                eof=-1;
        MPI_Bcast(&eof,1,MPI_INT,0,world);
        if(eof==-1)
             error->abort("file %s ended immaturely",file_name);
        
        
        if(atoms->my_p_no==0)
        {
            fgets(line,MAXCHAR,file);
            lngth=static_cast<int>(strlen(line))+1;
        }
        MPI_Bcast(&lngth,1,MPI_INT,0,world);
        MPI_Bcast(line,lngth,MPI_CHAR,0,world);
        
        char* srch=NULL;
        srch=strchr(line,'#');
        if(srch!=NULL)
            *srch='\0';
        
        lngth=static_cast<int>(strlen(line));
        icurs=0;
        
        while(icurs<lngth)
        {
            if(isspace(line[icurs])&&icurs<lngth)
            {
                while(isspace(line[icurs])&&icurs<lngth)
                    icurs++;
            }
            
            if(isspace(line[icurs])==0&&icurs<lngth)
            {
                narg++;
                while(isspace(line[icurs])==0&&icurs<lngth)
                    icurs++;
            }
        }
    }
    
    
    if(narg!=1)
        error->abort("wrong format in id file %s",file_name);
    int tmp;
    sscanf(line,"%d",&tmp);
        
    return tmp;
}
