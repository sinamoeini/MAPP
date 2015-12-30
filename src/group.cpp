#include "group.h"
#include "atoms.h"
#include "atom_types.h"
#include "error.h"
#include "memory.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 read line and broadcast
 --------------------------------------------*/
Group::Group(MAPP* mapp,const char* _group_name)
:InitPtrs(mapp)
{
    grp_sz=0;
    int lngth=static_cast<int>(strlen(_group_name))+1;
    name=new char[lngth];
    memcpy(name,_group_name,lngth*sizeof(char));
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Group::~Group()
{
    if (grp_sz)
        delete [] grp_idx;
    int lngth=static_cast<int>(strlen(name))+1;
    if(lngth)
        delete [] name;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void Group::get_idx(int& grp_sz_,int*& grp_idx_)
{
    update();
    grp_idx_=grp_idx;
    grp_sz_=grp_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Group_general::Group_general(MAPP* mapp,
const char* group_name,int nfiles,char** files)
:Group(mapp,group_name)
{
    
    int* id=mapp->id->begin();
    int natms=atoms->natms;
    int tot_natms=atoms->tot_natms;
    int nargs,iatm;
    char* line;
    CREATE1D(line,MAXCHAR);
    
    int* rank=new int[natms];
    for(int i=0;i<natms;i++) rank[i]=i;
    
    
    XMath* xmath=new XMath();
    
    xmath->quicksort(rank,rank+natms
    ,[&id](int* rank_i,int* rank_j){return (id[*rank_i]<id[*rank_j]);}
    ,[](int* rank_i,int* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    delete xmath;
   
    auto binary_srch=
    [] (int i,int* id,int* rank,int size)->int
    {
        if(i>id[rank[size-1]] || i<id[rank[0]])
            return -1;
        int mid;
        while(size)
        {
            mid=(size/2);
            
            if(i<id[rank[mid]])
                size=mid;
            else if(i>id[rank[mid]])
            {
                rank+=mid+1;
                size-=mid+1;
            }
            else
                return rank[mid];
        }
        
        return -1;
    };
    
    
    int* list=NULL;
    int list_size=0;
    int list_cpcty=0;
    int list_grow=10;
    int list_size_tot=0;
    int found;
    
    FILE* fp;
    
    for(int ifile=0;ifile<nfiles;ifile++)
    {
        mapp->open_file(fp,files[ifile],"r");

        while(mapp->read_line(fp,line)!=-1)
        {
            nargs=mapp->hash_remover(line);
            if(nargs==0)
                continue;
            
            if(sscanf(line,"%d",&iatm)!=1)
                error->abort("every line of file %s can contain only 1 argument",files[ifile]);
            if(iatm<0 || iatm>=tot_natms)
                error->abort("the id of atoms in file %s should be between 0 & %d",files[ifile],tot_natms);

            found=binary_srch(iatm,id,rank,natms);

            if(found!=-1)
            {
                if(list_size+1>list_cpcty)
                {
                    GROW(list,list_size,list_size+1+list_grow);
                    list_cpcty=list_size+1+list_grow;
                }
                
                list[list_size++]=found;
            }
            list_size_tot++;
        }

        if(atoms->my_p==0)
            fclose(fp);
    }
    
    delete [] line;

    grp_sz_tot=list_size_tot;
    grp_sz=list_size;
    grp_id=new int[grp_sz];
    grp_idx=new int[grp_sz];
    memcpy(grp_idx,list,grp_sz*sizeof(int));
    delete [] list;
    for(int i=0;i<grp_sz;i++)
        grp_id[i]=id[grp_idx[i]];
    atoms->xchng_chk(xchng_id);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Group_general::~Group_general()
{
    if (grp_sz)
        delete [] grp_id;
}
/*--------------------------------------------

 --------------------------------------------*/
void Group_general::update()
{
    if(atoms->xchng_chk(xchng_id))
        return;
    
    int natms=atoms->natms;
    int tot_p=atoms->tot_p;
    int* grp_sz_per_p=new int[tot_p];
    MPI_Allgather(&grp_sz,1,MPI_INT,grp_sz_per_p,1,MPI_INT,world);
    
    int* disp=new int[tot_p];
    *disp=0;
    for(int i=1;i<tot_p;i++)
        disp[i]=disp[i-1]+grp_sz_per_p[i-1];

    int* grp_id_tot=new int[grp_sz_tot];
    MPI_Allgatherv(grp_id,grp_sz,MPI_INT,grp_id_tot,grp_sz_per_p,disp,MPI_INT,world);
    delete [] disp;
    delete [] grp_sz_per_p;
    
    if(grp_sz)
    {
        delete [] grp_idx;
        delete [] grp_id;
    }
    grp_idx=new int[natms];
    grp_id=new int[natms];
    grp_sz=0;
    int* id=mapp->id->begin();
    
    class act
    {
    public:
        virtual void eq(int)=0;
        virtual void neq(int)=0;
    };
    
    class act_orig: public act
    {
    public:
        act_orig(){}
        void eq(int){}
        void neq(int){}
    };
    
    class act_grp: public act
    {
    public:
        int* id;
        
        int*& grp_id;
        int*& grp_idx;
        int& grp_sz;
        
        act_grp(
        int*& _grp_id,
        int*& _grp_idx,
        int& _grp_sz,
        int* _id):
        grp_id(_grp_id),
        grp_idx(_grp_idx),
        grp_sz(_grp_sz),
        id(_id){}
        void eq(int idx)
        {
            grp_idx[grp_sz]=idx;
            grp_id[grp_sz++]=id[idx];
        }
        void neq(int){}
    };
    
    
    act_orig* _act_orig_=new act_orig();
    act_grp* _act_grp_=new act_grp(grp_id,grp_idx,grp_sz,id);
    
    act* _act_orig=_act_orig_;
    act* _act_grp=_act_grp_;
    XMath* xmath=new XMath();
    xmath->srch_lst_lst(id,natms,_act_grp,grp_id_tot,grp_sz_tot,_act_orig);
    if(grp_sz_tot) delete [] grp_id_tot;
    delete _act_orig_;
    delete _act_grp_;
    delete xmath;
    
    int* _grp_idx=new int[grp_sz];
    int* _grp_id=new int[grp_sz];
    memcpy(_grp_id,grp_id,grp_sz*sizeof(int));
    memcpy(_grp_idx,grp_idx,grp_sz*sizeof(int));
    if(natms)
    {
        delete [] grp_id;
        delete [] grp_idx;
    }
    
    grp_id=_grp_id;
    grp_idx=_grp_idx;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
Group_type::Group_type(MAPP* mapp,const char*
group_name,int ntypes_,char** types)
:Group(mapp,group_name)
{
    ntypes=atom_types->no_types;
    type=new bool[ntypes];
    for(int ityp=0;ityp<ntypes;ityp++)
        type[ityp]=false;
    
    for(int ityp=0;ityp<ntypes;ityp++)
        type[atom_types->find_type(types[ityp])]=true;
    
    update();
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
Group_type::~Group_type()
{
    if(ntypes)
        delete [] type;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Group_type::update()
{
    int natms=atoms->natms;
    int* _grp_idx=new int[natms];
    dmd_type* _type=mapp->type->begin();
    grp_sz=0;
    for(int i=0;i<natms;i++)
        if(type[_type[i]])
            _grp_idx[grp_sz++]=i;
    grp_idx=new int[grp_sz];
    memcpy(grp_idx,_grp_idx,grp_sz*sizeof(int));
    if(natms)
        delete [] _grp_idx;
    MPI_Allreduce(&grp_sz,&grp_sz_tot,1,MPI_INT,MPI_SUM,world);
}
/*--------------------------------------------
 read line and broadcast
 --------------------------------------------*/
GroupCollection::GroupCollection(MAPP* mapp)
:InitPtrs(mapp)
{
    ngrps=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
GroupCollection::~GroupCollection()
{
    for(int igrp=0;igrp<ngrps;igrp++)
        delete grps[igrp];
    if(ngrps)
        delete [] grps;
}
/*--------------------------------------------
 find a group
 --------------------------------------------*/
Group* GroupCollection::find_grp(const char* name)
{
    for(int igrp=0;igrp<ngrps;igrp++)
        if(strcmp(name,grps[igrp]->name)==0)
            return grps[igrp];
    
    error->abort("group %s does not exist",name);
    
    return NULL;
}
/*--------------------------------------------
 add a new new group
 --------------------------------------------*/
Group* GroupCollection::add_grp_file(const char* name,int nfiles,char** files)
{
    for(int igrp=0;igrp<ngrps;igrp++)
        if(strcmp(name,grps[igrp]->name)==0)
        {
            error->abort("group %s already exists",name);
            return NULL;
        }
    
    Group* grp=new Group_general(mapp,name,nfiles,files);
    
    Group** _grps=new Group*[ngrps+1];
    memcpy(_grps,grps,ngrps*sizeof(Group*));
    if(ngrps)
        delete [] grps;
    grps=_grps;
    grps[ngrps]=grp;
    ngrps++;
    return grp;
}
/*--------------------------------------------
 add a new new group
 --------------------------------------------*/
Group* GroupCollection::add_grp_type(const char* name,int ntypes,char** types)
{
    for(int igrp=0;igrp<ngrps;igrp++)
        if(strcmp(name,grps[igrp]->name)==0)
        {
            error->abort("group %s already exists",name);
            return NULL;
        }
    
    Group* grp=new Group_type(mapp,name,ntypes,types);
    
    Group** _grps=new Group*[ngrps+1];
    memcpy(_grps,grps,ngrps*sizeof(Group*));
    if(ngrps)
        delete [] grps;
    grps=_grps;
    grps[ngrps]=grp;
    ngrps++;
    return grp;
    
    return NULL;
}
/*--------------------------------------------
 remove a vector
 --------------------------------------------*/
void GroupCollection::del_grp(Group* grp)
{
    int igrp=0;
    while(grps[igrp]!=grp && igrp<ngrps)
        igrp++;
    if(igrp==ngrps)
        error->abort("group was not found");
    
    Group** _grps=new Group*[ngrps-1];
    int j=0;
    for(int i=0;i<ngrps;i++)
        if(grps[i]!=grp)
            _grps[j++]=grps[i];
    delete [] grps;
    grps=_grps;
    ngrps--;

}





