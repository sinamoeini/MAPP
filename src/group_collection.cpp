#include "group_collection.h"
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
GroupCollection::GroupCollection(MAPP* mapp)
:InitPtrs(mapp)
{
    ngroups=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
GroupCollection::~GroupCollection()
{
    if(ngroups)
    {
        for(int igroup=0;igroup<ngroups;igroup++)
            delete groups[igroup];
        
        delete [] groups;
    }
}
/*--------------------------------------------
 add a group
 --------------------------------------------*/
int GroupCollection::add(char* group_name
,char* file_name)
{
    Group** new_groups;
    CREATE1D(new_groups,ngroups+1);
    
    for(int igroup=0;igroup<ngroups;igroup++)
        new_groups[igroup]=groups[igroup];
    
    new_groups[ngroups]=new Group(mapp,group_name,file_name);
    
    if(ngroups)
        delete [] groups;
    groups=new_groups;
    ngroups++;
    return ngroups-1;
}
/*--------------------------------------------
 add a group
 --------------------------------------------*/
int GroupCollection::add(char* group_name
,Region* region,int inc)
{
    Group** new_groups;
    CREATE1D(new_groups,ngroups+1);
    
    for(int igroup=0;igroup<ngroups;igroup++)
        new_groups[igroup]=groups[igroup];
    
    new_groups[ngroups]=new Group(mapp,group_name,region,inc);
    
    if(ngroups)
        delete [] groups;
    groups=new_groups;
    ngroups++;
    return ngroups-1;
}
/*--------------------------------------------
 find the group id
 --------------------------------------------*/
int GroupCollection::find(char* group_name)
{
    for(int igroup=0;igroup<ngroups;igroup++)
        if(strcmp(group_name,groups[igroup]->group_name)==0)
            return igroup;
    return -1;
}
/*--------------------------------------------
 remove a group
 --------------------------------------------*/
void GroupCollection::del(int igroup)
{
    Group** new_groups;
    CREATE1D(new_groups,ngroups-1);
    
    for(int i=0;i<igroup;i++)
        new_groups[i]=groups[i];
    
    for(int i=igroup+1;i<ngroups;i++)
        new_groups[i-1]=groups[i];
    
    if(ngroups)
        delete [] groups;
    groups=new_groups;
    ngroups--;
}




