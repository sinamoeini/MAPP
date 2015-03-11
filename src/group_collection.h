#ifndef __MAPP__group_collection__
#define __MAPP__group_collection__
#include "init.h"
#include "group.h"

namespace MAPP_NS
{
    class GroupCollection:protected InitPtrs
    {
    private:
    protected:
    public:
        GroupCollection(MAPP*);
        ~GroupCollection();
        int ngroups;
        Group** groups;
        int find(char*);
        int add(char*,char*);
        int add(char*,Region*,int);
        void del(int);
    };
}
#endif 
