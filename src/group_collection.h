#ifndef __MAPP__group_collection__
#define __MAPP__group_collection__
#include <iostream>
#include "init.h"
#include "group.h"
#include "region.h"
using namespace std;
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
#endif /* defined(__MAPP__group_collection__) */
