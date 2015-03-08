#ifndef __MAPP__group__
#define __MAPP__group__
#include <iostream>
#include "init.h"
#include "region.h"
using namespace std;
namespace MAPP_NS
{
    class Group:protected InitPtrs
    {
    private:
        int inline line_read(FILE*,char*&,char*);
    protected:
    public:
        Group(MAPP*,char*,char*);
        Group(MAPP*,char*,Region*,int);
        ~Group();
        
        char* group_name;
        int* group_id_list;
        int group_natms;
    };
}
#endif
