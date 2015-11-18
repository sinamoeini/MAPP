#ifndef __MAPP__group__
#define __MAPP__group__
#include "init.h"
namespace MAPP_NS
{
    class Group:protected InitPtrs
    {
    private:
        int inline line_read(FILE*,char*&,char*);
    protected:
    public:
        Group(MAPP*);
        ~Group();
        
        void update_ids();
        int* group_list;
        int group_list_size;
    };
}
#endif
