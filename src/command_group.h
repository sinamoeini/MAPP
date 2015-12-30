#ifdef Command_Style
    CommandStyle(Command_group,group)
#else
#ifndef __MAPP__command_group__
#define __MAPP__command_group__
#include "init.h"
namespace MAPP_NS
{
    class Command_group: protected InitPtrs
    {
    private:
    protected:
    public:
        Command_group(MAPP*,int,char**);
        ~Command_group();
    };
}

#endif 
#endif
