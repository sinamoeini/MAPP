#ifdef Command_Style
    CommandStyle(Command_del_write,del_write)
#else
#ifndef __MAPP__command_del_write__
#define __MAPP__command_del_write__

#include <iostream>
#include "init.h"
namespace MAPP_NS
{
    class Command_del_write:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_del_write(MAPP*,int,char**);
        ~Command_del_write();
    };
}
#endif
#endif

