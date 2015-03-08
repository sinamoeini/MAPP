#ifdef Command_Style
    CommandStyle(Command_comm_mode,comm_mode)
#else
#ifndef __MAPP__command_comm_mode__
#define __MAPP__command_comm_mode__

#include <stdio.h>
#include "init.h"
namespace MAPP_NS
{
    class Command_comm_mode:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_comm_mode(MAPP*,int,char**);
        ~Command_comm_mode();
        
    };
}
#endif
#endif
