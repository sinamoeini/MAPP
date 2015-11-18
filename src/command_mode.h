
#ifdef Command_Style
    CommandStyle(Command_mode,mode)
#else
#ifndef __MAPP__command_mode__
#define __MAPP__command_mode__
#include "init.h"
namespace MAPP_NS
{
    class Command_mode:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_mode(MAPP*,int,char**);
        ~Command_mode();
    };
}
#endif
#endif
