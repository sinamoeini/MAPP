
#ifdef Command_Style
    CommandStyle(Command_temp,temp)
#else
#ifndef __MAPP__command_temp__
#define __MAPP__command_temp__
#include "init.h"
namespace MAPP_NS
{
    class Command_temp:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_temp(MAPP*,int,char**);
        ~Command_temp();
    };
}
#endif
#endif