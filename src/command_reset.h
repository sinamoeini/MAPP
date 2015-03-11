#ifdef Command_Style
    CommandStyle(Command_reset,reset)
#else
#ifndef __MAPP__command_reset__
#define __MAPP__command_reset__
#include "init.h"
namespace MAPP_NS
{
    class Command_reset:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_reset(MAPP*,int,char**);
        ~Command_reset();
    };
}
#endif
#endif

