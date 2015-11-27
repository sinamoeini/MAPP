#ifdef Command_Style
    CommandStyle(Command_ff,ff)
#else
#ifndef __MAPP__command_ff__
#define __MAPP__command_ff__
#include "init.h"
namespace MAPP_NS
{
    class Command_ff:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_ff(MAPP*,int,char**);
        ~Command_ff();
    };
}
#endif
#endif

