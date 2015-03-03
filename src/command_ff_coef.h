#ifdef Command_Style
    CommandStyle(Command_ff_coef,ff_coef)
#else
#ifndef __MAPP__command_ff_coef__
#define __MAPP__command_ff_coef__

#include <iostream>
#include "init.h"
namespace MAPP_NS
{
    class Command_ff_coef:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_ff_coef(MAPP*,int,char**);
        ~Command_ff_coef();
    };
}
#endif
#endif

