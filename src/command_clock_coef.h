#ifdef Command_Style
    CommandStyle(Command_clock_coef,clock_coef)
#else
#ifndef __MAPP__command_clock_coef__
#define __MAPP__command_clock_coef__
#include "init.h"
namespace MAPP_NS
{
    class Command_clock_coef:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_clock_coef(MAPP*,int,char**);
        ~Command_clock_coef();
    };
}
#endif
#endif

