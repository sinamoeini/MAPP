#ifdef Command_Style
    CommandStyle(Command_dmd_coef,dmd_coef)
#else
#ifndef __MAPP__command_dmd_coef__
#define __MAPP__command_dmd_coef__
#include "init.h"
namespace MAPP_NS
{
    class Command_dmd_coef:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_dmd_coef(MAPP*,int,char**);
        ~Command_dmd_coef();
    };
}
#endif
#endif

