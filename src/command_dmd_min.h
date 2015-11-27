#ifdef Command_Style
    CommandStyle(Command_dmd_min,dmd_min)
#else
#ifndef __MAPP__command_dmd_min__
#define __MAPP__command_dmd_min__
#include "init.h"
namespace MAPP_NS
{
    class Command_dmd_min:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_dmd_min(MAPP*,int,char**);
        ~Command_dmd_min();
    };
}
#endif
#endif

