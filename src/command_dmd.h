#ifdef Command_Style
    CommandStyle(Command_dmd,dmd)
#else
#ifndef __MAPP__command_dmd__
#define __MAPP__command_dmd__
#include "init.h"
namespace MAPP_NS
{
    class Command_dmd:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_dmd(MAPP*,int,char**);
        ~Command_dmd();
    };
}
#endif
#endif

