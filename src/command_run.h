#ifdef Command_Style
    CommandStyle(Command_run,run)
#else
#ifndef __MAPP__command_run__
#define __MAPP__command_run__
#include "init.h"
namespace MAPP_NS
{
    class Command_run:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_run(MAPP*,int,char**);
        ~Command_run();
    };
}
#endif
#endif

