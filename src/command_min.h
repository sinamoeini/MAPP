#ifdef Command_Style
    CommandStyle(Command_min,min)
#else
#ifndef __MAPP__command_min__
#define __MAPP__command_min__
#include "init.h"
namespace MAPP_NS
{
    class Command_min:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_min(MAPP*,int,char**);
        ~Command_min();
    };
}
#endif
#endif

