#ifdef Command_Style
    CommandStyle(Command_precision,precision)
#else
#ifndef __MAPP__command_precision__
#define __MAPP__command_precision__
#include "init.h"
namespace MAPP_NS
{
    class Command_precision:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_precision(MAPP*,int,char**);
        ~Command_precision();
    };
}
#endif
#endif
