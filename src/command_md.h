#ifdef Command_Style
    CommandStyle(Command_md,md)
#else
#ifndef __MAPP__command_md__
#define __MAPP__command_md__
#include "init.h"
namespace MAPP_NS
{
    class Command_md:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_md(MAPP*,int,char**);
        ~Command_md();
    };
}
#endif
#endif

