#ifdef Command_Style
    CommandStyle(Command_write,write)
#else
#ifndef __MAPP__command_write__
#define __MAPP__command_write__
#include "init.h"
namespace MAPP_NS
{
    class Command_write:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_write(MAPP*,int,char**);
        ~Command_write();
    };
}
#endif
#endif

