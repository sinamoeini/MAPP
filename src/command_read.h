#ifdef Command_Style
    CommandStyle(Command_read,read)
#else
#ifndef __MAPP__command_read__
#define __MAPP__command_read__
#include "init.h"
namespace MAPP_NS
{
    class Command_read:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_read(MAPP*,int,char**);
        ~Command_read();
    };
}
#endif
#endif

