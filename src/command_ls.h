#ifdef Command_Style
    CommandStyle(Command_ls,ls)
#else
#ifndef __MAPP__command_ls__
#define __MAPP__command_ls__
#include "init.h"
namespace MAPP_NS
{
    class Command_ls:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_ls(MAPP*,int,char**);
        ~Command_ls();
    };
}
#endif
#endif

