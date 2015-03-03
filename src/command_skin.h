#ifdef Command_Style
    CommandStyle(Command_skin,skin)
#else
#ifndef __MAPP__command_skin__
#define __MAPP__command_skin__

#include <iostream>
#include "init.h"
namespace MAPP_NS
{
    class Command_skin:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_skin(MAPP*,int,char**);
        ~Command_skin();
    };
}
#endif
#endif

