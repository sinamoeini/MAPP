#ifdef Command_Style
    CommandStyle(Command_boltzmann,boltzmann)
#else
#ifndef __MAPP__command_boltzmann__
#define __MAPP__command_boltzmann__

#include <stdlib.h>
#include "init.h"
namespace MAPP_NS
{
    class Command_boltzmann:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_boltzmann(MAPP*,int,char**);
        ~Command_boltzmann();
    };
}
#endif
#endif

