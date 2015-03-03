#ifdef Command_Style
    CommandStyle(Command_grid,grid)
#else
#ifndef __MAPP__command_grid__
#define __MAPP__command_grid__

#include <iostream>
#include "init.h"
namespace MAPP_NS
{
    class Command_grid:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_grid(MAPP*,int,char**);
        ~Command_grid();
    };
}
#endif
#endif

