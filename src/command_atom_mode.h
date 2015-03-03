#ifdef Command_Style
    CommandStyle(Command_atom_mode,atom_mode)
#else
#ifndef __MAPP__command_atom_mode__
#define __MAPP__command_atom_mode__

#include <iostream>
#include "init.h"
namespace MAPP_NS
{
    class Command_atom_mode:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_atom_mode(MAPP*,int,char**);
        ~Command_atom_mode();
        
    };
}
#endif
#endif

