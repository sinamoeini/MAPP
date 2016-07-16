#ifdef Command_Style
    CommandStyle(Command_planck,planck)
#else
#ifndef __MAPP__command_planck__
#define __MAPP__command_planck__
#include "init.h"
namespace MAPP_NS
{
    class Command_planck:protected InitPtrs
    {
    private:
    protected:
    public:
        Command_planck(MAPP*,int,char**);
        ~Command_planck();
    };
}
#endif
#endif

