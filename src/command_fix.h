#ifdef Command_Style
    CommandStyle(Command_fix,fix)
#else
#ifndef __MAPP__command_fix__
#define __MAPP__command_fix__
#include "init.h"
namespace MAPP_NS {
    class Command_fix: protected InitPtrs
    {
    private:
    protected:
        void md(int,char**);
        void dmd(int,char**);
    public:
        Command_fix(MAPP*,int,char**);
        ~Command_fix();
    };
}

#endif 
#endif
