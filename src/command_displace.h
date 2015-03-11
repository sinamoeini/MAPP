#ifdef Command_Style
    CommandStyle(Command_displace,displace)
#else
#ifndef __MAPP__command_displace__
#define __MAPP__command_displace__
#include "init.h"
#include "atoms.h"
namespace MAPP_NS {
    class Command_displace: protected InitPtrs
    {
    private:
        void move(int*,int,type0*);
        int id_n;
    protected:
    public:
        Command_displace(MAPP*,int,char**);
        ~Command_displace();
        
    };
}

#endif 
#endif
