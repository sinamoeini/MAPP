#include "command_dmd_min.h"
#include "error.h"
#include "global.h"
#include "dmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_dmd_min::Command_dmd_min
(int nargs,char** args)

{
    if(nargs<2)
        Error::abort("wrong command: %s"
        ,args[0]);
    if(dmd==NULL)
        Error::abort("cannot add the min"
        " before the dmd is initiated");
    dmd->dmd_min(nargs,args);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_dmd_min::~Command_dmd_min()
{
}

