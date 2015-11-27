#include "command_dmd_min.h"
#include "error.h"
#include "dmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_dmd_min::Command_dmd_min
(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(nargs<2)
        error->abort("wrong command: %s"
        ,args[0]);
    if(mapp->dmd==NULL)
        error->abort("cannot add the min"
        " before the dmd is initiated");
    mapp->dmd->dmd_min(nargs,args);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_dmd_min::~Command_dmd_min()
{
}

