#include "command_dmd_coef.h"
#include "error.h"
#include "dmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_dmd_coef::Command_dmd_coef
(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(mapp->dmd==NULL)
        error->abort("cannot add the coefficients"
        " before the dmd is initiated");
    mapp->dmd->coef(nargs,args);

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_dmd_coef::~Command_dmd_coef()
{
}

