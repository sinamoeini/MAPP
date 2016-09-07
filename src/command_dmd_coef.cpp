#include "command_dmd_coef.h"
#include "init.h"
#include "error.h"
#include "dmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_dmd_coef::Command_dmd_coef
(int nargs,char** args)

{
    if(dmd==NULL)
        error->abort("cannot add the coefficients"
        " before the dmd is initiated");
    dmd->coef(nargs,args);

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_dmd_coef::~Command_dmd_coef()
{
}

