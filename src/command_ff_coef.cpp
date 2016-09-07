#include "command_ff_coef.h"
#include "error.h"
#include "ff.h"
#include "init.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_ff_coef::Command_ff_coef(int nargs,char** args)

{
    if(forcefield==NULL)
        error->abort("cannot add the coefficients"
        " before the forcefield is initiated");
    forcefield->coef(nargs,args);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_ff_coef::~Command_ff_coef()
{

}

