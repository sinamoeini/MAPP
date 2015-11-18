#include "command_ff_coef.h"
#include "error.h"
#include "ff.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_ff_coef::Command_ff_coef(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
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

