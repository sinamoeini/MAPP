#include "command_clock_coef.h"
#include "error.h"
#include "clock.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_clock_coef::Command_clock_coef
(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(mapp->clock==NULL)
        error->abort("cannot add the coefficients"
        " before the forcefield is initiated");
    mapp->clock->coef(nargs,args);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_clock_coef::~Command_clock_coef()
{

}

