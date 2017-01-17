#include <stdlib.h>
#include "command_step_tally.h"
#include "global.h"
#include "error.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_step_tally::Command_step_tally(
int nargs,char** args)
{
    if(nargs!=2)
        Error::abort("step_tally command "
        "needs 2arguments");
    
    step_tally=atoi(args[1]);
    
    if(step_tally<=0)
        Error::abort("step_tally cannot "
        "be equal or less than zero");
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_step_tally::~Command_step_tally()
{

}

