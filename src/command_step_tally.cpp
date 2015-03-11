#include <stdlib.h>
#include "command_step_tally.h"
#include "error.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_step_tally::Command_step_tally(MAPP* mapp,
int narg,char** args):InitPtrs(mapp)
{
    if(narg!=2)
        error->abort("step_tally command "
        "needs 2 arguements");
    
    step_tally=atoi(args[1]);
    
    if(step_tally<=0)
        error->abort("step_tally cannot "
        "be equal or less than zero");
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_step_tally::~Command_step_tally()
{

}

