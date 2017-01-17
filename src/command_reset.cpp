#include "command_reset.h"
#include "global.h"
#include "error.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_reset::Command_reset(int nargs,char** args)

{
    if(nargs!=1)
        Error::abort("reset needs 0 argument");
    step_no=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_reset::~Command_reset()
{

}

