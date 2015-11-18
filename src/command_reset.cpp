#include "command_reset.h"
#include "error.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_reset::Command_reset(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(nargs!=1)
        error->abort("reset needs 0 argument");
    step_no=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_reset::~Command_reset()
{

}

