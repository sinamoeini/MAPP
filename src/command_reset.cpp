#include "command_reset.h"
#include "error.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_reset::Command_reset(MAPP* mapp,int narg,char** args)
:InitPtrs(mapp)
{
    if(narg!=1)
        error->abort("reset needs 0 argument");
    step_no=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_reset::~Command_reset()
{

}

