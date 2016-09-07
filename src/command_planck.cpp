#include <stdlib.h>
#include "command_planck.h"
#include "md.h"
#include "error.h"
#include "init.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_planck::Command_planck(int nargs,char** args)

{
    if(md==NULL)
        error->abort("before assigning the planck"
        " constant, md ensemble should be initialized");
    
    if(nargs!=2)
        error->abort("planck command should have only 1 argument");
    
    type0 planck=atof(args[1]);
    
    if(planck<=0.0)
        error->abort("planck should be greater than 0.0");
    md->hplanck=planck;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_planck::~Command_planck()
{

}

