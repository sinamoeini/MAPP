#include <stdlib.h>
#include "command_time_step.h"
#include "md.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_time_step::Command_time_step(MAPP* mapp,
int narg,char** args):InitPtrs(mapp)
{
    if(mapp->md==NULL)
        error->abort("before adjusting the "
        "time_step, ensemble should be initialized");
    
    if(narg!=2)
        error->abort("time_step should have 1 arguement");
    
    type0 dt=atof(args[1]);
    if(dt<=0.0)
        error->abort("time_step should be"
        " greater than 0.0");
    
    mapp->md->dt=dt;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_time_step::~Command_time_step()
{

}

