#include <stdlib.h>
#include "command_run.h"
#include "error.h"
#include "md.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_run::Command_run(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(mapp->md==NULL)
        error->abort("before run, ensemble should be initialized");
    
    if(nargs!=2)
        error->abort("run should have 1argument");
    
    int steps=atoi(args[1]);
    if(steps<=0)
        error->abort("run should be "
                     "greater than 0");
    
    mapp->md->init();
    mapp->md->run(steps);
    mapp->md->fin();
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_run::~Command_run()
{

}

