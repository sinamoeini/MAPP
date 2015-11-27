#include <stdlib.h>
#include "command_run.h"
#include "error.h"
#include "md.h"
#include "dmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_run::Command_run(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(mapp->mode==MD_mode)
    {
        if(mapp->md==NULL)
            error->abort("before run, ensemble should be initialized");
        
        if(nargs!=2)
            error->abort("run should have 1 argument");
        
        int steps=atoi(args[1]);
        if(steps<=0)
            error->abort("run should be "
                         "greater than 0");
        
        mapp->md->init();
        mapp->md->run(steps);
        mapp->md->fin();
    }
    else if(mapp->mode==DMD_mode)
    {
        if(mapp->dmd==NULL)
            error->abort("before run, dmd should be initialized");
        if(nargs!=2)
            error->abort("run should have 1 argument");
        
        mapp->dmd->init();
        mapp->dmd->run(atof(args[1]));
        mapp->dmd->fin();
        
    }

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_run::~Command_run()
{

}

