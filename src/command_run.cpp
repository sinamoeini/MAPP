#include <stdlib.h>
#include "command_run.h"
#include "error.h"
#include "md.h"
#include "dmd.h"
#include "global.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_run::Command_run(int nargs,char** args)

{
    if(mode==MD_mode)
    {
        if(md==NULL)
            error->abort("before run, ensemble should be initialized");
        
        if(nargs!=2)
            error->abort("run should have 1 argument");
        
        int steps=atoi(args[1]);
        if(steps<0)
            error->abort("run should be "
                         "greater than 0");
        
        md->init();
        md->run(steps);
        md->fin();
    }
    else if(mode==DMD_mode)
    {
        if(dmd==NULL)
            error->abort("before run, dmd should be initialized");
        if(nargs!=2)
            error->abort("run should have 1 argument");
        
        dmd->init();
        dmd->run(atof(args[1]));
        dmd->fin();
    }

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_run::~Command_run()
{

}

