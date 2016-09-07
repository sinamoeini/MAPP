#include <stdlib.h>
#include "command_temp.h"
#include "init.h"
#include "error.h"
#include "ff.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_temp::Command_temp(int nargs,char** args)

{
    
    if(mode!=DMD_mode)
        error->abort("temp works only in dmd mode");
            
    if(nargs!=2)
        error->abort("temp should have 1 argument");
    
    type0 T=atof(args[1]);
    
    if(T<=0.0)
        error->abort("temperature must be greater than 0.0");
    
    ForceFieldDMD* forcefield_dmd=static_cast<ForceFieldDMD*>(forcefield);
    
    forcefield_dmd->set_temp(T);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_temp::~Command_temp()
{
}

