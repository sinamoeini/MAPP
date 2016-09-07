#include <stdlib.h>
#include "command_boltzmann.h"
#include "md.h"
#include "error.h"
#include "init.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_boltzmann::Command_boltzmann(int nargs,char** args)

{
    if(md==NULL)
        error->abort("before assigning the boltzmann"
        " constant, md ensemble should be initialized");
    
    if(nargs!=2)
        error->abort("boltzmann command should have only 1argument");
    
    type0 boltz=atof(args[1]);
    
    if(boltz<=0.0)
        error->abort("boltzmann should be greater than 0.0");
    md->boltz=boltz;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_boltzmann::~Command_boltzmann()
{

}

