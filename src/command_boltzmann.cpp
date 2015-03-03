#include "command_boltzmann.h"
#include "error.h"
#include "md.h"
using namespace MAPP_NS;
using namespace std;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_boltzmann::Command_boltzmann(MAPP* mapp,int narg,char** args)
:InitPtrs(mapp)
{
    if(mapp->md==NULL)
        error->abort("before assigning the boltzmann"
        " constant, md ensemble should be initialized");
    
    if(narg!=2)
        error->abort("boltzmann command should have only 1 arguement");
    
    type0 boltz=atof(args[1]);
    
    if(boltz<=0.0)
        error->abort("boltzmann should be greater than 0.0");
    mapp->md->boltz=boltz;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_boltzmann::~Command_boltzmann()
{

}

