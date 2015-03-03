#include "command_ff_coef.h"
#include "ff.h"
#include "error.h"
using namespace MAPP_NS;
using namespace std;
/*--------------------------------------------
 you might need to include ff_styles.h
 
 
 constructor
 --------------------------------------------*/
Command_ff_coef::Command_ff_coef(MAPP* mapp,int narg,char** args)
:InitPtrs(mapp)
{
    if(forcefield==NULL)
        error->abort("cannot add the coefficients"
        " before the forcefield is initiated");
    forcefield->coef(narg,args);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_ff_coef::~Command_ff_coef()
{

}

