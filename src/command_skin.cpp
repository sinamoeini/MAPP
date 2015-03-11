#include <stdlib.h>
#include "command_skin.h"
#include "atoms.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_skin::Command_skin(MAPP* mapp,int narg
,char** args):InitPtrs(mapp)
{
    if(narg!=2)
        error->abort("skin command needs only 1 argument");
    type0 s=atof(args[1]);
    if(s<=0.0)
        error->abort("skin cannot be equal or less than 0.0");

    atoms->chng_skin(s);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_skin::~Command_skin()
{

}

