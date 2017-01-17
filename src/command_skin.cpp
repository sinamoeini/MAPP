#include <stdlib.h>
#include "command_skin.h"
#include "comm.h"
#include "error.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_skin::Command_skin(int nargs
,char** args)
{
    if(nargs!=2)
        Error::abort("skin command needs only 1 argument");
    type0 s=atof(args[1]);
    if(s<=0.0)
        Error::abort("skin cannot be equal or less than 0.0");

    comm->skin=s;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_skin::~Command_skin()
{

}

