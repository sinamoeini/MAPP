#include "command_precision.h"
#include "global.h"
#include "error.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_precision::Command_precision(int nargs,char** args)

{
    if(nargs!=2)
        Error::abort("precision needs 1 argument");
    int p=atoi(args[1]);
    if(p<=0)
        Error::abort("precision must be greater than 0");
    precision=p;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_precision::~Command_precision()
{

}

