#include "command_precision.h"
#include "error.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_precision::Command_precision(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(nargs!=2)
        error->abort("precision needs 1 argument");
    int p=atoi(args[1]);
    if(p<=0)
        error->abort("precision must be greater than 0");
    mapp->precision=p;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_precision::~Command_precision()
{

}

