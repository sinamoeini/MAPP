#include "command_del_write.h"
#include "error.h"
#include "write.h"
#include "init.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_del_write::Command_del_write(int nargs,char** args)

{
    if(nargs!=1)
        error->abort("del_write needs 0 argument");
    
    if(write!=NULL)
        delete write;
    write=NULL;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_del_write::~Command_del_write()
{

}

