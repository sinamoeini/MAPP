#include "command_group.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "atom_types.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_group::Command_group(int nargs,char** args)
{
    if(nargs<3)
        error->abort("group command needs at least 2 arguments");
    
    if(strcmp(args[1],"add")==0)
    {
        if(nargs<5)
            error->abort("group add command needs at least 3 arguments");
        
        if(strcmp(args[3],"file")==0)
        {
            groups->add_grp_file(args[2],nargs-4,args+4);
        }
        else if(strcmp(args[3],"type")==0)
        {
            
            groups->add_grp_type(args[2],nargs-4,args+4);
        }
        else
            error->abort("unknown keyword: %s",args[1]);
    }
    else if(strcmp(args[1],"del")==0)
    {
        int iarg=2;
        while(iarg<nargs)
        {
            Group* grp=groups->find_grp(args[iarg]);
            groups->del_grp(grp);
            iarg++;
        }
    
    }
    else
        error->abort("unknown keyword: %s",args[1]);

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_group::~Command_group()
{
}



