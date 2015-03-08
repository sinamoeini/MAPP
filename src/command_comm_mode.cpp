#include "command_comm_mode.h"
#include "atoms.h"
#include "error.h"
using namespace MAPP_NS;
using namespace std;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_comm_mode::Command_comm_mode(MAPP* mapp
,int narg,char** args):InitPtrs(mapp)
{
    if(narg!=2)
        error->abort("atom_mode needs 1 argument");
    
    if(strcmp(args[1],"0")==0)
        atoms->comm_mode=COMM_MODE_0;
    else if(strcmp(args[1],"1")==0)
        atoms->comm_mode=COMM_MODE_1;
    else if(strcmp(args[1],"2")==0)
        atoms->comm_mode=COMM_MODE_2;
    else if(strcmp(args[1],"3")==0)
        atoms->comm_mode=COMM_MODE_3;
    else if(strcmp(args[1],"4")==0)
        atoms->comm_mode=COMM_MODE_4;
    else if(strcmp(args[1],"5")==0)
        atoms->comm_mode=COMM_MODE_5;
    else
        error->abort("unknown keyword %s "
        "for atom_mode command",args[1]);
    
    if(atoms->my_p_no==0)
        printf("communication mode is set to %s\n",args[1]);
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_comm_mode::~Command_comm_mode()
{
    
}

