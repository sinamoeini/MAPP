#include "command_atom_mode.h"
#include "atoms.h"
#include "error.h"
using namespace MAPP_NS;
using namespace std;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_atom_mode::Command_atom_mode(MAPP* mapp
,int narg,char** args):InitPtrs(mapp)
{
    if(narg!=2)
        error->abort("atom_mode needs 1 argument");
    
    if(strcmp(args[1],"basic")==0)
        atoms->atom_mode=BASIC_mode;
    else if(strcmp(args[1],"type")==0)
    {
        if(mapp->mode!=MD_mode)
        {
            error->abort("keyword %s "
            "for atom_mode command can "
            "only be used when mode is "
            "set to md",args[1]);
        }
        atoms->atom_mode=TYPE_mode;
    }
    else
        error->abort("unknown keyword %s "
        "for atom_mode command",args[1]);
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_atom_mode::~Command_atom_mode()
{
    
}

