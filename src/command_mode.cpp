#include "command_mode.h"
#include "atoms.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_mode::Command_mode(MAPP* mapp,int narg,char** args)
:InitPtrs(mapp)
{
    if(mapp->no_commands!=1)
        error->abort("mode command should be the first command");
    
    if(narg!=2)
        error->abort("wrong command: %s",args[0]);
    int new_mode;
    if(strcmp(args[1],"md")==0)
    {
        new_mode=MD_mode;
        if(atoms->my_p_no==0)
            printf("mode set to md\n");
    }
    else if(strcmp(args[1],"dmd")==0)
    {
        new_mode=DMD_mode;
        if(atoms->my_p_no==0)
            printf("mode set to dmd\n");
    }
    else
        error->abort("unknown mode: %s",args[1]);
    if(new_mode==mapp->mode)
        return;
    mapp->mode=new_mode;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_mode::~Command_mode()
{

}

