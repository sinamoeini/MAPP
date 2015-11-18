#include "command_mode.h"
#include "atoms.h"
#include "error.h"
#include "ff.h"
#include "neighbor_md.h"
#include "neighbor_dmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_mode::Command_mode(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(mapp->no_commands!=0)
        error->abort("mode command should be the first command");
    
    if(nargs!=2)
        error->abort("wrong command: %s",args[0]);
    int new_mode;
    if(strcmp(args[1],"md")==0)
    {
        new_mode=MD_mode;
        delete mapp->neighbor;
        mapp->neighbor=new Neighbor_md(mapp);
        if(atoms->my_p==0)
            printf("mode set to md\n");
    }
    else if(strcmp(args[1],"dmd")==0)
    {
        new_mode=DMD_mode;
        delete mapp->neighbor;
        mapp->neighbor=new Neighbor_dmd(mapp);
        delete mapp->forcefield;
        if(atoms->my_p==0)
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

