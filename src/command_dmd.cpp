#include "command_dmd.h"
#include "error.h"
#include "dmd_styles.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_dmd::Command_dmd
(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(nargs<2)
        error->abort("wrong command: %s"
        ,args[0]);
    
    if(mapp->dmd!=NULL)
    {
        delete mapp->dmd;
        mapp->dmd=NULL;
    }
    
    #define DMD_Style
    #define DMDStyle(class_name,style_name)\
    else if(strcmp(args[1],#style_name)==0)\
    {mapp->dmd=new class_name(mapp,nargs,args);}

    
    if(0){}
    #include "dmd_styles.h"
    else
        error->abort("wrong style of dmd"
                     ": %s",args[1]);
    #undef DMDStyle    
    #undef DMD_Style

    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_dmd::~Command_dmd()
{

}

