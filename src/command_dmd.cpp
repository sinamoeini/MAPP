#include "command_dmd.h"
#include "global.h"
#include "error.h"
#include "dmd_styles.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_dmd::Command_dmd
(int nargs,char** args)

{
    if(nargs<2)
        error->abort("wrong command: %s"
        ,args[0]);
    
    if(dmd)
    {
        delete dmd;
        dmd=NULL;
    }
    
    #define DMD_Style
    #define DMDStyle(class_name,style_name)\
    else if(strcmp(args[1],#style_name)==0)\
    {dmd=new class_name(nargs,args);}

    
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

