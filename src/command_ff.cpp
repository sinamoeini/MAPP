#include "command_ff.h"
#include "global.h"
#include "error.h"
#include "ff_styles.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_ff::Command_ff
(int nargs,char** args)

{
    if(nargs!=2)
        Error::abort("wrong command: %s",args[0]);
    
    if(forcefield)
        delete forcefield;
    
    #define FF_Style
    #define FFStyle(class_name,style_name)     \
    else if(strcmp(args[1],#style_name)==0)    \
    forcefield=new class_name();
    
    //different forcefileds
    if(0){}
    #include "ff_styles.h"
    else
        Error::abort("unknown forcefield: %s"
                     ,args[1]);
    #undef FFStyle
    #undef FF_Style


    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_ff::~Command_ff()
{

}

