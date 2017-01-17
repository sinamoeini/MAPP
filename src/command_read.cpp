#include "command_read.h"
#include "global.h"
#include "error.h"
#include "read_styles.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_read::Command_read
(int nargs,char** args)

{
    if(nargs<2)
        Error::abort("wrong command: %s",args[0]);
    
    Read* read;
    #define Read_Style
    #define ReadStyle(class_name,style_name)    \
    else if(strcmp(args[1],#style_name)==0)     \
    read= new class_name(nargs,args);
    
    if(0){}
    #include "read_styles.h"
    else
        Error::abort("wrong style of read: %s",args[1]);
    #undef ReadStyle
    #undef Read_Style
    delete read;


    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_read::~Command_read()
{

}

