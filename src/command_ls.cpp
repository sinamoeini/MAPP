#include "command_ls.h"
#include "error.h"
#include "ls_styles.h"
#include "min.h"
using namespace MAPP_NS;
LineSearch<Min>* Min::ls;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_ls::Command_ls(int nargs,char** args)
{
    if(nargs<2)
        Error::abort("wrong command: %s",args[0]);
    
    if(Min::ls)
        delete Min::ls;
    
    #define LS_Style
    #define LSStyle(class_name,style_name)     \
    else if(strcmp(args[1],#style_name)==0)    \
    Min::ls=new class_name<Min>(nargs,args);
    
    if(0){}
    #include "ls_styles.h"
    else
        Error::abort("unknown line search: %s"
                     ,args[1]);
    #undef LSStyle
    #undef LS_Style
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_ls::~Command_ls()
{

}

