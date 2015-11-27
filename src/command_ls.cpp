#include "command_ls.h"
#include "error.h"
#include "ls_styles.h"
#include "min.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_ls::Command_ls
(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(nargs<2)
        error->abort("wrong command: %s",args[0]);
    
    if(ls!=NULL)
        delete ls;
    
    #define LS_Style
    #define LSStyle(class_name,style_name)     \
    else if(strcmp(args[1],#style_name)==0)    \
    ls=new class_name<Min>(mapp,nargs,args);
    
    if(0){}
    #include "ls_styles.h"
    else
        error->abort("unknown line search: %s"
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

