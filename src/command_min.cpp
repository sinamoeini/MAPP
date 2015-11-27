#include "command_min.h"
#include "error.h"
#include "min_styles.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_min::Command_min
(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(nargs<2)
        error->abort("wrong command: %s"
        ,args[0]);
    Min*& min=mapp->min;
    if(min!=NULL)
        delete min;
    
    #define Min_Style
    #define MinStyle(class_name,style_name)     \
    else if(strcmp(args[1],#style_name)==0)     \
        {if(min!=NULL)delete min;               \
        min= new class_name(mapp,nargs,args);   \
        min->init();min->run();min->fin();      \
        delete min;min=NULL;}
    
    if(0){}
    #include "min_styles.h"
    else
        error->abort("wrong style of minimization"
        ": %s",args[1]);
    
    #undef MinStyle
    #undef Min_Style
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_min::~Command_min()
{

}

