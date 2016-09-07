#include "command_min.h"
#include "error.h"
#include "min_styles.h"
#include "MAPP.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_min::Command_min
(int nargs,char** args)

{
    if(nargs<2)
        error->abort("wrong command: %s"
        ,args[0]);
    
    /*
    Min*& min=mapp->min;
    if(min!=NULL)
        delete min;
    */
    Min* min;
    
    #define Min_Style
    #define MinStyle(class_name,style_name)\
    else if(strcmp(args[1],#style_name)==0)\
        min=new class_name(nargs,args);
    
    if(0){}
    #include "min_styles.h"
    else
        error->abort("wrong style of minimization"
        ": %s",args[1]);
    
    #undef Min_Style
    #undef MinStyle
    
    min->init();
    min->run();
    min->fin();
    delete min;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_min::~Command_min()
{

}

