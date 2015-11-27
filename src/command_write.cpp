#include "command_write.h"
#include "error.h"
#include "write_styles.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_write::Command_write
(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(nargs<2)
        error->abort("wrong command: %s",args[0]);
    
    if(write!=NULL)
        delete write;

    #define Write_Style
    #define WriteStyle(class_name,style_name)   \
    else if(strcmp(args[1],#style_name)==0)     \
    write= new class_name(mapp,nargs,args);
    
    if(0){}
    #include "write_styles.h"
    else
        error->abort("wrong style of write:"
            " %s",args[1]);
    #undef Write_Style
    #undef WriteStyle


    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_write::~Command_write()
{

}

