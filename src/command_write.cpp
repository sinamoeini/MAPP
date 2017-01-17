#include "command_write.h"
#include "error.h"
#include "write_styles.h"
#include "global.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_write::Command_write
(int nargs,char** args)

{
    if(nargs<2)
        Error::abort("wrong command: %s",args[0]);
    
    if(write)
        delete write;

    #define Write_Style
    #define WriteStyle(class_name,style_name)   \
    else if(strcmp(args[1],#style_name)==0)     \
    write= new class_name(nargs,args);
    
    if(0){}
    #include "write_styles.h"
    else
        Error::abort("wrong style of write:"
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

