#include "command_md.h"
#include "error.h"
#include "md_styles.h"
#include "global.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_md::Command_md
(int nargs,char** args)

{
    if(nargs<2)
        Error::abort("wrong command: %s"
        ,args[0]);
    
    int nh_xist=0;
    type0 t_step = 0.0,boltz=0.0;
    if(md)
    {
        nh_xist=1;
        t_step=md->dt;
        boltz=md->boltz;
        delete md;
    }
    #define MD_Style
    #define MDStyle(class_name,style_name)      \
    else if(strcmp(args[1],#style_name)==0)     \
    md=new class_name(nargs,args);

    if(0){}
    #include "md_styles.h"
    else
        Error::abort("wrong style of md: %s"
        ,args[1]);
    #undef MDStyle    
    #undef MD_Style
    if(nh_xist)
    {
        md->dt=t_step;
        md->boltz=boltz;
    }

    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_md::~Command_md()
{

}

