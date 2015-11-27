#include "command_md.h"
#include "error.h"
#include "md_styles.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_md::Command_md
(MAPP* mapp,int nargs,char** args)
:InitPtrs(mapp)
{
    if(nargs<2)
        error->abort("wrong command: %s"
        ,args[0]);
    
    int nh_xist=0;
    type0 t_step = 0.0,boltz=0.0;
    if(mapp->md!=NULL)
    {
        nh_xist=1;
        t_step=mapp->md->dt;
        boltz=mapp->md->boltz;
        delete mapp->md;
    }
    #define MD_Style
    #define MDStyle(class_name,style_name)      \
    else if(strcmp(args[1],#style_name)==0)     \
    mapp->md= new class_name(mapp,nargs,args);

    if(0){}
    #include "md_styles.h"
    else
        error->abort("wrong style of md: %s"
        ,args[1]);
    #undef MDStyle    
    #undef MD_Style
    if(nh_xist)
    {
        mapp->md->dt=t_step;
        mapp->md->boltz=boltz;
    }

    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_md::~Command_md()
{

}

