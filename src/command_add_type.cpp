#include <stdlib.h>
#include "command_add_type.h"
#include "atom_types.h"
#include "MAPP.h"
#include "global.h"
#include "error.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_add_type::Command_add_type(int nargs,char** args)

{

    if(nargs!=3)
        Error::abort("add_type command should have 2 arguments");
    
    type0 mass=atof(args[2]);
    
    if(mass<=0.0)
        Error::abort("mass should be greater than 0.0");
    
    int itype=mapp->atom_types->find_type_exist(args[1]);
    if(itype!=-1)
        mapp->atom_types->mass[itype]=mass;
    else
    {
        mapp->atom_types->add_type(mass,args[1]);
        if(forcefield)
            Error::warning("ff && ff_coef commands should be redeclared");
    }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_add_type::~Command_add_type()
{

}

