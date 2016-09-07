#include <stdlib.h>
#include "command_add_type.h"
#include "atom_types.h"
#include "init.h"
#include "error.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Command_add_type::Command_add_type(int nargs,char** args)

{

    if(nargs!=3)
        error->abort("add_type command should have 2 arguments");
    
    type0 mass=atof(args[2]);
    
    if(mass<=0.0)
        error->abort("mass should be greater than 0.0");
    
    int itype=atom_types->find_type_exist(args[1]);
    if(itype!=-1)
        atom_types->mass[itype]=mass;
    else
    {
        atom_types->add_type(mass,args[1]);
        if(forcefield!=NULL)
            error->warning("ff && ff_coef commands should be redeclared");
    }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Command_add_type::~Command_add_type()
{

}

