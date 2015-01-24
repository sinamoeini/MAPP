/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "ff.h"
#include "atom_types.h"
using namespace MAPP_NS;
using namespace std;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField::
ForceField(MAPP* mapp) : InitPtrs(mapp)
{
    if(atoms->tot_natms==0 ||
       atom_types->no_types==0)
        error->abort("system configuration "
        "should be loaded before initiating ff");

    if(atoms->dimension!=3)
        error->abort("the dimension of the box for ff");
        
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField::~ForceField()
{

}

