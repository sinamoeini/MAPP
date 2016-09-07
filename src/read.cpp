/*--------------------------------------------
 Created by Sina on 07/23/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "read.h"
#include "ff.h"
#include "atom_types.h"
#include "error.h"
#include "atoms.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Read::Read()
{
    
    if(forcefield!=NULL)
    {
        delete forcefield;
        forcefield=NULL;
        error->warning("removing the old ff while "
        "loading new configuration, you should "
        "reinitiate ff & ff_coef");
    }
    
    atoms->restart();
    delete atom_types;
    atom_types=new AtomTypes();
     
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Read::~Read()
{
    
}

