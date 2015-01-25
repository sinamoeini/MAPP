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
    cut_alloc=0;
    ns_alloc=0;
    if(atoms->tot_natms==0 ||
       atom_types->no_types==0)
        error->abort("system configuration "
        "should be loaded before initiating ff");

    if(atoms->dimension!=3)
        error->abort("the dimension of the box for ff");
    
    int no_types=atom_types->no_types;
    if(no_types)
    {
        CREATE1D(cut_sk_sq,no_types*(no_types+1)/2);
        CREATE1D(cut_sq,no_types*(no_types+1)/2);
        cut_alloc=1;
    }
    int dim=atoms->dimension;
    if(dim)
    {
        CREATE1D(nrgy_strss,dim*(dim+1)/2+1);
        ns_alloc=1;
    }

    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField::~ForceField()
{
    if(cut_alloc)
    {
        delete [] cut_sk_sq;
        delete [] cut_sq;        
    }
    if(ns_alloc)
        delete [] nrgy_strss;
}

