/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "ff.h"
#include "atom_types.h"
#include "timer.h"
using namespace MAPP_NS;
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
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField::force_calc_timer(int chk,type0* enst)
{
    timer->start(FORCE_TIME_mode);
    force_calc(chk,enst);
    timer->stop(FORCE_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceField::energy_calc_timer()
{
    type0 en;
    timer->start(FORCE_TIME_mode);
    en=energy_calc();
    timer->stop(FORCE_TIME_mode);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceField::g_calc_timer(int chk
,type0 alpha,type0* a,type0* g)
{
    type0 en;
    timer->start(FORCE_TIME_mode);
    en=g_calc_timer(chk,alpha,a,g);
    timer->stop(FORCE_TIME_mode);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField::c_d_calc_timer()
{
    timer->start(FORCE_TIME_mode);
    c_d_calc();
    timer->stop(FORCE_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField::create_2nd_neigh_lst_timer()
{
    timer->start(FORCE_TIME_mode);
    create_2nd_neigh_lst();
    timer->stop(FORCE_TIME_mode);
}


