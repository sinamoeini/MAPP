/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "ff.h"
#include "atom_types.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField::
ForceField(MAPP* mapp) : InitPtrs(mapp)
{
    gcmc_tag_enabled=false;
    gcmc_n_vars=gcmc_n_cutoff=1;
    cut_sz=0;
    ns_alloc=0;
    cut=NULL;
    cut_sk_sq=NULL;
    cut_sq=NULL;
    rsq_crd=NULL;
    
    if(atoms->tot_natms==0 || atom_types->no_types==0)
        error->abort("system configuration "
        "should be loaded before initiating ff");

    if(atoms->dimension!=3)
        error->abort("the dimension of the box for ff");
    
    cut_off_alloc();
    int dim=atoms->dimension;
    if(dim)
    {
        f=new Vec<type0>(atoms,mapp->x->dim);
        CREATE1D(nrgy_strss,dim*(dim+1)/2+1);
        CREATE1D(nrgy_strss_lcl,dim*(dim+1)/2+1);
        ns_alloc=1;
    }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField::~ForceField()
{
    cut_off_dealloc();
    if(ns_alloc)
    {
        delete [] nrgy_strss_lcl;
        delete [] nrgy_strss;
    }
}
/*--------------------------------------------
 allocate cutoff
 --------------------------------------------*/
void ForceField::cut_off_alloc()
{
    int no_types=atom_types->no_types;
    int cut_sz_=no_types;
    if(cut_sz_==cut_sz)
        return;
    cut_off_dealloc();
    cut_sz=cut_sz_;
    
    CREATE_2D(cut,no_types,no_types);
    
    
    CREATE_2D(cut_sq,no_types,no_types);
    
    
    CREATE_2D(cut_sk_sq,no_types,no_types);
    
    CREATE1D(rsq_crd,no_types);
}
/*--------------------------------------------
 deallocate cutoff
 --------------------------------------------*/
void ForceField::cut_off_dealloc()
{
    DEL_2D(cut);
    
    DEL_2D(cut_sq);
    
    
    DEL_2D(cut_sk_sq);
    
    
    delete [] rsq_crd;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
type0 ForceField::max_cut()
{
    type0 skin=atoms->get_skin();
    type0 tmp;
    type0 max_cut=0.0;
    for(int i=0;i<cut_sz;i++)
        for(int j=0;j<cut_sz;j++)
        {
            tmp=sqrt(cut_sq[i][j])+skin;
            cut_sk_sq[i][j]=tmp*tmp;
            max_cut=MAX(max_cut,tmp);
        }
    
    return max_cut;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField::force_calc_timer(bool flag)
{
    timer->start(FORCE_TIME_mode);
    force_calc(flag);
    if(flag)
    {
        type0** H=atoms->H;
        type0 vol=1.0;
        for(int idim=0;idim<atoms->dimension;idim++)
            vol*=H[idim][idim];
        for(int i=1;i<7;i++)
            nrgy_strss[i]/=vol;
    }
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
void ForceFieldMD::pre_xchng_energy_timer(GCMC* gcmc)
{
    timer->start(FORCE_TIME_mode);
    pre_xchng_energy(gcmc);
    timer->stop(FORCE_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldMD::xchng_energy_timer(GCMC* gcmc)
{
    timer->start(FORCE_TIME_mode);
    type0 en=xchng_energy(gcmc);
    timer->stop(FORCE_TIME_mode);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldMD::post_xchng_energy_timer(GCMC* gcmc)
{
    timer->start(FORCE_TIME_mode);
    post_xchng_energy(gcmc);
    timer->stop(FORCE_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceFieldDMD::ddc_norm_timer()
{
    timer->start(FORCE_TIME_mode);
    type0 en=ddc_norm();
    timer->stop(FORCE_TIME_mode);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::ddc_timer(type0* ddc_)
{
    timer->start(FORCE_TIME_mode);
    ddc(ddc_);
    timer->stop(FORCE_TIME_mode);

}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceFieldDMD::dc_timer()
{
    timer->start(FORCE_TIME_mode);
    dc();
    timer->stop(FORCE_TIME_mode);
}

