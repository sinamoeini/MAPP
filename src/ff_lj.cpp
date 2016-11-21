/*--------------------------------------------
 Created by Sina on 07/15/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <stdlib.h>
#include "ff_lj.h"
#include "neighbor.h"
#include "atom_types.h"
#include "memory.h"
#include "error.h"
#include "atoms.h"
#include "cmd.h"
#include "MAPP.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_lj::
ForceField_lj():ForceFieldMD()
{
    if(mode!=MD_mode)
        error->abort("ff lj works only "
        "for md mode");
    
    shift=0;
    no_types=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_lj::~ForceField_lj()
{
    deallocate();
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ForceField_lj::coef(int nargs,char** args)
{
    
    if(nargs!=2)
        error->abort("wrong coeff command "
        "for lj Force Field");
    
    cut_off_alloc();
    allocate();
    
    read_file(args[1]);
}
/*--------------------------------------------
 allocation
 --------------------------------------------*/
void ForceField_lj::allocate()
{
    if(no_types==atom_types->no_types)
        return;
    deallocate();
    no_types=atom_types->no_types;
    
    sigma=new type0*[no_types];
    *sigma=new type0[no_types*no_types];
    epsilon=new type0*[no_types];
    *epsilon=new type0[no_types*no_types];
    offset=new type0*[no_types];
    *offset=new type0[no_types*no_types];
    
    for(int i=1;i<no_types;i++)
    {
        sigma[i]=sigma[i-1]+no_types;
        epsilon[i]=epsilon[i-1]+no_types;
        offset[i]=offset[i-1]+no_types;
    }
}
/*--------------------------------------------
 allocation
 --------------------------------------------*/
void ForceField_lj::deallocate()
{
    if(!no_types)
        return;
    
    delete [] *sigma;
    delete [] sigma;
    delete [] *epsilon;
    delete [] epsilon;
    delete [] *offset;
    delete [] offset;
    no_types=0;
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_lj::init()
{
    neighbor->pair_wise=true;
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceField_lj::fin()
{
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceField_lj::init_xchng()
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceField_lj::fin_xchng()
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceField_lj::pre_xchng_energy(GCMC*)
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 xchng energy
 --------------------------------------------*/
type0 ForceField_lj::xchng_energy(GCMC*)
{
    error->abort("exchange has not been set for this forcefield");
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceField_lj::post_xchng_energy(GCMC*)
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceField_lj::
force_calc(bool st_clc)
{


    type0* x=atoms->x->begin();
    type0* fvec=f->begin();
    md_type* type=mapp->type->begin();
    
    int natms=atoms->natms;
    int iatm,jatm;
    int itype,jtype,icomp,jcomp;
    type0 dx0,dx1,dx2,rsq,sig,eps;
    type0 sig2,sig6,sig12,fpair,en;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    nrgy_strss_lcl[0]=0.0;
    if(st_clc)
        for(int i=1;i<7;i++)
            nrgy_strss_lcl[i]=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            icomp=3*iatm;
            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            if(rsq<cut_sq[itype][jtype])
            {
                sig=sigma[itype][jtype];
                eps=epsilon[itype][jtype];
                sig2=sig*sig/rsq;
                sig6=sig2*sig2*sig2;
                sig12=sig6*sig6;
                
                fpair=24.0*eps*(2.0*sig12-sig6)/rsq;
                en=4.0*eps*(sig12-sig6)
                +offset[itype][jtype];
                
                fvec[icomp]+=fpair*dx0;
                fvec[icomp+1]+=fpair*dx1;
                fvec[icomp+2]+=fpair*dx2;
                if(jatm<natms)
                {
                    fvec[jcomp]-=fpair*dx0;
                    fvec[jcomp+1]-=fpair*dx1;
                    fvec[jcomp+2]-=fpair*dx2;
                }
                
                
                if(jatm>=natms)
                {
                    fpair*=0.5;
                    en*=0.5;
                }
                nrgy_strss_lcl[0]+=en;
                
                if(st_clc)
                {
                    nrgy_strss_lcl[1]-=fpair*dx0*dx0;
                    nrgy_strss_lcl[2]-=fpair*dx1*dx1;
                    nrgy_strss_lcl[3]-=fpair*dx2*dx2;
                    nrgy_strss_lcl[4]-=fpair*dx1*dx2;
                    nrgy_strss_lcl[5]-=fpair*dx2*dx0;
                    nrgy_strss_lcl[6]-=fpair*dx0*dx1;
                }
            }
        }
    }

    
    if(st_clc)
    {
        for(int i=0;i<7;i++)
            nrgy_strss[i]=0.0;
        
        MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,7,MPI_TYPE0,MPI_SUM,world);
    }
    else
    {
        nrgy_strss[0]=0.0;
        MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,1,MPI_TYPE0,MPI_SUM,world);
    }
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
type0 ForceField_lj::energy_calc()
{
    type0* x=atoms->x->begin();
    md_type* type=mapp->type->begin();
    
    int natms=atoms->natms;
    int iatm,jatm;
    
    int itype,jtype,icomp,jcomp;
    type0 dx0,dx1,dx2,rsq;
    type0 eps,sig,sig2,sig6,sig12;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    type0 en=0.0;
    type0 en_tot=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            icomp=3*iatm;
            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            if(rsq<cut_sq[itype][jtype])
            {
                sig=sigma[itype][jtype];
                eps=epsilon[itype][jtype];
                sig2=sig*sig/rsq;
                sig6=sig2*sig2*sig2;
                sig12=sig6*sig6;
                
                if(jatm<natms)
                    en+=4.0*eps*(sig12-sig6)
                    +offset[itype][jtype];
                else
                    en+=2.0*eps*(sig12-sig6)
                    +offset[itype][jtype]*0.5;
                
            }
        }
    }
    
    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    return en_tot;
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_lj::read_file(char* file_name)
{
    
    FileReader fr;
    
    fr.add_2D("sigma",sigma);
    fr.symmetric();
    fr.add_vlog()=VLogics("gt",0.0);
    
    
    fr.add_2D("epsilon",epsilon);
    fr.symmetric();
    fr.add_vlog()=VLogics("gt",0.0);
    
    fr.add_2D("r_c",cut);
    fr.symmetric();
    fr.add_vlog()=VLogics("gt",0.0);
    

    

    fr.read_file(file_name);
    
    for(int itype=0;itype<no_types;itype++)
        for(int jtype=0;jtype<no_types;jtype++)
            cut_sq[itype][jtype]=cut[itype][jtype]*cut[itype][jtype];
    
    if(shift)
    {
        type0 sig2,sig6,sig12;
        for(int itype=0;itype<no_types;itype++)
            for(int jtype=0;jtype<no_types;jtype++)
            {
                sig2=sigma[itype][jtype]*sigma[itype][jtype]/cut_sq[itype][jtype];
                sig6=sig2*sig2*sig2;
                sig12=sig6*sig6;
                offset[itype][jtype]=-4.0*epsilon[itype][jtype]*(sig12-sig6);
            }
    }
    else
    {
        for(int itype=0;itype<no_types;itype++)
            for(int jtype=0;jtype<no_types;jtype++)
                offset[itype][jtype]=0.0;
    }
}



