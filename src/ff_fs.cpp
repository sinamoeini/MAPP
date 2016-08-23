/*--------------------------------------------
 Created by Sina on 01/29/14.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <stdlib.h>
#include "ff_fs.h"
#include "neighbor.h"
#include "atom_types.h"
#include "error.h"
#include "memory.h"
#include "cmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 Finnis-Sinclair (FS) potential
 ref:
 T. T. Lau, C. J. Forst, X. Lin, J. D. Gale,
 S. Yip, & K. J. Van Vliet
 Many-Body Potential for Point Defect Clusters
 in Fe-C Alloys
 Phys. Rev. Lett. Vol. 98, pp. 215501, 2007
 --------------------------------------------*/

/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_fs::
ForceField_fs(MAPP* mapp):ForceFieldMD(mapp)
{
    if(mapp->mode!=MD_mode)
        error->abort("ff fs works only "
        "for md mode");
    max_pairs=0;
    no_types=0;
    
    mat_t_1=NULL;
    mat_t_2=NULL;
    cut_phi=NULL;
    cut_rho=NULL;
    mat_k_1=NULL;
    mat_k_2=NULL;
    mat_k_3=NULL;
    mat_A=NULL;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_fs::~ForceField_fs()
{
    deallocate();
}
/*--------------------------------------------
 read the force field file 
 --------------------------------------------*/
void ForceField_fs::coef(int nargs,char** args)
{
    if(nargs!=2)
        error->abort("wrong coeff command "
        "for Finnis-Sinclair Force Field");
    
    cut_off_alloc();
    allocate();

    read_file(args[1]);
    
    for(int i=0;i<no_types;i++)
        for(int j=0;j<no_types;j++)
        {
            cut[i][j]=MAX(cut_phi[i][j],cut_rho[i][j]);
            cut_sq[i][j]=cut[i][j]*cut[i][j];
        }
}
/*--------------------------------------------
 allocation
 --------------------------------------------*/
void ForceField_fs::allocate()
{
    if(no_types==atom_types->no_types)
        return;
    deallocate();
    no_types=atom_types->no_types;
    
    CREATE_2D(mat_t_1,no_types,no_types);
    CREATE_2D(mat_t_2,no_types,no_types);
    CREATE_2D(cut_phi,no_types,no_types);
    CREATE_2D(cut_rho,no_types,no_types);
    CREATE_2D(mat_k_1,no_types,no_types);
    CREATE_2D(mat_k_2,no_types,no_types);
    CREATE_2D(mat_k_3,no_types,no_types);
    CREATE1D(mat_A,no_types);
    
}
/*--------------------------------------------
 allocation
 --------------------------------------------*/
void ForceField_fs::deallocate()
{
    if(!no_types)
        return;
    
    DEL_2D(mat_t_1);
    DEL_2D(mat_t_2);
    DEL_2D(mat_k_1);
    DEL_2D(mat_k_2);
    DEL_2D(mat_k_3);
    DEL_2D(cut_phi);
    DEL_2D(cut_rho);
    delete [] mat_A;
    
    mat_t_1=NULL;
    mat_t_2=NULL;
    cut_phi=NULL;
    cut_rho=NULL;
    mat_k_1=NULL;
    mat_k_2=NULL;
    mat_k_3=NULL;
    mat_A=NULL;
    
    no_types=0;
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_fs::read_file(char* file_name)
{    
    FileReader fr(mapp);
    
    fr.add_1D("A",mat_A);
    
    
    fr.add_2D("r_c_rho",cut_rho);
    fr.symmetric();
    fr.add_vlog()=vlogic("gt",0.0);
    
    fr.add_2D("r_c_phi",cut_phi);
    fr.symmetric();
    fr.add_vlog()=vlogic("gt",0.0);
    
    fr.add_2D("k1",mat_k_1);
    fr.symmetric();
    
    fr.add_2D("k2",mat_k_2);
    fr.symmetric();
    
    fr.add_2D("k3",mat_k_3);
    fr.symmetric();
    
    
    fr.add_2D("t1",mat_t_1);
    
    fr.add_2D("t2",mat_t_2);
    
    
    fr.read_file(file_name);
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_fs::init()
{
    neighbor->pair_wise=true;
    
    rho_ptr=new Vec<type0>(atoms,1);
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceField_fs::fin()
{
    if(max_pairs)
    {
        delete [] drhoi_dr;
        delete [] drhoj_dr;
        max_pairs=0;
    }
    
    delete rho_ptr;
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceField_fs::init_xchng()
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceField_fs::fin_xchng()
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceField_fs::pre_xchng_energy(GCMC*)
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 xchng energy
 --------------------------------------------*/
type0 ForceField_fs::xchng_energy(GCMC*)
{
    error->abort("exchange has not been set for this forcefield");
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceField_fs::post_xchng_energy(GCMC*)
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceField_fs::
force_calc(bool st_clc)
{
    if(max_pairs<neighbor->no_pairs)
    {
        if(max_pairs)
        {
            delete [] drhoi_dr;
            delete [] drhoj_dr;
        }
        
        max_pairs=neighbor->no_pairs;
        CREATE1D(drhoi_dr,max_pairs);
        CREATE1D(drhoj_dr,max_pairs);
        for(int i=0;i<max_pairs;i++)
            drhoi_dr[i]=drhoj_dr[i]=0.0;
    }

    
    type0* x=mapp->x->begin();
    type0* fvec=f->begin();
    type0* rho=rho_ptr->begin();
    md_type* type=mapp->type->begin();
    
    
    int iatm,jatm;
    
    int itype,jtype,icomp,jcomp;
    type0 dx0,dx1,dx2,rsq,en;
    type0 dr_rho,dr_phi,r,rho_coef,phi_coef,r_inv,rho_sqd;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    nrgy_strss_lcl[0]=0.0;
    if(st_clc)
        for(int i=1;i<7;i++)
            nrgy_strss_lcl[i]=0.0;
    
    int istart;
    int natms=atoms->natms;
    
    for(iatm=0;iatm<natms;iatm++) rho[iatm]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=3*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            
            drhoi_dr[istart]=drhoj_dr[istart]=0.0;
            
                
            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            if(rsq<cut_sq[itype][jtype])
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                
                if(r < cut_rho[itype][jtype])
                {
                    dr_rho=r-cut_rho[itype][jtype];
                    rho[iatm]+=dr_rho*dr_rho*(mat_t_1[jtype][itype]
                    +mat_t_2[jtype][itype]*dr_rho);
                    
                    
                    drhoi_dr[istart]=-r_inv*dr_rho*(mat_t_1[jtype][itype]
                    +1.5*mat_t_2[jtype][itype]*dr_rho);
                    
                    drhoj_dr[istart]=-r_inv*dr_rho*(mat_t_1[itype][jtype]
                    +1.5*mat_t_2[itype][jtype]*dr_rho);
                    
                    if(jatm<natms)
                        rho[jatm]+=dr_rho*dr_rho*(mat_t_1[itype][jtype]
                        +mat_t_2[itype][jtype]*dr_rho);
                    
                }
                
                if(r < cut_phi[itype][jtype])
                {
                    dr_phi=r-cut_phi[itype][jtype];
                    phi_coef=2.0*dr_phi*(mat_k_1[itype][jtype]
                    +mat_k_2[itype][jtype]*r+mat_k_3[itype][jtype]*rsq)
                    +dr_phi*dr_phi*(mat_k_2[itype][jtype]+2.0*mat_k_3[itype][jtype]*r);
                    phi_coef*=-r_inv;
                    
                    en=dr_phi*dr_phi*
                    (mat_k_1[itype][jtype]+mat_k_2[itype][jtype]*r+mat_k_3[itype][jtype]*rsq);
                    
                    fvec[icomp]+=dx0*phi_coef;
                    fvec[icomp+1]+=dx1*phi_coef;
                    fvec[icomp+2]+=dx2*phi_coef;
                    
                    if(jatm<natms)
                    {
                        fvec[jcomp]-=dx0*phi_coef;
                        fvec[jcomp+1]-=dx1*phi_coef;
                        fvec[jcomp+2]-=dx2*phi_coef;
                    }
                    
                    if(jatm>=natms)
                    {
                        phi_coef*=0.5;
                        en*=0.5;
                    }
                    
                    nrgy_strss_lcl[0]+=en;
                    
                    if(st_clc)
                    {
                        nrgy_strss_lcl[1]-=phi_coef*dx0*dx0;
                        nrgy_strss_lcl[2]-=phi_coef*dx1*dx1;
                        nrgy_strss_lcl[3]-=phi_coef*dx2*dx2;
                        nrgy_strss_lcl[4]-=phi_coef*dx1*dx2;
                        nrgy_strss_lcl[5]-=phi_coef*dx2*dx0;
                        nrgy_strss_lcl[6]-=phi_coef*dx0*dx1;
                    }
                }
            }
        
            istart++;
        }
        
        if(rho[iatm]>0.0)
        {
            rho_sqd=sqrt(rho[iatm]);
            nrgy_strss_lcl[0]-=mat_A[itype]*rho_sqd;
            rho[iatm]=-mat_A[itype]/rho_sqd;
        }
        else
        {
            rho[iatm]=0.0;
        }
    }
    
    atoms->update(rho_ptr);
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=3*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            
            if(drhoi_dr[istart]!=0.0 || drhoj_dr[istart]!=0.0)
            {
                jatm=neighbor_list[iatm][j];
                jcomp=3*jatm;
                dx0=x[icomp]-x[jcomp];
                dx1=x[icomp+1]-x[jcomp+1];
                dx2=x[icomp+2]-x[jcomp+2];
                
                rho_coef=drhoi_dr[istart]*rho[iatm]
                +drhoj_dr[istart]*rho[jatm];
                
                fvec[icomp]+=dx0*rho_coef;
                fvec[icomp+1]+=dx1*rho_coef;
                fvec[icomp+2]+=dx2*rho_coef;
                
                if(jatm<natms)
                {
                    fvec[jcomp]-=dx0*rho_coef;
                    fvec[jcomp+1]-=dx1*rho_coef;
                    fvec[jcomp+2]-=dx2*rho_coef;
                }
                
                if(jatm>=natms)
                    rho_coef*=0.5;
                
                if(st_clc)
                {
                    nrgy_strss_lcl[1]-=rho_coef*dx0*dx0;
                    nrgy_strss_lcl[2]-=rho_coef*dx1*dx1;
                    nrgy_strss_lcl[3]-=rho_coef*dx2*dx2;
                    nrgy_strss_lcl[4]-=rho_coef*dx1*dx2;
                    nrgy_strss_lcl[5]-=rho_coef*dx2*dx0;
                    nrgy_strss_lcl[6]-=rho_coef*dx0*dx1;
                }
            }
            istart++;

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
        MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,1,MPI_TYPE0,MPI_SUM,world);
    }
    
    
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
type0 ForceField_fs::energy_calc()
{
    type0* x=mapp->x->begin();
    type0* rho=rho_ptr->begin();
    md_type* type=mapp->type->begin();
    
    int natms=atoms->natms;
    int iatm,jatm;
    
    int itype,jtype,icomp,jcomp;
    type0 dx0,dx1,dx2,rsq;
    type0 dr_rho,dr_phi,r;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    for(int i=0;i<natms;i++)
        rho[i]=0.0;
    
    type0 en=0.0;
    type0 en_tot=0.0;
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=3*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];

            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;

            
            if(rsq <cut_sq[itype][jtype])
            {
                r=sqrt(rsq);
                if(r < cut_rho[itype][jtype])
                {
                    dr_rho=r-cut_rho[itype][jtype];
                    rho[iatm]+=dr_rho*dr_rho*(mat_t_1[jtype][itype]
                    +mat_t_2[jtype][itype]*dr_rho);
                    if(jatm<natms)
                        rho[jatm]+=dr_rho*dr_rho*(mat_t_1[itype][jtype]
                        +mat_t_2[itype][jtype]*dr_rho);
                }
                
                if(r < cut_phi[itype][jtype])
                {
                    dr_phi=r-cut_phi[itype][jtype];
                    if(jatm<natms)
                    {
                        en+=dr_phi*dr_phi*(mat_k_1[itype][jtype]
                        +mat_k_2[itype][jtype]*r+mat_k_3[itype][jtype]*rsq);
                    }
                    else
                    {
                        en+=0.5*dr_phi*dr_phi*(mat_k_1[itype][jtype]
                        +mat_k_2[itype][jtype]*r+mat_k_3[itype][jtype]*rsq);
                    }
                }
            }
            
        }
        
        if(rho[iatm]>0.0)
            en-=mat_A[itype]*sqrt(rho[iatm]);
    }

    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    return en_tot;
}

