/*--------------------------------------------
 Created by Sina on 01/29/14.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <stdlib.h>
#include "ff_fsm.h"
#include "neighbor.h"
#include "atom_types.h"
#include "error.h"
#include "memory.h"
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
ForceField_fsm::
ForceField_fsm(MAPP* mapp):ForceFieldMD(mapp)
{
    if(mapp->mode!=MD_mode)
        error->abort("ff fsm works only "
        "for md mode");
    no_types=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_fsm::~ForceField_fsm()
{
    deallocate();
}
/*--------------------------------------------
 read the force field file
 --------------------------------------------*/
void ForceField_fsm::coef(int nargs,char** args)
{
    if (nargs!=2)
        error->abort("wrong coeff command "
        "for Finnis-Sinclair Force Field");
    
    cut_off_alloc();
    allocate();
    
    read_file(args[1]);
    
    for(int i=0;i<no_types*(no_types+1)/2;i++)
        cut_sq[i]=MAX(cut_sq_phi[i],cut_sq_rho[i]);
    
}
/*--------------------------------------------
 allocate
 --------------------------------------------*/
void ForceField_fsm::allocate()
{
    if(no_types==atom_types->no_types)
        return;
    deallocate();
    no_types=atom_types->no_types;
    CREATE2D(mat_t_1,no_types,no_types);
    CREATE2D(mat_t_2,no_types,no_types);
    CREATE1D(mat_A,no_types);
    for(int i=0;i<no_types;i++)
    {
        for(int j=0;j<no_types;j++)
            mat_t_1[i][j]=mat_t_2[i][j]=0.0;
        mat_A[i]=0.0;
    }
    
    int arr_size=no_types*(no_types+1)/2;
    CREATE1D(cut_sq_phi,arr_size);
    CREATE1D(cut_sq_rho,arr_size);
    CREATE1D(mat_k_1,arr_size);
    CREATE1D(mat_k_2,arr_size);
    CREATE1D(mat_k_3,arr_size);
    
    for(int i=0;i<arr_size;i++)
        cut_sq_phi[i]=cut_sq_rho[i]
        =mat_k_1[i]=mat_k_2[i]=mat_k_3[i]=0.0;
}
/*--------------------------------------------
 allocate
 --------------------------------------------*/
void ForceField_fsm::deallocate()
{
    if(!no_types)
        return;

    for(int i=0;i<no_types;i++)
    {
        delete [] mat_t_1[i];
        delete [] mat_t_2[i];
    }
    
    delete [] mat_t_1;
    delete [] mat_t_2;
    delete [] mat_A;
    delete [] cut_sq_phi;
    delete [] cut_sq_rho;
    delete [] mat_k_1;
    delete [] mat_k_2;
    delete [] mat_k_3;
    
    no_types=0;
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_fsm::read_file(char* file_name)
{
    
    int* type_ref;
    
    int* A_chk;
    
    int** t1_chk;
    int** t2_chk;
    
    int* k1_chk;
    int* k2_chk;
    int* k3_chk;
    int* r_c_phi_chk;
    int* r_c_rho_chk;
    
    CREATE1D(A_chk,no_types);
    
    CREATE2D(t1_chk,no_types,no_types);
    CREATE2D(t2_chk,no_types,no_types);
    
    CREATE1D(k1_chk,no_types*(no_types+1));
    CREATE1D(k2_chk,no_types*(no_types+1));
    CREATE1D(k3_chk,no_types*(no_types+1));
    CREATE1D(r_c_phi_chk,no_types*(no_types+1));
    CREATE1D(r_c_rho_chk,no_types*(no_types+1));
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    
    char** args=NULL;
    int args_cpcty=0;
    int nargs;
    int no_types_file;
    
    for(int i=0;i<no_types;i++)
    {
        A_chk[i]=0;
        for(int j=0;j<no_types;j++)
            t1_chk[i][j]=t2_chk[i][j]=0;
    }
    
    for(int i=0;i<no_types*(no_types+1);i++)
        k1_chk[i]=k2_chk[i]=k3_chk[i]=r_c_phi_chk[i]=r_c_rho_chk[i]=0;
    
    
    if(atoms->my_p==0)
    {
        fp=fopen(file_name,"r");
        if(fp==NULL)
            error->abort("ff fsm file %s not found",file_name);
    }
    
    /*
     reading the header of the file
     find the first line and read the
     atomic types in the file
     */
    
    nargs=0;
    while(nargs==0 && mapp->read_line(fp,line) !=-1)
        nargs=mapp->parse_line(line,args,args_cpcty);
    
    if(nargs==0)
        error->abort("%s file ended immaturely",file_name);
    
    if(nargs<no_types)
        error->abort("the number of atoms in %s file"
                     " is less than the number of atom types present in the system",file_name);
    
    no_types_file=nargs;
    
    
    CREATE1D(type_ref,no_types);
    
    for(int i=0;i<no_types_file;i++)
    {
        type_ref[i]=atom_types->find_type_exist(args[i]);
    }
    
    if(args_cpcty)
        delete [] args;
    
    
    //lngth=read_line(fp,line);
    int icmp,jcmp,curs;
    type0 tmp;
    while(mapp->read_line(fp,line)!=-1)
    {
        if(mapp->hash_remover(line)==0)
            continue;

        if(sscanf(line,"A(%d) = %lf",&icmp,&tmp)==2)
        {
            if(type_ref[icmp]!=-1)
            {
                A_chk[type_ref[icmp]]=1;
                mat_A[type_ref[icmp]]=tmp;
            }
        }
        else if(sscanf(line,"t1(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<0 || icmp>no_types_file-1)
                error->abort("wrong component in %s file for t1(%d,%d)",file_name,icmp,jcmp);
            if(jcmp<0 || jcmp>no_types_file-1)
                error->abort("wrong component in %s file for t1(%d,%d)",file_name,icmp,jcmp);
            if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
            {
                t1_chk[type_ref[icmp]][type_ref[jcmp]]=1;
                mat_t_1[type_ref[icmp]][type_ref[jcmp]]=tmp;
            }
        }
        else if(sscanf(line,"t2(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<0 || icmp>no_types_file-1)
                error->abort("wrong component in %s file for t2(%d,%d)",file_name,icmp,jcmp);
            if(jcmp<0 || jcmp>no_types_file-1)
                error->abort("wrong component in %s file for t2(%d,%d)",file_name,icmp,jcmp);
            if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
            {
                t2_chk[type_ref[icmp]][type_ref[jcmp]]=1;
                mat_t_2[type_ref[icmp]][type_ref[jcmp]]=tmp;
            }
        }
        else if(sscanf(line,"k1(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<0 || icmp>no_types_file-1)
                error->abort("wrong component in %s file for k1(%d,%d)",file_name,icmp,jcmp);
            if(jcmp<0 || jcmp>no_types_file-1)
                error->abort("wrong component in %s file for k1(%d,%d)",file_name,icmp,jcmp);
            if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
            {
                curs=COMP(type_ref[icmp],type_ref[jcmp]);
                k1_chk[curs]=1;
                mat_k_1[curs]=tmp;
            }
        }
        else if(sscanf(line,"k2(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<0 || icmp>no_types_file-1)
                error->abort("wrong component in %s file for k2(%d,%d)",file_name,icmp,jcmp);
            if(jcmp<0 || jcmp>no_types_file-1)
                error->abort("wrong component in %s file for k2(%d,%d)",file_name,icmp,jcmp);
            if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
            {
                curs=COMP(type_ref[icmp],type_ref[jcmp]);
                k2_chk[curs]=1;
                mat_k_2[curs]=tmp;
            }
        }
        else if(sscanf(line,"k3(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<0 || icmp>no_types_file-1)
                error->abort("wrong component in %s file for k3(%d,%d)",file_name,icmp,jcmp);
            if(jcmp<0 || jcmp>no_types_file-1)
                error->abort("wrong component in %s file for k3(%d,%d)",file_name,icmp,jcmp);
            if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
            {
                curs=COMP(type_ref[icmp],type_ref[jcmp]);
                k3_chk[curs]=1;
                mat_k_3[curs]=tmp;
            }
        }
        else if(sscanf(line,"r_c_phi(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<0 || icmp>no_types_file-1)
                error->abort("wrong component in %s file for r_c_phi(%d,%d)",file_name,icmp,jcmp);
            if(jcmp<0 || jcmp>no_types_file-1)
                error->abort("wrong component in %s file for r_c_phi(%d,%d)",file_name,icmp,jcmp);
            if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
            {
                if(tmp<=0.0)
                    error->abort("r_c_phi(%d,%d) in %s "
                                 "file should be greater than 0.0",file_name,icmp,jcmp);
                
                curs=COMP(type_ref[icmp],type_ref[jcmp]);
                r_c_phi_chk[curs]=1;
                cut_sq_phi[curs]=tmp*tmp;
            }
        }
        else if(sscanf(line,"r_c_rho(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<0 || icmp>no_types_file-1)
                error->abort("wrong component in %s file for r_c_rho(%d,%d)",file_name,icmp,jcmp);
            if(jcmp<0 || jcmp>no_types_file-1)
                error->abort("wrong component in %s file for r_c_rho(%d,%d)",file_name,icmp,jcmp);
            if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
            {
                if(tmp<=0.0)
                    error->abort("r_c_rho(%d,%d) in %s "
                                 "file should be greater than 0.0",file_name,icmp,jcmp);
                
                curs=COMP(type_ref[icmp],type_ref[jcmp]);
                r_c_rho_chk[curs]=1;
                cut_sq_rho[curs]=tmp*tmp;
            }
        }
        else
            error->abort("invalid line in %s file: %s",file_name,line);
        
    }
    
    if(atoms->my_p==0)
        fclose(fp);
    
    delete [] line;
    
    /*
     check wether all the values are set or not
     */
    
    for(int i=0;i<no_types;i++)
        if(A_chk[i]==0)
            error->abort("A(%s) was not set by %s file "
                         ,file_name,atom_types->atom_names[i]);
    
    for(int i=0;i<no_types;i++)
        for(int j=0;j<no_types;j++)
        {
            if(t1_chk[i][j]==0)
                error->abort("t1(%s,%s) was not set by %s file "
                ,atom_types->atom_names[i],atom_types->atom_names[j],file_name);
            if(t2_chk[i][j]==0)
                error->abort("t2(%s,%s) was not set by %s file "
                ,atom_types->atom_names[i],atom_types->atom_names[j],file_name);
        }
    
    for(int i=0;i<no_types;i++)
        for(int j=i;j<no_types;j++)
        {
            curs=COMP(i,j);
            if(k1_chk[curs]==0)
                error->abort("k1(%s,%s) was not set by %s file "
                ,atom_types->atom_names[i],atom_types->atom_names[j],file_name);
            if(k2_chk[curs]==0)
                error->abort("k2(%s,%s) was not set by %s file "
                ,atom_types->atom_names[i],atom_types->atom_names[j],file_name);
            if(k3_chk[curs]==0)
                error->abort("k3(%s,%s) was not set by %s file "
                ,atom_types->atom_names[i],atom_types->atom_names[j],file_name);
            if(r_c_phi_chk[curs]==0)
                error->abort("r_c_phi(%s,%s) was not set by %s file "
                ,atom_types->atom_names[i],atom_types->atom_names[j],file_name);
            if(r_c_rho_chk[curs]==0)
                error->abort("r_c_rho(%s,%s) was not set by %s file "
                ,atom_types->atom_names[i],atom_types->atom_names[j],file_name);
            
        }
    
    //clean up
    
    if(no_types)
    {
        delete [] A_chk;
        
        for(int i=0; i<no_types;i++)
        {
            delete [] t1_chk[i];
            delete [] t2_chk[i];
        }
        
        delete [] t1_chk;
        delete [] t2_chk;
        
        
        
        delete [] k1_chk;
        delete [] k2_chk;
        delete [] k3_chk;
        delete [] r_c_phi_chk;
        delete [] r_c_rho_chk;
    }
    
    if(no_types_file)
        delete [] type_ref;
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_fsm::init()
{
    neighbor->pair_wise=true;
    
    rho_ptr=new Vec<type0>(atoms,1);
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceField_fsm::fin()
{
    delete rho_ptr;
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceField_fsm::
force_calc(bool st_clc)
{
    type0* x=mapp->x->begin();
    type0* fvec=f->begin();
    type0* rho=rho_ptr->begin();
    md_type* type=mapp->type->begin();
    
    int iatm,jatm;
    
    int itype,jtype,curs,icomp,jcomp;
    type0 dx0,dx1,dx2,rsq,csq,en;
    type0 dr_rho,dr_phi,r,rho_coef,phi_coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    nrgy_strss_lcl[0]=0.0;
    if (st_clc)
        for(int i=1;i<7;i++)
            nrgy_strss_lcl[i]=0.0;
    
    
    int natms=atoms->natms;
    
    for(int i=0;i<natms;i++)
        rho[i]=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=3*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            curs=COMP(itype,jtype);
            
            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            if(rsq < cut_sq_rho[curs])
            {
                r=sqrt(rsq);
                dr_rho=r-sqrt(cut_sq_rho[curs]);
                rho[iatm]+=dr_rho*dr_rho*(mat_t_1[jtype][itype]
                +mat_t_2[jtype][itype]*dr_rho);
                
                if(jatm<natms)
                    rho[jatm]+=dr_rho*dr_rho*(mat_t_1[itype][jtype]
                    +mat_t_2[itype][jtype]*dr_rho);
                
            }
        }
    }
    
    atoms->update(rho_ptr);
    
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=3*iatm;
        
        if(rho[iatm]>0.0)
            nrgy_strss_lcl[0]+=-mat_A[itype]*sqrt(rho[iatm]);
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            curs=COMP(itype,jtype);
            
            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            csq=cut_sq[curs];
            
            if(rsq < csq)
            {
                r=sqrt(rsq);
                
                if(rsq < cut_sq_rho[curs])
                {
                    dr_rho=r-sqrt(cut_sq_rho[curs]);
                    
                    rho_coef=0.0;
                    
                    if(rho[iatm]>0.0)
                        rho_coef+=mat_A[itype]*(dr_rho*(1.0*mat_t_1[jtype][itype]
                        +1.5*mat_t_2[jtype][itype]*dr_rho))/sqrt(rho[iatm]);

                    if(rho[jatm]>0.0)
                        rho_coef+=mat_A[jtype]*(dr_rho*(1.0*mat_t_1[itype][jtype]
                        +1.5*mat_t_2[itype][jtype]*dr_rho))/sqrt(rho[jatm]);
                    
                    rho_coef*=1.0/r;
                    
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
                
                if(rsq < cut_sq_phi[curs])
                {
                    dr_phi=r-sqrt(cut_sq_phi[curs]);
                    phi_coef=2.0*dr_phi*(mat_k_1[curs]+mat_k_2[curs]*r+mat_k_3[curs]*rsq)
                    +dr_phi*dr_phi*(mat_k_2[curs]+2.0*mat_k_3[curs]*r);
                    phi_coef*=-1.0/r;
                    en=dr_phi*dr_phi*(mat_k_1[curs]+mat_k_2[curs]*r+mat_k_3[curs]*rsq);
                    
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
                    
                    if (st_clc)
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
type0 ForceField_fsm::energy_calc()
{
    type0* x=mapp->x->begin();
    type0* rho=rho_ptr->begin();
    md_type* type=mapp->type->begin();
    
    int natms=atoms->natms;
    int iatm,jatm;
    
    int itype,jtype,curs,icomp,jcomp;
    type0 dx0,dx1,dx2,rsq,csq;
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
            curs=COMP(itype,jtype);
            
            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            csq=cut_sq[curs];
            
            if(rsq < csq)
            {
                r=sqrt(rsq);
                if(rsq < cut_sq_rho[curs])
                {
                    dr_rho=r-sqrt(cut_sq_rho[curs]);
                    rho[iatm]+=dr_rho*dr_rho*(mat_t_1[jtype][itype]
                                              +mat_t_2[jtype][itype]*dr_rho);
                    if(jatm<natms)
                        rho[jatm]+=dr_rho*dr_rho*(mat_t_1[itype][jtype]
                                                  +mat_t_2[itype][jtype]*dr_rho);
                }
                
                if(rsq < cut_sq_phi[curs])
                {
                    dr_phi=r-sqrt(cut_sq_phi[curs]);
                    if(jatm<natms)
                    {
                        en+=dr_phi*dr_phi*(mat_k_1[curs]+mat_k_2[curs]*r+mat_k_3[curs]*rsq);
                    }
                    else
                    {
                        en+=0.5*dr_phi*dr_phi*(mat_k_1[curs]+mat_k_2[curs]*r+mat_k_3[curs]*rsq);
                    }
                }
            }
        }
        
        if(rho[iatm]>0.0)
            en+=-mat_A[itype]*sqrt(rho[iatm]);
    }
    
    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    return en_tot;
}

