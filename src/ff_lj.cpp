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
using namespace MAPP_NS;
/*--------------------------------------------
 Fisher-Sinclair (FS) potential
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
ForceField_lj::
ForceField_lj(MAPP* mapp):ForceFieldMD(mapp)
{
    if(mapp->mode!=MD_mode)
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
    
    if (nargs!=2)
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
    int arr_size=no_types*(no_types+1)/2;
    CREATE1D(sigma,arr_size);
    CREATE1D(epsilon,arr_size);
    CREATE1D(offset,arr_size);
}
/*--------------------------------------------
 allocation
 --------------------------------------------*/
void ForceField_lj::deallocate()
{
    if(!no_types)
        return;
    
    delete [] sigma;
    delete [] epsilon;
    delete [] offset;
    
    no_types=0;
    
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_lj::init()
{
    neighbor->pair_wise=1;
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceField_lj::fin()
{
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceField_lj::
force_calc(bool st_clc,type0* en_st)
{
    type0* x=mapp->x->begin();
    type0* f=mapp->f->begin();
    md_type* type=mapp->type->begin();

    int natms=atoms->natms;
    int iatm,jatm;
    int itype,jtype,curs,icomp,jcomp;
    type0 dx0,dx1,dx2,rsq,sig,eps,csq;
    type0 sig2,sig6,sig12,fpair,en;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    nrgy_strss[0]=0.0;
    if (st_clc)
        for(int i=1;i<7;i++)
            nrgy_strss[i]=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            curs=COMP(itype,jtype);
            icomp=3*iatm;
            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            csq=cut_sq[curs];
            if (rsq<csq)
            {
                sig=sigma[curs];
                eps=epsilon[curs];
                sig2=sig*sig/rsq;
                sig6=sig2*sig2*sig2;
                sig12=sig6*sig6;
                
                fpair=24.0*eps*(2.0*sig12-sig6)/rsq;
                en=4.0*eps*(sig12-sig6)
                +offset[curs];
                
                f[icomp]+=fpair*dx0;
                f[icomp+1]+=fpair*dx1;
                f[icomp+2]+=fpair*dx2;
                if(jatm<natms)
                {
                    f[jcomp]-=fpair*dx0;
                    f[jcomp+1]-=fpair*dx1;
                    f[jcomp+2]-=fpair*dx2;
                }

               
                if(jatm>=natms)
                {
                    fpair*=0.5;
                    en*=0.5;
                }
                nrgy_strss[0]+=en;
                
                if (st_clc)
                {
                    nrgy_strss[1]-=fpair*dx0*dx0;
                    nrgy_strss[2]-=fpair*dx1*dx1;
                    nrgy_strss[3]-=fpair*dx2*dx2;
                    nrgy_strss[4]-=fpair*dx1*dx2;
                    nrgy_strss[5]-=fpair*dx2*dx0;
                    nrgy_strss[6]-=fpair*dx0*dx1;
                }
            }
        }
    }
    
    if(st_clc)
    {
        for(int i=0;i<7;i++)
            en_st[i]=0.0;
        
        MPI_Allreduce(nrgy_strss,en_st,7,MPI_TYPE0,MPI_SUM,world);
    }
    else
    {
        en_st[0]=0.0;
        MPI_Allreduce(nrgy_strss,en_st,1,MPI_TYPE0,MPI_SUM,world);
    }
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
type0 ForceField_lj::energy_calc()
{
    type0* x=mapp->x->begin();
    md_type* type=mapp->type->begin();
    
    int natms=atoms->natms;
    int iatm,jatm;
    
    int itype,jtype,curs,icomp,jcomp;
    type0 dx0,dx1,dx2,rsq,csq;
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
            curs=COMP(itype,jtype);
            icomp=3*iatm;
            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            csq=cut_sq[curs];
            if (rsq<csq)
            {
                sig=sigma[curs];
                eps=epsilon[curs];
                sig2=sig*sig/rsq;
                sig6=sig2*sig2*sig2;
                sig12=sig6*sig6;
                
                if(jatm<natms)
                    en+=4.0*eps*(sig12-sig6)
                    +offset[curs];
                else
                    en+=2.0*eps*(sig12-sig6)
                    +offset[curs]*0.5;
                
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
    int no_types=atom_types->no_types;
    int* type_ref;
    int* eps_chk;
    int* sigma_chk;
    int* r_c_chk;

    CREATE1D(eps_chk,no_types*(no_types+1));
    CREATE1D(sigma_chk,no_types*(no_types+1));
    CREATE1D(r_c_chk,no_types*(no_types+1));

    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    
    char** args=NULL;
    int nargs;
    int args_cpcty=0;
    int no_types_file;
    

    for(int i=0;i<no_types*(no_types+1);i++)
        eps_chk[i]=sigma_chk[i]=eps_chk[i]=0;
    
    
    if(atoms->my_p==0)
    {
        fp=fopen(file_name,"r");
        if(fp==NULL)
            error->abort("ff lj file %s not found",file_name);
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
        type_ref[i]=atom_types->find_type_exist(args[i]);
    
    if(args_cpcty)
        delete [] args;
    
    
    int icmp,jcmp,curs;
    type0 tmp;
    while(mapp->read_line(fp,line)!=-1)
    {
        if(mapp->hash_remover(line)==0)
            continue;
        if(sscanf(line,"r_c(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<0 || icmp>no_types_file-1)
                error->abort("wrong component in %s file for r_c(%i,%i)",file_name,icmp,jcmp);
            if(jcmp<0 || jcmp>no_types_file-1)
                error->abort("wrong component in %s file for r_c(%i,%i)",file_name,icmp,jcmp);
            
            if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
            {
                curs=COMP(type_ref[icmp],type_ref[jcmp]);
                if(tmp<=0.0)
                    error->abort("r_c(%d,%d) in %s "
                                 "file should be greater than 0.0",icmp,jcmp,file_name);
                r_c_chk[curs]=1;
                cut_sq[curs]=tmp*tmp;
                type0 skin=atoms->get_skin();
                cut_sk_sq[curs]=(tmp+skin)*(tmp+skin);
            }
        }
        else if(sscanf(line,"epsilon(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<0 || icmp>no_types_file-1)
                error->abort("wrong component in %s file for epsilon(%i,%i)",file_name,icmp,jcmp);
            if(jcmp<0 || jcmp>no_types_file-1)
                error->abort("wrong component in %s file for epsilon(%i,%i)",file_name,icmp,jcmp);
            if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
            {
                if(tmp<=0.0)
                    error->abort("epsilon(%d,%d) in %s "
                                 "file should be greater than 0.0",icmp,jcmp,file_name);
                curs=COMP(type_ref[icmp],type_ref[jcmp]);
                eps_chk[curs]=1;
                epsilon[curs]=tmp;
            }
        }
        else if(sscanf(line,"sigma(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            if(icmp<0 || icmp>no_types_file-1)
                error->abort("wrong component in %s file for sigma(%i,%i)",file_name,icmp,jcmp);
            if(jcmp<0 || jcmp>no_types_file-1)
                error->abort("wrong component in %s file for sigma(%i,%i)",file_name,icmp,jcmp);
            if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
            {
                if(tmp<=0.0)
                    error->abort("sigma(%d,%d) in %s "
                                 "file should be greater than 0.0",icmp,jcmp,file_name);
                curs=COMP(type_ref[icmp],type_ref[jcmp]);
                sigma_chk[curs]=1;
                sigma[curs]=tmp;
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
        for(int j=i;j<no_types;j++)
        {
            curs=COMP(i,j);
            if(eps_chk[curs]==0)
                error->abort("epsilon(%s,%s) was not set by %s file "
                ,atom_types->atom_names[i],atom_types->atom_names[j],file_name);
            if(sigma_chk[curs]==0)
                error->abort("sigma(%s,%s) was not set by %s file "
                ,atom_types->atom_names[i],atom_types->atom_names[j],file_name);
            if(r_c_chk[curs]==0)
                error->abort("r_c(%s,%s) was not set by %s file "
                ,atom_types->atom_names[i],atom_types->atom_names[j],file_name);
        }
    
    
    //clean up
    
    if(no_types)
    {
        delete [] eps_chk;
        delete [] sigma_chk;
        delete [] r_c_chk;

    }
    
    if(no_types_file)
        delete [] type_ref;
    
    
    if(shift)
    {
        type0 sig2,sig6,sig12;
        for(int icurs=0;icurs<no_types*(no_types+1)/2;icurs++)
        {
            sig2=sigma[icurs]*sigma[icurs]/cut_sq[icurs];
            sig6=sig2*sig2*sig2;
            sig12=sig6*sig6;
            offset[icurs]=-4.0*epsilon[icurs]*(sig12-sig6);
        }
    }
    else
    {
        for(int icurs=0;icurs<no_types*(no_types+1)/2;icurs++)
            offset[icurs]=0.0;
    }
}



