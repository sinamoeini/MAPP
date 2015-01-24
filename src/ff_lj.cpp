/*--------------------------------------------
 Created by Sina on 07/15/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <cmath>
#include "ff_lj.h"
#include "memory.h"
#include "neighbor.h"
#include "atom_types.h"
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
ForceField_lj(MAPP* mapp) : ForceField(mapp)
{
    if(mapp->mode!=MD_mode)
        error->abort("ff lj works only "
        "for md mode");
    
    arr_size=shift=0;
    int no_types=atom_types->no_types;
    
    int size=static_cast<int>((no_types+2)*(no_types+1)/2);
    CREATE1D(cut_sq,size);
    CREATE1D(cut_sk_sq,size);
    CREATE1D(sigma,size);
    CREATE1D(epsilon,size);
    CREATE1D(offset,size);
    arr_size=size;
    
    CREATE1D(nrgy_strss,7);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_lj::~ForceField_lj()
{
    if(arr_size)
    {
        delete [] cut_sq;
        delete [] cut_sk_sq;
        delete [] sigma;
        delete [] epsilon;
        delete [] offset;
    }
    
    delete [] nrgy_strss;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ForceField_lj::coef(int narg,char** arg)
{
    //itype jtype epsilon sigma cuttoff
    
    if (narg!=6)
        error->abort("wrong coeff command for LJ FF");
    
    int ityp=atom_types->find_type(arg[1]);
    int jtyp=atom_types->find_type(arg[2]);
    int curs=COMP(ityp,jtyp);
    TYPE0 skin=atoms->skin;
    
    
    TYPE0 eps=epsilon[curs]=atof(arg[3]);
    TYPE0 sig=sigma[curs]=atof(arg[4]);
    TYPE0 cut=atof(arg[5]);

    if (eps<0.0)
    {
        error->abort("LJ Epsilon cannot be smaller than zero");
    }
    if (sig<0.0)
    {
        error->abort("LJ Sigma cannot be smaller than zero");
    }
    if (cut<0.0)
    {
        error->abort("LJ cutoff cannot be smaller than zero");
    }
    
    cut_sq[curs]=cut*cut;
    cut_sk_sq[curs]=(cut+skin)*(cut+skin);
    if(shift)
    {
        TYPE0 sig2=sig*sig/(cut*cut);
        TYPE0 sig6=sig2*sig2*sig2;
        TYPE0 sig12=sig6*sig6;
        
        offset[curs]=-4.0*eps*(sig12-sig6);
    }
    else
    {
        offset[curs]=0.0;
    }
    

}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_lj::init()
{
    TYPE0 skin=atoms->skin;
    TYPE0 ph_cut=0.0;
    
    for (int i=0;i<arr_size;i++)
    {
        cut_sk_sq[i]=cut_sq[i]+skin*skin
        +2*sqrt(cut_sq[i])*skin;
        ph_cut=MAX(ph_cut,sqrt(cut_sq[i]));
    }

    
    atoms->set_ph(ph_cut);
    
    x_n=atoms->find("x");
    f_n=atoms->find("f");
    type_n=atoms->find("type");
    
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
force_calc(int st_clc,TYPE0* en_st)
{

    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    TYPE0* f;
    atoms->vectors[f_n].ret(f);
    int* type;
    atoms->vectors[type_n].ret(type);

    int natms=atoms->natms;
    int iatm,jatm;
    int itype,jtype,curs,icomp,jcomp;
    TYPE0 dx0,dx1,dx2,rsq,sig,eps,csq;
    TYPE0 sig2,sig6,sig12,ft;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    nrgy_strss[0]=0.0;
    if (st_clc)
        for (int i=1;i<7;i++)
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
                
                ft=24.0*eps*(2.0*sig12-sig6)/rsq;
                
                
                f[icomp]+=ft*dx0;
                f[icomp+1]+=ft*dx1;
                f[icomp+2]+=ft*dx2;
                if(jatm<natms)
                {
                    nrgy_strss[0]+=4.0*eps*(sig12-sig6)
                    +offset[curs];
                    f[jcomp]-=ft*dx0;
                    f[jcomp+1]-=ft*dx1;
                    f[jcomp+2]-=ft*dx2;
                    if (st_clc)
                    {
                        nrgy_strss[1]-=ft*dx0*dx0;
                        nrgy_strss[2]-=ft*dx1*dx1;
                        nrgy_strss[3]-=ft*dx2*dx2;
                        nrgy_strss[4]-=ft*dx1*dx2;
                        nrgy_strss[5]-=ft*dx2*dx0;
                        nrgy_strss[6]-=ft*dx0*dx1;
                    }
                    
                }
                else
                {
                    nrgy_strss[0]+=2.0*eps*(sig12-sig6)
                    +offset[curs]*0.5;
                    if (st_clc)
                    {
                        nrgy_strss[1]-=0.5*ft*dx0*dx0;
                        nrgy_strss[2]-=0.5*ft*dx1*dx1;
                        nrgy_strss[3]-=0.5*ft*dx2*dx2;
                        nrgy_strss[4]-=0.5*ft*dx1*dx2;
                        nrgy_strss[5]-=0.5*ft*dx2*dx0;
                        nrgy_strss[6]-=0.5*ft*dx0*dx1;
                    }
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
TYPE0 ForceField_lj::energy_calc()
{
    /*
    TYPE0* x=(TYPE0*)atoms->vectors[x_n].ret_vec();
    int* type=(int*)atoms->vectors[type_n].ret_vec();
     */
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    int* type;
    atoms->vectors[type_n].ret(type);
    
    int natms=atoms->natms;
    int iatm,jatm;
    
    int itype,jtype,curs,icomp,jcomp;
    TYPE0 dx0,dx1,dx2,rsq,csq;
    TYPE0 eps,sig,sig2,sig6,sig12;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    TYPE0 en=0.0;
    TYPE0 en_tot=0.0;
    
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
    
    int lngth;
    char** args=NULL;
    int narg;
    int no_types_file;
    

    for(int i=0;i<no_types*(no_types+1);i++)
        eps_chk[i]=sigma_chk[i]=eps_chk[i]=0;
    
    
    if(atoms->my_p_no==0)
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
    
    lngth=read_line(fp,line);
    narg=0;
    while(narg==0)
    {
        if(lngth!=-1)
            narg=mapp->parse_line(line,args);
        else
            error->abort("%s file ended immaturely",file_name);
        
        lngth=read_line(fp,line);
    }
    
    if(narg<no_types)
        error->abort("the number of atoms in %s file"
        " is less than the number of atom types present in the system",file_name);
    
    no_types_file=narg;
    
    
    
    CREATE1D(type_ref,no_types);
    
    for(int i=0;i<no_types_file;i++)
    {
        type_ref[i]=atom_types->find_type_exist(args[i]);
    }
    
    for(int i=0;i<narg;i++)
        delete [] args[i];
    if(narg)
        delete [] args;
    
    
    //lngth=read_line(fp,line);
    int icmp,jcmp,curs;
    TYPE0 tmp;
    while(lngth!=-1)
    {
        if(lngth>1)
        {
            if(sscanf(line,"r_c(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
            {
                if(icmp<0 || icmp>no_types_file-1)
                    error->abort("wrong component in ff lj file for r_c(%i,%i)",icmp,jcmp);
                if(jcmp<0 || jcmp>no_types_file-1)
                    error->abort("wrong component in ff lj file for r_c(%i,%i)",icmp,jcmp);
                
                if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
                {
                    curs=COMP(type_ref[icmp],type_ref[jcmp]);
                    if(tmp<=0.0)
                        error->abort("r_c(%d,%d) in %s "
                        "file should be greater than 0.0",icmp,jcmp,file_name);
                    r_c_chk[curs]=1;
                    cut_sq[curs]=tmp*tmp;
                    TYPE0 skin=atoms->skin;
                    cut_sk_sq[curs]=(tmp+skin)*(tmp+skin);
                }
            }
            else if(sscanf(line,"epsilon(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
            {
                if(icmp<0 || icmp>no_types_file-1)
                    error->abort("wrong component in ff lj file for epsilon(%i,%i)",icmp,jcmp);
                if(jcmp<0 || jcmp>no_types_file-1)
                    error->abort("wrong component in ff lj file for epsilon(%i,%i)",icmp,jcmp);
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
                    error->abort("wrong component in ff lj file for sigma(%i,%i)",icmp,jcmp);
                if(jcmp<0 || jcmp>no_types_file-1)
                    error->abort("wrong component in ff lj file for sigma(%i,%i)",icmp,jcmp);
                if(type_ref[icmp]!=-1 && type_ref[jcmp]!=-1)
                {
                    if(tmp<=0.0)
                        error->abort("sigma(%d,%d) in %s "
                        "file should be greater than 0.0",icmp,jcmp,file_name);
                    curs=COMP(type_ref[icmp],type_ref[jcmp]);
                    eps_chk[curs]=1;
                    sigma[curs]=tmp;
                }
            }
            else
                error->abort("invalid line in %s file: %s",file_name,line);
        }
        
        lngth=read_line(fp,line);
    }
    
    if(atoms->my_p_no==0)
        fclose(fp);
    
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
    
    delete [] line;
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
int ForceField_lj::read_line(FILE* file,char*& line)
{
    int lenght;
    int eof=0;
    if(atoms->my_p_no==0)
        if(feof(file))
            eof=-1;
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if(eof==-1)
        return -1;
    
    if(atoms->my_p_no==0)
    {
        char* n_line;
        fgets(line,MAXCHAR,file);
        mapp->hash_remover(line,n_line);
        delete [] line;
        line =n_line;
        lenght=static_cast<int>(strlen(line))+1;
    }
    MPI_Bcast(&lenght,1,MPI_INT,0,world);
    MPI_Bcast(line,lenght,MPI_CHAR,0,world);
    
    return lenght;
}


