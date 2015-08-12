#include <stdlib.h>
#include "neighbor.h"
#include "ff_eam.h"
#include "atom_types.h"
using namespace MAPP_NS;
enum{NOT_SET,FUNC_FL,SET_FL,FINNIS_FL};

/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam::
ForceField_eam(MAPP* mapp) : ForceField(mapp)
{
    if(mapp->mode!=MD_mode)
        error->abort("ff eam works only "
        "for md mode");

    no_types=atom_types->no_types;
    
    allocated=0;
    eam_mode=NOT_SET;
    max_pairs=0;


    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_eam::~ForceField_eam()
{

    
    if(allocated) clean_up();
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceField_eam::
force_calc(int st_clc,type0* en_st)
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
    }
    
    type0* x;
    atoms->vectors[x_n]->ret(x);
    type0* f;
    atoms->vectors[f_n]->ret(f);
    type0* rho;
    atoms->vectors[rho_n]->ret(rho);
    /*
    type0* dF;
    atoms->vectors[dF_n]->ret(dF);
     */
    int* type;
    atoms->vectors[type_n]->ret(type);
    
    int iatm,jatm;
    
    int itype,jtype,icomp,jcomp;
    type0 dx0,dx1,dx2,rsq,z2p,z2;
    type0 r,p,r_inv,fpair,tmp0,tmp1;
    type0 drho_i_dr,drho_j_dr,dphi_dr;
    type0 rho_i,rho_j,phi;
    int m,istart;
    type0* coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    nrgy_strss[0]=0.0;
    if (st_clc)
        for (int i=1;i<7;i++)
            nrgy_strss[i]=0.0;
    
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
            
            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            drhoi_dr[istart]=drhoj_dr[istart]=0.0;
            if(rsq < cut_sq_0)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                p=r*dr_inv;
                m=static_cast<int>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);
                
                coef=rho_arr[type2rho[jtype][itype]][m];
                rho_i=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                drho_i_dr=(coef[6]*p+coef[5])*p+coef[4];
                coef=rho_arr[type2rho[itype][jtype]][m];
                rho_j=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                drho_j_dr=(coef[6]*p+coef[5])*p+coef[4];
                
                coef=phi_r_arr[type2phi[itype][jtype]][m];
                z2p=(coef[6]*p + coef[5])*p+coef[4];
                z2=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                
                phi=z2*r_inv;
                dphi_dr=z2p*r_inv-phi*r_inv;
            
                rho[iatm]+=rho_i;
                fpair=-dphi_dr*r_inv;
                f[icomp]+=fpair*dx0;
                f[icomp+1]+=fpair*dx1;
                f[icomp+2]+=fpair*dx2;
                if(jatm<natms)
                {
                    rho[jatm]+=rho_j;
                    f[jcomp]-=fpair*dx0;
                    f[jcomp+1]-=fpair*dx1;
                    f[jcomp+2]-=fpair*dx2;
                    
                    nrgy_strss[0]+=phi;
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
                else
                {
                    nrgy_strss[0]+=0.5*phi;
                    if (st_clc)
                    {
                        nrgy_strss[1]-=0.5*fpair*dx0*dx0;
                        nrgy_strss[2]-=0.5*fpair*dx1*dx1;
                        nrgy_strss[3]-=0.5*fpair*dx2*dx2;
                        nrgy_strss[4]-=0.5*fpair*dx1*dx2;
                        nrgy_strss[5]-=0.5*fpair*dx2*dx0;
                        nrgy_strss[6]-=0.5*fpair*dx0*dx1;
                    }
                }
                
                drhoi_dr[istart]=-drho_i_dr*r_inv;
                drhoj_dr[istart]=-drho_j_dr*r_inv;
            }
            
            istart++;
        }
        p=rho[iatm]*drho_inv;
        m=static_cast<int> (p);
        m=MIN(m,nr-2);
        p-=m;
        p=MIN(p,1.0);
        itype=type[iatm];
        coef=F_arr[itype][m];
        tmp1=(coef[6]*p+coef[5])*p+coef[4];
        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho[iatm]>rho_max)
            tmp0+=tmp1*(rho[iatm]-rho_max);
        nrgy_strss[0]+=tmp0;
        rho[iatm]=tmp1;
        
    }
    
    atoms->update_ph(rho_n);
    
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
                jtype=type[jatm];
                jcomp=3*jatm;
                
                fpair=rho[iatm]*drhoi_dr[istart]+rho[jatm]*drhoj_dr[istart];
                
                dx0=x[icomp]-x[jcomp];
                dx1=x[icomp+1]-x[jcomp+1];
                dx2=x[icomp+2]-x[jcomp+2];
                
                f[icomp]+=dx0*fpair;
                f[icomp+1]+=dx1*fpair;
                f[icomp+2]+=dx2*fpair;
                
                if(jatm<natms)
                {
                    f[jcomp]-=dx0*fpair;
                    f[jcomp+1]-=dx1*fpair;
                    f[jcomp+2]-=dx2*fpair;
                    
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
                else
                {
                    if (st_clc)
                    {
                        nrgy_strss[1]-=0.5*fpair*dx0*dx0;
                        nrgy_strss[2]-=0.5*fpair*dx1*dx1;
                        nrgy_strss[3]-=0.5*fpair*dx2*dx2;
                        nrgy_strss[4]-=0.5*fpair*dx1*dx2;
                        nrgy_strss[5]-=0.5*fpair*dx2*dx0;
                        nrgy_strss[6]-=0.5*fpair*dx0*dx1;
                    }
                }
            }
            istart++;
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
        MPI_Allreduce(nrgy_strss,en_st,1,MPI_TYPE0,MPI_SUM,world);
    }
}
/*--------------------------------------------
 energy calculation
 --------------------------------------------*/
type0 ForceField_eam::energy_calc()
{
    type0* x;
    atoms->vectors[x_n]->ret(x);
    type0* f;
    atoms->vectors[f_n]->ret(f);
    type0* rho;
    atoms->vectors[rho_n]->ret(rho);
    int* type;
    atoms->vectors[type_n]->ret(type);
    
    int iatm,jatm;
    
    int itype,jtype,icomp,jcomp;
    type0 dx0,dx1,dx2,rsq;
    type0 r,p,phi,tmp0;
    int m;
    type0* coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    type0 en=0.0;
    type0 en_tot=0.0;
    
    int natms=atoms->natms;
    
    for(iatm=0;iatm<natms;iatm++) rho[iatm]=0.0;
    
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
            
            if(rsq<cut_sq_0)
            {
                r=sqrt(rsq);
                
                p=r*dr_inv;
                m=static_cast<int>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);

                coef=phi_r_arr[type2phi[itype][jtype]][m];
                phi=(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0])/r;
                
                coef=rho_arr[type2rho[jtype][itype]][m];
                rho[iatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(jatm<natms)
                {
                    coef=rho_arr[type2rho[itype][jtype]][m];
                    rho[jatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                    en+=phi;
                }
                else
                    en+=0.5*phi;
                

            }
        }
        
        p=rho[iatm]*drho_inv;
        m=static_cast<int>(p);
        m=MIN(m,nr-2);
        p-=m;
        p=MIN(p,1.0);
        itype=type[iatm];
        coef=F_arr[itype][m];
        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho[iatm]>rho_max)
            tmp0+=((coef[6]*p+coef[5])*p+coef[4])*(rho[iatm]-rho_max);
        en+=tmp0;

    }
    
    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    return en_tot;
}
/*--------------------------------------------
 init before running
 --------------------------------------------*/
void ForceField_eam::init()
{
    for (int i=0;i<no_types*(no_types+1)/2;i++)
        cut_sq[i]=cut_sq_0;
    
    x_n=atoms->find("x");
    f_n=atoms->find("f");
    type_n=atoms->find("type");
    
    neighbor->pair_wise=1;
    
    rho_n=atoms->add<type0>(1,1,"rho");
}
/*--------------------------------------------
 fin after running
 --------------------------------------------*/
void ForceField_eam::fin()
{
    if(max_pairs)
    {
        delete [] drhoi_dr;
        delete [] drhoj_dr;
        max_pairs=0;
    }
    
    atoms->del(rho_n);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ForceField_eam::coef(int narg,char** arg)
{
    if (narg<3)
        error->abort("ff_coef for ff eam "
        "should at least have 2 arguments");
    
    clean_up();
    if(strcmp(arg[1],"FS")==0)
    {
        eam_mode=FINNIS_FL;
        set_fs(narg-2,&arg[2]);
    }
    else if(strcmp(arg[1],"SetFL")==0)
    {
        eam_mode=SET_FL;
        set_setfl(narg-2,&arg[2]);
    }
    else if(strcmp(arg[1],"FuncFL")==0)
    {
        eam_mode=FUNC_FL;
        set_funcfl(narg-2,&arg[2]);
    }
    else
    {
        error->abort("unknown format "
        "%s for ff eam file",arg[1]);
    }
    
    cut_sq_0=static_cast<type0>(nr*nr)*dr*dr;
    rho_max=static_cast<type0>(nrho)*drho;
    
    set_arrays();

}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam::set_funcfl(int no_files
,char** file_names)
{
    if(no_files!=no_types)
        error->abort("for FuncFL mode number of"
        " ff eam files should be equal to the number"
        " of atom types present in the system");
    
    type0* drs;
    type0* drhos;
    int* nrs;
    int* nrhos;
    CREATE1D(drs,no_types);
    CREATE1D(drhos,no_types);
    CREATE1D(nrs,no_types);
    CREATE1D(nrhos,no_types);
    
    type0** tmp_F;
    type0** tmp_rho;
    type0** tmp_zi;
    CREATE1D(tmp_F,no_types);
    CREATE1D(tmp_rho,no_types);
    CREATE1D(tmp_zi,no_types);
    type0* tmp;
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    int narg;
    char** arg;
    type0 mass;
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        fp=NULL;
        if(atoms->my_p_no==0)
        {
            fp=fopen(file_names[ityp],"r");
            if(fp==NULL)
                error->abort("ff eam file %s not found",file_names[ityp]);
        }
        
        for(int i=0;i<2;i++)
            if(line_read(fp,line)==-1)
                error->abort("%s file ended immaturely",file_names[ityp]);
        
        narg=mapp->parse_line(line,arg);
        if(narg!=4)
            error->abort("invalid line in %s file: %s",file_names[ityp],line);
        
        mass=atof(arg[1]);
        
        if(mass!=atom_types->mass[ityp])
            error->abort("mass of element %s in %s file (%lf) does not match "
            "the mass that is already assigned to the element (%lf), make "
            "sure that the right sequence of ff eam files are "
            "used",atom_types->atom_names[ityp],file_names[ityp],mass
            ,atom_types->mass[ityp]);
        
        for(int i=0;i<narg;i++)
            delete arg[i];
        if(narg)
            delete [] arg;
        
        
        
        
        if(line_read(fp,line)==-1)
            error->abort("%s file ended immaturely",file_names[ityp]);
        
        narg=mapp->parse_line(line,arg);
        if(narg!=5)
            error->abort("invalid line in %s file: %s",file_names[ityp],line);
        
        nrhos[ityp]=atoi(arg[0]);
        nrs[ityp]=atoi(arg[2]);
        drhos[ityp]=atof(arg[1]);
        drs[ityp]=atof(arg[3]);
        
        if(nrhos[ityp]<5)
            error->abort("nrho in %s file must be larger than 5",file_names[ityp]);
        if(nrs[ityp]<5)
            error->abort("nr in %s file must be larger than 5",file_names[ityp]);
        if(drhos[ityp]<=0.0)
            error->abort("drho in %s file must be larger than 0.0",file_names[ityp]);
        if(drs[ityp]<=0.0)
            error->abort("dr in %s file must be larger than 0.0",file_names[ityp]);
        
        for(int i=0;i<narg;i++)
            delete arg[i];
        if(narg)
            delete [] arg;
        
        
        
        int tot=nrhos[ityp]+2*nrs[ityp];
        
        CREATE1D(tmp,tot);
        
        int ipos=0;
        while (ipos<tot)
        {
            if(line_read(fp,line)==-1)
                error->abort("%s file ended immaturely",file_names[ityp]);
            
            narg=mapp->parse_line(line,arg);
            
            if(ipos+narg>tot)
                error->abort("%s file ended immaturely",file_names[ityp]);
            
            for(int i=0;i<narg;i++)
            {
                tmp[ipos]=atof(arg[i]);
                ipos++;
                delete [] arg[i];
            }
            if(narg)
                delete [] arg;
            
        }
        
        if(atoms->my_p_no==0)
            fclose(fp);
        
        CREATE1D(tmp_F[ityp],nrhos[ityp]);
        CREATE1D(tmp_zi[ityp],nrs[ityp]);
        CREATE1D(tmp_rho[ityp],nrs[ityp]);
        

        memcpy(tmp_F[ityp],&tmp[0],nrhos[ityp]*sizeof(type0));
        memcpy(tmp_zi[ityp],&tmp[nrhos[ityp]],nrs[ityp]*sizeof(type0));
        memcpy(tmp_rho[ityp],&tmp[nrhos[ityp]+nrs[ityp]],nrs[ityp]*sizeof(type0));
        
        delete [] tmp;
        
    }
    delete [] line;
    
    
    type0 maxr=0.0;
    type0 maxrho=0.0;
    type0 maxdr=0.0;
    type0 maxdrho=0.0;
    
    for(int i=0;i<no_types;i++)
    {
        maxr=MAX(maxr,static_cast<type0>(nrs[i]-1)*drs[i]);
        maxrho=MAX(maxrho,static_cast<type0>(nrhos[i]-1)*drhos[i]);
        maxdr=MAX(maxdr,drs[i]);
        maxdrho=MAX(maxdrho,drhos[i]);
    }
    
    nr=static_cast<int>(maxr/maxdr+0.5);
    nrho=static_cast<int>(maxrho/maxdrho+0.5);
    dr=maxdr;
    drho=maxdrho;
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    allocate();
    
    type0 r,p,tmp0,tmp1;
    int k;
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        for(int i=0;i<nrho;i++)
        {
            r=static_cast<type0>(i)*drho;
            p=r/drhos[ityp];
            k=static_cast<int> (p);
            k=MIN(k,nrhos[ityp]-2);
            p-=k;
            p=MIN(p,1.0);

            F_arr[ityp][i][0]=interpolate(tmp_F[ityp],nrs[ityp],p,k);
        }
    }
    
    
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        for(int i=0;i<nr;i++)
        {
            r=static_cast<type0>(i)*dr;
            p=r/drs[ityp];
            k=static_cast<int> (p);
            k=MIN(k,nrs[ityp]-2);
            p-=k;
            p=MIN(p,1.0);
            rho_arr[type2rho[ityp][0]][i][0]=interpolate(tmp_rho[ityp],nrs[ityp],p,k);
        }
    }
    
    
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        for(int jtyp=0;jtyp<ityp+1;jtyp++)
        {
            for(int i=0;i<nr;i++)
            {
                
                r=static_cast<type0>(i)*dr;
                
                p=r/drs[ityp];
                k=static_cast<int> (p);
                k=MIN(k,nrs[ityp]-2);
                p-=k;
                p=MIN(p,1.0);
                tmp0=interpolate(tmp_zi[ityp],nrs[ityp],p,k);

                
                p=r/drs[jtyp];
                k=static_cast<int> (p);
                k=MIN(k,nrs[jtyp]-2);
                p-=k;
                p=MIN(p,1.0);
                tmp1=interpolate(tmp_zi[jtyp],nrs[jtyp],p,k);
                
                phi_r_arr[type2phi[ityp][jtyp]][i][0]=27.2*0.529*tmp0*tmp1;
            }
        }
    }
    
    
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        if(nrhos[ityp])
            delete [] tmp_F[ityp];
        
        if(nrs[ityp])
        {
            delete [] tmp_rho[ityp];
            delete [] tmp_zi[ityp];
        }
    }
    if(no_types)
    {
        delete [] tmp_F;
        delete [] tmp_rho;
        delete [] tmp_zi;
        delete [] drhos;
        delete [] drs;
        delete [] nrhos;
        delete [] nrs;
    }

}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam::set_setfl(int no_files
,char** file_names)
{
    if(no_files!=1)
        error->abort("for SetFL mode number of"
        " ff eam files should be one");
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    
    fp=NULL;
    if(atoms->my_p_no==0)
    {
        fp=fopen(file_names[0],"r");
        if(fp==NULL)
            error->abort("%s file not found",file_names[0]);
    }
    
    for(int i=0;i<4;i++)
        if(line_read(fp,line)==-1)
            error->abort("%s file ended immaturely",file_names[0]);
    
    char** arg;
    int narg;
    
    narg=mapp->parse_line(line,arg);
    if(narg<2 || narg-1 < no_types)
        error->abort("invalid line in %s file: %s",file_names[0],line);
    
    int tot_no_types=narg-1;
    
    int* type_ref;
    CREATE1D(type_ref,tot_no_types);
    for(int i=0;i<tot_no_types;i++)
        type_ref[i]=atom_types->find_type_exist(arg[i+1]);

    for(int i=0;i<narg;i++)
        delete [] arg[i];
    if(narg)
        delete [] arg;
   
    int* tmp_type_ref;
    CREATE1D(tmp_type_ref,no_types);
    for(int i=0;i<no_types;i++)
        tmp_type_ref[i]=-1;
    for(int ityp=0;ityp<tot_no_types;ityp++)
        if(type_ref[ityp]!=-1)
            tmp_type_ref[type_ref[ityp]]=0;
    
    for(int i=0;i<no_types;i++)
        if(tmp_type_ref[i]==-1)
            error->abort("%s file does not contain "
            "parameters for element %s",file_names[0]
            ,atom_types->atom_names[i]);
    if(no_types)
        delete [] tmp_type_ref;
    
    if(line_read(fp,line)==-1)
        error->abort("%s file ended immaturely",file_names[0]);
    narg=mapp->parse_line(line,arg);
    
    
    if(narg!=5)
        error->abort("invalid line in %s file: %s",file_names[0],line);
    
    nrho=atoi(arg[0]);
    nr=atoi(arg[2]);
    drho=atof(arg[1]);
    dr=atof(arg[3]);
    
    
    for(int j=0;j<narg;j++)
        delete [] arg[j];
    if(narg)
        delete [] arg;
    
    if(nrho<5)
        error->abort("nrho in %s file must be larger than 5",file_names[0]);
    if(nr<5)
        error->abort("nr in %s file must be larger than 5",file_names[0]);
    if(drho<=0.0)
        error->abort("drho in %s file must be larger than 0.0",file_names[0]);
    if(dr<=0.0)
        error->abort("dr in %s file must be larger than 0.0",file_names[0]);
    
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    
    if(line_read(fp,line)==-1)
        error->abort("%s file ended immaturely",file_names[0]);
    
    
    int ipos=0;
    int tot=tot_no_types*(nrho+nr)+tot_no_types*(tot_no_types+1)*nr/2;
    
    type0* tmp;
    CREATE1D(tmp,tot);
    while (ipos<tot)
    {
        if(line_read(fp,line)==-1)
            error->abort("%s file ended immaturely",file_names[0]);
        
        narg=mapp->parse_line(line,arg);
        
        if(ipos+narg>tot)
            error->abort("%s file ended immaturely",file_names[0]);
        
        for(int i=0;i<narg;i++)
        {
            tmp[ipos]=atof(arg[i]);
            ipos++;
            delete [] arg[i];
        }
        if(narg)
            delete [] arg;
        
    }
    
    delete [] line;
    if(atoms->my_p_no==0)
        fclose(fp);
    
    
    allocate();
    
    
    int component;

    ipos=0;
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(type_ref[ityp]==-1)
        {
            ipos+=nrho+nr;
        }
        else
        {
            component=type_ref[ityp];
            
            for(int i=0;i<nrho;i++)
                F_arr[component][i][0]=tmp[ipos+i];
            
            ipos+=nrho;
            
            component=type2rho[type_ref[ityp]][type_ref[0]];
            for(int i=0;i<nr;i++)
                rho_arr[component][i][0]=tmp[ipos+i];
            
            ipos+=nr;
        }
        
    }
    

    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(type_ref[ityp]==-1)
        {
            ipos+=nr*(ityp+1);
        }
        else
        {
            for(int jtyp=0;jtyp<ityp+1;jtyp++)
            {
                if(type_ref[jtyp]==-1)
                {
                    ipos+=nr;
                }
                else
                {
                    component=type2phi[type_ref[ityp]][type_ref[jtyp]];
                    
                    for(int i=0;i<nr;i++)
                        phi_r_arr[component][i][0]=tmp[ipos+i];
                    
                    ipos+=nr;
                }
                
            }
        }
    }
    
    
    delete [] tmp;
    
    if(tot_no_types)
        delete [] type_ref;

}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam::set_fs(int no_files
,char** file_names)
{
    if(no_files!=1)
        error->abort("one file is needed for eam");
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    
    fp=NULL;
    if(atoms->my_p_no==0)
    {
        fp=fopen(file_names[0],"r");
        if(fp==NULL)
            error->abort("%s file not found",file_names[0]);
    }
    
    for(int i=0;i<4;i++)
        if(line_read(fp,line)==-1)
            error->abort("%s file ended immaturely",file_names[0]);
    
    char** arg;
    int narg;
    
    narg=mapp->parse_line(line,arg);
    if(narg<2 || narg-1 < no_types)
        error->abort("invalid line in %s file: %s",file_names[0],line);
    
    int tot_no_types=narg-1;
    
    int* type_ref;
    CREATE1D(type_ref,tot_no_types);
    for(int i=0;i<tot_no_types;i++)
        type_ref[i]=atom_types->find_type_exist(arg[i+1]);
    
    for(int i=0;i<narg;i++)
        delete [] arg[i];
    if(narg)
        delete [] arg;
    
    int* tmp_type_ref;
    CREATE1D(tmp_type_ref,no_types);
    for(int i=0;i<no_types;i++)
        tmp_type_ref[i]=-1;
    for(int ityp=0;ityp<tot_no_types;ityp++)
        if(type_ref[ityp]!=-1)
            tmp_type_ref[type_ref[ityp]]=0;
    
    for(int i=0;i<no_types;i++)
        if(tmp_type_ref[i]==-1)
            error->abort("%s file does not contain "
            "parameters for element %s",file_names[0]
            ,atom_types->atom_names[i]);
    if(no_types)
        delete [] tmp_type_ref;
    
    if(line_read(fp,line)==-1)
        error->abort("%s file ended immaturely",file_names[0]);
    narg=mapp->parse_line(line,arg);
    
    
    if(narg!=5)
        error->abort("invalid line in %s file: %s",file_names[0],line);
    
    nrho=atoi(arg[0]);
    nr=atoi(arg[2]);
    drho=atof(arg[1]);
    dr=atof(arg[3]);
    
    for(int j=0;j<narg;j++)
        delete [] arg[j];
    if(narg)
        delete [] arg;
    
    if(nrho<5)
        error->abort("nrho in %s file must be larger than 5",file_names[0]);
    if(nr<5)
        error->abort("nr in %s file must be larger than 5",file_names[0]);
    if(drho<=0.0)
        error->abort("drho in %s file must be larger than 0.0",file_names[0]);
    if(dr<=0.0)
        error->abort("dr in %s file must be larger than 0.0",file_names[0]);
    
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    
    
    
    
    
    int ipos;
    int tot=tot_no_types*nrho+tot_no_types*tot_no_types*nr
    +tot_no_types*(tot_no_types+1)*nr/2;
    int tot_loc=nrho+tot_no_types*nr;
    int ipos_loc;
    
    type0* tmp;
    CREATE1D(tmp,tot);
    
    
    ipos=0;
    
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(line_read(fp,line)==-1)
            error->abort("%s file ended immaturely",file_names[0]);
        ipos_loc=0;
        while (ipos_loc<tot_loc)
        {
            if(line_read(fp,line)==-1)
                error->abort("%s file ended immaturely",file_names[0]);
            
            narg=mapp->parse_line(line,arg);
            
            if(ipos_loc+narg>tot_loc)
                error->abort("%s file ended immaturely",file_names[0]);
            
            for(int i=0;i<narg;i++)
            {
                tmp[ipos+ipos_loc]=atof(arg[i]);
                ipos_loc++;
                delete [] arg[i];
            }
            if(narg)
                delete [] arg;
        }
        ipos+=ipos_loc;
    }
    
    tot_loc=tot_no_types*(tot_no_types+1)*nr/2;
    ipos_loc=0;
    while (ipos_loc<tot_loc)
    {
        if(line_read(fp,line)==-1)
            error->abort("%s file ended immaturely",file_names[0]);
        
        narg=mapp->parse_line(line,arg);
        
        if(ipos_loc+narg>tot_loc)
            error->abort("%s file ended immaturely",file_names[0]);
        
        for(int i=0;i<narg;i++)
        {
            tmp[ipos+ipos_loc]=atof(arg[i]);
            ipos_loc++;
            delete [] arg[i];
        }
        if(narg)
            delete [] arg;
    }
    
    
    
    delete [] line;
    if(atoms->my_p_no==0)
        fclose(fp);
    
    allocate();
    
    
    int component;
    
    ipos=0;
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(type_ref[ityp]==-1)
        {
            ipos+=nr*tot_no_types+nrho;
        }
        else
        {
            component=type_ref[ityp];
            
            for(int i=0;i<nrho;i++)
                F_arr[component][i][0]=tmp[ipos+i];
            
            ipos+=nrho;
            
            for(int jtyp=0;jtyp<tot_no_types;jtyp++)
            {
                if(type_ref[jtyp]==-1)
                {
                    ipos+=nr;
                }
                else
                {
                    component=type2rho[type_ref[ityp]][type_ref[jtyp]];

                    for(int i=0;i<nr;i++)
                        rho_arr[component][i][0]=tmp[ipos+i];
                    
                    ipos+=nr;
                }

            }
        }
    }
    
    
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(type_ref[ityp]==-1)
        {
            ipos+=nr*(ityp+1);
        }
        else
        {
            for(int jtyp=0;jtyp<ityp+1;jtyp++)
            {
                if(type_ref[jtyp]==-1)
                {
                    ipos+=nr;
                }
                else
                {
                    component=type2phi[type_ref[ityp]][type_ref[jtyp]];
                    
                    for(int i=0;i<nr;i++)
                        phi_r_arr[component][i][0]=tmp[ipos+i];
                    
                    ipos+=nr;
                }
            }
        }
    }
    
    delete [] tmp;
    if(tot_no_types)
        delete [] type_ref;
}
/*--------------------------------------------
 read line and broadcast
 --------------------------------------------*/
int ForceField_eam::line_read(FILE* file
,char*& line)
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
        fgets(line,MAXCHAR,file);
        lenght=static_cast<int>(strlen(line))+1;
    }
    MPI_Bcast(&lenght,1,MPI_INT,0,world);
    MPI_Bcast(line,lenght,MPI_CHAR,0,world);
    
    return lenght;
}
/*--------------------------------------------
 clean up the arrays
 --------------------------------------------*/
void ForceField_eam::clean_up()
{
    if(allocated==0)
        return;
    
    int no_types=atom_types->no_types;
    
    for(int i=0;i<no_types;i++)
    {
        delete [] type2phi[i];
        delete [] type2rho[i];
    }
    if(no_types)
    {
        delete [] type2phi;
        delete [] type2rho;
    }
    
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
        {
            for(int i=0;i<nr;i++)
                delete [] phi_r_arr[ityp][i];
            delete [] phi_r_arr[ityp];
        }
        delete [] phi_r_arr;
        
        for(int ityp=0;ityp<no_types;ityp++)
        {
            for(int i=0;i<nr;i++)
                delete [] rho_arr[ityp][i];
            delete [] rho_arr[ityp];
        }
        delete [] rho_arr;
        
        
        for(int ityp=0;ityp<no_types;ityp++)
        {
            for(int i=0;i<nrho;i++)
                delete [] F_arr[ityp][i];
            delete [] F_arr[ityp];
        }
        delete [] F_arr;
    }
    else if(eam_mode==FINNIS_FL)
    {
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
        {
            for(int i=0;i<nr;i++)
                delete [] phi_r_arr[ityp][i];
            delete [] phi_r_arr[ityp];
        }
        delete [] phi_r_arr;
        
        for(int ityp=0;ityp<no_types*no_types;ityp++)
        {
            for(int i=0;i<nr;i++)
                delete [] rho_arr[ityp][i];
            delete [] rho_arr[ityp];
        }
        delete [] rho_arr;
        
        for(int ityp=0;ityp<no_types;ityp++)
        {
            for(int i=0;i<nrho;i++)
                delete [] F_arr[ityp][i];
            delete [] F_arr[ityp];
        }
        delete [] F_arr;
    }
    
    allocated=0;
}
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void ForceField_eam::allocate()
{
    int no_types=atom_types->no_types;
    
    CREATE2D(type2phi,no_types,no_types);
    CREATE2D(type2rho,no_types,no_types);
    
    for(int ityp=0;ityp<no_types;ityp++)
        for(int jtyp=0;jtyp<no_types;jtyp++)
        {
            type2phi[ityp][jtyp]=COMP(ityp,jtyp);
            type2rho[ityp][jtyp]=jtyp*no_types+ityp;
        }
    
    
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        
        CREATE1D(phi_r_arr,no_types*(no_types+1)/2);
        
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
            CREATE1D(phi_r_arr[ityp],nr);
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(phi_r_arr[ityp][i],7);
        
        CREATE1D(rho_arr,no_types*no_types);
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(rho_arr[ityp],nr);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(rho_arr[type2rho[ityp][0]][i],7);
         
        
        for(int ityp=0;ityp<no_types;ityp++)
            for(int jtyp=0;jtyp<no_types;jtyp++)
                rho_arr[type2rho[ityp][jtyp]]
                =rho_arr[type2rho[ityp][0]];
        
        CREATE1D(F_arr,no_types);
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(F_arr[ityp],nrho);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<nrho;i++)
                CREATE1D(F_arr[ityp][i],7);
        
    }
    else if(eam_mode==FINNIS_FL)
    {
        
        CREATE1D(phi_r_arr,no_types*(no_types+1)/2);
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
            CREATE1D(phi_r_arr[ityp],nr);
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(phi_r_arr[ityp][i],7);

        CREATE1D(rho_arr,no_types*no_types);
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            CREATE1D(rho_arr[ityp],nr);
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(rho_arr[ityp][i],7);
        
        CREATE1D(F_arr,no_types);
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(F_arr[ityp],nrho);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<nrho;i++)
                CREATE1D(F_arr[ityp][i],7);
    }
    
    allocated=1;
}
/*--------------------------------------------
 set the arrays
 --------------------------------------------*/
void ForceField_eam::set_arrays()
{
    int no_types=atom_types->no_types;
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        for(int i=0;i<no_types;i++)
            interpolate(nrho,drho,F_arr[i]);
        
        for(int i=0;i<no_types;i++)
            interpolate(nr,dr,rho_arr[i]);
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<i+1;j++)
                interpolate(nr,dr,phi_r_arr[type2phi[i][j]]);
    }
    else if(eam_mode==FINNIS_FL)
    {
        for(int i=0;i<no_types;i++)
            interpolate(nrho,drho,F_arr[i]);
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<no_types;j++)
                interpolate(nr,dr,rho_arr[type2rho[i][j]]);
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<i+1;j++)
                interpolate(nr,dr,phi_r_arr[type2phi[i][j]]);
    }
}
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void ForceField_eam::interpolate(int n,type0 delta
,type0** spline)
{
    spline[0][1]=spline[1][0]-spline[0][0];
    spline[1][1]=0.5*(spline[2][0]-spline[0][0]);
    spline[n-2][1]=0.5*(spline[n-1][0]-spline[n-3][0]);
    spline[n-1][1]=spline[n-1][0]-spline[n-2][0];
    
    for (int i=2;i<n-2;i++)
        spline[i][1]=((spline[i-2][0]-spline[i+2][0])+
            8.0*(spline[i+1][0]-spline[i-1][0]))/12.0;
    
    for (int i=0;i<n-1;i++)
    {
        spline[i][2]=3.0*(spline[i+1][0]-spline[i][0])-
        2.0*spline[i][1]-spline[i+1][1];
        spline[i][3]=spline[i][1]+spline[i+1][1]-
        2.0*(spline[i+1][0]-spline[i][0]);
    }
    
    spline[n-1][2]=0.0;
    spline[n-1][3]=0.0;
    for(int i=0;i<n;i++)
    {
        spline[i][4]=spline[i][1]/delta;
        spline[i][5]=2.0*spline[i][2]/delta;
        spline[i][6]=3.0*spline[i][3]/delta;
    }
        
}
/*--------------------------------------------
 interpolate for a single array
 --------------------------------------------*/
type0 ForceField_eam::interpolate(type0* arr
,int n,type0 p,int k)
{
    
    // k_min=0 k_max=n-1
    type0 coef0,coef1,coef2,coef3,tmp;
    
    coef0=arr[k];
    
    if(k==0)
    {
        coef1=arr[1]-arr[0];
        tmp=0.5*(arr[2]-arr[0]);
        coef2=(-arr[0]+2.0*arr[1]-arr[2])/2.0;
        coef3=(arr[0]-2.0*arr[1]+arr[2])/2.0;
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }
    else if(k==1)
    {
        coef1=0.5*(arr[2]-arr[0]);
        tmp=((arr[0]-arr[4])+8.0*(arr[3]-arr[1]))/12.0;
        coef2=(11.0*arr[0]-28.0*arr[1]+24.0*arr[2]-8.0*arr[3]+arr[4])/12.0;
        coef3=(-5.0*arr[0]+16.0*arr[1]-18.0*arr[2]+8.0*arr[3]-arr[4])/12.0;
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }
    else if(k==n-2)
    {
        coef1=0.5*(arr[n-1]-arr[n-3]);
        tmp=arr[n-1]-arr[n-2];
        coef2=arr[n-3]-3.0*arr[n-2]+2.0*arr[n-1];
        coef3=0.5*(-arr[n-3]+3.0*arr[n-2]-2.0*arr[n-1]);
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }
    else if (k==n-1)
    {
        coef1=arr[n-1]-arr[n-2];
        return coef0+coef1*p;
    }
    else
    {
        coef1=((arr[k-2]-arr[k+2])+
               8.0*(arr[k+1]-arr[k-1]))/12.0;
        tmp=((arr[k-1]-arr[k+3])+
               8.0*(arr[k+2]-arr[k]))/12.0;
        
        coef2=(-2.0*arr[k-2]+15.0*arr[k-1]-28.0*arr[k]+21.0*arr[k+1]-6.0*arr[k+2])/12.0;
        coef3=(arr[k-2]-7.0*arr[k-1]+16.0*arr[k]-17.0*arr[k+1]+7.0*arr[k+2])/12.0;
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }

    
}
