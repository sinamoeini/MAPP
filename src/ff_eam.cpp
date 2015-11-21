#include <stdlib.h>
#include "neighbor.h"
#include "ff_eam.h"
#include "atom_types.h"
#include "error.h"
#include "memory.h"
using namespace MAPP_NS;
enum{NOT_SET,FUNC_FL,SET_FL,FINNIS_FL};

/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam::
ForceField_eam(MAPP* mapp):ForceFieldMD(mapp)
{
    if(mapp->mode!=MD_mode)
        error->abort("ff eam works only "
        "for md mode");

    max_pairs=0;
    eam_reader=NULL;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_eam::~ForceField_eam()
{
    delete eam_reader;
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
    
    type0* x=mapp->x->begin();
    type0* f=mapp->f->begin();
    type0* rho=rho_ptr->begin();
    md_type* type=mapp->type->begin();
    
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
        for(int i=1;i<7;i++)
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
            if(rsq < cut_sq[COMP(itype,jtype)])
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
                if(jatm<natms)
                    rho[jatm]+=rho_j;
                    
                fpair=-dphi_dr*r_inv;
                
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
                    phi*=0.5;
                }
                
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
                
                drhoi_dr[istart]=-drho_i_dr*r_inv;
                drhoj_dr[istart]=-drho_j_dr*r_inv;
            }
            
            istart++;
        }
        p=rho[iatm]*drho_inv;
        m=static_cast<int> (p);
        m=MIN(m,nrho-2);
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
                }
                
                if(jatm>=natms)
                    fpair*=0.5;
                
                if(st_clc)
                {
                    nrgy_strss[1]-=fpair*dx0*dx0;
                    nrgy_strss[2]-=fpair*dx1*dx1;
                    nrgy_strss[3]-=fpair*dx2*dx2;
                    nrgy_strss[4]-=fpair*dx1*dx2;
                    nrgy_strss[5]-=fpair*dx2*dx0;
                    nrgy_strss[6]-=fpair*dx0*dx1;
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
    type0* x=mapp->x->begin();
    type0* rho=rho_ptr->begin();
    md_type* type=mapp->type->begin();
    
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
            
            if(rsq<cut_sq[COMP(itype,jtype)])
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
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        itype=type[iatm];
        coef=F_arr[itype][m];
        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho[iatm]>rho_max)
        {
            //tmp0+=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv*(rho[iatm]-rho_max);
            tmp0+=((coef[6]*p+coef[5])*p+coef[4])*(rho[iatm]-rho_max);
        }
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
    neighbor->pair_wise=1;
    
    rho_ptr=new Vec<type0>(atoms,1);
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
    
    delete rho_ptr;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ForceField_eam::coef(int nargs,char** args)
{
    if (nargs<3)
        error->abort("ff_coef for ff eam "
        "should at least have 2 arguments");
    
    cut_off_alloc();
    
    delete eam_reader;
    eam_reader=new EAMFileReader(mapp);
    eam_reader->file_format(args[1]);
    int iarg=2;
    while(iarg<nargs)
    {
        eam_reader->add_file(args[iarg],iarg-2);
        iarg++;
    }
    eam_reader->setup();
    
    nr=eam_reader->nr;
    nrho=eam_reader->nrho;
    dr=eam_reader->dr;
    drho=eam_reader->drho;
    dr_inv=eam_reader->dr_inv;
    drho_inv=eam_reader->drho_inv;
    rho_max=eam_reader->rho_max;
    
    F_arr=eam_reader->F_arr;
    phi_r_arr=eam_reader->phi_r_arr;
    rho_arr=eam_reader->rho_arr;
    
    type2rho=eam_reader->type2rho;
    type2phi=eam_reader->type2phi;
    
    int no_types=atom_types->no_types;
    memcpy(cut_sq,eam_reader->cut_sq,(no_types*(no_types+1)/2)*sizeof(type0));
}
 