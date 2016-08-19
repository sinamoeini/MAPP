#include <stdlib.h>
#include "neighbor.h"
#include "ff_eam.h"
#include "atom_types.h"
#include "error.h"
#include "memory.h"
#include "cmd.h"
using namespace MAPP_NS;
enum{NOT_SET,FUNC_FL,SET_FL,FINNIS_FL};

/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam::
ForceField_eam(MAPP* mapp):ForceFieldMD(mapp)
{
    gcmc_n_cutoff=2;
    gcmc_n_vars=2;
    gcmc_tag_enabled=true;

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
    }
    
    type0* x=mapp->x->begin();
    type0* fvec=f->begin();
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
    
    nrgy_strss_lcl[0]=0.0;
    if(st_clc)
        for(int i=1;i<7;i++)
            nrgy_strss_lcl[i]=0.0;
    
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
            if(rsq < cut_sq[itype][jtype])
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
                    phi*=0.5;
                }
                
                nrgy_strss_lcl[0]+=phi;
                
                if(st_clc)
                {
                    nrgy_strss_lcl[1]-=fpair*dx0*dx0;
                    nrgy_strss_lcl[2]-=fpair*dx1*dx1;
                    nrgy_strss_lcl[3]-=fpair*dx2*dx2;
                    nrgy_strss_lcl[4]-=fpair*dx1*dx2;
                    nrgy_strss_lcl[5]-=fpair*dx2*dx0;
                    nrgy_strss_lcl[6]-=fpair*dx0*dx1;
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
        nrgy_strss_lcl[0]+=tmp0;
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
                
                fvec[icomp]+=dx0*fpair;
                fvec[icomp+1]+=dx1*fpair;
                fvec[icomp+2]+=dx2*fpair;
                
                if(jatm<natms)
                {
                    fvec[jcomp]-=dx0*fpair;
                    fvec[jcomp+1]-=dx1*fpair;
                    fvec[jcomp+2]-=dx2*fpair;
                }
                
                if(jatm>=natms)
                    fpair*=0.5;
                
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
            
            if(rsq<cut_sq[itype][jtype])
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
            tmp0+=((coef[6]*p+coef[5])*p+coef[4])*(rho[iatm]-rho_max);
        en+=tmp0;

    }

    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    return en_tot;
}
/*--------------------------------------------
 calculate the exchange energy due to deletion
 or insertion
 --------------------------------------------*/
void ForceField_eam::xchng_energy(GCMC* gcmc)
{
    
    type0 en=0.0,en_tot=0.0;
    
    int& iatm=gcmc->iatm;
    md_type& itype=gcmc->itype;
    
    int& jatm=gcmc->jatm;
    md_type& jtype=gcmc->jtype;
    
    type0&rsq=gcmc->rsq;
    
    type0* rho=rho_ptr->begin();
    type0* F=F_ptr->begin();
    type0* rho_xchng=rho_xchng_ptr->begin();
    type0* F_xchng=F_xchng_ptr->begin();
    type0 rho_iatm=0.0,F_iatm=0.0,rho_iatm_lcl=0.0;
    type0 r,p;
    int m;
    type0* coef;
    int natms=atoms->natms;

    for(gcmc->reset_iatm();iatm!=-1;gcmc->next_iatm())
        for(gcmc->reset_jatm();jatm!=-1;gcmc->next_jatm())
        {
            r=sqrt(rsq);
            
            p=r*dr_inv;
            m=static_cast<int>(p);
            m=MIN(m,nr-2);
            p-=m;
            p=MIN(p,1.0);
            
            coef=phi_r_arr[type2phi[itype][jtype]][m];
            
            if(jatm<natms)
            {
                en+=(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0])/r;
                coef=rho_arr[type2rho[itype][jtype]][m];
                rho_xchng[jatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
            }
            else
                en+=0.5*(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0])/r;
            
            coef=rho_arr[type2rho[jtype][itype]][m];
            rho_iatm_lcl+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        }
    
    MPI_Allreduce(&rho_iatm_lcl,&rho_iatm,1,MPI_TYPE0,MPI_SUM,world);
    
    if(gcmc->im_root)
    {
        p=rho_iatm*drho_inv;
        m=static_cast<int>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[itype][m];
        F_iatm=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho_iatm>rho_max)
            F_iatm+=((coef[6]*p+coef[5])*p+coef[4])*(rho_iatm-rho_max);
        
        en+=F_iatm;
    }
    
    md_type* type=mapp->type->begin();
    type0 tmp0,c0=1.0,en0=0.0;
    if(gcmc->xchng_mode==DEL_MODE) c0=-1.0;
    
    for(int i=0;i<natms;i++)
        if(rho_xchng[i])
        {
            tmp0=rho[i]+c0*rho_xchng[i];
            p=tmp0*drho_inv;
            m=static_cast<int>(p);
            m=MIN(m,nrho-2);
            p-=m;
            p=MIN(p,1.0);
            coef=F_arr[type[i]][m];
            F_xchng[i]=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
            if(tmp0>rho_max)
                F_xchng[i]+=((coef[6]*p+coef[5])*p+coef[4])*(tmp0-rho_max);
            en0+=F_xchng[i]-F[i];
        }
    
    if(gcmc->xchng_mode==DEL_MODE)
        en0*=-1.0;

    en+=en0;

    
    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    
    bool succ=gcmc->decide(en_tot);
    
    if(succ)
    {
        natms=atoms->natms;
        rho=rho_ptr->begin();
        F=F_ptr->begin();
        
        rho_xchng=rho_xchng_ptr->begin();
        F_xchng=F_xchng_ptr->begin();
        for(int i=0;i<natms;i++)
            if(rho_xchng[i])
            {
                rho[i]+=c0*rho_xchng[i];
                F[i]=F_xchng[i];
            }
        
        if(gcmc->im_root && gcmc->xchng_mode==INS_MODE)
        {
            rho[natms-1]=rho_iatm;
            F[natms-1]=F_iatm;
            rho_xchng[natms-1]=0.0;
            F_xchng[natms-1]=0.0;
        }
    }
    else
    {
        natms=atoms->natms;
        rho_xchng=rho_xchng_ptr->begin();
        F_xchng=F_xchng_ptr->begin();
        for(int i=0;i<natms;i++)
            rho_xchng[i]=F_xchng[i]=0.0;
    }
    
}
/*--------------------------------------------
 pre gcmc energy claculate the increase or
 decrease in electron density
 --------------------------------------------*/
void ForceField_eam::pre_gcmc_energy(GCMC* gcmc)
{
    int& icomm=gcmc->icomm;
    type0 en;
    type0 rho_iatm_lcl;
    
    int& iatm=gcmc->iatm;
    md_type& itype=gcmc->itype;
    
    int& jatm=gcmc->jatm;
    md_type& jtype=gcmc->jtype;
    
    type0&rsq=gcmc->rsq;
    
    type0 r,p;
    int m;
    type0* coef;
    type0 c0=1.0,en0=0.0;
    
    type0* rho=rho_ptr->begin();
    type0* rho_xchng=rho_xchng_ptr->begin();
    int natms=atoms->natms;
    int* tag=gcmc->tag_vec_p->begin();
    for(int i=0;i<natms;i++)
        rho_xchng[i]=rho[i];
    
    for(gcmc->reset_icomm();icomm!=-1;gcmc->next_icomm())
    {
        c0=1.0;
        if(gcmc->xchng_mode==DEL_MODE)
            c0=-1.0;
        
        en=rho_iatm_lcl=0.0;
        for(gcmc->reset_iatm();iatm!=-1;gcmc->next_iatm())
            for(gcmc->reset_jatm();jatm!=-1;gcmc->next_jatm())
            {
                r=sqrt(rsq);
                
                p=r*dr_inv;
                m=static_cast<int>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);
                
                coef=phi_r_arr[type2phi[itype][jtype]][m];
                
                if(jatm<natms)
                {
                    en+=(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0])/r;
                    coef=rho_arr[type2rho[itype][jtype]][m];
                    rho_xchng[jatm]+=c0*(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0]);
                }
                else
                    en+=0.5*(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0])/r;
                
                
                coef=rho_arr[type2rho[jtype][itype]][m];
                rho_iatm_lcl+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
            }
        
        
        
        type0* F=F_ptr->begin();
        type0* F_xchng=F_xchng_ptr->begin();
        md_type* type=mapp->type->begin();
        
        type0 tmp0;
        en0=0.0;
        for(int i=0;i<natms;i++)
            if(tag[i]==icomm)
            {
                tmp0=rho_xchng[i];
                p=tmp0*drho_inv;
                m=static_cast<int>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                coef=F_arr[type[i]][m];
                F_xchng[i]=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(tmp0>rho_max)
                    F_xchng[i]+=((coef[6]*p+coef[5])*p+coef[4])*(tmp0-rho_max);
                en0+=F_xchng[i]-F[i];
            }
        
        en+=en0*c0;
        
        gcmc->lcl_vars[0]=en;
        gcmc->lcl_vars[1]=rho_iatm_lcl;
        
    }

}
/*--------------------------------------------
 calculate the energy if I am root if not
 pass the lcl variables
 --------------------------------------------*/
type0 ForceField_eam::gcmc_energy(GCMC* gcmc)
{
    
    int& icomm=gcmc->icomm;
    for(gcmc->reset_icomm();icomm!=-1;gcmc->next_icomm())
        MPI_Reduce(gcmc->lcl_vars,gcmc->vars,2,MPI_TYPE0,MPI_SUM,gcmc->curr_root,*gcmc->curr_comm);
 
    
    if(gcmc->im_root)
    {
        //restart the comms
        gcmc->reset_icomm();
        type0 rho_iatm,en;
        int m;
        type0* coef;
        
        en=gcmc->vars[0];
        rho_iatm=gcmc->vars[1];
        type0 p=rho_iatm*drho_inv;
        m=static_cast<int>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        coef=F_arr[gcmc->itype][m];
        type0 F_iatm=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho_iatm>rho_max)
            F_iatm+=((coef[6]*p+coef[5])*p+coef[4])*(rho_iatm-rho_max);
        
        en+=F_iatm;
        return en;
    }
    
    return 0.0;
}
/*--------------------------------------------
 calculate the energy if I am root if not
 pass the lcl variables
 --------------------------------------------*/
void ForceField_eam::post_gcmc_energy(GCMC* gcmc)
{
    int* tag=gcmc->tag_vec_p->begin();
    type0* rho=rho_ptr->begin();
    type0* F=F_ptr->begin();
    
    type0* rho_xchng=rho_xchng_ptr->begin();
    type0* F_xchng=F_xchng_ptr->begin();
    int natms=atoms->natms;
    for(int i=0;i<natms;i++)
        if(tag[i]==0)
        {
            rho[i]=rho_xchng[i];
            F[i]=F_xchng[i];
        }
    
    if(gcmc->im_root && gcmc->xchng_mode==INS_MODE && gcmc->root_succ)
    {
        rho[natms-1]=gcmc->vars[1];
        type0 p=rho[natms-1]*drho_inv;
        int m=static_cast<int>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        type0* coef=F_arr[gcmc->itype][m];
        F[natms-1]=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho[natms-1]>rho_max)
            F[natms-1]+=((coef[6]*p+coef[5])*p+coef[4])*(rho[natms-1]-rho_max);
    }

}
/*--------------------------------------------
 init before running
 --------------------------------------------*/
void ForceField_eam::init()
{
    neighbor->pair_wise=true;
    
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
 init xchng
 --------------------------------------------*/
void ForceField_eam::init_xchng()
{
    F_ptr=new Vec<type0>(atoms,1);
    rho_xchng_ptr=new Vec<type0>(atoms,1);
    F_xchng_ptr=new Vec<type0>(atoms,1);
    
    type0* x=mapp->x->begin();
    type0* rho=rho_ptr->begin();
    type0* F=F_ptr->begin();
    type0* rho_xchng=rho_xchng_ptr->begin();
    type0* F_xchng=F_xchng_ptr->begin();
    md_type* type=mapp->type->begin();
    
    int iatm,jatm;
    
    int itype,jtype,icomp,jcomp;
    type0 dx0,dx1,dx2,rsq;
    type0 r,p;
    int m;
    type0* coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    int natms=atoms->natms;
    
    for(iatm=0;iatm<natms;iatm++) rho[iatm]=rho_xchng[iatm]=F_xchng[iatm]=0.0;
    
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
            
            if(rsq<cut_sq[itype][jtype])
            {
                r=sqrt(rsq);
                
                p=r*dr_inv;
                m=static_cast<int>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);
                
                coef=phi_r_arr[type2phi[itype][jtype]][m];
                
                coef=rho_arr[type2rho[jtype][itype]][m];
                rho[iatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(jatm<natms)
                {
                    coef=rho_arr[type2rho[itype][jtype]][m];
                    rho[jatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                }
            }
        }
        
        p=rho[iatm]*drho_inv;
        m=static_cast<int>(p);
        m=MIN(m,nrho-2);
        p-=m;
        p=MIN(p,1.0);
        itype=type[iatm];
        coef=F_arr[itype][m];
        F[iatm]=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
    }
}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceField_eam::fin_xchng()
{
    delete F_xchng_ptr;
    delete rho_xchng_ptr;
    delete F_ptr;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ForceField_eam::coef(int nargs,char** args)
{
    
    char* file_format;
    char** files;
    int nfiles;
    
    Pattern cmd(error);
    add_var(atom_types->no_types,"no. of elements in simulation");
    
    /*----------------------------*/
    cmd.cmd("ff_coef");
    cmd.add_var(file_format,"file_format");
    cmd.add_vdesc(0,"defines the format of EAM files");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=vlogic("eq","FS")
    +vlogic("eq","SetFL")+vlogic("eq","FuncFL");
    /*------------------------------------------------------------------------------------*/

    /*----------------------------*/
    cmd.cmd_dy("files",files,nfiles);
    cmd.add_vdesc(0,"defines the path to the potential file/s");
    /*--------------------------------------------------------*/
    cmd.add_clog()=logic(file_format,"eq","FuncFL")/logic(nfiles,"eq",atom_types->no_types);
    cmd.add_clog().link()=(logic(file_format,"eq","FS")+logic(file_format,"eq","SetFL"))
    /logic(nfiles,"eq",1);
    /*------------------------------------------------------------------------------------*/
    
    cmd.scan(args,nargs);

    cut_off_alloc();
    
    delete eam_reader;
    eam_reader=new EAMFileReader(mapp);
    eam_reader->file_format(file_format);
    for(int i=0;i<nfiles;i++)
        eam_reader->add_file(files[i],i);
    
    setup();
}
/*--------------------------------------------
 setup
 --------------------------------------------*/
void ForceField_eam::setup()
{
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
    for(int i=0;i<no_types;i++)
        for(int j=0;j<no_types;j++)
        {
            cut[i][j]=eam_reader->cut[COMP(i,j)];
            cut_sq[i][j]=eam_reader->cut_sq[COMP(i,j)];
        }
}

 