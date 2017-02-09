#include <stdlib.h>
#include "eam_file_reader.h"
#include "ff_eam_dmd.h"
#include "neighbor_dmd.h"
#include "error.h"
#include "memory.h"
#include "atom_types.h"
#include "xmath.h"
#include "cmd.h"
#include "atoms.h"
#include "MAPP.h"
#include "dynamic.h"
#include <limits>
#define PI_IN_SQ 0.564189583547756286948079451561
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam_dmd::
ForceField_eam_dmd():ForceFieldDMD()
{    
    if(mode!=DMD_mode)
        Error::abort("ff eam_dmd works only "
        "for dmd mode");

    eam_reader=NULL;
    no_types=0;
    no_i=0;
    max_pairs=0;
    phi_psi_sz_sz=0;
    n_phi_psi=0;
    M_N_sz_sz=0;
    kbT=beta=-1.0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_eam_dmd::~ForceField_eam_dmd()
{
    deallocate();
    delete eam_reader;

    if(no_i)
    {
        delete [] xi;
        delete [] wi_0;
        delete [] wi_1;
        delete [] wi_2;
    }
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceField_eam_dmd::
force_calc(bool st_clc)
{
    if(!dynamic_flag)
        return force_calc_static(st_clc);
    
    type0* xvec=x->begin();
    atom_type* typevec=type->begin();
    type0* fvec=f->begin();
    type0* E=E_ptr->begin();
    type0* dE=dE_ptr->begin();
    type0* mu=mu_ptr->begin();;
    type0* cvec=c->begin();
    
    int iatm,jatm;
    
    int icomp,jcomp,istart;
    type0 dx0,dx1,dx2,rsq;
    type0 r=0.0,r_inv=0.0;
    int m;
    type0* coef;
    type0 tmp0,tmp1;
    type0 fpair,apair;
    
    type0 p,c_iv;

    
    if(st_clc)
        for(int i=1;i<7;i++)
            nrgy_strss_lcl[i]=0.0;
    
    for(int i=0;i<natms*c_dim;i++) E[i]=0.0;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    if(max_pairs<neighbor->no_pairs)
    {
        if(max_pairs)
        {
            delete [] rho_phi;
            delete [] drho_phi_dr;
            delete [] drho_phi_dalpha;
        }
        max_pairs=neighbor->no_pairs;
        int no_0=max_pairs*stride;
        CREATE1D(rho_phi,no_0);
        CREATE1D(drho_phi_dr,no_0);
        CREATE1D(drho_phi_dalpha,no_0);
    }
    
    for(int i=0;i<max_pairs*stride;i++)
        rho_phi[i]=drho_phi_dr[i]=drho_phi_dalpha[i]=0.0;
    
    nrgy_strss_lcl[0]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+c_dim)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+c_dim)*jatm;
            dx0=xvec[icomp]-xvec[jcomp];
            dx1=xvec[icomp+1]-xvec[jcomp+1];
            dx2=xvec[icomp+2]-xvec[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            if(rsq<cut_sq_mod_0)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                
                for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && cvec[ic_dim]>=0.0;ic_dim++)
                {
                    itype=typevec[ic_dim];
                    for(int jc_dim=jatm*c_dim,jtype;jc_dim<(jatm+1)*c_dim && cvec[jc_dim]>=0.0;jc_dim++)
                    {
                        jtype=typevec[jc_dim];
                        
                        if(rsq>cut_sq[itype][jtype])
                            continue;
                        
                        type0 alpha_ij=sqrt(xvec[3*iatm+ic_dim+3]*xvec[3*iatm+ic_dim+3]+xvec[3*jatm+jc_dim+3]*xvec[3*jatm+jc_dim+3]);
                        
                        if(alpha_ij>=alpha_max)
                            continue;
                        
                        type0 upper=(r+cut[itype][jtype])/alpha_ij;
                        type0 lower=(r-cut[itype][jtype])/alpha_ij;
                        
                        if(lower>=xi[no_i-1])
                            continue;
                        
                        int rho_ij=istart+type2rho_pair_ij[itype][jtype];
                        int rho_ji=istart+type2rho_pair_ji[jtype][itype];
                        int phi_ij=istart+type2phi_pair_ij[itype][jtype];
                        
                        rho_phi[rho_ij]=0.0;
                        drho_phi_dr[rho_ij]=0.0;
                        drho_phi_dalpha[rho_ij]=0.0;
                        
                        rho_phi[rho_ji]=0.0;
                        drho_phi_dr[rho_ji]=0.0;
                        drho_phi_dalpha[rho_ji]=0.0;
                        
                        rho_phi[phi_ij]=0.0;
                        drho_phi_dr[phi_ij]=0.0;
                        drho_phi_dalpha[phi_ij]=0.0;
                        
                        int m;
                        type0 rtmp;
                        type0 p,tmp0,tmp1;
                        type0* coef;
                        for(int i=0;i<no_i;i++)
                        {
                            if(xi[i]<=lower && xi[i]>=upper)
                                continue;
                            rtmp=r-xi[i]*alpha_ij;
                            
                            p=fabs(rtmp)*dr_inv;
                            m=static_cast<int>(p);
                            m=MIN(m,nr-2);
                            p-=m;
                            p=MIN(p,1.0);
                            
                            coef=rho_r_arr[type2rho[itype][jtype]][m];
                            tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                            tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                            if(rtmp<0.0)
                                tmp0*=-1.0;
                            
                            rho_phi[rho_ij]+=wi_0[i]*tmp0;
                            drho_phi_dr[rho_ij]+=wi_0[i]*tmp1;
                            drho_phi_dalpha[rho_ij]+=wi_1[i]*tmp1;
                            
                            if(itype!=jtype)
                            {
                                coef=rho_r_arr[type2rho[jtype][itype]][m];
                                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                                tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                                if(rtmp<0.0)
                                    tmp0*=-1.0;
                                
                                rho_phi[rho_ji]+=wi_0[i]*tmp0;
                                drho_phi_dr[rho_ji]+=wi_0[i]*tmp1;
                                drho_phi_dalpha[rho_ji]+=wi_1[i]*tmp1;
                            }
                            
                            
                            
                            coef=phi_r_arr[type2phi[itype][jtype]][m];
                            tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                            tmp1=((3.0*coef[3]*p+2.0*coef[2])*p+coef[1])*dr_inv;
                            if(rtmp<0.0)
                                tmp0*=-1.0;
                            
                            rho_phi[phi_ij]+=wi_0[i]*tmp0;
                            drho_phi_dr[phi_ij]+=wi_0[i]*tmp1;
                            drho_phi_dalpha[phi_ij]+=wi_1[i]*tmp1;
                        }
                        
                        tmp0=PI_IN_SQ*r_inv;
                        
                        rho_phi[rho_ij]*=tmp0;
                        drho_phi_dr[rho_ij]*=tmp0;
                        drho_phi_dr[rho_ij]-=rho_phi[rho_ij]*r_inv;
                        drho_phi_dr[rho_ij]*=r_inv;
                        drho_phi_dalpha[rho_ij]*=-tmp0/alpha_ij;
                        
                        if(itype!=jtype)
                        {
                            rho_phi[rho_ji]*=tmp0;
                            drho_phi_dr[rho_ji]*=tmp0;
                            drho_phi_dr[rho_ji]-=rho_phi[rho_ji]*r_inv;
                            drho_phi_dr[rho_ji]*=r_inv;
                            drho_phi_dalpha[rho_ji]*=-tmp0/alpha_ij;
                        }
                        
                        rho_phi[phi_ij]*=tmp0;
                        drho_phi_dr[phi_ij]*=tmp0;
                        drho_phi_dr[phi_ij]-=rho_phi[phi_ij]*r_inv;
                        drho_phi_dr[phi_ij]*=r_inv;
                        drho_phi_dalpha[phi_ij]*=-tmp0/alpha_ij;
                        
                        E[ic_dim]+=cvec[jc_dim]*rho_phi[rho_ji];
                        
                        if(jatm<natms)
                        {
                            E[jc_dim]+=cvec[ic_dim]*rho_phi[rho_ij];
                            nrgy_strss_lcl[0]+=cvec[ic_dim]*cvec[jc_dim]*rho_phi[phi_ij];
                        }
                        else
                            nrgy_strss_lcl[0]+=0.5*cvec[ic_dim]*cvec[jc_dim]*rho_phi[phi_ij];
                    }
                }
            }
            
            
            istart+=stride;
        }
        //printf("%e %e %d\n",E[0],nrgy_strss_lcl[0],neighbor_list_size[0]);
        
        c_iv=1.0;
        for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && cvec[ic_dim]>=0.0;ic_dim++)
        {
            
            itype=typevec[ic_dim];
            c_iv-=cvec[ic_dim];
            p=E[ic_dim]*drho_inv;
            
            m=static_cast<int>(p);
            
            
            m=MIN(m,nrho-2);
            p-=m;
            p=MIN(p,1.0);
            coef=F_arr[itype][m];
            
            tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
            tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
            
            
            if(E[ic_dim]>rho_max)
                tmp0+=tmp1*(E[ic_dim]-rho_max);
            
            E[ic_dim]=tmp0;
            dE[ic_dim]=tmp1;
            
            mu[ic_dim]=tmp0;
            if(cvec[ic_dim]!=0.0)
            {
                nrgy_strss_lcl[0]+=cvec[ic_dim]*(tmp0+c_0[itype]-3.0*kbT*log(xvec[3*iatm+ic_dim+3]));
            }
            
            
            nrgy_strss_lcl[0]+=kbT*calc_ent(cvec[ic_dim]);
        }
        nrgy_strss_lcl[0]+=kbT*calc_ent(c_iv);
    }
    
    dynamic->update(dE_ptr);

    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+c_dim)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+c_dim)*jatm;
            
            dx0=xvec[icomp]-xvec[jcomp];
            dx1=xvec[icomp+1]-xvec[jcomp+1];
            dx2=xvec[icomp+2]-xvec[jcomp+2];
            r_inv=1.0/sqrt(dx0*dx0+dx1*dx1+dx2*dx2);
            
            for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && cvec[ic_dim]>=0.0;ic_dim++)
            {
                itype=typevec[ic_dim];
                for(int jc_dim=jatm*c_dim,jtype;jc_dim<(jatm+1)*c_dim && cvec[jc_dim]>=0.0;jc_dim++)
                {
                    jtype=typevec[jc_dim];
                    
                    fpair=-(drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]*dE[ic_dim]
                    +drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]*dE[jc_dim]
                    +drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]])*cvec[ic_dim]*cvec[jc_dim];
                    
                    apair=-(drho_phi_dalpha[istart+type2rho_pair_ji[jtype][itype]]*dE[ic_dim]
                    +drho_phi_dalpha[istart+type2rho_pair_ij[itype][jtype]]*dE[jc_dim]
                    +drho_phi_dalpha[istart+type2phi_pair_ij[itype][jtype]])*cvec[ic_dim]*cvec[jc_dim];
                    
                    mu[ic_dim]+=cvec[jc_dim]*(rho_phi[istart+type2phi_pair_ij[itype][jtype]]+
                    rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jc_dim]);
                    
                    if(jatm<natms)
                        mu[jc_dim]+=cvec[ic_dim]*(rho_phi[istart+type2phi_pair_ji[jtype][itype]]+
                        rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[ic_dim]);
                    
                    if(apair==0.0 && fpair==0.0)
                        continue;
                
                    
                    fvec[icomp]+=dx0*fpair;
                    fvec[icomp+1]+=dx1*fpair;
                    fvec[icomp+2]+=dx2*fpair;
                    fvec[3*iatm+ic_dim+3]+=apair*xvec[3*iatm+ic_dim+3];
                    if(jatm<natms)
                    {
                        fvec[jcomp]-=dx0*fpair;
                        fvec[jcomp+1]-=dx1*fpair;
                        fvec[jcomp+2]-=dx2*fpair;
                        fvec[3*jatm+jc_dim+3]+=apair*xvec[3*jatm+jc_dim+3];
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
            }
            istart+=stride;
        }
    }
    
    type0* alpha=x->begin()+__dim__;
    type0* fvec_alpha=f->begin()+__dim__;
    cvec=c->begin();
    for(iatm=0;iatm<natms;iatm++)
    {
        for(int i=0;i<c_dim;i++)
            if(cvec[i]>0.0)
                fvec_alpha[i]+=3.0*kbT*cvec[i]/alpha[i];
        
        fvec_alpha+=x_dim;
        alpha+=x_dim;
        cvec+=c_dim;
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
type0 ForceField_eam_dmd::energy_calc()
{
    type0 en=0.0,en_tot;
    type0* xvec=x->begin();
    atom_type* typevec=type->begin();
    type0* E=E_ptr->begin();
    type0* cvec=c->begin();
    
    int iatm,jatm;
    
    int icomp,jcomp;
    type0 dx0,dx1,dx2,rsq;
    type0 alpha;
    type0 r=0.0,r_inv=0.0;
    int m;
    type0* coef;
    type0 tmp0,tmp1,rtmp,c_iv;
    type0 upper,lower;
    
    type0 p;
    type0 phi_it_jt_0;
    type0 rho_it_jt_0;
    type0 rho_jt_it_0;
    type0 a0;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    for(int i=0;i<natms*c_dim;i++) E[i]=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+c_dim)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+c_dim)*jatm;
            dx0=xvec[icomp]-xvec[jcomp];
            dx1=xvec[icomp+1]-xvec[jcomp+1];
            dx2=xvec[icomp+2]-xvec[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            if(rsq<cut_sq_mod_0)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                
                for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && cvec[ic_dim]>=0.0;ic_dim++)
                {
                    itype=typevec[ic_dim];
                    for(int jc_dim=jatm*c_dim,jtype;jc_dim<(jatm+1)*c_dim && cvec[jc_dim]>=0.0;jc_dim++)
                    {
                        jtype=typevec[jc_dim];
                        
                        if(rsq>cut_sq[itype][jtype])
                            continue;
                        
                        alpha=sqrt(xvec[3*iatm+ic_dim+3]*xvec[3*iatm+ic_dim+3]+xvec[3*jatm+jc_dim+3]*xvec[3*jatm+jc_dim+3]);
                        
                        if(alpha>=alpha_max)
                            continue;
                        
                        upper=(r+cut[itype][jtype])/alpha;
                        lower=(r-cut[itype][jtype])/alpha;
                        
                        if(lower>=xi[no_i-1])
                            continue;

                        phi_it_jt_0=rho_it_jt_0=rho_jt_it_0=0.0;
                        
                        for(int i=0;i<no_i;i++)
                        {
                            if(xi[i]<=lower && xi[i]>=upper)
                                continue;
                            
                            rtmp=r-xi[i]*alpha;
                            p=fabs(rtmp)*dr_inv;
                            m=static_cast<int>(p);
                            m=MIN(m,nr-2);
                            p-=m;
                            p=MIN(p,1.0);
                            
                            coef=rho_r_arr[type2rho[itype][jtype]][m];
                            tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                            if(rtmp<0.0)
                                tmp0*=-1.0;
                            
                            rho_it_jt_0+=wi_0[i]*tmp0;
                            
                            if(itype!=jtype)
                            {
                                coef=rho_r_arr[type2rho[jtype][itype]][m];
                                tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                                if(rtmp<0.0)
                                    tmp0*=-1.0;
                            }
                            
                            rho_jt_it_0+=wi_0[i]*tmp0;
                            
                            coef=phi_r_arr[type2phi[itype][jtype]][m];
                            tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                            if(rtmp<0.0)
                                tmp0*=-1.0;
                            
                            phi_it_jt_0+=wi_0[i]*tmp0;
                        }
                        
                        a0=PI_IN_SQ*r_inv;
                        rho_it_jt_0*=a0;
                        rho_jt_it_0*=a0;
                        phi_it_jt_0*=a0;
                        
                        E[ic_dim]+=cvec[jc_dim]*rho_jt_it_0;
                        
                        if(jatm<natms)
                        {
                            E[jc_dim]+=cvec[ic_dim]*rho_it_jt_0;
                            en+=cvec[ic_dim]*cvec[jc_dim]*phi_it_jt_0;
                        }
                        else
                            en+=0.5*cvec[ic_dim]*cvec[jc_dim]*phi_it_jt_0;
                    }
                }
            }
        }
        
        c_iv=1.0;
        for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && cvec[ic_dim]>=0.0;ic_dim++)
        {
            itype=typevec[ic_dim];
            c_iv-=cvec[ic_dim];
            p=E[ic_dim]*drho_inv;

            m=static_cast<int>(p);
            m=MIN(m,nrho-2);
            p-=m;
            p=MIN(p,1.0);
            coef=F_arr[itype][m];
            
            tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
            if(E[ic_dim]>rho_max)
            {
                tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
                tmp0+=tmp1*(E[ic_dim]-rho_max);
            }
            
            E[ic_dim]=tmp0;

            if(cvec[ic_dim]!=0.0)
                en+=cvec[ic_dim]*(tmp0+c_0[itype])-3.0*kbT*cvec[ic_dim]*log(xvec[3*iatm+ic_dim+3]);
            
            en+=kbT*calc_ent(cvec[ic_dim]);
        }
        en+=kbT*calc_ent(c_iv);
    }
    
    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    return en_tot;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void ForceField_eam_dmd::init()
{    
    if(kbT==-1.0)
        Error::abort("temperature for the siumulation has not been set");
    
    
    cv_ptr=new Vec<type0>(atoms,1);
    E_ptr=new Vec<type0>(atoms,c_dim);
    dE_ptr=new Vec<type0>(atoms,c_dim);
    mu_ptr=new Vec<type0>(atoms,no_types,"mu");
    mu_ptr->change_dimension(0.0,no_types,no_types-c_dim);

    crd_ptr=new Vec<type0>(atoms,c_dim);
    s_ptr=new Vec<type0>(atoms,c_dim);
    
    x_tmp_ptr=new Vec<type0>(atoms,1);

    t_ptr=new Vec<type0>(atoms,c_dim);
    
    type0* mu=mu_ptr->begin();
    for(int i=0;i<c_dim*natms;i++)
        mu[i]=0.0;

}
/*--------------------------------------------
 fin
 --------------------------------------------*/
void ForceField_eam_dmd::fin()
{
    if(max_pairs)
    {
        delete [] rho_phi;
        delete [] drho_phi_dr;
        delete [] drho_phi_dalpha;
        max_pairs=0;
    }
    
    if(phi_psi_sz_sz)
        delete [] phi_psi_sz;
    phi_psi_sz_sz=0;
    
    if(n_phi_psi)
    {
        delete [] psi_IJ;
        delete [] psi_JI;
        delete [] phi_IJ;
        delete [] psi_r_IJ;
        delete [] psi_r_JI;
        delete [] phi_r_IJ;
        delete [] psi_alpha_IJ;
        delete [] psi_alpha_JI;
        delete [] phi_alpha_IJ;
        n_phi_psi=0;
    }
    
    if(M_N_sz_sz)
    {
        delete [] M_IJ;
        M_N_sz_sz=0;
    }
    
    delete t_ptr;
    delete s_ptr;
    delete crd_ptr;
    delete mu_ptr;
    delete dE_ptr;
    delete E_ptr;
    delete cv_ptr;
    delete x_tmp_ptr;
    
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceField_eam_dmd::init_xchng()
{
    Error::abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceField_eam_dmd::fin_xchng()
{
    Error::abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 ff_coef
 --------------------------------------------*/
void ForceField_eam_dmd::coef(int nargs,char** args)
{
    if(nargs!=2)
        Error::abort("wrong coeff command "
        "for ff dmd");
    
    cut_off_alloc();
    allocate();

    no_types=atom_types->no_types;
    c_dim=c->dim;
    x_dim=x->dim;
    
    read_file(args[1]);
}
/*-------------------------------------------
 
 --------------------------------------------*/
void ForceField_eam_dmd::allocate()
{
    if(no_types==atom_types->no_types)
        return;
    
    deallocate();
    no_types=atom_types->no_types;
    CREATE1D(c_0,no_types);
    CREATE1D(c_1,no_types);
    CREATE1D(g_fac,no_types);
}
/*-------------------------------------------
 
 --------------------------------------------*/
void ForceField_eam_dmd::deallocate()
{
    if(!no_types)
        return;
    
    delete [] c_0;
    delete [] c_1;
    delete [] g_fac;
    no_types=0;

}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_eam_dmd::read_file(char* file_name)
{
    int degree;
    type0 T0;
    char* file_format=NULL;
    char** files;
    files=new char*[no_types];
    for(int i=0;i<no_types;i++)
        files[i]=NULL;
    
    FileReader fr;
    
    fr.add_scl("degree",degree);
    fr.add_vlog()=VLogics("gt",0);
    
    fr.add_scl("T",T0);
    fr.add_vlog()=VLogics("gt",0.0);
    
    fr.add_scl("alpha_max",alpha_max);
    fr.add_vlog()=VLogics("gt",0.0);
        
    fr.add_1D("gamma",g_fac);
    fr.add_vlog()=VLogics("gt",0.0)*VLogics("le",1.0);

    fr.add_1D("r_crd",rsq_crd);
    fr.add_vlog()=VLogics("gt",0.0);
    
    fr.add_scl("file_format",file_format);
    fr.add_vlog()=VLogics("eq","FS")
    +VLogics("eq","SetFunc")+VLogics("eq","SetFL");

    
    fr.add_1D("files",files);
    fr.add_clog()=Logics(file_format,"eq","FuncFL")-VLogics("set");
    
    fr.add_scl("file",files[0]);
    fr.add_clog()=(Logics(file_format,"eq","SetFL")+Logics(file_format,"eq","FS"))-VLogics("set");

    fr.read_file(file_name);
    
    

    
    set_temp(T0);
    set_weight_abs(degree);
    
    for(int i=0;i<no_types;i++)
        rsq_crd[i]*=rsq_crd[i];
    
    delete eam_reader;
    eam_reader=new EAMFileReader();
    
    eam_reader->file_format(file_format);
    eam_reader->add_file(files[0],0);
    if(!strcmp(file_format,"FuncFL"))
        for(int i=1;i<no_types;i++)
            eam_reader->add_file(files[i],i);
    delete [] files;
    
    
    type0* xvec=x->begin();
    int* id=atoms->id->begin();
    atom_type* typevec=type->begin();
    type0 alpha_bound=alpha_max/sqrt(2.0);
    for(int iatm=0;iatm<natms;iatm++)
    {
        for(int ic_dim=iatm*c_dim;ic_dim<(iatm+1)*c_dim;ic_dim++)
            if(xvec[3*iatm+ic_dim+3]>=alpha_bound)
                Error::abort("alpha_%s of atom %d is greater than "
                "the maximum allowed %lf (alpha_max/sqrt(2.0)) set "
                "by %s file",atom_types->atom_names[typevec[ic_dim]]
                ,id[iatm],alpha_bound,file_name);
    }
    
    eam_reader->setup();
    
    stride=eam_reader->stride;
    nr=eam_reader->nr;
    nrho=eam_reader->nrho;
    dr=eam_reader->dr;
    drho=eam_reader->drho;
    dr_inv=eam_reader->dr_inv;
    drho_inv=eam_reader->drho_inv;
    rho_max=eam_reader->rho_max;
    
    F_arr=eam_reader->F_arr;
    phi_r_arr=eam_reader->phi_r_arr;
    rho_r_arr=eam_reader->rho_r_arr;
    type2rho=eam_reader->type2rho;
    type2phi=eam_reader->type2phi;
    type2rho_pair_ij=eam_reader->type2rho_pair_ij;
    type2rho_pair_ji=eam_reader->type2rho_pair_ji;
    type2phi_pair_ij=eam_reader->type2phi_pair_ij;
    type2phi_pair_ji=eam_reader->type2phi_pair_ji;
        
    type0 cut_max=0.0;
    for(int itype=0;itype<no_types;itype++)
        for(int jtype=0;jtype<no_types;jtype++)
        {
            type0 tmp=sqrt(eam_reader->cut_sq[COMP(itype,jtype)]);
            cut_max=MAX(cut_max,tmp);
            cut[itype][jtype]=tmp;
            cut_sq[itype][jtype]=(tmp+xi[no_i-1]*alpha_max)*(tmp+xi[no_i-1]*alpha_max);
        }
    
    type0 mod_rc=cut_max+xi[no_i-1]*alpha_max;
    cut_sq_mod_0=mod_rc*mod_rc;
    
    for(int itype=0;itype<no_types;itype++)
        if(rsq_crd[itype]>cut_sq[itype][itype])
            Error::abort("r_crd(%s) set by file %s should be less than %lf"
            ,file_name,atom_types->atom_names[itype],sqrt(cut_sq[itype][itype]));
}
/*--------------------------------------------
 set the temperature in the simulation
 --------------------------------------------*/
void ForceField_eam_dmd::set_temp(type0 T)
{
    type0 kb=8.617332478e-5;
    type0 hbar=6.5821192815e-16;
    type0 mass;
    type0 deb_l;
    
    for(int i=0;i<no_types;i++)
    {
        mass=atom_types->mass[i];
        c_1[i]=sqrt(0.5*kb*T/mass)/M_PI;
        mass*=1.0364269184093291236e-28;
        deb_l=hbar*hbar*2.0/(mass*kb*T);
        c_0[i]=1.5*kb*T*(log(deb_l)-1.0);
    }
    
    kbT=kb*T;
    beta=1.0/kbT;
}
/*--------------------------------------------
 Gaussian-Hermite quadrature weights and
 abscissas for 1 to 14 points
 --------------------------------------------*/
type0 ForceField_eam_dmd::set_weight_abs(int n)
{
    if(n<1)
        Error::abort("number of gaussian "
        "points for ff eam_dmd should be "
        "greater than 1");
    
    if(no_i)
    {
        delete [] xi;
        delete [] wi_0;
        delete [] wi_1;
        delete [] wi_2;
    }
    
    no_i=n;
    
    CREATE1D(xi,n);
    CREATE1D(wi_0,n);
    CREATE1D(wi_1,n);
    CREATE1D(wi_2,n);
    
    XMath::quadrature_hg(n,xi,wi_0);
    
    for(int i=0;i<no_i;i++)
    {
        wi_1[i]=wi_0[i]*xi[i];
        wi_2[i]=wi_0[i]*xi[i]*xi[i];
    }
    
    if(no_i>1)
    {
        type0 min_del=xi[no_i-1]-xi[0];
        for(int i=1;i<no_i;i++)
            min_del=MIN(min_del,xi[i]-xi[i-1]);
        return 2.0/min_del;
    }
    else
    {
        return -1.0;
    }
}
/*--------------------------------------------
 return M_{ij}^{\alpha}
 --------------------------------------------*/
type0 ForceField_eam_dmd::calc_ent(type0 x)
{
    type0 ans=x*log(x);
    if(isnan(ans))
        return 0.0;
    return ans;
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd::dc()
{
    calc_mu();

    type0* xvec=x->begin();
    type0* cvec=c->begin();
    type0* c_dvec=c_d->begin();
    atom_type* typevec=type->begin();
    type0* cv=cv_ptr->begin();
    type0* mu=mu_ptr->begin();
    
    for(int i=0;i<natms*c_dim;i++)
        c_dvec[i]=0.0;
    
    type0 gamma_i,gamma_j,x_ij_sq,alpha_Q_sq;
    type0* xi;
    type0* xj;
    type0 mu_ji,alpha_i,alpha_j,Qi;
    
    int iatm,jatm,icmp,jcmp;
    type0 d_i,d_j,dc_ij;
    
    int** neighbor_list_2nd=neighbor->neighbor_list_2nd;
    int* neighbor_list_size_2nd=neighbor->neighbor_list_size_2nd;
    
    for(int ic_dim=0,itype;ic_dim<natms*c_dim;ic_dim++)
    {
        itype=typevec[ic_dim];
        icmp=ic_dim%c_dim;
        iatm=ic_dim/c_dim;
        xi=xvec+x_dim*iatm;
        alpha_i=xi[3+icmp];
        
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            
            jc_dim=neighbor_list_2nd[ic_dim][j];
            jcmp=jc_dim%c_dim;
            jatm=jc_dim/c_dim;
            xj=xvec+x_dim*jatm;
            alpha_j=xj[3+icmp];
            mu_ji=beta*(mu[jc_dim]-mu[ic_dim]);
            
            x_ij_sq=(xi[0]-xj[0])*(xi[0]-xj[0])
            +(xi[1]-xj[1])*(xi[1]-xj[1])
            +(xi[2]-xj[2])*(xi[2]-xj[2]);
            gamma_i=x_ij_sq/(alpha_i*alpha_i);
            gamma_j=x_ij_sq/(alpha_j*alpha_j);
            
            calc_Q(itype,gamma_i,gamma_j,mu_ji,Qi,alpha_Q_sq);
            
            d_i=d_j=-c_1[itype]*x_ij_sq*alpha_Q_sq;
            d_i/=alpha_i*alpha_i*alpha_i;
            d_j/=alpha_j*alpha_j*alpha_j;
            
            dc_ij=d_i*exp(-Qi)*cvec[ic_dim]*cv[jatm]-d_j*exp(-Qi+mu_ji)*cvec[jc_dim]*cv[iatm];

            c_dvec[ic_dim]+=dc_ij;
            if(jc_dim<natms*c_dim)
            {
                c_dvec[jc_dim]-=dc_ij;
            }
        }
    }
}
/*--------------------------------------------
 calculate norm of d^2c/dt^2
 --------------------------------------------*/
type0 ForceField_eam_dmd::ddc_norm()
{
    int n=natms*c_dim;
    type0* tmp=NULL;
    CREATE1D(tmp,n);
    for(int i=0;i<n;i++)
        tmp[i]=0;
    update_J(1.0,tmp,tmp);
    delete [] tmp;
    
    
    type0* c_dvec=c_d->begin();
    type0* s=s_ptr->begin();
    memcpy(s,c_dvec,sizeof(type0)*n);
    operator()(s_ptr,t_ptr);
    type0* t=t_ptr->begin();
    type0* cvec=c->begin();
    type0 ans_lcl=0.0,ans,tmp0;
    
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        if(cvec[ic_dim]>=0.0)
        {
            tmp0=t[ic_dim]+s[ic_dim];
            ans_lcl+=tmp0*tmp0;
        }
    
    MPI_Allreduce(&ans_lcl,&ans,1,MPI_TYPE0,MPI_SUM,world);
    return sqrt(ans);
    
}
/*--------------------------------------------
 calculate norm of d^2c/dt^2
 --------------------------------------------*/
void ForceField_eam_dmd::ddc(type0* ddc_)
{
    
    
    
    Error::abort("not done");
}
/*--------------------------------------------
 calculate mu crd dE ddE and local energy
 --------------------------------------------*/
void ForceField_eam_dmd::calc_mu()
{
    type0* cvec=c->begin();
    atom_type* typevec=type->begin();
    type0* xvec=x->begin();
    type0* E=E_ptr->begin();
    type0* mu=mu_ptr->begin();
    type0* dE=dE_ptr->begin();
    type0* cv=cv_ptr->begin();
    type0* ddE=E;
    
    for(int i=0;i<natms*c_dim;i++)
        mu[i]=E[i]=0.0;
    int itype,m;
    type0 p,tmp0,tmp1,tmp2;
    type0* coef;
    
    nrgy_strss_lcl[0]=0.0;
    
    for(int i=0;i<natms+natms_ph;i++)
    {
        cv[i]=1.0;
        for(int j=0;j<c_dim;j++)
            if(cvec[i*c_dim+j]>=0.0)
                cv[i]-=cvec[i*c_dim+j];
    }
    
    int istart=0;
    for(int ic_dim=0;ic_dim<c_dim*natms;ic_dim++)
    {
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            E[ic_dim]+=cvec[jc_dim]*psi_JI[istart];
            mu[ic_dim]+=cvec[jc_dim]*phi_IJ[istart];
            
            
            if(jc_dim<phi_psi_sz_sz)
            {
                E[jc_dim]+=cvec[ic_dim]*psi_IJ[istart];
                mu[jc_dim]+=cvec[ic_dim]*phi_IJ[istart];
                nrgy_strss_lcl[0]+=cvec[ic_dim]*cvec[jc_dim]*phi_IJ[istart];
            }
            else
                nrgy_strss_lcl[0]+=0.5*cvec[ic_dim]*cvec[jc_dim]*phi_IJ[istart];
                
            istart++;
        }
        if(cvec[ic_dim]>=0.0)
        {
            itype=typevec[ic_dim];
            p=E[ic_dim]*drho_inv;
            m=static_cast<int>(p);
            m=MIN(m,nrho-2);
            p-=m;
            p=MIN(p,1.0);
            coef=F_arr[itype][m];
            tmp0=(((coef[4]*p+coef[3])*p+coef[2])*p+coef[1])*p+coef[0];
            tmp1=(((4.0*coef[4]*p+3.0*coef[3])*p+2.0*coef[2])*p+coef[1])*drho_inv;
            tmp2=(((12.0*coef[4]*p+6.0*coef[3])*p+2.0*coef[2]))*drho_inv*drho_inv;
            if(E[ic_dim]>rho_max)
            {
                tmp0+=tmp1*(E[ic_dim]-rho_max);
                tmp2=0.0;
            }
            
            dE[ic_dim]=tmp1;
            ddE[ic_dim]=tmp2;
            mu[ic_dim]+=tmp0;
            tmp2=tmp0+c_0[itype];
            if(cvec[ic_dim]!=0.0)
                tmp2-=3.0*kbT*log(xvec[3*(ic_dim/c_dim +1)+ic_dim]);
            nrgy_strss_lcl[0]+=cvec[ic_dim]*tmp2+kbT*calc_ent(cvec[ic_dim]);
        }
        if(ic_dim%c_dim==c_dim-1)
            nrgy_strss_lcl[0]+=kbT*calc_ent(cv[ic_dim/c_dim]);
    }
    
    dynamic->update(dE_ptr);
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
    {
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            mu[ic_dim]+=psi_IJ[istart]*cvec[jc_dim]*dE[jc_dim];
            
            if(jc_dim<phi_psi_sz_sz)
            {
                mu[jc_dim]+=psi_JI[istart]*cvec[ic_dim]*dE[ic_dim];
            }
            istart++;
        }
    }
    
    dynamic->update(mu_ptr);
}
/*--------------------------------------------
 create the sparse matrices
 --------------------------------------------*/
void ForceField_eam_dmd::init_static()
{

    if(phi_psi_sz_sz)
        delete [] phi_psi_sz;
        
    if(n_phi_psi)
    {
        delete [] psi_IJ;
        delete [] psi_JI;
        delete [] phi_IJ;
        delete [] psi_r_IJ;
        delete [] psi_r_JI;
        delete [] phi_r_IJ;
        delete [] psi_alpha_IJ;
        delete [] psi_alpha_JI;
        delete [] phi_alpha_IJ;
    }
    
    if(M_N_sz_sz)
    {
        delete [] M_IJ;
    }
    
    
    int no_pairs=neighbor->no_pairs;
    phi_psi_sz_sz=natms*c_dim;
    
    CREATE1D(phi_psi_sz,phi_psi_sz_sz);
    CREATE1D(phi_psi_cmp,c_dim*c_dim*no_pairs);
    CREATE1D(psi_IJ,c_dim*c_dim*no_pairs);
    CREATE1D(psi_JI,c_dim*c_dim*no_pairs);
    CREATE1D(phi_IJ,c_dim*c_dim*no_pairs);
    CREATE1D(psi_r_IJ,c_dim*c_dim*no_pairs);
    CREATE1D(psi_r_JI,c_dim*c_dim*no_pairs);
    CREATE1D(phi_r_IJ,c_dim*c_dim*no_pairs);
    CREATE1D(psi_alpha_IJ,c_dim*c_dim*no_pairs);
    CREATE1D(psi_alpha_JI,c_dim*c_dim*no_pairs);
    CREATE1D(phi_alpha_IJ,c_dim*c_dim*no_pairs);
    
    for(int i=0;i<phi_psi_sz_sz;i++)
        phi_psi_sz[i]=0;
    
    
    type0* cvec=c->begin();
    atom_type* typevec=type->begin();
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    int istart=0;
    int iistart=0;
    int isps=0;
    type0 psi_ij,psi_ji,phi_ij;
    type0 psi_r_ij,psi_r_ji,phi_r_ij;
    type0 psi_alpha_ij,psi_alpha_ji,phi_alpha_ij;
    for(int iatm=0;iatm<natms;iatm++)
    {
        
        for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && cvec[ic_dim]>=0.0;ic_dim++)
        {
            itype=typevec[ic_dim];
            
            iistart=istart;
            for(int j=0,jatm;j<neighbor_list_size[iatm];j++)
            {
                jatm=neighbor_list[iatm][j];

                for(int jc_dim=jatm*c_dim,jtype;jc_dim<(jatm+1)*c_dim && cvec[jc_dim]>=0.0;jc_dim++)
                {
                    jtype=typevec[jc_dim];
                    
                    psi_ij=rho_phi[iistart+type2rho_pair_ij[itype][jtype]];
                    psi_ji=rho_phi[iistart+type2rho_pair_ji[jtype][itype]];
                    phi_ij=rho_phi[iistart+type2phi_pair_ij[itype][jtype]];
                    
                    psi_r_ij=drho_phi_dr[iistart+type2rho_pair_ij[itype][jtype]];
                    psi_r_ji=drho_phi_dr[iistart+type2rho_pair_ji[jtype][itype]];
                    phi_r_ij=drho_phi_dr[iistart+type2phi_pair_ij[itype][jtype]];
                    
                    psi_alpha_ij=drho_phi_dalpha[iistart+type2rho_pair_ij[itype][jtype]];
                    psi_alpha_ji=drho_phi_dalpha[iistart+type2rho_pair_ji[jtype][itype]];
                    phi_alpha_ij=drho_phi_dalpha[iistart+type2phi_pair_ij[itype][jtype]];
                    
                    if(psi_ij!=0.0 || psi_ji!=0.0 || phi_ij!=0.0 ||
                       psi_r_ij!=0.0 || psi_r_ji!=0.0 || phi_r_ij!=0.0 ||
                       psi_alpha_ij!=0.0 || psi_alpha_ji!=0.0 || phi_alpha_ij!=0.0)
                    {
                        phi_psi_sz[ic_dim]++;
                        phi_psi_cmp[isps]=jc_dim;
                        
                        psi_IJ[isps]=psi_ij;
                        psi_JI[isps]=psi_ji;
                        phi_IJ[isps]=phi_ij;
                        
                        psi_r_IJ[isps]=psi_r_ij;
                        psi_r_JI[isps]=psi_r_ji;
                        phi_r_IJ[isps]=phi_r_ij;
                        
                        psi_alpha_IJ[isps]=psi_alpha_ij;
                        psi_alpha_JI[isps]=psi_alpha_ji;
                        phi_alpha_IJ[isps]=phi_alpha_ij;
                        
                        isps++;
                    }
                }
                iistart+=stride;
            }
        }
        istart=iistart;
    }
    if(max_pairs)
    {
        delete [] rho_phi;
        delete [] drho_phi_dr;
        delete [] drho_phi_dalpha;
        max_pairs=0;
    }
    
    GROW(phi_psi_cmp,c_dim*c_dim*no_pairs,isps);
    GROW(psi_IJ,c_dim*c_dim*no_pairs,isps);
    GROW(psi_JI,c_dim*c_dim*no_pairs,isps);
    GROW(phi_IJ,c_dim*c_dim*no_pairs,isps);
    GROW(psi_r_IJ,c_dim*c_dim*no_pairs,isps);
    GROW(psi_r_JI,c_dim*c_dim*no_pairs,isps);
    GROW(phi_r_IJ,c_dim*c_dim*no_pairs,isps);
    GROW(psi_alpha_IJ,c_dim*c_dim*no_pairs,isps);
    GROW(psi_alpha_JI,c_dim*c_dim*no_pairs,isps);
    GROW(phi_alpha_IJ,c_dim*c_dim*no_pairs,isps);
    n_phi_psi=isps;

    
    int* neighbor_list_size_2nd=neighbor->neighbor_list_size_2nd;
    M_N_sz_sz=0;
    for(int i=0;i<natms*c_dim;i++)
        M_N_sz_sz+=neighbor_list_size_2nd[i];
    
    CREATE1D(M_IJ,3*M_N_sz_sz);
    
}
/*--------------------------------------------
 create the sparse matrices
 --------------------------------------------*/
type0 ForceField_eam_dmd::update_J(type0 alpha,type0* a,type0* g)
{
    calc_mu();
    
    type0* cvec=c->begin();
    type0* xvec=x->begin();
    atom_type* typevec=type->begin();
    type0* mu=mu_ptr->begin();
    type0* cv=cv_ptr->begin();
    
    int istart;
    
    
    type0* c_dvec=c_d->begin();
    for(int i=0;i<phi_psi_sz_sz;i++)
    {
        c_dvec[i]=0.0;
        if(cvec[i]>=0.0)
            g[i]=cvec[i]+a[i];
        else
            g[i]=0.0;
    }
    
    if(alpha<0.0)
        Error::abort("alpha less than 0.0");
    
    type0 iota=log(alpha);
    alpha_tmp=alpha;
    
    type0 gamma_i,gamma_j,x_ij_sq,alpha_Q_sq;
    type0 theta_i,theta_j;

    type0* xi;
    type0* xj;
    type0 ans_lcl=0.0,dc_ij,dg_ij,mu_ji,alpha_i,alpha_j;
    type0 d_i,d_j,Qi;
    type0 exp_mod_Qi,exp_mod_Qj;
    int iatm,jatm,icmp,jcmp;
    int** neighbor_list_2nd=neighbor->neighbor_list_2nd;
    int* neighbor_list_size_2nd=neighbor->neighbor_list_size_2nd;
    
    istart=0;
    for(int ic_dim=0,itype;ic_dim<phi_psi_sz_sz;ic_dim++)
    {
        itype=typevec[ic_dim];
        icmp=ic_dim%c_dim;
        iatm=ic_dim/c_dim;
        xi=xvec+x_dim*iatm;
        alpha_i=xi[3+icmp];
        
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            jc_dim=neighbor_list_2nd[ic_dim][j];
            jcmp=jc_dim%c_dim;
            jatm=jc_dim/c_dim;
            xj=xvec+x_dim*jatm;
            alpha_j=xj[3+icmp];
            mu_ji=beta*(mu[jc_dim]-mu[ic_dim]);
            
            x_ij_sq=(xi[0]-xj[0])*(xi[0]-xj[0])
            +(xi[1]-xj[1])*(xi[1]-xj[1])
            +(xi[2]-xj[2])*(xi[2]-xj[2]);
            gamma_i=x_ij_sq/(alpha_i*alpha_i);
            gamma_j=x_ij_sq/(alpha_j*alpha_j);
            
            calc_Q(itype,gamma_i,gamma_j,mu_ji,Qi,alpha_Q_sq,theta_i);
            
            exp_mod_Qi=exp(-Qi+iota);
            exp_mod_Qj=exp(-Qi+mu_ji+iota);
            
            d_i=d_j=-c_1[itype]*x_ij_sq*alpha_Q_sq;
            d_i/=alpha_i*alpha_i*alpha_i;
            d_j/=alpha_j*alpha_j*alpha_j;
            theta_j=theta_i+1.0;
            
            dc_ij=d_i*exp(-Qi)*cvec[ic_dim]*cv[jatm]-d_j*exp(-Qi+mu_ji)*cvec[jc_dim]*cv[iatm];
            dg_ij=d_i*exp_mod_Qi*cvec[ic_dim]*cv[jatm]-d_j*exp_mod_Qj*cvec[jc_dim]*cv[iatm];

            M_IJ[istart++]=kbT*(d_i*theta_i*exp_mod_Qi*cvec[ic_dim]*cv[jatm]-d_j*theta_j*exp_mod_Qj*cvec[jc_dim]*cv[iatm]);
            M_IJ[istart++]=d_i*exp_mod_Qi;
            M_IJ[istart++]=d_j*exp_mod_Qj;
            g[ic_dim]-=dg_ij;
            c_dvec[ic_dim]+=dc_ij;
            
            if(jc_dim<phi_psi_sz_sz)
            {
                g[jc_dim]+=dg_ij;
                c_dvec[jc_dim]-=dc_ij;
            }
        }
        if(cvec[ic_dim]>=0.0)
        {
            
            ans_lcl+=g[ic_dim]*g[ic_dim];
        }
    }
    type0 ans;
    MPI_Allreduce(&ans_lcl,&ans,1,MPI_TYPE0,MPI_SUM,world);
    return sqrt(ans);
    
}
/*--------------------------------------------
 create the sparse matrices
 --------------------------------------------*/
void ForceField_eam_dmd::operator()(Vec<type0>* x_ptr,Vec<type0>* Ax_ptr)
{
    dynamic->update(x_ptr);
    
    type0 tmp0;
    type0* cvec=c->begin();
    type0* ddE=E_ptr->begin();
    type0* dE=dE_ptr->begin();
    type0* x=x_ptr->begin();
    type0* Ax=Ax_ptr->begin();
    type0* b0=crd_ptr->begin();
    
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        b0[ic_dim]=Ax[ic_dim]=0.0;
    
    int istart=0;
    
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            Ax[ic_dim]+=cvec[ic_dim]*ddE[ic_dim]*psi_JI[istart]*x[jc_dim];
            if(jc_dim<phi_psi_sz_sz)
                Ax[jc_dim]+=cvec[jc_dim]*ddE[jc_dim]*psi_IJ[istart]*x[ic_dim];
            istart++;
        }
    
    dynamic->update(Ax_ptr);
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            tmp0=psi_IJ[istart]*dE[ic_dim]+psi_JI[istart]*dE[jc_dim]+phi_IJ[istart];
            
            b0[ic_dim]+=psi_IJ[istart]*Ax[jc_dim]+tmp0*x[jc_dim];
            if(jc_dim<phi_psi_sz_sz)
                b0[jc_dim]+=psi_JI[istart]*Ax[ic_dim]+tmp0*x[ic_dim];
            
            istart++;
        }
    
    dynamic->update(crd_ptr);
    
    type0* x_tmp=x_tmp_ptr->begin();
    type0* cv=cv_ptr->begin();
    for(int i=0;i<natms+natms_ph;i++)
    {
        x_tmp[i]=0.0;
        for(int j=0;j<c_dim;j++)
            if(cvec[i*c_dim+j]>=0.0)
                x_tmp[i]+=x[i*c_dim+j];
    }

    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        Ax[ic_dim]=-x[ic_dim];
    
    
    int** neighbor_list_2nd=neighbor->neighbor_list_2nd;
    int* neighbor_list_size_2nd=neighbor->neighbor_list_size_2nd;
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            jc_dim=neighbor_list_2nd[ic_dim][j];
            
            tmp0=M_IJ[istart++]*(b0[jc_dim]-b0[ic_dim]);
            tmp0+=M_IJ[istart++]*(cv[jc_dim/c_dim]*x[ic_dim]-cvec[ic_dim]*x_tmp[jc_dim/c_dim]);
            tmp0-=M_IJ[istart++]*(cv[ic_dim/c_dim]*x[jc_dim]-cvec[jc_dim]*x_tmp[ic_dim/c_dim]);
            Ax[ic_dim]+=tmp0;
            
            if(jc_dim<phi_psi_sz_sz)
                Ax[jc_dim]-=tmp0;
        }
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceField_eam_dmd::
force_calc_static(bool st_clc)
{
    type0* xvec=x->begin();
    type0* fvec=f->begin();
    type0* dE=dE_ptr->begin();
    type0* cvec=c->begin();
    
    int i_comp,I_comp;
    int j_comp,J_comp;
    type0 dx0,dx1,dx2;
    type0 fpair,apair;
    
    if(st_clc)
        for(int i=1;i<7;i++)
            nrgy_strss_lcl[i]=0.0;
   
    int istart=0;
    
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
    {
        i_comp=(3+c_dim)*(ic_dim/c_dim);
        I_comp=3*(ic_dim/c_dim)+ic_dim;
    
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            fpair=-(psi_r_JI[istart]*dE[ic_dim]+psi_r_IJ[istart]*dE[jc_dim]+phi_r_IJ[istart])*cvec[ic_dim]*cvec[jc_dim];
            apair=-(psi_alpha_JI[istart]*dE[ic_dim]+psi_alpha_IJ[istart]*dE[jc_dim]+phi_alpha_IJ[istart])*cvec[ic_dim]*cvec[jc_dim];
            
            j_comp=(3+c_dim)*(jc_dim/c_dim);
            J_comp=3*(jc_dim/c_dim)+jc_dim;
            
            dx0=xvec[i_comp]-xvec[j_comp];
            dx1=xvec[i_comp+1]-xvec[j_comp+1];
            dx2=xvec[i_comp+2]-xvec[j_comp+2];
            
            fvec[i_comp]+=dx0*fpair;
            fvec[i_comp+1]+=dx1*fpair;
            fvec[i_comp+2]+=dx2*fpair;
            fvec[I_comp]+=apair*xvec[I_comp];
            
            if(jc_dim<phi_psi_sz_sz)
            {
                fvec[j_comp]-=dx0*fpair;
                fvec[j_comp+1]-=dx1*fpair;
                fvec[j_comp+2]-=dx2*fpair;
                fvec[J_comp]+=apair*xvec[J_comp];
            }
            
            if(jc_dim>=phi_psi_sz_sz)
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
            
            istart++;
        }
        
        if(cvec[ic_dim]>0.0)
            fvec[I_comp]+=3.0*kbT*cvec[ic_dim]/xvec[I_comp];
    }
    
    if(st_clc)
    {
        for(int i=0;i<7;i++)
            nrgy_strss[i]=0.0;
        
        MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,7,MPI_TYPE0,MPI_SUM,world);
    }
    else
        MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,1,MPI_TYPE0,MPI_SUM,world);
}

/*--------------------------------------------
 gamma_i = x_ij*x_ij/(alpha_i*alpha_i)
 gamma_j = x_ij*x_ij/(alpha_j*alpha_j)
 mu_ji = beta*(mu_j-mu_i)+3*log(alpha_j/alpha_i)
 
 Q=beta*U
 
 alpha_Q_sq=alpha_U_sq/(x_ij*x_ij)
 --------------------------------------------*/
inline void ForceField_eam_dmd::
calc_Q(int& itype,type0& gamma_i,type0& gamma_j,type0& mu_ji
,type0& Q,type0& alpha_Q_sq)
{
    type0 z0;
    type0 a=5.0*(gamma_i-gamma_j-6.0*mu_ji);
    type0 xi=-2.0*(gamma_i+gamma_j);
    type0 d=sqrt((xi-a)*(xi-a)-8.0*a*gamma_i);
    z0=-4.0*gamma_i/((xi-a)-d);
    
    
    if(z0<0.0 || z0>1.0)
    {
        Error::abort("could not find z0");
    }
    
    
    Q=z0*z0*(1.0-z0)*(1.0-z0)*((1.0-z0)*gamma_i+z0*gamma_j);
    Q+=z0*z0*z0*(6.0*z0*z0-15.0*z0+10.0)*mu_ji;
    Q*=g_fac[itype];
    Q+=0.5*(1.0-g_fac[itype])*mu_ji;
    alpha_Q_sq=1.0/((1.0-z0)*gamma_i+z0*gamma_j);
    
}
/*--------------------------------------------
 gamma_i = x_ij*x_ij/(alpha_i*alpha_i)
 gamma_j = x_ij*x_ij/(alpha_j*alpha_j)
 mu_ji = beta*(mu_j-mu_i)+3*log(alpha_j/alpha_i)
 
 Q=beta*U
 dQ=dU/d((mu_j-mu_i)+3*log(alpha_j/alpha_i)/beta)
 
 alpha_Q_sq=alpha_U_sq/(x_ij*x_ij)
 dalpha_Q_sq=1/(beta*x_ij*x_ij)
 dalpha_Q_sq/d((mu_j-mu_i)+3*log(alpha_j/alpha_i)/beta)
 --------------------------------------------*/
inline void ForceField_eam_dmd::
calc_Q(int& itype,type0& gamma_i,type0& gamma_j,type0& mu_ji
,type0& Q,type0& alpha_Q_sq,type0& theta)
{
    type0 z0,dz0,dQ;
    type0 a=5.0*(gamma_i-gamma_j-6.0*mu_ji);
    type0 xi=-2.0*(gamma_i+gamma_j);
    type0 d=sqrt((xi-a)*(xi-a)-8.0*a*gamma_i);
    z0=-4.0*gamma_i/((xi-a)-d);
    dz0=30.0*z0*(1.0-z0)/d;
    
    if(z0<0.0 || z0>1.0)
    {
        Error::abort("could not find z0");
    }
    
    dQ=z0*z0*z0*(6.0*z0*z0-15.0*z0+10.0);
    Q=z0*z0*(1.0-z0)*(1.0-z0)*((1.0-z0)*gamma_i+z0*gamma_j);
    Q+=dQ*mu_ji;
    Q*=g_fac[itype];
    Q+=0.5*(1.0-g_fac[itype])*mu_ji;
    dQ*=g_fac[itype];
    dQ+=0.5*(1.0-g_fac[itype]);
    
    alpha_Q_sq=1.0/((1.0-z0)*gamma_i+z0*gamma_j);
    theta=alpha_Q_sq*(gamma_i-gamma_j)*dz0-dQ;
}









