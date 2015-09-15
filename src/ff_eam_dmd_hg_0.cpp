#include <stdlib.h>
#include "ff_eam_dmd_hg_0.h"
#include "neighbor.h"
#include "atom_types.h"
#include "xmath.h"
using namespace MAPP_NS;
#define PI_IN_SQ 0.564189583547756286948079451561
#define CRD_ENBL 0
#define OLD_NEIGH 0
#define DIFF_ENBL 1
enum{NOT_SET,FUNC_FL,SET_FL,FINNIS_FL};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam_dmd_hg_0::
ForceField_eam_dmd_hg_0(MAPP* mapp) : ForceField(mapp)
{
    if(mapp->mode!=DMD_mode)
        error->abort("ff eam_dmd works only "
                     "for dmd mode");
    
    no_i=0;
    allocated=0;
    eam_mode=NOT_SET;
    max_pairs=0;
    neigh_lst_sz_sz=0;
    
    
    no_types=atom_types->no_types;
    if(atoms->vectors[0]->dim!=3+no_types)
        error->abort("the dimension of x"
        " vector should be 3 + no of types");
    
    
    
    CREATE1D(c_0,no_types);
    CREATE1D(f_t,no_types);
    CREATE1D(rsq_crd,no_types);
    CREATE1D(Q_nrm,no_types);
    
    c_d_n=atoms->find_exist("c_d");
    
    if(c_d_n<0)
    {
        int tmp0=atoms->find("c");
        tmp0=atoms->vectors[tmp0]->dim;
        c_d_n=atoms->add<type0>(0,tmp0,"c_d");
    }
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_eam_dmd_hg_0::~ForceField_eam_dmd_hg_0()
{
    if(no_types)
    {
        delete [] c_0;
        delete [] rsq_crd;
        delete [] f_t;
        delete [] Q_nrm;
    }
    
    
    if(allocated) clean_up();
    
    if(no_i)
    {
        delete [] xi;
        delete [] wi_0;
        delete [] wi_1;
        delete [] wi_2;
    }
    
    if(max_pairs)
    {
        delete [] rho_phi;
        delete [] drho_phi_dr;
        delete [] drho_phi_dalpha;
    }
    
    if(neigh_lst_sz_sz)
    {
        for(int i=0;i<neigh_lst_sz_sz;i++)
            if(neigh_lst_sz[i])
                delete [] neigh_lst[i];
        
        delete [] neigh_lst_sz;
        delete [] neigh_lst;
    }
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::
force_calc(int st_clc,type0* en_st)
{
    if(DIFF_ENBL)
        force_calc_diff(st_clc,en_st);
    else
        force_calc_ndiff(st_clc,en_st);
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::
force_calc_ndiff(int st_clc,type0* en_st)
{
    
    type0* x;
    atoms->vectors[x_n]->ret(x);
    type0* f;
    atoms->vectors[f_n]->ret(f);
    type0* E;
    atoms->vectors[E_n]->ret(E);
    type0* dE;
    atoms->vectors[dE_n]->ret(dE);
    type0* mu;
    atoms->vectors[mu_n]->ret(mu);
    type0* c;
    atoms->vectors[c_n]->ret(c);
    
    int iatm,jatm;
    
    int icomp,jcomp,istart;
    type0 dx0,dx1,dx2,rsq;
    type0 alpha;
    type0 r=0.0,r_inv=0.0;
    int m;
    type0* coef;
    type0 alpha_inv,tmp0,tmp1,rtmp;
    type0 upper,lower;
    type0 fpair,apair;
    
    type0 p,p2,p3,p4,c_iv;
    type0 phi_it_jt_0,phi_it_jt_1,phi_it_jt_2;
    type0 rho_it_jt_0,rho_it_jt_1,rho_it_jt_2;
    type0 rho_jt_it_0,rho_jt_it_1,rho_jt_it_2;
    type0 a0,a1,a2;
    
    
    
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
    
    nrgy_strss[0]=0.0;
    if (st_clc)
        for (int i=1;i<7;i++)
            nrgy_strss[i]=0.0;
    int natms=atoms->natms;
    for(int i=0;i<natms*no_types;i++) E[i]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+no_types)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+no_types)*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            if(rsq<cut_sq_mod_0)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                
                for(int itype=0;itype<no_types; itype++)
                {
                    for(int jtype=0;jtype<no_types; jtype++)
                    {
                        if(c[no_types*iatm+itype]>=0.0 && c[no_types*jatm+jtype]>=0.0)
                        {
                            alpha=sqrt(x[icomp+3+itype]*x[icomp+3+itype]+x[jcomp+3+jtype]*x[jcomp+3+jtype]);
                            
                            if(alpha<alpha_max)
                            {
                                alpha_inv=1.0/alpha;
                                
                                upper=(r+rc)*alpha_inv;
                                lower=(r-rc)*alpha_inv;
                                
                                if(lower<xi[no_i-1])
                                {
                                    phi_it_jt_0=phi_it_jt_1=phi_it_jt_2=0.0;
                                    rho_it_jt_0=rho_it_jt_1=rho_it_jt_2=0.0;
                                    rho_jt_it_0=rho_jt_it_1=rho_jt_it_2=0.0;
                                    
                                    for(int i=0;i<no_i;i++)
                                    {
                                        if(xi[i]>lower && xi[i]<upper)
                                        {
                                            rtmp=r-xi[i]*alpha;
                                            
                                            p=fabs(rtmp)*dr_inv;
                                            m=static_cast<int>(p);
                                            m=MIN(m,nr-2);
                                            p-=m;
                                            p=MIN(p,1.0);
                                            p2=p*p;
                                            p3=p2*p;
                                            coef=rho_r_arr[type2rho[itype][jtype]][m];
                                            tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                            if(rtmp<0.0)
                                                tmp0*=-1.0;
                                            rho_it_jt_0+=wi_0[i]*tmp0;
                                            rho_it_jt_1+=wi_1[i]*tmp0;
                                            rho_it_jt_2+=wi_2[i]*tmp0;
                                            
                                            if(itype!=jtype)
                                            {
                                                coef=rho_r_arr[type2rho[jtype][itype]][m];
                                                tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                                if(rtmp<0.0)
                                                    tmp0*=-1.0;
                                            }
                                            rho_jt_it_0+=wi_0[i]*tmp0;
                                            rho_jt_it_1+=wi_1[i]*tmp0;
                                            rho_jt_it_2+=wi_2[i]*tmp0;
                                            
                                            
                                            coef=phi_r_arr[type2phi[itype][jtype]][m];
                                            tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                            if(rtmp<0.0)
                                                tmp0*=-1.0;
                                            
                                            phi_it_jt_0+=wi_0[i]*tmp0;
                                            phi_it_jt_1+=wi_1[i]*tmp0;
                                            phi_it_jt_2+=wi_2[i]*tmp0;
                                        }
                                    }
                                    a0=PI_IN_SQ*r_inv;
                                    a1=a0*r_inv;
                                    a2=a0*alpha_inv;
                                    
                                    rho_it_jt_2*=2.0*a2;
                                    rho_it_jt_2-=a2*rho_it_jt_0;
                                    rho_it_jt_1*=-2.0*a2;
                                    rho_it_jt_1-=a1*rho_it_jt_0;
                                    rho_it_jt_0*=a0;
                                    
                                    rho_jt_it_2*=2.0*a2;
                                    rho_jt_it_2-=a2*rho_jt_it_0;
                                    rho_jt_it_1*=-2.0*a2;
                                    rho_jt_it_1-=a1*rho_jt_it_0;
                                    rho_jt_it_0*=a0;
                                    
                                    phi_it_jt_2*=2.0*a2;
                                    phi_it_jt_2-=a2*phi_it_jt_0;
                                    phi_it_jt_1*=-2.0*a2;
                                    phi_it_jt_1-=a1*phi_it_jt_0;
                                    phi_it_jt_0*=a0;
                                    
                                    
                                    rho_phi[istart+type2rho_pair_ij[itype][jtype]]=rho_it_jt_0;
                                    rho_phi[istart+type2rho_pair_ji[jtype][itype]]=rho_jt_it_0;
                                    rho_phi[istart+type2phi_pair_ij[itype][jtype]]=phi_it_jt_0;
                                    
                                    drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]=rho_it_jt_1;
                                    drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]=rho_jt_it_1;
                                    drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]]=phi_it_jt_1;
                                    
                                    drho_phi_dalpha[istart+type2rho_pair_ij[itype][jtype]]=rho_it_jt_2;
                                    drho_phi_dalpha[istart+type2rho_pair_ji[jtype][itype]]=rho_jt_it_2;
                                    drho_phi_dalpha[istart+type2phi_pair_ij[itype][jtype]]=phi_it_jt_2;
                                    
                                    E[iatm*no_types+itype]+=c[jatm*no_types+jtype]*rho_jt_it_0;
                                    
                                    if (jatm<natms)
                                    {
                                        E[jatm*no_types+jtype]+=c[iatm*no_types+itype]*rho_it_jt_0;
                                        nrgy_strss[0]+=c[iatm*no_types+itype]*c[jatm*no_types+jtype]*phi_it_jt_0;
                                    }
                                    else
                                    {
                                        nrgy_strss[0]+=0.5*c[iatm*no_types+itype]*c[jatm*no_types+jtype]*phi_it_jt_0;
                                    }
                                }
                            }
                            
                        }
                    }
                    
                }
            }
            
            
            
            istart+=stride;
        }
        
        
        c_iv=1.0;
        for(int itype=0;itype<no_types; itype++)
        {
            if(c[no_types*iatm+itype]>=0.0)
            {
                c_iv-=c[iatm*no_types+itype];
                p=E[iatm*no_types+itype]*drho_inv;
                m=static_cast<int>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                p2=p*p;
                p3=p2*p;
                p4=p3*p;
                coef=F_arr[itype][m];
                
                tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                tmp1=(4.0*coef[4]*p3+3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*drho_inv;
                
                
                if(E[iatm*no_types+itype]>rho_max)
                    tmp0+=tmp1*(E[iatm*no_types+itype]-rho_max);
                
                E[iatm*no_types+itype]=tmp0;
                dE[iatm*no_types+itype]=tmp1;
                
                
                nrgy_strss[0]+=c[iatm*no_types+itype]*(tmp0+c_0[itype]
                -3.0*kbT*log(x[icomp+3+itype]));
                
                mu[iatm*no_types+itype]=tmp0-3.0*kbT*log(x[icomp+3+itype]);
                
                nrgy_strss[0]+=kbT*calc_ent(c[iatm*no_types+itype]);
                
                f[icomp+3+itype]+=3.0*kbT*c[iatm*no_types+itype]/x[icomp+3+itype];
            }
        }
        nrgy_strss[0]+=kbT*calc_ent(c_iv);
        
    }
    
    
    atoms->update_ph(dE_n);
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+no_types)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+no_types)*jatm;
            
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            r_inv=1.0/sqrt(dx0*dx0+dx1*dx1+dx2*dx2);
            
            for(int itype=0;itype<no_types; itype++)
            {
                for(int jtype=0;jtype<no_types; jtype++)
                {
                    if(c[no_types*iatm+itype]>=0.0 && c[no_types*jatm+jtype]>=0.0)
                    {
                        fpair=-(drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm*no_types+itype]
                                +drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm*no_types+jtype]
                                +drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]])*c[iatm*no_types+itype]*c[jatm*no_types+jtype]*r_inv;
                        
                        apair=-(
                                drho_phi_dalpha[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm*no_types+itype]
                                +drho_phi_dalpha[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm*no_types+jtype]
                                +drho_phi_dalpha[istart+type2phi_pair_ij[itype][jtype]])*c[iatm*no_types+itype]*c[jatm*no_types+jtype]/
                        sqrt(x[icomp+3+itype]*x[icomp+3+itype]+x[jcomp+3+jtype]*x[jcomp+3+jtype]);
                        
                        mu[iatm*no_types+itype]+=c[jatm*no_types+jtype]*(drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]]
                        +drho_phi_dalpha[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm*no_types+jtype]);
                        
                        if(jatm<natms)
                            mu[jatm*no_types+jtype]+=c[iatm*no_types+itype]*(drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]]
                            +drho_phi_dalpha[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm*no_types+itype]);
                        
                        
                        if(apair!=0.0 || fpair!=0.0)
                        {
                            f[icomp]+=dx0*fpair;
                            f[icomp+1]+=dx1*fpair;
                            f[icomp+2]+=dx2*fpair;
                            f[icomp+3+itype]+=apair*x[icomp+3+itype];
                            
                            if(jatm<natms)
                            {
                                f[jcomp]-=dx0*fpair;
                                f[jcomp+1]-=dx1*fpair;
                                f[jcomp+2]-=dx2*fpair;
                                f[jcomp+3+jtype]+=apair*x[jcomp+3+jtype];
                                
                                if(st_clc==1)
                                {
                                    nrgy_strss[1]-=fpair*dx0*dx0;
                                    nrgy_strss[2]-=fpair*dx1*dx1;
                                    nrgy_strss[3]-=fpair*dx2*dx2;
                                    nrgy_strss[4]-=fpair*dx1*dx2;
                                    nrgy_strss[5]-=fpair*dx2*dx0;
                                    nrgy_strss[6]-=fpair*dx0*dx1;
                                }
                                else if(st_clc==2)
                                {
                                    nrgy_strss[1]-=fpair*dx0*(image[3*jatm]-image[3*iatm]);
                                    nrgy_strss[2]-=fpair*dx1*(image[3*jatm+1]-image[3*iatm+1]);
                                    nrgy_strss[3]-=fpair*dx2*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[4]-=fpair*dx1*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[5]-=fpair*dx0*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[6]-=fpair*dx0*(image[3*jatm+1]-image[3*iatm+1]);
                                }
                            }
                            else
                            {
                                if(st_clc==1)
                                {
                                    nrgy_strss[1]-=0.5*fpair*dx0*dx0;
                                    nrgy_strss[2]-=0.5*fpair*dx1*dx1;
                                    nrgy_strss[3]-=0.5*fpair*dx2*dx2;
                                    nrgy_strss[4]-=0.5*fpair*dx1*dx2;
                                    nrgy_strss[5]-=0.5*fpair*dx2*dx0;
                                    nrgy_strss[6]-=0.5*fpair*dx0*dx1;
                                }
                                else if(st_clc==2)
                                {
                                    nrgy_strss[1]-=0.5*fpair*dx0*(image[3*jatm]-image[3*iatm]);
                                    nrgy_strss[2]-=0.5*fpair*dx1*(image[3*jatm+1]-image[3*iatm+1]);
                                    nrgy_strss[3]-=0.5*fpair*dx2*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[4]-=0.5*fpair*dx1*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[5]-=0.5*fpair*dx0*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[6]-=0.5*fpair*dx0*(image[3*jatm+1]-image[3*iatm+1]);
                                }
                            }
                        }
                    }
                }
            }
            istart+=stride;
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
 force calculation
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::
force_calc_diff(int st_clc,type0* en_st)
{
    
    type0* x;
    atoms->vectors[x_n]->ret(x);
    type0* f;
    atoms->vectors[f_n]->ret(f);
    type0* E;
    atoms->vectors[E_n]->ret(E);
    type0* dE;
    atoms->vectors[dE_n]->ret(dE);
    type0* mu;
    atoms->vectors[mu_n]->ret(mu);
    type0* c;
    atoms->vectors[c_n]->ret(c);
    
    int iatm,jatm;
    
    int icomp,jcomp,istart;
    type0 dx0,dx1,dx2,rsq;
    type0 alpha;
    type0 r=0.0,r_inv=0.0;
    int m;
    type0* coef;
    type0 alpha_inv,tmp0,tmp1,rtmp;
    type0 upper,lower;
    type0 fpair,apair;
    
    type0 p,p2,p3,p4,c_iv;
    type0 phi_it_jt_0,phi_it_jt_1,phi_it_jt_2;
    type0 rho_it_jt_0,rho_it_jt_1,rho_it_jt_2;
    type0 rho_jt_it_0,rho_jt_it_1,rho_jt_it_2;
    type0 a0;
    
    
    
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
    
    nrgy_strss[0]=0.0;
    if (st_clc)
        for (int i=1;i<7;i++)
            nrgy_strss[i]=0.0;
    int natms=atoms->natms;
    for(int i=0;i<natms*no_types;i++) E[i]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+no_types)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+no_types)*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            if(rsq<cut_sq_mod_0)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                
                for(int itype=0;itype<no_types; itype++)
                {
                    for(int jtype=0;jtype<no_types; jtype++)
                    {
                        if(c[no_types*iatm+itype]>=0.0 && c[no_types*jatm+jtype]>=0.0)
                        {
                            alpha=sqrt(x[icomp+3+itype]*x[icomp+3+itype]+x[jcomp+3+jtype]*x[jcomp+3+jtype]);
                            
                            if(alpha<alpha_max)
                            {
                                alpha_inv=1.0/alpha;
                                
                                upper=(r+rc)*alpha_inv;
                                lower=(r-rc)*alpha_inv;
                                
                                if(lower<xi[no_i-1])
                                {
                                    phi_it_jt_0=phi_it_jt_1=phi_it_jt_2=0.0;
                                    rho_it_jt_0=rho_it_jt_1=rho_it_jt_2=0.0;
                                    rho_jt_it_0=rho_jt_it_1=rho_jt_it_2=0.0;
                                    
                                    for(int i=0;i<no_i;i++)
                                    {
                                        if(xi[i]>lower && xi[i]<upper)
                                        {
                                            rtmp=r-xi[i]*alpha;
                                            
                                            p=fabs(rtmp)*dr_inv;
                                            m=static_cast<int>(p);
                                            m=MIN(m,nr-2);
                                            p-=m;
                                            p=MIN(p,1.0);
                                            p2=p*p;
                                            p3=p2*p;
                                            coef=rho_r_arr[type2rho[itype][jtype]][m];
                                            tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                            tmp1=(3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*dr_inv;
                                            if(rtmp<0.0)
                                                tmp0*=-1.0;
                                            rho_it_jt_0+=wi_0[i]*tmp0;
                                            rho_it_jt_1+=wi_0[i]*tmp1;
                                            rho_it_jt_2+=wi_1[i]*tmp1;
                                            
                                            if(itype!=jtype)
                                            {
                                                coef=rho_r_arr[type2rho[jtype][itype]][m];
                                                tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                                tmp1=(3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*dr_inv;
                                                if(rtmp<0.0)
                                                    tmp0*=-1.0;
                                            }
                                            rho_jt_it_0+=wi_0[i]*tmp0;
                                            rho_jt_it_1+=wi_0[i]*tmp1;
                                            rho_jt_it_2+=wi_1[i]*tmp1;
                                            
                                            
                                            coef=phi_r_arr[type2phi[itype][jtype]][m];
                                            tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                            tmp1=(3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*dr_inv;
                                            if(rtmp<0.0)
                                                tmp0*=-1.0;
                                            
                                            phi_it_jt_0+=wi_0[i]*tmp0;
                                            phi_it_jt_1+=wi_0[i]*tmp1;
                                            phi_it_jt_2+=wi_1[i]*tmp1;
                                        }
                                    }
                                    a0=PI_IN_SQ*r_inv;
                                    rho_it_jt_0*=a0;
                                    rho_it_jt_1*=a0;
                                    rho_it_jt_1-=rho_it_jt_0*r_inv;
                                    rho_it_jt_2*=-a0;
                                    
                                    rho_jt_it_0*=a0;
                                    rho_jt_it_1*=a0;
                                    rho_jt_it_1-=rho_jt_it_0*r_inv;
                                    rho_jt_it_2*=-a0;
                                    
                                    phi_it_jt_0*=a0;
                                    phi_it_jt_1*=a0;
                                    phi_it_jt_1-=phi_it_jt_0*r_inv;
                                    phi_it_jt_2*=-a0;
                                    
                                    
                                    rho_phi[istart+type2rho_pair_ij[itype][jtype]]=rho_it_jt_0;
                                    rho_phi[istart+type2rho_pair_ji[jtype][itype]]=rho_jt_it_0;
                                    rho_phi[istart+type2phi_pair_ij[itype][jtype]]=phi_it_jt_0;
                                    
                                    drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]=rho_it_jt_1;
                                    drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]=rho_jt_it_1;
                                    drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]]=phi_it_jt_1;
                                    
                                    drho_phi_dalpha[istart+type2rho_pair_ij[itype][jtype]]=rho_it_jt_2;
                                    drho_phi_dalpha[istart+type2rho_pair_ji[jtype][itype]]=rho_jt_it_2;
                                    drho_phi_dalpha[istart+type2phi_pair_ij[itype][jtype]]=phi_it_jt_2;
                                    
                                    E[iatm*no_types+itype]+=c[jatm*no_types+jtype]*rho_jt_it_0;
                                    
                                    if (jatm<natms)
                                    {
                                        E[jatm*no_types+jtype]+=c[iatm*no_types+itype]*rho_it_jt_0;
                                        nrgy_strss[0]+=c[iatm*no_types+itype]*c[jatm*no_types+jtype]*phi_it_jt_0;
                                        
                                    }
                                    else
                                    {
                                        nrgy_strss[0]+=0.5*c[iatm*no_types+itype]*c[jatm*no_types+jtype]*phi_it_jt_0;
                                        
                                    }
                                }
                            }
                            
                        }
                    }
                    
                }
            }
            
            
            
            istart+=stride;
        }
        
        
        c_iv=1.0;
        for(int itype=0;itype<no_types; itype++)
        {
            if(c[no_types*iatm+itype]>=0.0)
            {
                c_iv-=c[iatm*no_types+itype];
                p=E[iatm*no_types+itype]*drho_inv;
                m=static_cast<int>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                p2=p*p;
                p3=p2*p;
                p4=p3*p;
                coef=F_arr[itype][m];
                
                tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                tmp1=(4.0*coef[4]*p3+3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*drho_inv;
                
                
                if(E[iatm*no_types+itype]>rho_max)
                    tmp0+=tmp1*(E[iatm*no_types+itype]-rho_max);
                
                E[iatm*no_types+itype]=tmp0;
                dE[iatm*no_types+itype]=tmp1;
                
                nrgy_strss[0]+=c[iatm*no_types+itype]*E[iatm*no_types+itype];
                
                nrgy_strss[0]+=kbT*calc_ent(c[iatm*no_types+itype]);
                
                nrgy_strss[0]-=3.0*kbT*c[iatm*no_types+itype]*log(x[icomp+3+itype])
                -c[iatm*no_types+itype]*c_0[itype];
                
                mu[iatm*no_types+itype]=tmp0-3.0*kbT*log(x[icomp+3+itype]);
                
                f[icomp+3+itype]+=3.0*kbT*c[iatm*no_types+itype]/x[icomp+3+itype];
            }
        }
        nrgy_strss[0]+=kbT*calc_ent(c_iv);
        
    }
    
    atoms->update_ph(dE_n);
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+no_types)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+no_types)*jatm;
            
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            r_inv=1.0/sqrt(dx0*dx0+dx1*dx1+dx2*dx2);
            
            for(int itype=0;itype<no_types; itype++)
            {
                for(int jtype=0;jtype<no_types; jtype++)
                {
                    if(c[no_types*iatm+itype]>=0.0 && c[no_types*jatm+jtype]>=0.0)
                    {
                        fpair=-(drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm*no_types+itype]
                                +drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm*no_types+jtype]
                                +drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]])
                                *c[iatm*no_types+itype]*c[jatm*no_types+jtype]*r_inv;
                        
                        apair=-(drho_phi_dalpha[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm*no_types+itype]
                                +drho_phi_dalpha[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm*no_types+jtype]
                                +drho_phi_dalpha[istart+type2phi_pair_ij[itype][jtype]])
                                *c[iatm*no_types+itype]*c[jatm*no_types+jtype]/
                                sqrt(x[icomp+3+itype]*x[icomp+3+itype]+x[jcomp+3+jtype]*x[jcomp+3+jtype]);
                        
                        mu[iatm*no_types+itype]+=c[jatm*no_types+jtype]*(drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]]
                        +drho_phi_dalpha[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm*no_types+jtype]);
                        
                        if(jatm<natms)
                            mu[jatm*no_types+jtype]+=c[iatm*no_types+itype]*(drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]]
                            +drho_phi_dalpha[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm*no_types+itype]);
                        
                        if(apair!=0.0 || fpair!=0.0)
                        {
                            f[icomp]+=dx0*fpair;
                            f[icomp+1]+=dx1*fpair;
                            f[icomp+2]+=dx2*fpair;
                            f[icomp+3+itype]+=apair*x[icomp+3+itype];
                            
                            if(jatm<natms)
                            {
                                f[jcomp]-=dx0*fpair;
                                f[jcomp+1]-=dx1*fpair;
                                f[jcomp+2]-=dx2*fpair;
                                f[jcomp+3+jtype]+=apair*x[jcomp+3+jtype];
                                
                                if(st_clc==1)
                                {
                                    nrgy_strss[1]-=fpair*dx0*dx0;
                                    nrgy_strss[2]-=fpair*dx1*dx1;
                                    nrgy_strss[3]-=fpair*dx2*dx2;
                                    nrgy_strss[4]-=fpair*dx1*dx2;
                                    nrgy_strss[5]-=fpair*dx2*dx0;
                                    nrgy_strss[6]-=fpair*dx0*dx1;
                                }
                                else if(st_clc==2)
                                {
                                    nrgy_strss[1]-=fpair*dx0*(image[3*jatm]-image[3*iatm]);
                                    nrgy_strss[2]-=fpair*dx1*(image[3*jatm+1]-image[3*iatm+1]);
                                    nrgy_strss[3]-=fpair*dx2*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[4]-=fpair*dx1*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[5]-=fpair*dx0*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[6]-=fpair*dx0*(image[3*jatm+1]-image[3*iatm+1]);
                                }
                            }
                            else
                            {
                                if(st_clc==1)
                                {
                                    nrgy_strss[1]-=0.5*fpair*dx0*dx0;
                                    nrgy_strss[2]-=0.5*fpair*dx1*dx1;
                                    nrgy_strss[3]-=0.5*fpair*dx2*dx2;
                                    nrgy_strss[4]-=0.5*fpair*dx1*dx2;
                                    nrgy_strss[5]-=0.5*fpair*dx2*dx0;
                                    nrgy_strss[6]-=0.5*fpair*dx0*dx1;
                                }
                                else if(st_clc==2)
                                {
                                    nrgy_strss[1]-=0.5*fpair*dx0*(image[3*jatm]-image[3*iatm]);
                                    nrgy_strss[2]-=0.5*fpair*dx1*(image[3*jatm+1]-image[3*iatm+1]);
                                    nrgy_strss[3]-=0.5*fpair*dx2*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[4]-=0.5*fpair*dx1*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[5]-=0.5*fpair*dx0*(image[3*jatm+2]-image[3*iatm+2]);
                                    nrgy_strss[6]-=0.5*fpair*dx0*(image[3*jatm+1]-image[3*iatm+1]);
                                }
                            }
                        }
                    }
                }
            }
            istart+=stride;
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
type0 ForceField_eam_dmd_hg_0::energy_calc()
{
    
    type0 en=0.0,en_tot;
    
    type0* x;
    atoms->vectors[x_n]->ret(x);
    type0* E;
    atoms->vectors[E_n]->ret(E);
    type0* c;
    atoms->vectors[c_n]->ret(c);
    
    int iatm,jatm;
    
    int icomp,jcomp;
    type0 dx0,dx1,dx2,rsq;
    type0 alpha;
    type0 r=0.0,r_inv=0.0;
    int m;
    type0* coef;
    type0 alpha_inv,tmp0,tmp1,rtmp,c_iv;
    type0 upper,lower;
    
    type0 p,p2,p3,p4;
    type0 phi_it_jt_0;
    type0 rho_it_jt_0;
    type0 rho_jt_it_0;
    type0 a0;
    
    
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    int natms=atoms->natms;
    for(int i=0;i<natms*no_types;i++) E[i]=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        
        icomp=(3+no_types)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+no_types)*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            if(rsq<cut_sq_mod_0)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                
                for(int itype=0;itype<no_types; itype++)
                {
                    for(int jtype=0;jtype<no_types; jtype++)
                    {
                        if(c[no_types*iatm+itype]>=0.0 && c[no_types*jatm+jtype]>=0.0)
                        {
                            alpha=sqrt(x[icomp+3+itype]*x[icomp+3+itype]+x[jcomp+3+jtype]*x[jcomp+3+jtype]);
                            
                            if(alpha<alpha_max)
                            {
                                alpha_inv=1.0/alpha;
                                
                                upper=(r+rc)*alpha_inv;
                                lower=(r-rc)*alpha_inv;
                                
                                if(lower<xi[no_i-1])
                                {
                                    phi_it_jt_0=0.0;
                                    rho_it_jt_0=0.0;
                                    rho_jt_it_0=0.0;
                                    
                                    for(int i=0;i<no_i;i++)
                                    {
                                        if(xi[i]>lower && xi[i]<upper)
                                        {
                                            rtmp=r-xi[i]*alpha;
                                            
                                            p=fabs(rtmp)*dr_inv;
                                            m=static_cast<int>(p);
                                            m=MIN(m,nr-2);
                                            p-=m;
                                            p=MIN(p,1.0);
                                            p2=p*p;
                                            p3=p2*p;
                                            
                                            coef=rho_r_arr[type2rho[itype][jtype]][m];
                                            tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                            if(rtmp<0.0)
                                                tmp0*=-1.0;
                                            
                                            rho_it_jt_0+=wi_0[i]*tmp0;
                                            
                                            
                                            if(itype!=jtype)
                                            {
                                                coef=rho_r_arr[type2rho[jtype][itype]][m];
                                                tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                                if(rtmp<0.0)
                                                    tmp0*=-1.0;
                                            }
                                            
                                            rho_jt_it_0+=wi_0[i]*tmp0;
                                            
                                            
                                            
                                            coef=phi_r_arr[type2phi[itype][jtype]][m];
                                            tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                            if(rtmp<0.0)
                                                tmp0*=-1.0;
                                            
                                            phi_it_jt_0+=wi_0[i]*tmp0;
                                            
                                        }
                                    }
                                    
                                    a0=PI_IN_SQ*r_inv;
                                    
                                    rho_it_jt_0*=a0;
                                    rho_jt_it_0*=a0;
                                    phi_it_jt_0*=a0;
                                    
                                    
                                    
                                    E[iatm*no_types+itype]+=c[jatm*no_types+jtype]*rho_jt_it_0;
                                    
                                    if (jatm<natms)
                                    {
                                        E[jatm*no_types+jtype]+=c[iatm*no_types+itype]*rho_it_jt_0;
                                        en+=c[iatm*no_types+itype]*c[jatm*no_types+jtype]*phi_it_jt_0;
                                        
                                    }
                                    else
                                    {
                                        en+=0.5*c[iatm*no_types+itype]*c[jatm*no_types+jtype]*phi_it_jt_0;
                                        
                                    }
                                }
                            }
                            
                        }
                        
                    }
                    
                }
            }
            
            
        }
        
        c_iv=1.0;
        for(int itype=0;itype<no_types; itype++)
        {
            if(c[no_types*iatm+itype]>=0.0)
            {
                c_iv-=c[no_types*iatm+itype];
                p=E[iatm*no_types+itype]*drho_inv;
                m=static_cast<int>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                p2=p*p;
                p3=p2*p;
                p4=p3*p;
                coef=F_arr[itype][m];
                
                tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                
                
                if(E[iatm*no_types+itype]>rho_max)
                {
                    tmp1=(4.0*coef[4]*p3+3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*drho_inv;
                    tmp0+=tmp1*(E[iatm*no_types+itype]-rho_max);
                }
                
                E[iatm*no_types+itype]=tmp0;
                
                en+=c[iatm*no_types+itype]*E[iatm*no_types+itype];
                
                en+=kbT*calc_ent(c[iatm*no_types+itype]);
                
                en-=3.0*kbT*c[iatm*no_types+itype]*log(x[icomp+3+itype])
                -c[iatm*no_types+itype]*c_0[itype];
            }
        }
        en+=kbT*calc_ent(c_iv);
        
    }
    
    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    
    return en_tot;
}
/*--------------------------------------------
 init before running
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::init()
{
    
    for (int i=0;i<no_types*(no_types+1)/2;i++)
        cut_sq[i]=cut_sq_mod_0;
    
    x_n=atoms->find("x");
    f_n=atoms->find("f");
    c_n=atoms->find("c");
    
    neighbor->pair_wise=1;
    
    E_n=atoms->add<type0>(1,no_types,"E");
    dE_n=atoms->add<type0>(1,no_types,"dE");
    mu_n=atoms->add<type0>(1,no_types,"mu");
    s_n=atoms->add<type0>(1,no_types,"s");
    crd_n=atoms->add<type0>(1,no_types,"crd");
    if(CRD_ENBL)
        t_n=atoms->add<type0>(1,2*no_types,"t");
    else
        t_n=atoms->add<type0>(1,no_types,"t");
    
    type0* mu;
    atoms->vectors[mu_n]->ret(mu);
    for(int i=0;i<no_types*atoms->natms;i++)
        mu[i]=0.0;
    
}
/*--------------------------------------------
 fin after running
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::fin()
{
    atoms->del(t_n);
    atoms->del(crd_n);
    atoms->del(s_n);
    atoms->del(mu_n);
    atoms->del(dE_n);
    atoms->del(E_n);
    
    if(max_pairs)
    {
        delete [] rho_phi;
        delete [] drho_phi_dr;
        delete [] drho_phi_dalpha;
        max_pairs=0;
    }
    
}
/*--------------------------------------------
 ff_coef no_i alpha_min alpha_max T kB hbar
 SetFL Ni_u3.eam
 hbar 6.5821192815e-16 eVs
 kB 8.617332478e-5 eV/K
 T 300.0 K
 mass conversion from amu to eVs^2/A^2:
 1.0364269184093291236e-28
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::coef(int narg,char** arg)
{
    
    if (narg!=2)
        error->abort("wrong coeff command "
        "for ff dmd");

    read_file(arg[1]);
    
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::read_file(char* file_name)
{
    
    int narg;
    int tot_no_types,lngth;
    int icmp;
    int str_lngth;
    type0 tmp,alpha_coef=0.0;
    int tmp_i;
    
    
    FILE* fp=NULL;
    char** arg=NULL;
    char** file_list=NULL;
    int* file_chk;
    int* type_ref;
    char* line;
    CREATE1D(line,MAXCHAR);
    char* strtmp;
    CREATE1D(strtmp,MAXCHAR);
    
    
    
    if(atoms->my_p_no==0)
    {
        fp=fopen(file_name,"r");
        if(fp==NULL)
            error->abort("ff eam_dmd file %s not found",file_name);
    }
    
    narg=0;
    while(narg==0)
    {
        if(read_line(fp,line)==-1)
            error->abort("%s file ended immaturely",file_name);
        
        narg=mapp->parse_line(line,arg);
    }
    
    tot_no_types=narg;
    
    if(tot_no_types<no_types)
        error->abort("the number of atoms in %s file"
                     " is less than the number of atom types"
                     " present in the system",file_name);
    
    
    
    CREATE1D(type_ref,tot_no_types);
    for(int ityp=0;ityp<tot_no_types;ityp++)
        type_ref[ityp]=atom_types->find_type_exist(arg[ityp]);
    
    for(int i=0;i<narg;i++)
        delete [] arg[i];
    if(narg)
        delete [] arg;
    
    CREATE1D(file_list,no_types);
    CREATE1D(file_chk,no_types);
    for(int i=0;i<no_types;i++)
        file_chk[i]=0;
    for(int i=0;i<no_types;i++)
        rsq_crd[i]=f_t[i]=Q_nrm[i]=-1.0;
    
    
    alpha_min=alpha_max=0.0;
    kbT=-1.0;
    
    lngth=read_line(fp,line);
    while(lngth!=-1)
    {
        if(sscanf(line,"degree = %d",&tmp_i)==1)
        {
            alpha_coef=set_weight_abs(tmp_i);
        }
        else if(sscanf(line,"alpha_min = %lf",&tmp)==1)
        {
            if(tmp<=0.0)
                error->abort("minimum alpha for ff eam_dmd "
                             "should be greater than 0.0");
            alpha_min=tmp;
        }
        else if(sscanf(line,"alpha_max = %lf",&tmp)==1)
        {
            if(tmp<=0.0)
                error->abort("maximum alpha for ff eam_dmd "
                             "should be greater than 0.0");
            alpha_max=tmp;
        }
        else if(sscanf(line,"T = %lf",&tmp)==1)
        {
            if(tmp<=0.0)
                error->abort("temperature for ff eam_dmd "
                             "should be greater than 0.0");
            
            type0 kb=8.617332478e-5;
            type0 hbar=6.5821192815e-16;
            type0 mass;
            type0 deb_l;
            
            for(int i=0;i<no_types;i++)
            {
                mass=atom_types->mass[i];
                mass*=1.0364269184093291236e-28;
                deb_l=hbar*hbar*2.0/(mass*kb*tmp);
                c_0[i]=1.5*kb*tmp*(log(deb_l)-1.0);
            }
            
            kbT=kb*tmp;
            beta=1.0/kbT;
            
        }
        else if(sscanf(line,"r_crd(%d) = %lf",&icmp,&tmp)==2)
        {
            if(icmp<0 || icmp>tot_no_types-1)
                error->abort("wrong component in %s "
                "file for r_crd(%d)",file_name,icmp);
            
            if(tmp<=0.0)
                error->abort("r_crd(%d) for ff eam_dmd "
                "in %s file should be greater than 0.0",icmp,file_name);
            
            if(type_ref[icmp]!=-1)
            {
                rsq_crd[type_ref[icmp]]=tmp*tmp;
            }
        }
        else if(sscanf(line,"f_t(%d) = %lf",&icmp,&tmp)==2)
        {
            if(icmp<0 || icmp>tot_no_types-1)
                error->abort("wrong component in %s "
                             "file for f_t(%d)",file_name,icmp);
            
            if(tmp<=0.0)
                error->abort("f_t(%d) for ff eam_dmd "
                             "in %s file should be greater than 0.0",icmp,file_name);
            
            if(type_ref[icmp]!=-1)
            {
                f_t[type_ref[icmp]]=tmp;
            }
        }
        else if(sscanf(line,"Q(%d) = %lf",&icmp,&tmp)==2)
        {
            if(icmp<0 || icmp>tot_no_types-1)
                error->abort("wrong component in %s "
                             "file for Q(%d)",file_name,icmp);
            
            if(tmp<=0.0)
                error->abort("Q(%d) for ff eam_dmd "
                             "in %s file should be greater than 0.0",icmp,file_name);
            
            if(type_ref[icmp]!=-1)
            {
                Q_nrm[type_ref[icmp]]=tmp;
            }
        }
        else if(sscanf(line,"file_format = %s",strtmp)==1)
        {
            if(strcmp(strtmp,"FS")==0)
                eam_mode=FINNIS_FL;
            else if(strcmp(strtmp,"SetFL")==0)
                eam_mode=SET_FL;
            else if(strcmp(strtmp,"FuncFL")==0)
                eam_mode=FUNC_FL;
            else
                error->abort("unknown file format "
                             "for ff eam_dmd: %s",strtmp);
        }
        else if(sscanf(line,"file(%d) = %s",&icmp,strtmp)==2)
        {
            if(icmp<0 || icmp>tot_no_types-1)
                error->abort("wrong component in %s "
                             "file for file(%d)",file_name,icmp);
            if(type_ref[icmp]!=-1)
            {
                str_lngth=static_cast<int>(strlen(strtmp))+1;
                CREATE1D(file_list[type_ref[icmp]],str_lngth);
                memcpy(file_list[type_ref[icmp]],strtmp,str_lngth*sizeof(char));
            }
        }
        
        
        lngth=read_line(fp,line);
    }
    
    
    
    if(no_i==0)
    {
        error->abort("degree was not set by %s file ",file_name);
    }
    if(kbT==-1.0)
    {
        error->abort("T was not set by %s file ",file_name);
    }
    /*
     if(alpha_min==0.0)
     {
     error->abort("minimum alpha was not set by %s file ",file_name);
     }
     */
    if(alpha_max==0.0)
    {
        error->abort("miximum alpha was not set by %s file ",file_name);
    }
    for(int i=0;i<no_types;i++)
    {
        if(f_t[i]==-1.0)
            error->abort("f_t(%s) was not set by %s file "
            ,atom_types->atom_names[i],file_name);
        if(rsq_crd[i]==-1.0)
            error->abort("r_crd(%s) was not set by %s file "
            ,atom_types->atom_names[i],file_name);
        if(Q_nrm[i]==-1.0)
            error->abort("Q(%s) was not set by %s file "
                         ,atom_types->atom_names[i],file_name);
    }
    max_rsq_crd=0.0;
    for(int i=0;i<no_types;i++)
        max_rsq_crd=MAX(max_rsq_crd,rsq_crd[i]);
    
    for(int i=0;i<no_types;i++)
        Q_nrm[i]/=f_t[i]*sqrt(f_t[i]);
    
    
    if(alpha_max<=alpha_min)
        error->abort("maximum alpha should be less "
                     "than minimum alpha for ff eam_dmd");
    
    
    int icomp;
    type0* x;
    atoms->vectors[0]->ret(x);
    type0 alpha_bound=alpha_max/sqrt(2.0);
    int* id;
    int id_n=atoms->find("id");
    atoms->vectors[id_n]->ret(id);
    for(int iatm=0;iatm<atoms->natms;iatm++)
    {
        icomp=(3+no_types)*iatm+3;
        for(int ityp=0;ityp<no_types;ityp++)
            if(x[icomp+ityp]>=alpha_bound)
                error->abort("alpha_%d of atomn %d is greater than "
                             "the maximum allowed %lf (alpha_max/sqrt(2.0)) set "
                             "by %s file",ityp,id[iatm],alpha_bound,file_name);
    }
    
    
    clean_up();
    
    if(eam_mode==FINNIS_FL)
    {
        if(file_list[0]==0)
            error->abort("file(%s) was not set by %s file "
                         ,atom_types->atom_names[0],file_name);
        
        set_fs(1,file_list);
    }
    else if(eam_mode==SET_FL)
    {
        if(file_list[0]==0)
            error->abort("file(%s) was not set by %s file "
                         ,atom_types->atom_names[0],file_name);
        
        set_setfl(1,file_list);
    }
    else if(eam_mode==FUNC_FL)
    {
        for(int i=0;i<no_types;i++)
            if(file_chk[i]==0)
                error->abort("file(%s) was not set by %s file "
                             ,atom_types->atom_names[i],file_name);
        
        set_funcfl(no_types,file_list);
    }
    set_arrays();
    
    rc=(static_cast<type0>(nr)-1.0)*dr;
    rho_max=(static_cast<type0>(nrho)-1.0)*drho;
    
    alpha_min=0.0;
    if(no_i!=1)
    {
        if(alpha_max==0.0)
            alpha_max=alpha_coef*rc;
        else
            alpha_max=MIN(alpha_max,alpha_coef*rc);
    }
    
    cut_sq_0=rc*rc;
    mod_rc=rc+xi[no_i-1]*alpha_max;
    cut_sq_mod_0=mod_rc*mod_rc;
    if(max_rsq_crd>cut_sq_mod_0)
        error->abort("r_crd set by %s file must be less than %lf",file_name,mod_rc);
    
    if(no_types)
    {
        for(int i=0;i<no_types;i++)
            if(file_chk[i])
                delete [] file_list[i];
        delete [] file_list;
        
        delete [] file_chk;
        
    }
    
    delete [] line;
    delete [] strtmp;
    if(tot_no_types)
        delete [] type_ref;
}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::set_funcfl(int no_files
,char** file_names)
{
    if(no_files!=no_types)
        error->abort("for FuncFL mode number of"
                     " ff eam_dmd files should be equal to the number"
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
                error->abort("ff eam_dmd file %s not found",file_names[ityp]);
        }
        
        for(int i=0;i<2;i++)
            if(read_line(fp,line)==-1)
                error->abort("%s file ended immaturely",file_names[ityp]);
        
        narg=mapp->parse_line(line,arg);
        if(narg!=4)
            error->abort("invalid line in %s file: %s",file_names[ityp],line);
        
        mass=atof(arg[1]);
        
        if(mass!=atom_types->mass[ityp])
            error->abort("mass of element %s in %s file (%lf) does not match "
                         "the mass that is already assigned to the element (%lf), make "
                         "sure that the right sequence of ff eam_dmd files are "
                         "used",atom_types->atom_names[ityp],file_names[ityp],mass
                         ,atom_types->mass[ityp]);
        
        for(int i=0;i<narg;i++)
            delete arg[i];
        if(narg)
            delete [] arg;
        
        
        
        
        if(read_line(fp,line)==-1)
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
            if(read_line(fp,line)==-1)
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
            rho_r_arr[type2rho[ityp][0]][i][0]=interpolate(tmp_rho[ityp],nrs[ityp],p,k)*static_cast<type0>(i)*dr;
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
void ForceField_eam_dmd_hg_0::set_setfl(int no_files
                                        ,char** file_names)
{
    if(no_files!=1)
        error->abort("for SetFL mode number of"
                     " ff eam_dmd files should be one");
    
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
        if(read_line(fp,line)==-1)
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
    
    
    if(read_line(fp,line)==-1)
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
    
    
    if(read_line(fp,line)==-1)
        error->abort("%s file ended immaturely",file_names[0]);
    
    
    int ipos=0;
    int tot=tot_no_types*(nrho+nr)+tot_no_types*(tot_no_types+1)*nr/2;
    
    type0* tmp;
    CREATE1D(tmp,tot);
    while (ipos<tot)
    {
        if(read_line(fp,line)==-1)
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
                rho_r_arr[component][i][0]=tmp[ipos+i]*static_cast<type0>(i)*dr;
            
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
void ForceField_eam_dmd_hg_0::set_fs(int no_files
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
        if(read_line(fp,line)==-1)
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
    
    
    if(read_line(fp,line)==-1)
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
        if(read_line(fp,line)==-1)
            error->abort("%s file ended immaturely",file_names[0]);
        ipos_loc=0;
        while (ipos_loc<tot_loc)
        {
            if(read_line(fp,line)==-1)
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
        if(read_line(fp,line)==-1)
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
                        rho_r_arr[component][i][0]=tmp[ipos+i]*static_cast<type0>(i)*dr;
                    
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
int ForceField_eam_dmd_hg_0::read_line(FILE* file
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
void ForceField_eam_dmd_hg_0::clean_up()
{
    if(allocated==0)
        return;
    
    for(int i=0;i<no_types;i++)
    {
        delete [] type2phi[i];
        delete [] type2rho[i];
        delete [] type2rho_pair_ij[i];
        delete [] type2rho_pair_ji[i];
        delete [] type2phi_pair_ij[i];
        delete [] type2phi_pair_ji[i];
    }
    if(no_types)
    {
        delete [] type2phi;
        delete [] type2rho;
        delete [] type2rho_pair_ij;
        delete [] type2rho_pair_ji;
        delete [] type2phi_pair_ij;
        delete [] type2phi_pair_ji;
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
                delete [] rho_r_arr[ityp][i];
            delete [] rho_r_arr[ityp];
        }
        delete [] rho_r_arr;
        
        
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
                delete [] rho_r_arr[ityp][i];
            delete [] rho_r_arr[ityp];
        }
        delete [] rho_r_arr;
        
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
void ForceField_eam_dmd_hg_0::allocate()
{
    int no_types=atom_types->no_types;
    
    CREATE2D(type2phi,no_types,no_types);
    CREATE2D(type2rho,no_types,no_types);
    CREATE2D(type2rho_pair_ij,no_types,no_types);
    CREATE2D(type2rho_pair_ji,no_types,no_types);
    CREATE2D(type2phi_pair_ij,no_types,no_types);
    CREATE2D(type2phi_pair_ji,no_types,no_types);
    
    for(int ityp=0;ityp<no_types;ityp++)
        for(int jtyp=0;jtyp<no_types;jtyp++)
        {
            type2phi[ityp][jtyp]=COMP(ityp,jtyp);
            type2rho[ityp][jtyp]=jtyp*no_types+ityp;
        }
    
    
    stride=0;
    for(int ityp=0;ityp<no_types;ityp++)
    {
        type2rho_pair_ij[ityp][ityp]=type2rho_pair_ji[ityp][ityp]=stride;
        stride++;
    }
    
    for(int ityp=0;ityp<no_types;ityp++)
        for(int jtyp=0;jtyp<no_types;jtyp++)
        {
            if(ityp!=jtyp)
            {
                type2rho_pair_ij[ityp][jtyp]=stride;
                stride++;
            }
        }
    
    for(int ityp=0;ityp<no_types;ityp++)
        for(int jtyp=0;jtyp<no_types;jtyp++)
        {
            if(ityp!=jtyp)
            {
                type2rho_pair_ji[ityp][jtyp]=stride;
                stride++;
            }
        }
    
    for(int ityp=0;ityp<no_types;ityp++)
        for(int jtyp=0;jtyp<no_types;jtyp++)
        {
            type2phi_pair_ij[ityp][jtyp]=type2phi_pair_ji[jtyp][ityp]=stride;
            stride++;
            
        }
    
    
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        
        CREATE1D(phi_r_arr,no_types*(no_types+1)/2);
        
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
            CREATE1D(phi_r_arr[ityp],nr);
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(phi_r_arr[ityp][i],4);
        
        CREATE1D(rho_r_arr,no_types*no_types);
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(rho_r_arr[ityp],nr);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(rho_r_arr[type2rho[ityp][0]][i],4);
        
        
        for(int ityp=0;ityp<no_types;ityp++)
            for(int jtyp=0;jtyp<no_types;jtyp++)
                rho_r_arr[type2rho[ityp][jtyp]]
                =rho_r_arr[type2rho[ityp][0]];
        
        CREATE1D(F_arr,no_types);
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(F_arr[ityp],nrho);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<nrho;i++)
                CREATE1D(F_arr[ityp][i],5);
        
    }
    else if(eam_mode==FINNIS_FL)
    {
        
        CREATE1D(phi_r_arr,no_types*(no_types+1)/2);
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
            CREATE1D(phi_r_arr[ityp],nr);
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(phi_r_arr[ityp][i],4);
        
        CREATE1D(rho_r_arr,no_types*no_types);
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            CREATE1D(rho_r_arr[ityp],nr);
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(rho_r_arr[ityp][i],4);
        
        CREATE1D(F_arr,no_types);
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(F_arr[ityp],nrho);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<nrho;i++)
                CREATE1D(F_arr[ityp][i],5);
    }
    
    allocated=1;
}
/*--------------------------------------------
 set the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::set_arrays()
{
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        for(int i=0;i<no_types;i++)
            interpolate_m(nrho,drho,F_arr[i]);
        
        for(int i=0;i<no_types;i++)
            interpolate(nr,dr,rho_r_arr[i]);
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<i+1;j++)
                interpolate(nr,dr,phi_r_arr[type2phi[i][j]]);
    }
    else if(eam_mode==FINNIS_FL)
    {
        for(int i=0;i<no_types;i++)
            interpolate_m(nrho,drho,F_arr[i]);
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<no_types;j++)
                interpolate(nr,dr,rho_r_arr[type2rho[i][j]]);
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<i+1;j++)
                interpolate(nr,dr,phi_r_arr[type2phi[i][j]]);
    }
}
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::interpolate(int n
,type0 delta,type0** spline)
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

}
/*--------------------------------------------
 interpolate for a single array
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::interpolate(type0* arr
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
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::interpolate_m(int n
,type0 delta,type0** spline)
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
    
    
    
    for (int i=0;i<n-1;i++)
    {
        spline[i][4]=(spline[i+1][2]-spline[i][2])/6.0-0.5*spline[i][3];
    }
    spline[n-1][4]=0.0;
    

}
/*--------------------------------------------
 Gaussian-Hermite quadrature weights and
 abscissas for 1 to 14 points
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::set_weight_abs(int n)
{
    
    if(n<1)
        error->abort("number of gaussian "
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
    
    XMath* xmath= new XMath(mapp);
    xmath->quadrature_hg(n,xi,wi_0);
    delete xmath;
    
    
    
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
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::create_2nd_neigh_lst()
{
    if(OLD_NEIGH)
    {
        type0 dx0,dx1,dx2,rsq;
        int iatm,jatm,icomp,jcomp;
        int natms=atoms->natms;
        int** neighbor_list=neighbor->neighbor_list;
        int* neighbor_list_size=neighbor->neighbor_list_size;
        
        type0* x;
        atoms->vectors[x_n]->ret(x);
        
        
        int* tmp_neigh_list;
        int tmp_neigh_list_size=1024;
        CREATE1D(tmp_neigh_list,1024);
        
        if(neigh_lst_sz_sz)
        {
            for(int i=0;i<neigh_lst_sz_sz;i++)
                if(neigh_lst_sz[i])
                    delete [] neigh_lst[i];
            
            delete [] neigh_lst_sz;
            delete [] neigh_lst;
        }
        neigh_lst_sz_sz=natms;
        
        CREATE1D(neigh_lst,natms);
        CREATE1D(neigh_lst_sz,natms);
        
        for(int i=0;i<natms;i++)
            neigh_lst_sz[i]=0;
        
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=(3+no_types)*iatm;
            
            for(int j=0;j<neighbor_list_size[iatm];j++)
            {
                jatm=neighbor_list[iatm][j];
                jcomp=(3+no_types)*jatm;
                dx0=x[icomp]-x[jcomp];
                dx1=x[icomp+1]-x[jcomp+1];
                dx2=x[icomp+2]-x[jcomp+2];
                rsq=dx0*dx0+dx1*dx1+dx2*dx2;
                
                if(rsq<max_rsq_crd)
                {
                    if(neigh_lst_sz[iatm]+1>tmp_neigh_list_size)
                    {
                        GROW(tmp_neigh_list,tmp_neigh_list_size,tmp_neigh_list_size+1);
                        tmp_neigh_list_size++;
                    }
                    tmp_neigh_list[neigh_lst_sz[iatm]]=jatm;
                    neigh_lst_sz[iatm]++;
                    
                }
            }
            
            CREATE1D(neigh_lst[iatm],neigh_lst_sz[iatm]);
            memcpy(neigh_lst[iatm],tmp_neigh_list,neigh_lst_sz[iatm]*sizeof(int));
        }
        
        if(tmp_neigh_list_size)
            delete [] tmp_neigh_list;
    }
    else
    {
        type0 dx0,dx1,dx2,rsq;
        int iatm,jatm,icomp,jcomp,iicomp,jjcomp;
        int natms=atoms->natms;
        int** neighbor_list=neighbor->neighbor_list;
        int* neighbor_list_size=neighbor->neighbor_list_size;
        
        type0* x;
        atoms->vectors[x_n]->ret(x);
        type0* c;
        atoms->vectors[c_n]->ret(c);
        
        
        int** tmp_neigh_list;
        int* tmp_neigh_list_size;
        CREATE1D(tmp_neigh_list_size,no_types);
        CREATE2D(tmp_neigh_list,no_types,1024);
        
        for(int itype=0;itype<no_types;itype++)
            tmp_neigh_list_size[itype]=1024;
        
        if(neigh_lst_sz_sz)
        {
            for(int i=0;i<neigh_lst_sz_sz;i++)
                if(neigh_lst_sz[i])
                    delete [] neigh_lst[i];
            
            delete [] neigh_lst_sz;
            delete [] neigh_lst;
        }
        
        neigh_lst_sz_sz=natms*no_types;
        CREATE1D(neigh_lst,neigh_lst_sz_sz);
        CREATE1D(neigh_lst_sz,neigh_lst_sz_sz);
        
        for(int i=0;i<neigh_lst_sz_sz;i++)
            neigh_lst_sz[i]=0;
        
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=(3+no_types)*iatm;
            iicomp=no_types*iatm;
            for(int j=0;j<neighbor_list_size[iatm];j++)
            {
                jatm=neighbor_list[iatm][j];
                jcomp=(3+no_types)*jatm;
                jjcomp=no_types*jatm;
                dx0=x[icomp]-x[jcomp];
                dx1=x[icomp+1]-x[jcomp+1];
                dx2=x[icomp+2]-x[jcomp+2];
                rsq=dx0*dx0+dx1*dx1+dx2*dx2;
                
                
                for(int itype=0;itype<no_types;itype++)
                {
                    if(c[iicomp+itype]>=0.0 && c[jjcomp+itype]>=0.0 && rsq<rsq_crd[itype])
                    {
                        if(neigh_lst_sz[iicomp+itype]+1>tmp_neigh_list_size[itype])
                        {
                            GROW(tmp_neigh_list[itype],tmp_neigh_list_size[itype],tmp_neigh_list_size[itype]+1);
                            tmp_neigh_list_size[itype]++;
                        }
                        tmp_neigh_list[itype][neigh_lst_sz[iicomp+itype]]=jjcomp+itype;
                        neigh_lst_sz[iicomp+itype]++;
                    }
                }
            }
            
            for(int itype=0;itype<no_types;itype++)
            {
                CREATE1D(neigh_lst[iicomp+itype],neigh_lst_sz[iicomp+itype]);
                memcpy(neigh_lst[iicomp+itype],tmp_neigh_list[itype],neigh_lst_sz[iicomp+itype]*sizeof(int));
            }
        }
        
        for(int itype=0;itype<no_types;itype++)
            if(tmp_neigh_list_size[itype])
                delete [] tmp_neigh_list[itype];
        if(no_types)
        {
            delete [] tmp_neigh_list_size;
            delete [] tmp_neigh_list;
        }
    }
    
}
/*--------------------------------------------
 return M_{ij}^{\alpha}
 this fucntion should be symmetric wrt fi &
 fj
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::
mat(type0 fi,type0 crdi,type0 fj,type0 crdj
    ,int itype)
{
    
    if(fi>fj)
    {
        if(fi-fj>f_t[itype])
        {
            return fi;
        }
        else
        {
            return fi+Q_nrm[itype]*(f_t[itype]-fi+fj)*sqrt(f_t[itype]-fi+fj);
        }
    }
    else
    {
        if(fj-fi>f_t[itype])
        {
            return fj;
        }
        else
        {
            return fj+Q_nrm[itype]*(f_t[itype]-fj+fi)*sqrt(f_t[itype]-fj+fi);
        }
    }
    
}
/*--------------------------------------------
 dmat(fi,fj)/dfi
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::
dmat0(type0 fi,type0 crdi,type0 fj,type0 crdj
      ,int itype)
{

    if(fi>fj)
    {
        if(fi-fj>f_t[itype])
        {
            return 1.0;
        }
        else
        {
            return 1.0-1.5*Q_nrm[itype]*sqrt(f_t[itype]-fi+fj);
        }
    }
    else
    {
        if(fj-fi>f_t[itype])
        {
            return 0.0;
        }
        else
        {
            return 1.5*Q_nrm[itype]*sqrt(f_t[itype]-fj+fi);
        }
    }
    
}
/*--------------------------------------------
 dmat(fi,fj)/dfi
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::
dmat1(type0 fi,type0 crdi,type0 fj,type0 crdj
,int itype)
{
    
    return 0.0;
}
/*--------------------------------------------
 return M_{ij}^{\alpha}
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::calc_ent(type0 x)
{
    type0 ans=x*log(x);
    if (isnan(ans))
        return 0.0;
    return ans;
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::g_calc(int chk
,type0 alpha,type0* a,type0* g,type0* en_st)
{
    if(CRD_ENBL)
        return g_calc_crd(chk,alpha,a,g,en_st);
    else
        return g_calc_ncrd(chk,alpha,a,g,en_st);
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::c_d_calc(int chk,type0* en_st)
{
    if(CRD_ENBL)
        return c_d_calc_crd(chk,en_st);
    else
        return c_d_calc_ncrd(chk,en_st);
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::c_dd_norm()
{
    if(CRD_ENBL)
        return c_dd_norm_crd();
    else
        return c_dd_norm_ncrd();
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::c_d_calc_crd(int chk,type0* en_st)
{
    
    type0* x;
    atoms->vectors[x_n]->ret(x);
    type0* dE;
    atoms->vectors[dE_n]->ret(dE);
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    type0* mu;
    atoms->vectors[mu_n]->ret(mu);
    type0* crd;
    atoms->vectors[crd_n]->ret(crd);
    
    type0 fi,fj,exp_fi,exp_fj;
    type0 p,p2,p3,p4,tmp0,tmp1,s_ij,crdi,crdj,w_ij,w_ji,c_iv,c_jv;
    type0* coef;
    type0 dx0,dx1,dx2,r_inv,fpair;
    int m;
    int iatm,jatm;
    int icomp,jcomp,istart;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    int natms=atoms->natms;
    
    for (int i=0;i<7;i++)
        nrgy_strss[i]=0.0;
    
    for(int i=0;i<natms*no_types;i++)
        mu[i]=c_d[i]=0.0;
    
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        crd[iatm+itype]+=c[jcomp+jtype];
                        mu[icomp+itype]+=c[jcomp+jtype]*rho_phi[istart+type2rho_pair_ji[jtype][itype]];
                        
                        if(jatm<natms)
                        {
                            crd[jatm+jtype]+=c[icomp+itype];
                            mu[jcomp+jtype]+=c[icomp+itype]*rho_phi[istart+type2rho_pair_ij[itype][jtype]];
                        }
                    }
                }
            }
            
            istart+=stride;
            
        }
        
        for(int itype=0;itype<no_types; itype++)
        {
            if(c[icomp+itype]>=0.0)
            {
                p=mu[iatm*no_types+itype]*drho_inv;
                m=static_cast<int>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                p2=p*p;
                p3=p2*p;
                p4=p3*p;
                coef=F_arr[itype][m];
                
                
                tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                tmp1=(4.0*coef[4]*p3+3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*drho_inv;
                
                if(mu[iatm*no_types+itype]>rho_max)
                    tmp0+=tmp1*(mu[icomp+itype]-rho_max);
                
                mu[icomp+itype]=tmp0+c_0[itype]-3.0*kbT*log(x[(3+no_types)*iatm+3+itype]);
                dE[icomp+itype]=tmp1;
            }
        }
    }
    
    atoms->update_ph(dE_n);
    
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        c_iv=1.0;
        for(int itype=0;itype<no_types;itype++)
            if(c[icomp+itype]>=0.0)
            {
                nrgy_strss[0]+=c[icomp+itype]*mu[icomp+itype]+kbT*calc_ent(c[icomp+itype]);
                c_iv-=c[icomp+itype];
            }
        nrgy_strss[0]+=kbT*calc_ent(c_iv);
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            dx0=x[iatm*(3+no_types)]-x[jatm*(3+no_types)];
            dx1=x[iatm*(3+no_types)+1]-x[jatm*(3+no_types)+1];
            dx2=x[iatm*(3+no_types)+2]-x[jatm*(3+no_types)+2];
            r_inv=1.0/sqrt(dx0*dx0+dx1*dx1+dx2*dx2);
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        mu[icomp+itype]+=c[jcomp+jtype]*(rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype]
                                                         +rho_phi[istart+type2phi_pair_ij[itype][jtype]]);
                        
                        if(jatm<natms)
                            mu[jcomp+jtype]+=c[icomp+itype]*(rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                                                             +rho_phi[istart+type2phi_pair_ji[jtype][itype]]);
                        
                        fpair=-(drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                                +drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype]
                                +drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]])
                        *c[icomp+itype]*c[jcomp+jtype]*r_inv;
                        
                        if(jatm<natms)
                        {
                            nrgy_strss[0]+=c[icomp+itype]*c[jcomp+jtype]
                            *rho_phi[istart+type2phi_pair_ij[itype][jtype]];
                            nrgy_strss[1]-=fpair*dx0*dx0;
                            nrgy_strss[2]-=fpair*dx1*dx1;
                            nrgy_strss[3]-=fpair*dx2*dx2;
                            nrgy_strss[4]-=fpair*dx1*dx2;
                            nrgy_strss[5]-=fpair*dx2*dx0;
                            nrgy_strss[6]-=fpair*dx0*dx1;
                        }
                        else
                        {
                            nrgy_strss[0]+=0.5*c[icomp+itype]*c[jcomp+jtype]
                            *rho_phi[istart+type2phi_pair_ij[itype][jtype]];
                            nrgy_strss[1]-=0.5*fpair*dx0*dx0;
                            nrgy_strss[2]-=0.5*fpair*dx1*dx1;
                            nrgy_strss[3]-=0.5*fpair*dx2*dx2;
                            nrgy_strss[4]-=0.5*fpair*dx1*dx2;
                            nrgy_strss[5]-=0.5*fpair*dx2*dx0;
                            nrgy_strss[6]-=0.5*fpair*dx0*dx1;
                        }
                    }
                }
            }
            istart+=stride;
        }
    }
    
    atoms->update_ph(mu_n);
    atoms->update_ph(crd_n);
    
    
   if(OLD_NEIGH)
   {
       for(iatm=0;iatm<natms;iatm++)
       {
           icomp=no_types*iatm;
           
           c_iv=1.0;
           for(int itype=0;itype<no_types;itype++)
               if(c[icomp+itype]>=0.0)
                   c_iv-=c[icomp+itype];
           
           for(int j=0;j<neigh_lst_sz[iatm];j++)
           {
               jatm=neigh_lst[iatm][j];
               jcomp=no_types*jatm;
               
               c_jv=1.0;
               for(int itype=0;itype<no_types;itype++)
                   if(c[jcomp+itype]>=0.0)
                       c_jv-=c[jcomp+itype];
               
               for(int itype=0;itype<no_types;itype++)
               {
                   if(c[icomp+itype]>=0.0 && c[jcomp+itype]>=0.0)
                   {
                       crdi=crd[icomp+itype];
                       crdj=crd[jcomp+itype];
                       fi=mu[icomp+itype];
                       fj=mu[jcomp+itype];
                       
                       
                       exp_fi=exp(beta*(fi-mat(fi,crdi,fj,crdj,itype)));
                       exp_fj=exp(beta*(fj-mat(fi,crdi,fj,crdj,itype)));
                       w_ij=-c[icomp+itype]*c_jv*exp_fi;
                       w_ji=-c[jcomp+itype]*c_iv*exp_fj;
                       
                       s_ij=w_ij-w_ji;
                       
                       
                       c_d[icomp+itype]+=s_ij;
                       if(jatm<natms)
                       {
                           c_d[jcomp+itype]-=s_ij;
                       }
                   }
               }
           }
           
       }
   }
   else
   {
       for(iatm=0;iatm<natms;iatm++)
       {
           icomp=no_types*iatm;
           c_iv=1.0;
           for(int itype=0;itype<no_types;itype++)
               if(c[icomp+itype]>=0.0)
                   c_iv-=c[icomp+itype];
           
           for(int itype=0;itype<no_types;itype++)
           {
               for(int j=0;j<neigh_lst_sz[icomp];j++)
               {
                   c_jv=1.0;
                   jcomp=neigh_lst[icomp][j]-itype;
                   for(int jtype=0;jtype<no_types;jtype++)
                   {
                       if(c[jcomp]>=0.0)
                           c_jv-=c[jcomp];
                       jcomp++;
                   }
                   jcomp=neigh_lst[icomp][j];
                   
                   crdi=crd[icomp];
                   crdj=crd[jcomp];
                   fi=mu[icomp];
                   fj=mu[jcomp];
                   
                   tmp0=mat(fi,crdi,fj,crdj,itype);
                   
                   exp_fi=exp(beta*(fi-tmp0));
                   exp_fj=exp(beta*(fj-tmp0));
                   w_ij=-c[icomp]*c_jv*exp_fi;
                   w_ji=-c[jcomp]*c_iv*exp_fj;
                   
                   s_ij=w_ij-w_ji;
                   
                   c_d[icomp]+=s_ij;
                   if(jcomp<natms*no_types)
                       c_d[jcomp]-=s_ij;
               }
               
               icomp++;
           }
       }
   }
    
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_hg_0::c_d_calc_ncrd(int chk,type0* en_st)
{
    
    type0* x;
    atoms->vectors[x_n]->ret(x);
    type0* dE;
    atoms->vectors[dE_n]->ret(dE);
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    type0* mu;
    atoms->vectors[mu_n]->ret(mu);
    
    
    type0 fi,fj,exp_fi,exp_fj;
    type0 p,p2,p3,p4,tmp0,tmp1,s_ij,w_ij,w_ji,c_iv,c_jv;
    type0* coef;
    type0 dx0,dx1,dx2,r_inv,fpair;
    int m;
    int iatm,jatm;
    int icomp,jcomp,istart;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    int natms=atoms->natms;
    
    if(chk)
        for (int i=0;i<7;i++)
            nrgy_strss[i]=0.0;
    
    for(int i=0;i<natms*no_types;i++)
        mu[i]=c_d[i]=0.0;
    
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        mu[icomp+itype]+=c[jcomp+jtype]*rho_phi[istart+type2rho_pair_ji[jtype][itype]];

                        
                        if(jatm<natms)
                            mu[jcomp+jtype]+=c[icomp+itype]*rho_phi[istart+type2rho_pair_ij[itype][jtype]];
                    }
                }
            }
            
            istart+=stride;
            
        }
        
        for(int itype=0;itype<no_types; itype++)
        {
            if(c[icomp+itype]>=0.0)
            {
                p=mu[iatm*no_types+itype]*drho_inv;
                m=static_cast<int>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                p2=p*p;
                p3=p2*p;
                p4=p3*p;
                coef=F_arr[itype][m];
                
                
                tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                tmp1=(4.0*coef[4]*p3+3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*drho_inv;

                if(mu[iatm*no_types+itype]>rho_max)
                    tmp0+=tmp1*(mu[icomp+itype]-rho_max);
                
                mu[icomp+itype]=tmp0+c_0[itype]-3.0*kbT*log(x[(3+no_types)*iatm+3+itype]);
                dE[icomp+itype]=tmp1;
            }
        }
    }
    
    atoms->update_ph(dE_n);
    

    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        if(chk)
        {
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                {
                    nrgy_strss[0]+=c[icomp+itype]*mu[icomp+itype]+kbT*calc_ent(c[icomp+itype]);
                    c_iv-=c[icomp+itype];
                }
            nrgy_strss[0]+=kbT*calc_ent(c_iv);
        }
        
        
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            dx0=x[iatm*(3+no_types)]-x[jatm*(3+no_types)];
            dx1=x[iatm*(3+no_types)+1]-x[jatm*(3+no_types)+1];
            dx2=x[iatm*(3+no_types)+2]-x[jatm*(3+no_types)+2];
            r_inv=1.0/sqrt(dx0*dx0+dx1*dx1+dx2*dx2);
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        mu[icomp+itype]+=c[jcomp+jtype]*(rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype]
                                                         +rho_phi[istart+type2phi_pair_ij[itype][jtype]]);
                        
                        if(jatm<natms)
                            mu[jcomp+jtype]+=c[icomp+itype]*(rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                                                             +rho_phi[istart+type2phi_pair_ji[jtype][itype]]);
                        
                        
                        if(chk)
                        {
                            fpair=-(drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                                    +drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype]
                                    +drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]])
                            *c[icomp+itype]*c[jcomp+jtype]*r_inv;
                            if(jatm<natms)
                            {
                                nrgy_strss[0]+=c[icomp+itype]*c[jcomp+jtype]
                                *rho_phi[istart+type2phi_pair_ij[itype][jtype]];
                                nrgy_strss[1]-=fpair*dx0*dx0;
                                nrgy_strss[2]-=fpair*dx1*dx1;
                                nrgy_strss[3]-=fpair*dx2*dx2;
                                nrgy_strss[4]-=fpair*dx1*dx2;
                                nrgy_strss[5]-=fpair*dx2*dx0;
                                nrgy_strss[6]-=fpair*dx0*dx1;
                            }
                            else
                            {
                                nrgy_strss[0]+=0.5*c[icomp+itype]*c[jcomp+jtype]
                                *rho_phi[istart+type2phi_pair_ij[itype][jtype]];
                                nrgy_strss[1]-=0.5*fpair*dx0*dx0;
                                nrgy_strss[2]-=0.5*fpair*dx1*dx1;
                                nrgy_strss[3]-=0.5*fpair*dx2*dx2;
                                nrgy_strss[4]-=0.5*fpair*dx1*dx2;
                                nrgy_strss[5]-=0.5*fpair*dx2*dx0;
                                nrgy_strss[6]-=0.5*fpair*dx0*dx1;
                            }
                        }
                        
                    }
                }
            }
            istart+=stride;
        }
    }
    
    atoms->update_ph(mu_n);
    
    if(OLD_NEIGH)
    {
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int j=0;j<neigh_lst_sz[iatm];j++)
            {
                jatm=neigh_lst[iatm][j];
                jcomp=no_types*jatm;
                
                c_jv=1.0;
                for(int itype=0;itype<no_types;itype++)
                    if(c[jcomp+itype]>=0.0)
                        c_jv-=c[jcomp+itype];
                
                for(int itype=0;itype<no_types;itype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+itype]>=0.0)
                    {
                        fi=mu[icomp+itype];
                        fj=mu[jcomp+itype];
                        
                        
                        exp_fi=exp(beta*(fi-mat(fi,0.0,fj,0.0,itype)));
                        exp_fj=exp(beta*(fj-mat(fj,0.0,fi,0.0,itype)));
                        w_ij=-c[icomp+itype]*c_jv*exp_fi;
                        w_ji=-c[jcomp+itype]*c_iv*exp_fj;
                        
                        s_ij=w_ij-w_ji;
                        c_d[icomp+itype]+=s_ij;
                        
                        if(jatm<natms)
                            c_d[jcomp+itype]-=s_ij;
                    }
                }
            }
            
        }
    }
    else
    {
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int j=0;j<neigh_lst_sz[icomp];j++)
                {
                    c_jv=1.0;
                    jcomp=neigh_lst[icomp][j]-itype;
                    for(int jtype=0;jtype<no_types;jtype++)
                    {
                        if(c[jcomp]>=0.0)
                            c_jv-=c[jcomp];
                        jcomp++;
                    }
                    jcomp=neigh_lst[icomp][j];
                    
                    fi=mu[icomp];
                    fj=mu[jcomp];
                    
                    tmp0=mat(fi,0.0,fj,0.0,itype);
                    
                    exp_fi=exp(beta*(fi-tmp0));
                    exp_fj=exp(beta*(fj-tmp0));
                    w_ij=-c[icomp]*c_jv*exp_fi;
                    w_ji=-c[jcomp]*c_iv*exp_fj;
                    
                    s_ij=w_ij-w_ji;
                    
                    c_d[icomp]+=s_ij;
                    if(jcomp<natms*no_types)
                        c_d[jcomp]-=s_ij;
                }
                
                icomp++;
            }
        }
    }
    
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::g_calc_crd(int chk
,type0 alpha,type0* a,type0* g,type0* en_st)
{
    type0* x;
    atoms->vectors[x_n]->ret(x);
    type0* E;
    atoms->vectors[E_n]->ret(E);
    type0* dE;
    atoms->vectors[dE_n]->ret(dE);
    type0* ddE=E;
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    type0* mu;
    atoms->vectors[mu_n]->ret(mu);
    type0* s;
    atoms->vectors[s_n]->ret(s);
    type0* t;
    atoms->vectors[t_n]->ret(t);
    type0* crd;
    atoms->vectors[crd_n]->ret(crd);
    int v_n=crd_n;
    type0* v=crd;
    
    type0 inner,ans;
    type0 fi,fj,exp_fi,exp_fj;
    type0 p,p2,p3,p4,tmp0,tmp1,tmp2,s_ij;
    type0 w_ij,w_ji;
    type0 crdi,crdj,c_iv,c_jv;
    type0* coef;
    int m;
    int iatm,jatm;
    int icomp,jcomp,istart;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    if (chk)
        for (int i=0;i<7;i++)
            en_st[i]=nrgy_strss[i]=0.0;
    
    int natms=atoms->natms;
    for(int i=0;i<natms*no_types;i++)
        crd[i]=c_d[i]=g[i]=E[i]=mu[i]=0.0;
    
    
    /*
     beginning of level 1
     */
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        crd[iatm+itype]+=c[jcomp+jtype];
                        E[icomp+itype]+=c[jcomp+jtype]*rho_phi[istart+type2rho_pair_ji[jtype][itype]];
                        
                        if(jatm<natms)
                        {
                            crd[jatm+jtype]+=c[icomp+itype];
                            E[jcomp+jtype]+=c[icomp+itype]*rho_phi[istart+type2rho_pair_ij[itype][jtype]];
                        }
                        
                    }
                }
            }
            
            istart+=stride;
        }
        
        c_iv=1.0;
        for(int itype=0;itype<no_types; itype++)
        {
            if(c[icomp+itype]>=0.0)
            {
                c_iv-=c[icomp+itype];
                p=E[iatm*no_types+itype]*drho_inv;
                m=static_cast<int>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                p2=p*p;
                p3=p2*p;
                p4=p3*p;
                coef=F_arr[itype][m];
                
                
                tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                tmp1=(4.0*coef[4]*p3+3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*drho_inv;
                tmp2=(12.0*coef[4]*p2+6.0*coef[3]*p+2.0*coef[2])*drho_inv*drho_inv;
                
                if(E[iatm*no_types+itype]>rho_max)
                {
                    tmp0+=tmp1*(E[icomp+itype]-rho_max);
                    tmp2=0.0;
                }
                
                mu[icomp+itype]=tmp0+c_0[itype]-3.0*kbT*log(x[(3+no_types)*iatm+3+itype]);
                nrgy_strss[0]+=c[icomp+itype]*mu[icomp+itype];
                nrgy_strss[0]+=kbT*calc_ent(c[icomp+itype]);
                dE[icomp+itype]=tmp1;
                ddE[icomp+itype]=tmp2;
            }
        }
        if(chk)
            nrgy_strss[0]+=kbT*calc_ent(c_iv);
    }
    
    
    /*
     beginning of level 2
     */
    
    atoms->update_ph(dE_n);
    
    type0 r_inv,dx0,dx1,dx2,fpair;
    
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            dx0=x[iatm*(3+no_types)]-x[jatm*(3+no_types)];
            dx1=x[iatm*(3+no_types)+1]-x[jatm*(3+no_types)+1];
            dx2=x[iatm*(3+no_types)+2]-x[jatm*(3+no_types)+2];
            r_inv=1.0/sqrt(dx0*dx0+dx1*dx1+dx2*dx2);
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[jcomp+jtype]>=0.0 && c[icomp+itype]>=0.0)
                    {
                        mu[icomp+itype]+=c[jcomp+jtype]*(rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype]
                        +rho_phi[istart+type2phi_pair_ij[itype][jtype]]);
                        
                        if(jatm<natms)
                        {
                            mu[jcomp+jtype]+=c[icomp+itype]*(rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                            +rho_phi[istart+type2phi_pair_ji[jtype][itype]]);
                        }
                        
                        if (chk)
                        {
                            fpair=-(drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                                    +drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype]
                                    +drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]])
                            *c[icomp+itype]*c[jcomp+jtype]*r_inv;
                            
                            if(jatm<natms)
                            {
                                nrgy_strss[0]+=c[icomp+itype]*c[jcomp+jtype]
                                *rho_phi[istart+type2phi_pair_ij[itype][jtype]];
                                nrgy_strss[1]-=fpair*dx0*dx0;
                                nrgy_strss[2]-=fpair*dx1*dx1;
                                nrgy_strss[3]-=fpair*dx2*dx2;
                                nrgy_strss[4]-=fpair*dx1*dx2;
                                nrgy_strss[5]-=fpair*dx2*dx0;
                                nrgy_strss[6]-=fpair*dx0*dx1;
                            }
                            else
                            {
                                nrgy_strss[0]+=0.5*c[icomp+itype]*c[jcomp+jtype]
                                *rho_phi[istart+type2phi_pair_ij[itype][jtype]];
                                nrgy_strss[1]-=0.5*fpair*dx0*dx0;
                                nrgy_strss[2]-=0.5*fpair*dx1*dx1;
                                nrgy_strss[3]-=0.5*fpair*dx2*dx2;
                                nrgy_strss[4]-=0.5*fpair*dx1*dx2;
                                nrgy_strss[5]-=0.5*fpair*dx2*dx0;
                                nrgy_strss[6]-=0.5*fpair*dx0*dx1;
                            }
                            
                        }
                    }
                }
            }
            istart+=stride;
        }
    }
    
    /*
     beginning of level 3
     */
    atoms->update_ph(mu_n);
    atoms->update_ph(crd_n);
    
    for(int i=0;i<natms*no_types;i++)
    {
        if(c[i]>=0.0)
        {
            s[i]=2.0*(c[i]+a[i]);
        }
    }
    
    if(OLD_NEIGH)
    {
        inner=0.0;
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int j=0;j<neigh_lst_sz[iatm];j++)
            {
                jatm=neigh_lst[iatm][j];
                jcomp=no_types*jatm;
                
                c_jv=1.0;
                for(int itype=0;itype<no_types;itype++)
                    if(c[jcomp+itype]>=0.0)
                        c_jv-=c[jcomp+itype];
                
                for(int itype=0;itype<no_types;itype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+itype]>=0.0)
                    {
                        crdi=crd[icomp+itype];
                        crdj=crd[jcomp+itype];
                        fi=mu[icomp+itype];
                        fj=mu[jcomp+itype];
                        
                        tmp0=mat(fi,crdi,fj,crdj,itype);
                        
                        exp_fi=exp(beta*(fi-tmp0));
                        exp_fj=exp(beta*(fj-tmp0));
                        w_ij=-c[icomp+itype]*c_jv*exp_fi;
                        w_ji=-c[jcomp+itype]*c_iv*exp_fj;
                        
                        s_ij=w_ij-w_ji;
                        
                        s[icomp+itype]-=2.0*alpha*s_ij;
                        c_d[icomp+itype]+=s_ij;
                        if(jatm<natms)
                        {
                            s[jcomp+itype]+=2.0*alpha*s_ij;
                            c_d[jcomp+itype]-=s_ij;
                        }
                    }
                }
            }
            for(int itype=0;itype<no_types;itype++)
            {
                if(c[icomp+itype]>=0.0)
                    inner+=s[icomp+itype]*s[icomp+itype];
            }
        }
        ans=0.0;
        MPI_Allreduce(&inner,&ans,1,MPI_TYPE0,MPI_SUM,world);
        ans*=0.25;
        ans=sqrt(ans);
    }
    else
    {
        inner=0.0;
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int j=0;j<neigh_lst_sz[icomp];j++)
                {
                    c_jv=1.0;
                    jcomp=neigh_lst[icomp][j]-itype;
                    for(int jtype=0;jtype<no_types;jtype++)
                    {
                        if(c[jcomp]>=0.0)
                            c_jv-=c[jcomp];
                        jcomp++;
                    }
                    jcomp=neigh_lst[icomp][j];
                    
                    crdi=crd[icomp];
                    crdj=crd[jcomp];
                    fi=mu[icomp];
                    fj=mu[jcomp];
                    
                    tmp0=mat(fi,crdi,fj,crdj,itype);
                    
                    exp_fi=exp(beta*(fi-tmp0));
                    exp_fj=exp(beta*(fj-tmp0));
                    w_ij=-c[icomp]*c_jv*exp_fi;
                    w_ji=-c[jcomp]*c_iv*exp_fj;
                    
                    s_ij=w_ij-w_ji;
                    
                    s[icomp]-=2.0*alpha*s_ij;
                    c_d[icomp]+=s_ij;
                    if(jcomp<natms*no_types)
                    {
                        s[jcomp]+=2.0*alpha*s_ij;
                        c_d[jcomp]-=s_ij;
                    }
                    
                    
                }
                inner+=s[icomp]*s[icomp];
                icomp++;
            }
        }
        ans=0.0;
        MPI_Allreduce(&inner,&ans,1,MPI_TYPE0,MPI_SUM,world);
        ans*=0.25;
        ans=sqrt(ans);
    }

    if(chk)
        MPI_Allreduce(nrgy_strss,en_st,7,MPI_TYPE0,MPI_SUM,world);
    
    
    if(chk==0)
        return ans;
    
    /*
     beginning of level 4
     */
    atoms->update_ph(s_n);
    
    for(int i=0;i<natms*no_types;i++)
    {
        if(c[i]>=0.0)
        {
            t[2*i]=t[2*i+1]=0.0;
            g[i]=-s[i];
        }
    }
    
    if(OLD_NEIGH)
    {
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int j=0;j<neigh_lst_sz[iatm];j++)
            {
                jatm=neigh_lst[iatm][j];
                jcomp=no_types*jatm;
                
                c_jv=1.0;
                for(int itype=0;itype<no_types;itype++)
                    if(c[jcomp+itype]>=0.0)
                        c_jv-=c[jcomp+itype];
                
                for(int itype=0;itype<no_types;itype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+itype]>=0.0)
                    {
                        s_ij=s[icomp+itype]-s[jcomp+itype];
                        fi=mu[icomp+itype];
                        fj=mu[jcomp+itype];
                        crdi=crd[icomp+itype];
                        crdj=crd[jcomp+itype];
                        tmp0=mat(fi,crdi,fj,crdj,itype);
                        exp_fi=exp(beta*(fi-tmp0));
                        exp_fj=exp(beta*(fj-tmp0));
                        w_ij=-c[icomp+itype]*c_jv*exp_fi;
                        w_ji=-c[jcomp+itype]*c_iv*exp_fj;
                        
                        g[icomp+itype]-=alpha*s_ij*(c_jv*exp_fi+c[jcomp+itype]*exp_fj);
                        t[2*(icomp+itype)]+=s_ij*(w_ij+(w_ji-w_ij)*dmat0(fi,crdi,fj,crdj,itype));
                        t[2*(icomp+itype)+1]+=s_ij*(w_ji-w_ij)*dmat1(fi,crdi,fj,crdj,itype);
                        
                        if(jatm<natms)
                        {
                            g[jcomp+itype]+=alpha*s_ij*(c_iv*exp_fj+c[icomp+itype]*exp_fi);
                            t[2*(jcomp+itype)]-=s_ij*(w_ji+(w_ij-w_ji)*dmat0(fj,crdj,fi,crdi,itype));
                            t[2*(jcomp+itype)+1]-=s_ij*(w_ij-w_ji)*dmat1(fj,crdj,fi,crdi,itype);
                        }
                    }
                }
            }
        }

    }
    else
    {
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int j=0;j<neigh_lst_sz[icomp];j++)
                {
                    c_jv=1.0;
                    jcomp=neigh_lst[icomp][j]-itype;
                    for(int jtype=0;jtype<no_types;jtype++)
                    {
                        if(c[jcomp]>=0.0)
                            c_jv-=c[jcomp];
                        jcomp++;
                    }
                    jcomp=neigh_lst[icomp][j];
                    
                    crdi=crd[icomp];
                    crdj=crd[jcomp];
                    fi=mu[icomp];
                    fj=mu[jcomp];
                    tmp0=mat(fi,crdi,fj,crdj,itype);
                    exp_fi=exp(beta*(fi-tmp0));
                    exp_fj=exp(beta*(fj-tmp0));
                    w_ij=-c[icomp]*c_jv*exp_fi;
                    w_ji=-c[jcomp]*c_iv*exp_fj;
                    s_ij=s[icomp]-s[jcomp];
                    
                    
                    g[icomp]-=alpha*s_ij*(c_jv*exp_fi+c[jcomp]*exp_fj);
                    t[2*icomp]+=s_ij*(w_ij+(w_ji-w_ij)*dmat0(fi,crdi,fj,crdj,itype));
                    t[2*icomp+1]+=s_ij*(w_ji-w_ij)*dmat1(fi,crdi,fj,crdj,itype);
                    
                    if(jcomp<natms*no_types)
                    {
                        g[jcomp]+=alpha*s_ij*(c_iv*exp_fj+c[icomp]*exp_fi);
                        t[2*jcomp]-=s_ij*(w_ji+(w_ij-w_ji)*dmat0(fj,crdj,fi,crdi,itype));
                        t[2*jcomp+1]-=s_ij*(w_ij-w_ji)*dmat1(fj,crdj,fi,crdi,itype);
                    }
                }
                icomp++;
            }
        }
    }

    
    /*
     beginning of level 5
     */
    atoms->update_ph(t_n);
    
    for(int i=0;i<no_types*natms;i++)
        v[i]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        tmp0=rho_phi[istart+type2phi_pair_ij[itype][jtype]]
                        +rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                        +rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype];
                        tmp0*=beta*alpha;
                        
                        v[icomp+itype]+=rho_phi[istart+type2rho_pair_ji[jtype][itype]]*t[2*(jcomp+jtype)];
                        g[icomp+itype]+=tmp0*t[2*(jcomp+jtype)]+alpha*beta*t[2*(jcomp+jtype)+1];
                        if(jatm<natms)
                        {
                            v[jcomp+jtype]+=rho_phi[istart+type2rho_pair_ij[itype][jtype]]*t[2*(icomp+itype)];
                            g[jcomp+jtype]+=tmp0*t[2*(icomp+itype)]+alpha*beta*t[2*(icomp+itype)+1];
                        }
                    }
                }
            }
            istart+=stride;
        }
        for(int itype=0;itype<no_types;itype++)
        {
            v[icomp+itype]*=c[icomp+itype]*ddE[icomp+itype];
        }
    }
    
    /*
     beginning of level 6
     */
    atoms->update_ph(v_n);
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[jcomp+jtype]>=0.0 && c[icomp+itype]>=0.0)
                    {
                        g[icomp+itype]+=alpha*beta*rho_phi[istart+type2rho_pair_ij[itype][jtype]]
                        *v[jcomp+jtype];
                        if(jatm<natms)
                            g[jcomp+jtype]+=alpha*beta*rho_phi[istart+type2rho_pair_ji[jtype][itype]]
                            *v[icomp+itype];
                    }
                }
            }
            istart+=stride;
        }
        
    }
    
    for(int i=0;i<no_types*natms;i++)
        if(c[i]>=0.0)
            g[i]/=2.0*ans;
    
    
    return ans;
    
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::g_calc_ncrd(int chk
,type0 alpha,type0* a,type0* g,type0* en_st)
{
    type0* x;
    atoms->vectors[x_n]->ret(x);
    type0* E;
    atoms->vectors[E_n]->ret(E);
    type0* dE;
    atoms->vectors[dE_n]->ret(dE);
    type0* ddE=E;
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* c_d;
    atoms->vectors[c_d_n]->ret(c_d);
    type0* mu;
    atoms->vectors[mu_n]->ret(mu);
    type0* s;
    atoms->vectors[s_n]->ret(s);
    type0* t;
    atoms->vectors[t_n]->ret(t);
    int v_n=crd_n;
    type0* v;
    atoms->vectors[v_n]->ret(v);
    
    type0 inner,ans;
    type0 fi,fj,exp_fi,exp_fj;
    type0 p,p2,p3,p4,tmp0,tmp1,tmp2,s_ij;
    type0 w_ij,w_ji;
    type0 c_iv,c_jv;
    type0* coef;
    int m;
    int iatm,jatm;
    int icomp,jcomp,istart;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    if (chk)
        for (int i=0;i<7;i++)
            en_st[i]=nrgy_strss[i]=0.0;
    
    int natms=atoms->natms;
    for(int i=0;i<natms*no_types;i++)
        c_d[i]=g[i]=E[i]=mu[i]=0.0;
    
    
    /*
     beginning of level 1
     */
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        E[icomp+itype]+=c[jcomp+jtype]*rho_phi[istart+type2rho_pair_ji[jtype][itype]];
                        
                        if(jatm<natms)
                        {
                            E[jcomp+jtype]+=c[icomp+itype]*rho_phi[istart+type2rho_pair_ij[itype][jtype]];
                        }
                        
                    }
                }
            }
            
            istart+=stride;
        }
        
        c_iv=1.0;
        for(int itype=0;itype<no_types; itype++)
        {
            if(c[icomp+itype]>=0.0)
            {
                c_iv-=c[icomp+itype];
                p=E[iatm*no_types+itype]*drho_inv;
                m=static_cast<int>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                p2=p*p;
                p3=p2*p;
                p4=p3*p;
                coef=F_arr[itype][m];
                
                
                tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                tmp1=(4.0*coef[4]*p3+3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*drho_inv;
                tmp2=(12.0*coef[4]*p2+6.0*coef[3]*p+2.0*coef[2])*drho_inv*drho_inv;
                
                if(E[iatm*no_types+itype]>rho_max)
                {
                    tmp0+=tmp1*(E[icomp+itype]-rho_max);
                    tmp2=0.0;
                }
                
                mu[icomp+itype]=tmp0+c_0[itype]-3.0*kbT*log(x[(3+no_types)*iatm+3+itype]);
                nrgy_strss[0]+=c[icomp+itype]*mu[icomp+itype];
                nrgy_strss[0]+=kbT*calc_ent(c[icomp+itype]);
                dE[icomp+itype]=tmp1;
                ddE[icomp+itype]=tmp2;
            }
        }
        if(chk)
            nrgy_strss[0]+=kbT*calc_ent(c_iv);
    }
    
    
    /*
     beginning of level 2
     */
    
    atoms->update_ph(dE_n);
    
    type0 r_inv,dx0,dx1,dx2,fpair;
    
    
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            dx0=x[iatm*(3+no_types)]-x[jatm*(3+no_types)];
            dx1=x[iatm*(3+no_types)+1]-x[jatm*(3+no_types)+1];
            dx2=x[iatm*(3+no_types)+2]-x[jatm*(3+no_types)+2];
            r_inv=1.0/sqrt(dx0*dx0+dx1*dx1+dx2*dx2);
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[jcomp+jtype]>=0.0 && c[icomp+itype]>=0.0)
                    {
                        
                        mu[icomp+itype]+=c[jcomp+jtype]*(rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype]
                                                         +rho_phi[istart+type2phi_pair_ij[itype][jtype]]);
                        
                        if(jatm<natms)
                        {
                            mu[jcomp+jtype]+=c[icomp+itype]*(rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                                                             +rho_phi[istart+type2phi_pair_ji[jtype][itype]]);
                        }
                        
                        if (chk)
                        {
                            fpair=-(drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                                    +drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype]
                                    +drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]])
                            *c[icomp+itype]*c[jcomp+jtype]*r_inv;
                            
                            if(jatm<natms)
                            {
                                
                                nrgy_strss[0]+=c[icomp+itype]*c[jcomp+jtype]
                                *rho_phi[istart+type2phi_pair_ij[itype][jtype]];
                                nrgy_strss[1]-=fpair*dx0*dx0;
                                nrgy_strss[2]-=fpair*dx1*dx1;
                                nrgy_strss[3]-=fpair*dx2*dx2;
                                nrgy_strss[4]-=fpair*dx1*dx2;
                                nrgy_strss[5]-=fpair*dx2*dx0;
                                nrgy_strss[6]-=fpair*dx0*dx1;
                            }
                            else
                            {
                                
                                nrgy_strss[0]+=0.5*c[icomp+itype]*c[jcomp+jtype]
                                *rho_phi[istart+type2phi_pair_ij[itype][jtype]];
                                nrgy_strss[1]-=0.5*fpair*dx0*dx0;
                                nrgy_strss[2]-=0.5*fpair*dx1*dx1;
                                nrgy_strss[3]-=0.5*fpair*dx2*dx2;
                                nrgy_strss[4]-=0.5*fpair*dx1*dx2;
                                nrgy_strss[5]-=0.5*fpair*dx2*dx0;
                                nrgy_strss[6]-=0.5*fpair*dx0*dx1;
                            }
                            
                        }
                    }
                }
            }
            istart+=stride;
        }
    }
    /*
     beginning of level 3
     */
    atoms->update_ph(mu_n);
    
    for(int i=0;i<natms*no_types;i++)
    {
        if(c[i]>=0.0)
        {
            s[i]=2.0*(c[i]+a[i]);
        }
    }
    if(OLD_NEIGH)
    {
        inner=0.0;
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int j=0;j<neigh_lst_sz[iatm];j++)
            {
                jatm=neigh_lst[iatm][j];
                jcomp=no_types*jatm;
                
                c_jv=1.0;
                for(int itype=0;itype<no_types;itype++)
                    if(c[jcomp+itype]>=0.0)
                        c_jv-=c[jcomp+itype];
                
                for(int itype=0;itype<no_types;itype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+itype]>=0.0)
                    {
                        fi=mu[icomp+itype];
                        fj=mu[jcomp+itype];
                        
                        tmp0=mat(fi,0.0,fj,0.0,itype);
                        
                        exp_fi=exp(beta*(fi-tmp0));
                        exp_fj=exp(beta*(fj-tmp0));
                        w_ij=-c[icomp+itype]*c_jv*exp_fi;
                        w_ji=-c[jcomp+itype]*c_iv*exp_fj;
                        
                        s_ij=w_ij-w_ji;
                        
                        s[icomp+itype]-=2.0*alpha*s_ij;
                        c_d[icomp+itype]+=s_ij;
                        if(jatm<natms)
                        {
                            s[jcomp+itype]+=2.0*alpha*s_ij;
                            c_d[jcomp+itype]-=s_ij;
                        }
                    }
                }
            }
            for(int itype=0;itype<no_types;itype++)
            {
                if(c[icomp+itype]>=0.0)
                    inner+=s[icomp+itype]*s[icomp+itype];
            }
        }
        ans=0.0;
        MPI_Allreduce(&inner,&ans,1,MPI_TYPE0,MPI_SUM,world);
        ans*=0.25;
        ans=sqrt(ans);

    }
    else
    {
        inner=0.0;
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int j=0;j<neigh_lst_sz[icomp];j++)
                {
                    c_jv=1.0;
                    jcomp=neigh_lst[icomp][j]-itype;
                    for(int jtype=0;jtype<no_types;jtype++)
                    {
                        if(c[jcomp]>=0.0)
                            c_jv-=c[jcomp];
                        jcomp++;
                    }
                    jcomp=neigh_lst[icomp][j];
                    
                    fi=mu[icomp];
                    fj=mu[jcomp];
                    
                    tmp0=mat(fi,0.0,fj,0.0,itype);
                    
                    exp_fi=exp(beta*(fi-tmp0));
                    exp_fj=exp(beta*(fj-tmp0));
                    w_ij=-c[icomp]*c_jv*exp_fi;
                    w_ji=-c[jcomp]*c_iv*exp_fj;
                    
                    s_ij=w_ij-w_ji;
                    
                    s[icomp]-=2.0*alpha*s_ij;
                    c_d[icomp]+=s_ij;
                    if(jcomp<natms*no_types)
                    {
                        s[jcomp]+=2.0*alpha*s_ij;
                        c_d[jcomp]-=s_ij;
                    }
                    
                    
                }
                inner+=s[icomp]*s[icomp];
                icomp++;
            }
        }
        ans=0.0;
        MPI_Allreduce(&inner,&ans,1,MPI_TYPE0,MPI_SUM,world);
        ans*=0.25;
        ans=sqrt(ans);

    }
    
    if(chk)
        MPI_Allreduce(nrgy_strss,en_st,7,MPI_TYPE0,MPI_SUM,world);
    
    
    
    if(chk==0)
        return ans;
    
    /*
     beginning of level 4
     */
    atoms->update_ph(s_n);
    
    for(int i=0;i<natms*no_types;i++)
    {
        if(c[i]>=0.0)
        {
            t[i]=0.0;
            g[i]=-s[i];
        }
    }
    
    
    if(OLD_NEIGH)
    {
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int j=0;j<neigh_lst_sz[iatm];j++)
            {
                jatm=neigh_lst[iatm][j];
                jcomp=no_types*jatm;
                
                c_jv=1.0;
                for(int itype=0;itype<no_types;itype++)
                    if(c[jcomp+itype]>=0.0)
                        c_jv-=c[jcomp+itype];
                
                for(int itype=0;itype<no_types;itype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+itype]>=0.0)
                    {
                        s_ij=s[icomp+itype]-s[jcomp+itype];
                        fi=mu[icomp+itype];
                        fj=mu[jcomp+itype];
                        tmp0=mat(fi,0.0,fj,0.0,itype);
                        exp_fi=exp(beta*(fi-tmp0));
                        exp_fj=exp(beta*(fj-tmp0));
                        w_ij=-c[icomp+itype]*c_jv*exp_fi;
                        w_ji=-c[jcomp+itype]*c_iv*exp_fj;
                        
                        g[icomp+itype]-=alpha*s_ij*(c_jv*exp_fi+c[jcomp+itype]*exp_fj);
                        t[icomp+itype]+=s_ij*(w_ij+(w_ji-w_ij)*dmat0(fi,0,fj,0,itype));
                        
                        
                        if(jatm<natms)
                        {
                            g[jcomp+itype]+=alpha*s_ij*(c_iv*exp_fj+c[icomp+itype]*exp_fi);
                            t[jcomp+itype]-=s_ij*(w_ji+(w_ij-w_ji)*dmat0(fj,0.0,fi,0.0,itype));
                        }
                    }
                }
            }
        }
    }
    else
    {
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int j=0;j<neigh_lst_sz[icomp];j++)
                {
                    c_jv=1.0;
                    jcomp=neigh_lst[icomp][j]-itype;
                    for(int jtype=0;jtype<no_types;jtype++)
                    {
                        if(c[jcomp]>=0.0)
                            c_jv-=c[jcomp];
                        jcomp++;
                    }
                    jcomp=neigh_lst[icomp][j];
                    
                    
                    s_ij=s[icomp]-s[jcomp];
                    fi=mu[icomp];
                    fj=mu[jcomp];
                    tmp0=mat(fi,0.0,fj,0.0,itype);
                    exp_fi=exp(beta*(fi-tmp0));
                    exp_fj=exp(beta*(fj-tmp0));
                    w_ij=-c[icomp]*c_jv*exp_fi;
                    w_ji=-c[jcomp]*c_iv*exp_fj;
                    
                    g[icomp]-=alpha*s_ij*(c_jv*exp_fi+c[jcomp]*exp_fj);
                    t[icomp]+=s_ij*(w_ij+(w_ji-w_ij)*dmat0(fi,0,fj,0,itype));
                    
                    
                    if(jcomp<natms*no_types)
                    {
                        g[jcomp]+=alpha*s_ij*(c_iv*exp_fj+c[icomp]*exp_fi);
                        t[jcomp]-=s_ij*(w_ji+(w_ij-w_ji)*dmat0(fj,0.0,fi,0.0,itype));
                    }
                }
                
                icomp++;
            }
        }
    }
     
    /*
     beginning of level 5
     */
    atoms->update_ph(t_n);
    
    for(int i=0;i<no_types*natms;i++)
        v[i]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        tmp0=rho_phi[istart+type2phi_pair_ij[itype][jtype]]
                        +rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                        +rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype];
                        tmp0*=beta*alpha;
                        
                        v[icomp+itype]+=rho_phi[istart+type2rho_pair_ji[jtype][itype]]*t[jcomp+jtype];
                        g[icomp+itype]+=tmp0*t[jcomp+jtype];
                        if(jatm<natms)
                        {
                            v[jcomp+jtype]+=rho_phi[istart+type2rho_pair_ij[itype][jtype]]*t[icomp+itype];
                            g[jcomp+jtype]+=tmp0*t[icomp+itype];
                        }
                    }
                }
            }
            istart+=stride;
        }
        for(int itype=0;itype<no_types;itype++)
        {
            v[icomp+itype]*=c[icomp+itype]*ddE[icomp+itype];
        }
    }
    
    /*
     beginning of level 6
     */
    atoms->update_ph(v_n);
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[jcomp+jtype]>=0.0 && c[icomp+itype]>=0.0)
                    {
                        g[icomp+itype]+=alpha*beta*rho_phi[istart+type2rho_pair_ij[itype][jtype]]
                        *v[jcomp+jtype];
                        if(jatm<natms)
                            g[jcomp+jtype]+=alpha*beta*rho_phi[istart+type2rho_pair_ji[jtype][itype]]
                            *v[icomp+itype];
                    }
                }
            }
            istart+=stride;
        }
        
    }
    
    for(int i=0;i<no_types*natms;i++)
        if(c[i]>=0.0)
            g[i]/=2.0*ans;
    
    
    return ans;
    
}
/*--------------------------------------------
 calculate norm of d^2c/dt^2
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::c_dd_norm_crd()
{
    type0* x;
    atoms->vectors[x_n]->ret(x);
    type0* dE;
    atoms->vectors[dE_n]->ret(dE);
    type0* ddE;
    atoms->vectors[E_n]->ret(ddE);
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* s;
    atoms->vectors[s_n]->ret(s);
    type0* mu;
    atoms->vectors[mu_n]->ret(mu);
    type0* crd;
    atoms->vectors[crd_n]->ret(crd);
    type0* t;
    atoms->vectors[t_n]->ret(t);
    
    type0 fi,fj,exp_fi,exp_fj,norm_lcl,norm;
    type0 p,p2,p3,p4,tmp0,tmp1,tmp2,s_ij,crdi,crdj,w_ij,w_ji,c_iv,c_jv;
    type0* coef;
    int m;
    int iatm,jatm;
    int icomp,jcomp,istart;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    int natms=atoms->natms;
    
    
    for(int i=0;i<natms*no_types;i++)
        mu[i]=s[i]=0.0;
    
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        crd[iatm+itype]+=c[jcomp+jtype];
                        mu[icomp+itype]+=c[jcomp+jtype]*rho_phi[istart+type2rho_pair_ji[jtype][itype]];
                        
                        if(jatm<natms)
                        {
                            crd[jatm+jtype]+=c[icomp+itype];
                            mu[jcomp+jtype]+=c[icomp+itype]*rho_phi[istart+type2rho_pair_ij[itype][jtype]];
                        }
                    }
                }
            }
            
            istart+=stride;
            
        }
        
        for(int itype=0;itype<no_types; itype++)
        {
            if(c[icomp+itype]>=0.0)
            {
                p=mu[iatm*no_types+itype]*drho_inv;
                m=static_cast<int>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                p2=p*p;
                p3=p2*p;
                p4=p3*p;
                coef=F_arr[itype][m];
                
                
                tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                tmp1=(4.0*coef[4]*p3+3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*drho_inv;
                tmp2=(12.0*coef[4]*p2+6.0*coef[3]*p+2.0*coef[2])*drho_inv*drho_inv;
                
                if(mu[iatm*no_types+itype]>rho_max)
                {
                    tmp0+=tmp1*(mu[icomp+itype]-rho_max);
                    tmp2=0.0;
                }
                
                mu[icomp+itype]=tmp0+c_0[itype]-3.0*kbT*log(x[(3+no_types)*iatm+3+itype]);
                dE[icomp+itype]=tmp1;
                ddE[icomp+itype]=tmp2;
            }
        }
    }
    
    atoms->update_ph(dE_n);
    
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        mu[icomp+itype]+=c[jcomp+jtype]*(rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype]
                            +rho_phi[istart+type2phi_pair_ij[itype][jtype]]);
                        
                        if(jatm<natms)
                            mu[jcomp+jtype]+=c[icomp+itype]*(rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                                +rho_phi[istart+type2phi_pair_ji[jtype][itype]]);

                    }
                }
            }
            istart+=stride;
        }
    }
    
    atoms->update_ph(mu_n);
    atoms->update_ph(crd_n);
    
    if(OLD_NEIGH)
    {
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int j=0;j<neigh_lst_sz[iatm];j++)
            {
                jatm=neigh_lst[iatm][j];
                jcomp=no_types*jatm;
                
                c_jv=1.0;
                for(int itype=0;itype<no_types;itype++)
                    if(c[jcomp+itype]>=0.0)
                        c_jv-=c[jcomp+itype];
                
                for(int itype=0;itype<no_types;itype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+itype]>=0.0)
                    {
                        crdi=crd[icomp+itype];
                        crdj=crd[jcomp+itype];
                        fi=mu[icomp+itype];
                        fj=mu[jcomp+itype];
                        
                        
                        exp_fi=exp(beta*(fi-mat(fi,crdi,fj,crdj,itype)));
                        exp_fj=exp(beta*(fj-mat(fi,crdi,fj,crdj,itype)));
                        w_ij=-c[icomp+itype]*c_jv*exp_fi;
                        w_ji=-c[jcomp+itype]*c_iv*exp_fj;
                        
                        s_ij=w_ij-w_ji;
                        
                        
                        s[icomp+itype]+=s_ij;
                        if(jatm<natms)
                        {
                            s[jcomp+itype]-=s_ij;
                        }
                    }
                }
            }
        }
    }
    else
    {
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int j=0;j<neigh_lst_sz[icomp];j++)
                {
                    c_jv=1.0;
                    jcomp=neigh_lst[icomp][j]-itype;
                    for(int jtype=0;jtype<no_types;jtype++)
                    {
                        if(c[jcomp]>=0.0)
                            c_jv-=c[jcomp];
                        jcomp++;
                    }
                    jcomp=neigh_lst[icomp][j];
                    
                    crdi=crd[icomp];
                    crdj=crd[jcomp];
                    fi=mu[icomp];
                    fj=mu[jcomp];
                    
                    tmp0=mat(fi,crdi,fj,crdj,itype);
                    
                    exp_fi=exp(beta*(fi-tmp0));
                    exp_fj=exp(beta*(fj-tmp0));
                    w_ij=-c[icomp]*c_jv*exp_fi;
                    w_ji=-c[jcomp]*c_iv*exp_fj;
                    
                    s_ij=w_ij-w_ji;
                    
                    s[icomp]+=s_ij;
                    if(jcomp<natms*no_types)
                    {
                        s[jcomp]-=s_ij;
                    }
                }
                icomp++;
            }
        }
    }
    

    atoms->update_ph(s_n);
    
    for(int i=0;i<2*natms*no_types;i++)
        t[i]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        
                        t[2*(icomp+itype)]+=rho_phi[istart+type2rho_pair_ji[jtype][itype]]*s[jcomp+jtype];
                        t[2*(icomp+itype)+1]+=s[jcomp+jtype];

                        if(jatm<natms)
                        {
                            t[2*(jcomp+jtype)]+=rho_phi[istart+type2rho_pair_ij[itype][jtype]]*s[icomp+itype];
                            t[2*(jcomp+jtype)+1]+=s[icomp+itype];
                        }
                    }
                }
            }
            istart+=stride;
        }
        
        for(int itype=0;itype<no_types;itype++)
            if(c[icomp+itype]>=0.0)
                t[2*(icomp+itype)]*=c[icomp+itype]*ddE[icomp+itype];
    }
    
    atoms->update_ph(t_n);
    
    
    for(int i=0;i<natms*no_types;i++)
        ddE[i]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        
                         tmp0=rho_phi[istart+type2phi_pair_ij[itype][jtype]]
                         +rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                         +rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype];

                        ddE[icomp+itype]+=tmp0*s[jcomp+jtype]+t[2*(jcomp+jtype)]*rho_phi[istart+type2rho_pair_ij[itype][jtype]];
                        
                        if(jatm<natms)
                        {
                            ddE[jcomp+jtype]+=tmp0*s[icomp+itype]+t[2*(icomp+itype)]*rho_phi[istart+type2rho_pair_ji[jtype][itype]];

                        }
                    }
                }
            }
            istart+=stride;
        }
        
    }
    atoms->update_ph(E_n);
    
    
    for(int i=0;i<natms*no_types;i++)
        dE[i]=0.0;
    
    if(OLD_NEIGH)
    {
        norm_lcl=0.0;
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int j=0;j<neigh_lst_sz[iatm];j++)
            {
                jatm=neigh_lst[iatm][j];
                jcomp=no_types*jatm;
                
                c_jv=1.0;
                for(int itype=0;itype<no_types;itype++)
                    if(c[jcomp+itype]>=0.0)
                        c_jv-=c[jcomp+itype];
                
                for(int itype=0;itype<no_types;itype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+itype]>=0.0)
                    {
                        crdi=crd[icomp+itype];
                        crdj=crd[jcomp+itype];
                        fi=mu[icomp+itype];
                        fj=mu[jcomp+itype];
                        
                        
                        exp_fi=exp(beta*(fi-mat(fi,crdi,fj,crdj,itype)));
                        exp_fj=exp(beta*(fj-mat(fi,crdi,fj,crdj,itype)));
                        w_ij=-c[icomp+itype]*c_jv*exp_fi;
                        w_ji=-c[jcomp+itype]*c_iv*exp_fj;
                        
                        s_ij=w_ij-w_ji;
                        
                        
                        tmp0=(s[jcomp+itype]*c[icomp+itype]-s[icomp+itype]*c_jv)*exp_fi
                        -(s[icomp+itype]*c[jcomp+itype]-s[jcomp+itype]*c_iv)*exp_fj;
                        
                        tmp0+=beta*(ddE[icomp+itype]*w_ij-ddE[jcomp+itype]*w_ji);
                        
                        tmp0+=-beta*(w_ij-w_ji)*(t[2*(icomp+itype)+1]*dmat1(fi,crdi,fj,crdj,itype)
                                                 +t[2*(jcomp+itype)+1]*dmat1(fj,crdj,fi,crdi,itype)
                                                 +ddE[icomp+itype]*dmat0(fi,crdi,fj,crdj,itype)
                                                 +ddE[jcomp+itype]*dmat0(fj,crdj,fi,crdi,itype));
                        
                        
                        dE[icomp+itype]+=tmp0;
                        if(jatm<natms)
                        {
                            dE[jcomp+itype]-=tmp0;
                        }
                    }
                }
            }
            
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    norm_lcl+=dE[icomp+itype]*dE[icomp+itype];
        }
        MPI_Allreduce(&norm_lcl,&tmp0,1,MPI_TYPE0,MPI_SUM,world);
        norm=sqrt(tmp0);
    }
    else
    {
        norm_lcl=0.0;
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int j=0;j<neigh_lst_sz[icomp];j++)
                {
                    c_jv=1.0;
                    jcomp=neigh_lst[icomp][j]-itype;
                    for(int jtype=0;jtype<no_types;jtype++)
                    {
                        if(c[jcomp]>=0.0)
                            c_jv-=c[jcomp];
                        jcomp++;
                    }
                    jcomp=neigh_lst[icomp][j];
                    
                    crdi=crd[icomp];
                    crdj=crd[jcomp];
                    fi=mu[icomp];
                    fj=mu[jcomp];
                    
                    tmp0=mat(fi,crdi,fj,crdj,itype);
                    
                    exp_fi=exp(beta*(fi-tmp0));
                    exp_fj=exp(beta*(fj-tmp0));
                    w_ij=-c[icomp]*c_jv*exp_fi;
                    w_ji=-c[jcomp]*c_iv*exp_fj;
                    
                    s_ij=w_ij-w_ji;
                    
                    tmp0=(s[jcomp]*c[icomp]-s[icomp]*c_jv)*exp_fi
                    -(s[icomp]*c[jcomp]-s[jcomp]*c_iv)*exp_fj;
                    
                    tmp0+=beta*(ddE[icomp]*w_ij-ddE[jcomp]*w_ji);
                    
                    tmp0+=-beta*(w_ij-w_ji)*(t[2*icomp+1]*dmat1(fi,crdi,fj,crdj,itype)
                                             +t[2*jcomp+1]*dmat1(fj,crdj,fi,crdi,itype)
                                             +ddE[icomp]*dmat0(fi,crdi,fj,crdj,itype)
                                             +ddE[jcomp]*dmat0(fj,crdj,fi,crdi,itype));
                    
                    
                    dE[icomp]+=tmp0;
                    if(jcomp<natms*no_types)
                    {
                        dE[jcomp]-=tmp0;
                    }
                }
                
                norm_lcl+=dE[icomp]*dE[icomp];
                icomp++;
            }
        }
        MPI_Allreduce(&norm_lcl,&tmp0,1,MPI_TYPE0,MPI_SUM,world);
        norm=sqrt(tmp0);
    }
    
    return norm;
    
}
/*--------------------------------------------
 calculate norm of d^2c/dt^2
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_0::c_dd_norm_ncrd()
{
    type0* x;
    atoms->vectors[x_n]->ret(x);
    type0* dE;
    atoms->vectors[dE_n]->ret(dE);
    type0* ddE;
    atoms->vectors[E_n]->ret(ddE);
    type0* c;
    atoms->vectors[c_n]->ret(c);
    type0* s;
    atoms->vectors[s_n]->ret(s);
    type0* mu;
    atoms->vectors[mu_n]->ret(mu);
    type0* t;
    atoms->vectors[t_n]->ret(t);
    
    type0 fi,fj,exp_fi,exp_fj,norm_lcl,norm;
    type0 p,p2,p3,p4,tmp0,tmp1,tmp2,s_ij,w_ij,w_ji,c_iv,c_jv;
    type0* coef;
    int m;
    int iatm,jatm;
    int icomp,jcomp,istart;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    int natms=atoms->natms;
    
    
    for(int i=0;i<natms*no_types;i++)
        mu[i]=s[i]=0.0;
    
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        mu[icomp+itype]+=c[jcomp+jtype]*rho_phi[istart+type2rho_pair_ji[jtype][itype]];
                        
                        if(jatm<natms)
                        {
                            mu[jcomp+jtype]+=c[icomp+itype]*rho_phi[istart+type2rho_pair_ij[itype][jtype]];
                        }
                    }
                }
            }
            
            istart+=stride;
            
        }
        
        for(int itype=0;itype<no_types; itype++)
        {
            if(c[icomp+itype]>=0.0)
            {
                p=mu[iatm*no_types+itype]*drho_inv;
                m=static_cast<int>(p);
                m=MIN(m,nrho-2);
                p-=m;
                p=MIN(p,1.0);
                p2=p*p;
                p3=p2*p;
                p4=p3*p;
                coef=F_arr[itype][m];
                
                
                tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                tmp1=(4.0*coef[4]*p3+3.0*coef[3]*p2+2.0*coef[2]*p+coef[1])*drho_inv;
                tmp2=(12.0*coef[4]*p2+6.0*coef[3]*p+2.0*coef[2])*drho_inv*drho_inv;
                
                if(mu[iatm*no_types+itype]>rho_max)
                {
                    tmp0+=tmp1*(mu[icomp+itype]-rho_max);
                    tmp2=0.0;
                }
                
                mu[icomp+itype]=tmp0+c_0[itype]-3.0*kbT*log(x[(3+no_types)*iatm+3+itype]);
                dE[icomp+itype]=tmp1;
                ddE[icomp+itype]=tmp2;
            }
        }
    }
    
    atoms->update_ph(dE_n);
    
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        mu[icomp+itype]+=c[jcomp+jtype]*(rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype]
                                +rho_phi[istart+type2phi_pair_ij[itype][jtype]]);
                        
                        if(jatm<natms)
                            mu[jcomp+jtype]+=c[icomp+itype]*(rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                                +rho_phi[istart+type2phi_pair_ji[jtype][itype]]);
                        
                    }
                }
            }
            istart+=stride;
        }
    }
    
    atoms->update_ph(mu_n);
    if(OLD_NEIGH)
    {
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int j=0;j<neigh_lst_sz[iatm];j++)
            {
                jatm=neigh_lst[iatm][j];
                jcomp=no_types*jatm;
                
                c_jv=1.0;
                for(int itype=0;itype<no_types;itype++)
                    if(c[jcomp+itype]>=0.0)
                        c_jv-=c[jcomp+itype];
                
                for(int itype=0;itype<no_types;itype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+itype]>=0.0)
                    {
                        fi=mu[icomp+itype];
                        fj=mu[jcomp+itype];
                        
                        
                        exp_fi=exp(beta*(fi-mat(fi,0,fj,0,itype)));
                        exp_fj=exp(beta*(fj-mat(fi,0,fj,0,itype)));
                        w_ij=-c[icomp+itype]*c_jv*exp_fi;
                        w_ji=-c[jcomp+itype]*c_iv*exp_fj;
                        
                        s_ij=w_ij-w_ji;
                        
                        
                        s[icomp+itype]+=s_ij;
                        if(jatm<natms)
                        {
                            s[jcomp+itype]-=s_ij;
                        }
                    }
                }
            }
        }
    }
    else
    {
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int j=0;j<neigh_lst_sz[icomp];j++)
                {
                    c_jv=1.0;
                    jcomp=neigh_lst[icomp][j]-itype;
                    for(int jtype=0;jtype<no_types;jtype++)
                    {
                        if(c[jcomp]>=0.0)
                            c_jv-=c[jcomp];
                        jcomp++;
                    }
                    jcomp=neigh_lst[icomp][j];
                    
                    fi=mu[icomp];
                    fj=mu[jcomp];
                    
                    tmp0=mat(fi,0.0,fj,0.0,itype);
                    
                    exp_fi=exp(beta*(fi-tmp0));
                    exp_fj=exp(beta*(fj-tmp0));
                    w_ij=-c[icomp]*c_jv*exp_fi;
                    w_ji=-c[jcomp]*c_iv*exp_fj;
                    
                    s_ij=w_ij-w_ji;
                    
                    s[icomp]+=s_ij;
                    if(jcomp<natms*no_types)
                    {
                        s[jcomp]-=s_ij;
                    }
                }
                icomp++;
            }
        }
    }

    atoms->update_ph(s_n);
    
    for(int i=0;i<natms*no_types;i++)
        t[i]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        
                        t[icomp+itype]+=rho_phi[istart+type2rho_pair_ji[jtype][itype]]*s[jcomp+jtype];
                        
                        if(jatm<natms)
                        {
                            t[jcomp+jtype]+=rho_phi[istart+type2rho_pair_ij[itype][jtype]]*s[icomp+itype];
                        }
                    }
                }
            }
            istart+=stride;
        }
        
        for(int itype=0;itype<no_types;itype++)
            if(c[icomp+itype]>=0.0)
                t[icomp+itype]*=c[icomp+itype]*ddE[icomp+itype];
    }
    
    atoms->update_ph(t_n);
    
    
    for(int i=0;i<natms*no_types;i++)
        ddE[i]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+jtype]>=0.0)
                    {
                        
                        tmp0=rho_phi[istart+type2phi_pair_ij[itype][jtype]]
                        +rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                        +rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype];
                        
                        ddE[icomp+itype]+=tmp0*s[jcomp+jtype]+t[jcomp+jtype]*rho_phi[istart+type2rho_pair_ij[itype][jtype]];
                        
                        if(jatm<natms)
                        {
                            ddE[jcomp+jtype]+=tmp0*s[icomp+itype]+t[icomp+itype]*rho_phi[istart+type2rho_pair_ji[jtype][itype]];
                            
                        }
                    }
                }
            }
            istart+=stride;
        }
        
    }
    atoms->update_ph(E_n);
    
    
    for(int i=0;i<natms*no_types;i++)
        dE[i]=0.0;
    
    if(OLD_NEIGH)
    {
        norm_lcl=0.0;
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int j=0;j<neigh_lst_sz[iatm];j++)
            {
                jatm=neigh_lst[iatm][j];
                jcomp=no_types*jatm;
                
                c_jv=1.0;
                for(int itype=0;itype<no_types;itype++)
                    if(c[jcomp+itype]>=0.0)
                        c_jv-=c[jcomp+itype];
                
                for(int itype=0;itype<no_types;itype++)
                {
                    if(c[icomp+itype]>=0.0 && c[jcomp+itype]>=0.0)
                    {
                        
                        fi=mu[icomp+itype];
                        fj=mu[jcomp+itype];
                        
                        
                        exp_fi=exp(beta*(fi-mat(fi,0,fj,0,itype)));
                        exp_fj=exp(beta*(fj-mat(fi,0,fj,0,itype)));
                        w_ij=-c[icomp+itype]*c_jv*exp_fi;
                        w_ji=-c[jcomp+itype]*c_iv*exp_fj;
                        
                        s_ij=w_ij-w_ji;
                        
                        
                        tmp0=(s[jcomp+itype]*c[icomp+itype]-s[icomp+itype]*c_jv)*exp_fi
                        -(s[icomp+itype]*c[jcomp+itype]-s[jcomp+itype]*c_iv)*exp_fj;
                        
                        tmp0+=beta*(ddE[icomp+itype]*w_ij-ddE[jcomp+itype]*w_ji);
                        
                        tmp0+=-beta*(w_ij-w_ji)*(
                                                 +ddE[icomp+itype]*dmat0(fi,0,fj,0,itype)
                                                 +ddE[jcomp+itype]*dmat0(fj,0,fi,0,itype));
                        
                        
                        dE[icomp+itype]+=tmp0;
                        if(jatm<natms)
                        {
                            dE[jcomp+itype]-=tmp0;
                        }
                    }
                }
            }
            
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    norm_lcl+=dE[icomp+itype]*dE[icomp+itype];
        }
        MPI_Allreduce(&norm_lcl,&tmp0,1,MPI_TYPE0,MPI_SUM,world);
        norm=sqrt(tmp0);
    }
    else
    {
        norm_lcl=0.0;
        for(iatm=0;iatm<natms;iatm++)
        {
            icomp=no_types*iatm;
            c_iv=1.0;
            for(int itype=0;itype<no_types;itype++)
                if(c[icomp+itype]>=0.0)
                    c_iv-=c[icomp+itype];
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int j=0;j<neigh_lst_sz[icomp];j++)
                {
                    c_jv=1.0;
                    jcomp=neigh_lst[icomp][j]-itype;
                    for(int jtype=0;jtype<no_types;jtype++)
                    {
                        if(c[jcomp]>=0.0)
                            c_jv-=c[jcomp];
                        jcomp++;
                    }
                    jcomp=neigh_lst[icomp][j];
                    
                    fi=mu[icomp];
                    fj=mu[jcomp];
                    
                    tmp0=mat(fi,0.0,fj,0.0,itype);
                    
                    exp_fi=exp(beta*(fi-tmp0));
                    exp_fj=exp(beta*(fj-tmp0));
                    w_ij=-c[icomp]*c_jv*exp_fi;
                    w_ji=-c[jcomp]*c_iv*exp_fj;
                    
                    s_ij=w_ij-w_ji;
                    
                    tmp0=(s[jcomp]*c[icomp]-s[icomp]*c_jv)*exp_fi
                    -(s[icomp]*c[jcomp]-s[jcomp]*c_iv)*exp_fj;
                    
                    tmp0+=beta*(ddE[icomp]*w_ij-ddE[jcomp]*w_ji);
                    
                    tmp0+=-beta*(w_ij-w_ji)*(
                                             ddE[icomp]*dmat0(fi,0.0,fj,0.0,itype)
                                             +ddE[jcomp]*dmat0(fj,0.0,fi,0.0,itype));
                    
                    
                    dE[icomp]+=tmp0;
                    if(jcomp<natms*no_types)
                    {
                        dE[jcomp]-=tmp0;
                    }
                }
                
                norm_lcl+=dE[icomp]*dE[icomp];
                icomp++;
            }
        }
        MPI_Allreduce(&norm_lcl,&tmp0,1,MPI_TYPE0,MPI_SUM,world);
        norm=sqrt(tmp0);
    }
    

    return norm;
    
}




