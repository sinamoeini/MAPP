#include "ff_eam_dmd_hg_1.h"
#include "neighbor.h"
#include "atom_types.h"
using namespace MAPP_NS;
#define PI_IN_SQ 0.564189583547756286948079451561
enum{NOT_SET,FUNC_FL,SET_FL,FINNIS_FL};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam_dmd_hg_1::
ForceField_eam_dmd_hg_1(MAPP* mapp) : ForceField(mapp)
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
    if(atoms->vectors[0].dim!=3+no_types)
        error->abort("the dimension of x"
                     " vector should be 3 + no of types");
    
    
    
    CREATE1D(c_0,no_types);
    CREATE1D(delta_e,no_types);
    
    c_d_n=atoms->find_exist("c_d");
    
    if(c_d_n<0)
    {
        int tmp0=atoms->find("c");
        tmp0=atoms->vectors[tmp0].dim;
        c_d_n=atoms->add<type0>(0,tmp0,"c_d");
    }
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_eam_dmd_hg_1::~ForceField_eam_dmd_hg_1()
{
    if(no_types)
    {
        delete [] c_0;
        delete [] delta_e;
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
void ForceField_eam_dmd_hg_1::
force_calc(int st_clc,type0* en_st)
{
    type0* x;
    atoms->vectors[x_n].ret(x);
    type0* f;
    atoms->vectors[f_n].ret(f);
    type0* E;
    atoms->vectors[E_n].ret(E);
    type0* dE;
    atoms->vectors[dE_n].ret(dE);
    type0* c;
    atoms->vectors[c_n].ret(c);
    
    int iatm,jatm;
    
    int icomp,jcomp,istart;
    type0 dx0,dx1,dx2,rsq,z2p,z2;
    type0 alpha;
    type0 r=0.0,r_inv=0.0;
    int m;
    type0* coef;
    type0 inv_alpha,tmp0,tmp1,rtmp;
    type0 upper,lower;
    type0 fpair,apair;
    
    type0 p,p2,p3,p4;
    type0 phi_it_jt_0,phi_it_jt_1,phi_it_jt_2;
    type0 rho_it_jt_0,rho_it_jt_1,rho_it_jt_2;
    type0 rho_jt_it_0,rho_jt_it_1,rho_jt_it_2;
    
    
    
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
        for(int i=0;i<no_0;i++)
            rho_phi[i]=drho_phi_dr[i]=drho_phi_dalpha[i]=0.0;
    }
    
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
                        alpha=sqrt(x[icomp+3+itype]*x[icomp+3+itype]+x[jcomp+3+jtype]*x[jcomp+3+jtype]);
                        
                        if(alpha_min<alpha && alpha<alpha_max)
                        {
                            inv_alpha=1.0/alpha;
                            upper=(r+rc)*inv_alpha;
                            lower=(r-rc)*inv_alpha;
                            
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
                                        tmp1=coef[6]*p2+coef[5]*p+coef[4];
                                        if(rtmp<0.0)
                                            tmp0*=-1.0;
                                        
                                        rho_it_jt_0+=wi_0[i]*tmp0;
                                        rho_it_jt_1+=wi_0[i]*tmp1;
                                        rho_it_jt_2+=wi_1[i]*tmp1;
                                        
                                        if(itype!=jtype)
                                        {
                                            coef=rho_r_arr[type2rho[jtype][itype]][m];
                                            tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                            tmp1=coef[6]*p2+coef[5]*p+coef[4];
                                            if(rtmp<0.0)
                                                tmp0*=-1.0;
                                        }
                                        
                                        rho_jt_it_0+=wi_0[i]*tmp0;
                                        rho_jt_it_1+=wi_0[i]*tmp1;
                                        rho_jt_it_2+=wi_1[i]*tmp1;
                                        
                                        
                                        coef=phi_r_arr[type2phi[itype][jtype]][m];
                                        tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                        tmp1=coef[6]*p2+coef[5]*p+coef[4];
                                        if(rtmp<0.0)
                                            tmp0*=-1.0;
                                        
                                        phi_it_jt_0+=wi_0[i]*tmp0;
                                        phi_it_jt_1+=wi_0[i]*tmp1;
                                        phi_it_jt_2+=wi_1[i]*tmp1;
                                    }
                                }
                                
                                
                                rho_it_jt_0*=PI_IN_SQ*r_inv;
                                rho_it_jt_1*=PI_IN_SQ*r_inv;
                                rho_it_jt_1-=rho_it_jt_0*r_inv;
                                rho_it_jt_2*=-PI_IN_SQ*r_inv;
                                
                                rho_jt_it_0*=PI_IN_SQ*r_inv;
                                rho_jt_it_1*=PI_IN_SQ*r_inv;
                                rho_jt_it_1-=rho_jt_it_0*r_inv;
                                rho_jt_it_2*=-PI_IN_SQ*r_inv;
                                
                                phi_it_jt_0*=PI_IN_SQ*r_inv;
                                phi_it_jt_1*=PI_IN_SQ*r_inv;
                                phi_it_jt_1-=phi_it_jt_0*r_inv;
                                phi_it_jt_2*=-PI_IN_SQ*r_inv;
                                
                                
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
                        else if (alpha<=alpha_min)
                        {
                            if(rsq<cut_sq_0)
                            {
                                r=sqrt(rsq);
                                p=r*dr_inv;
                                m=static_cast<int>(p);
                                m=MIN(m,nr-2);
                                p-=m;
                                p=MIN(p,1.0);
                                p2=p*p;
                                p3=p2*p;
                                
                                coef=phi_r_arr[type2phi[itype][jtype]][m];
                                z2p=coef[6]*p2+coef[5]*p+coef[4];
                                z2=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                phi_it_jt_0=z2*r_inv;
                                phi_it_jt_1=(z2p-z2*r_inv)*r_inv;
                                
                                coef=rho_r_arr[type2rho[itype][jtype]][m];
                                z2p=coef[6]*p2+coef[5]*p+coef[4];
                                z2=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                rho_it_jt_0=z2*r_inv;
                                rho_it_jt_1=(z2p-z2*r_inv)*r_inv;
                                
                                coef=rho_r_arr[type2rho[jtype][itype]][m];
                                z2p=coef[6]*p2+coef[5]*p+coef[4];
                                z2=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                rho_jt_it_0=z2*r_inv;
                                rho_jt_it_1=(z2p-z2*r_inv)*r_inv;
                                
                                
                                
                                rho_phi[istart+type2rho_pair_ij[itype][jtype]]=rho_it_jt_0;
                                rho_phi[istart+type2rho_pair_ji[jtype][itype]]=rho_jt_it_0;
                                rho_phi[istart+type2phi_pair_ij[itype][jtype]]=phi_it_jt_0;
                                
                                drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]=rho_it_jt_1;
                                drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]=rho_jt_it_1;
                                drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]]=phi_it_jt_1;
                                
                                drho_phi_dalpha[istart+type2rho_pair_ij[itype][jtype]]=0.0;
                                drho_phi_dalpha[istart+type2rho_pair_ji[jtype][itype]]=0.0;
                                drho_phi_dalpha[istart+type2phi_pair_ij[itype][jtype]]=0.0;
                                
                                
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
            
            
            
            istart+=stride;
        }
        
        
        for(int itype=0;itype<no_types; itype++)
        {
            p=E[iatm*no_types+itype]*drho_inv;
            m=static_cast<int>(p);
            m=MIN(m,nr-2);
            p-=m;
            p=MIN(p,1.0);
            p2=p*p;
            p3=p2*p;
            p4=p3*p;
            coef=F_arr[itype][m];
            
            tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
            tmp1=coef[8]*p3+coef[7]*p2+coef[6]*p+coef[5];
            
            
            if(E[iatm*no_types+itype]>rho_max)
                tmp0+=tmp1*(E[iatm*no_types+itype]-rho_max);
            
            E[iatm*no_types+itype]=tmp0;
            dE[iatm*no_types+itype]=tmp1;
            
            nrgy_strss[0]+=c[iatm*no_types+itype]*E[iatm*no_types+itype];
            
            nrgy_strss[0]+=kbT*calc_ent(c[iatm*no_types+itype]);
            
            nrgy_strss[0]-=3.0*kbT*c[iatm*no_types+itype]*log(x[icomp+3+itype])
            -c[iatm*no_types+itype]*c_0[itype];
            
            f[icomp+3+itype]+=3.0*kbT*c[iatm*no_types+itype]/x[icomp+3+itype];
        }
        
    }
    
    atoms->update(dE_n);
    
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
                    
                    
                    fpair=-(drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm*no_types+itype]
                            +drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm*no_types+jtype]
                            +drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]])*c[iatm*no_types+itype]*c[jatm*no_types+jtype]*r_inv;
                    
                    apair=-(
                            drho_phi_dalpha[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm*no_types+itype]
                            +drho_phi_dalpha[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm*no_types+jtype]
                            +drho_phi_dalpha[istart+type2phi_pair_ij[itype][jtype]])*c[iatm*no_types+itype]*c[jatm*no_types+jtype]/
                    sqrt(x[icomp+3+itype]*x[icomp+3+itype]+x[jcomp+3+jtype]*x[jcomp+3+jtype]);
                    
                    
                    
                    
                    if(apair!=0.0 || fpair!=0.0)
                    {
                        f[icomp]+=dx0*fpair;
                        f[icomp+1]+=dx1*fpair;
                        f[icomp+2]+=dx2*fpair;
                        f[icomp+3+itype]+=apair*x[icomp+3+itype];
                        
                        if (jatm<natms)
                        {
                            f[jcomp]-=dx0*fpair;
                            f[jcomp+1]-=dx1*fpair;
                            f[jcomp+2]-=dx2*fpair;
                            f[jcomp+3+jtype]+=apair*x[jcomp+3+jtype];
                            
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
type0 ForceField_eam_dmd_hg_1::energy_calc()
{
    type0 en=0.0,en_tot;
    
    type0* x;
    atoms->vectors[x_n].ret(x);
    type0* E;
    atoms->vectors[E_n].ret(E);
    type0* c;
    atoms->vectors[c_n].ret(c);
    
    int iatm,jatm;
    
    int icomp,jcomp;
    type0 dx0,dx1,dx2,rsq,z2;
    type0 alpha;
    type0 r=0.0,r_inv=0.0;
    int m;
    type0* coef;
    type0 inv_alpha,tmp0,tmp1,rtmp;
    type0 upper,lower;
    
    type0 p,p2,p3,p4;
    type0 phi_it_jt_0;
    type0 rho_it_jt_0;
    type0 rho_jt_it_0;
    
    
    
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
                        alpha=sqrt(x[icomp+3+itype]*x[icomp+3+itype]+x[jcomp+3+jtype]*x[jcomp+3+jtype]);
                        
                        if(alpha_min<alpha && alpha<alpha_max)
                        {
                            inv_alpha=1.0/alpha;
                            upper=(r+rc)*inv_alpha;
                            lower=(r-rc)*inv_alpha;
                            
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
                                
                                
                                rho_it_jt_0*=PI_IN_SQ*r_inv;
                                rho_jt_it_0*=PI_IN_SQ*r_inv;
                                phi_it_jt_0*=PI_IN_SQ*r_inv;
                                
                                
                                
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
                        else if (alpha<=alpha_min)
                        {
                            
                            if(rsq<cut_sq_0)
                            {
                                r=sqrt(rsq);
                                p=r*dr_inv;
                                m=static_cast<int>(p);
                                m=MIN(m,nr-2);
                                p-=m;
                                p=MIN(p,1.0);
                                p2=p*p;
                                p3=p2*p;
                                
                                coef=phi_r_arr[type2phi[itype][jtype]][m];
                                z2=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                phi_it_jt_0=z2*r_inv;
                                
                                
                                coef=rho_r_arr[type2rho[itype][jtype]][m];
                                z2=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                rho_it_jt_0=z2*r_inv;
                                
                                
                                coef=rho_r_arr[type2rho[jtype][itype]][m];
                                z2=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                rho_jt_it_0=z2*r_inv;
                                
                                
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
        
        
        for(int itype=0;itype<no_types; itype++)
        {
            p=E[iatm*no_types+itype]*drho_inv;
            m=static_cast<int>(p);
            m=MIN(m,nr-2);
            p-=m;
            p=MIN(p,1.0);
            p2=p*p;
            p3=p2*p;
            p4=p3*p;
            coef=F_arr[itype][m];
            
            tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
            tmp1=coef[8]*p3+coef[7]*p2+coef[6]*p+coef[5];
            
            
            if(E[iatm*no_types+itype]>rho_max)
                tmp0+=tmp1*(E[iatm*no_types+itype]-rho_max);
            
            E[iatm*no_types+itype]=tmp0;
            
            en+=c[iatm*no_types+itype]*E[iatm*no_types+itype];
            
            en+=kbT*calc_ent(c[iatm*no_types+itype]);
            
            en-=3.0*kbT*c[iatm*no_types+itype]*log(x[icomp+3+itype])
            -c[iatm*no_types+itype]*c_0[itype];
            
            
        }
        
    }
    
    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    
    return en_tot;
}
/*--------------------------------------------
 init before running
 --------------------------------------------*/
void ForceField_eam_dmd_hg_1::init()
{
    type0 skin=atoms->skin;
    type0 ph_cut=0.0;
    for (int i=0;i<no_types*(no_types+1)/2;i++)
        cut_sk_sq[i]=cut_sq_mod_0+(skin)*(skin)
        +2.0*sqrt(cut_sq_mod_0)*(skin);
    
    ph_cut=sqrt(cut_sq_mod_0);
    
    atoms->set_ph(ph_cut);
    
    x_n=atoms->find("x");
    f_n=atoms->find("f");
    c_n=atoms->find("c");
    
    neighbor->pair_wise=1;
    
    E_n=atoms->add<type0>(1,no_types,"E");
    dE_n=atoms->add<type0>(1,no_types,"dE");
    ddE_n=atoms->add<type0>(1,no_types,"ddE");
    n_n=atoms->add<type0>(1,no_types,"n");
    s_n=atoms->add<type0>(1,no_types,"s");
    t_n=atoms->add<type0>(1,2*no_types,"t");
    v_n=atoms->add<type0>(1,2*no_types,"v");
    crd_n=atoms->add<type0>(1,no_types,"crd");
    
}
/*--------------------------------------------
 fin after running
 --------------------------------------------*/
void ForceField_eam_dmd_hg_1::fin()
{
    atoms->del(crd_n);
    atoms->del(v_n);
    atoms->del(t_n);
    atoms->del(s_n);
    atoms->del(n_n);
    atoms->del(ddE_n);
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
void ForceField_eam_dmd_hg_1::coef(int narg,char** arg)
{
    type0 kb,T,hbar;
    if (narg<9)
        error->abort("ff_coef for ff eam_dmd "
                     "should at least have 8 arguments");
    
    
    rsq_crd=atof(arg[1]);
    if(rsq_crd<=0.0)
        error->abort("coordination cutoff for ff eam_dmd "
                     "should be greater than 0.0");
    
    rsq_crd*=rsq_crd;
    
    set_weight_abs(atoi(arg[2]));
    
    alpha_min=atof(arg[3]);
    if(alpha_min<=0.0)
        error->abort("minimum alpha for ff eam_dmd "
                     "should be greater than 0.0");
    
    alpha_max=atof(arg[4]);
    if(alpha_max<=0.0)
        error->abort("maximum alpha for ff eam_dmd "
                     "should be greater than 0.0");
    
    if(alpha_max<=alpha_min)
        error->abort("maximum alpha should be less "
                     "than minimum alpha for ff eam_dmd");
    
    
    T=atof(arg[5]);
    if(T<=0.0)
        error->abort("temperature for ff eam_dmd "
                     "should be greater than 0.0");
    
    
    kb=8.617332478e-5;
    hbar=6.5821192815e-16;
    
    
    setup_delta_e(arg[6]);
    
    clean_up();
    
    if(strcmp(arg[7],"FS")==0)
    {
        eam_mode=FINNIS_FL;
        set_fs(narg-8,&arg[8]);
    }
    else if(strcmp(arg[7],"SetFL")==0)
    {
        eam_mode=SET_FL;
        set_setfl(narg-8,&arg[8]);
    }
    else if(strcmp(arg[7],"FuncFL")==0)
    {
        eam_mode=FUNC_FL;
        set_funcfl(narg-8,&arg[8]);
    }
    else
        error->abort("unknown file format "
        "for ff eam_dmd: %s",arg[7]);
    
    set_arrays();
    
    type0 mass;
    type0 deb_l;
    
    for(int i=0;i<no_types;i++)
    {
        mass=atom_types->mass[i];
        mass*=1.0364269184093291236e-28;
        deb_l=hbar*hbar*2.0/(mass*kb*T);
        c_0[i]=1.5*kb*T*(log(deb_l)-1.0);
    }
    
    kbT=kb*T;
    beta=1.0/kbT;
    
    rc=(static_cast<type0>(nr)-1.0)*dr;
    rho_max=(static_cast<type0>(nrho)-1.0)*drho;
    cut_sq_0=rc*rc;
    mod_rc=rc+xi[no_i-1]*alpha_max;
    cut_sq_mod_0=mod_rc*mod_rc;
    
}
/*--------------------------------------------
 read from file and assign DeltaE
 --------------------------------------------*/
void ForceField_eam_dmd_hg_1::setup_delta_e(char* file_name)
{
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    char** arg=NULL;
    int narg;
    int tot_no_types;
    int* delta_e_chk;
    CREATE1D(delta_e_chk,no_types);
    
    for(int i=0;i<no_types;i++)
        delta_e_chk[i]=0;
    
    if(atoms->my_p_no==0)
    {
        fp=fopen(file_name,"r");
        if(fp==NULL)
            error->abort("ff eam_dmd file %s not found",file_name);
    }
    
    narg=0;
    while(narg==0)
    {
        if(line_read(fp,line)==-1)
            error->abort("%s file ended immaturely",file_name);
        
        narg=mapp->parse_line(line,arg);
    }
    
    tot_no_types=narg;
    
    if(tot_no_types<no_types)
        error->abort("the number of atoms in %s file"
                     " is less than the number of atom types"
                     " present in the system",file_name);
    
    
    int* type_ref;
    CREATE1D(type_ref,tot_no_types);
    for(int ityp=0;ityp<tot_no_types;ityp++)
        type_ref[ityp]=atom_types->find_type_exist(arg[ityp]);
    
    for(int i=0;i<narg;i++)
        delete [] arg[i];
    if(narg)
        delete [] arg;
    
    
    for(int itype=0;itype<tot_no_types;itype++)
    {
        narg=0;
        while(narg==0)
        {
            if(line_read(fp,line)==-1)
                error->abort("%s file ended immaturely",file_name);
            narg=mapp->parse_line(line,arg);
        }
        
        if(narg!=1)
            error->abort("invalid line in %s file: %s",file_name,line);
        if(type_ref[itype]!=-1)
        {
            delta_e[type_ref[itype]]=atof(arg[0]);
            delta_e_chk[type_ref[itype]]=1;
        }
        
        for(int i=0;i<narg;i++)
            delete [] arg[i];
        if(narg)
            delete [] arg;
    }
    
    
    for(int i=0;i<no_types;i++)
        if(delta_e_chk[i]==0)
            error->abort("delta_e(%s) was not set by %s file "
                         ,file_name,atom_types->atom_names[i]);
    
    delete [] line;
    if(tot_no_types)
        delete [] type_ref;
    if(no_types)
        delete [] delta_e_chk;
    
    
}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam_dmd_hg_1::set_funcfl(int no_files
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
            if(line_read(fp,line)==-1)
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
void ForceField_eam_dmd_hg_1::set_setfl(int no_files
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
void ForceField_eam_dmd_hg_1::set_fs(int no_files
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
    
    
    if(line_read(fp,line)==-1)
        error->abort("%s file ended immaturely",file_names[0]);
    
    
    int ipos=0;
    int tot=tot_no_types*nrho+tot_no_types*tot_no_types*nr
    +tot_no_types*(tot_no_types+1)*nr/2;
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
int ForceField_eam_dmd_hg_1::line_read(FILE* file
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
void ForceField_eam_dmd_hg_1::clean_up()
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
void ForceField_eam_dmd_hg_1::allocate()
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
                CREATE1D(phi_r_arr[ityp][i],7);
        
        CREATE1D(rho_r_arr,no_types*no_types);
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(rho_r_arr[ityp],nr);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(rho_r_arr[type2rho[ityp][0]][i],7);
        
        
        for(int ityp=0;ityp<no_types;ityp++)
            for(int jtyp=0;jtyp<no_types;jtyp++)
                rho_r_arr[type2rho[ityp][jtyp]]
                =rho_r_arr[type2rho[ityp][0]];
        
        CREATE1D(F_arr,no_types);
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(F_arr[ityp],nrho);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<nrho;i++)
                CREATE1D(F_arr[ityp][i],9);
        
    }
    else if(eam_mode==FINNIS_FL)
    {
        
        CREATE1D(phi_r_arr,no_types*(no_types+1)/2);
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
            CREATE1D(phi_r_arr[ityp],nr);
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(phi_r_arr[ityp][i],7);
        
        CREATE1D(rho_r_arr,no_types*no_types);
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            CREATE1D(rho_r_arr[ityp],nr);
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(rho_r_arr[ityp][i],7);
        
        CREATE1D(F_arr,no_types);
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(F_arr[ityp],nrho);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<nrho;i++)
                CREATE1D(F_arr[ityp][i],9);
    }
    
    allocated=1;
}
/*--------------------------------------------
 set the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_hg_1::set_arrays()
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
void ForceField_eam_dmd_hg_1::interpolate(int n
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
type0 ForceField_eam_dmd_hg_1::interpolate(type0* arr
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
void ForceField_eam_dmd_hg_1::interpolate_m(int n
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
    
    
    
    for(int i=0;i<n;i++)
    {
        spline[i][5]=spline[i][1]/delta;
        spline[i][6]=2.0*spline[i][2]/delta;
        spline[i][7]=3.0*spline[i][3]/delta;
        spline[i][8]=4.0*spline[i][4]/delta;
    }
    
}
/*--------------------------------------------
 Gaussian-Hermite quadrature weights and
 abscissas for 1 to 14 points
 --------------------------------------------*/
void ForceField_eam_dmd_hg_1::set_weight_abs(int n)
{
    
    if(n<1 || n>14)
        error->abort("number of gaussian "
                     "points for ff eam_dmd should be "
                     "between 1 & 14");
    
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
    
    if(no_i==1)
    {
        wi_0[0]=1.7724538509055160273;
        
        xi[0]=0.0;
    }
    else if(no_i==2)
    {
        wi_0[0]=0.8862269254527580136;
        wi_0[1]=0.8862269254527580136;
        
        xi[0]=-0.70710678118654752440;
        xi[1]=0.70710678118654752440;
    }
    else if(no_i==3)
    {
        wi_0[0]=0.2954089751509193379;
        wi_0[1]=1.1816359006036773515;
        wi_0[2]=0.2954089751509193379;
        
        xi[0]=-1.2247448713915890491;
        xi[1]=0.0;
        xi[2]=1.2247448713915890491;
    }
    else if(no_i==4)
    {
        wi_0[0]=0.0813128354472451771;
        wi_0[1]=0.8049140900055128365;
        wi_0[2]=0.8049140900055128365;
        wi_0[3]=0.0813128354472451771;
        
        xi[0]=-1.6506801238857845559;
        xi[1]=-0.52464762327529031788;
        xi[2]=0.52464762327529031788;
        xi[3]=1.6506801238857845559;
    }
    else if(no_i==5)
    {
        wi_0[0]=0.01995324205904591321;
        wi_0[1]=0.3936193231522411598;
        wi_0[2]=0.94530872048294188123;
        wi_0[3]=0.3936193231522411598;
        wi_0[4]=0.01995324205904591321;
        
        xi[0]=-2.0201828704560856329;
        xi[1]=-0.95857246461381850711;
        xi[2]=0.0;
        xi[3]=0.95857246461381850711;
        xi[4]=2.0201828704560856329;
    }
    else if(no_i==6)
    {
        wi_0[0]=0.00453000990550884564;
        wi_0[1]=0.1570673203228566439;
        wi_0[2]=0.7246295952243925241;
        wi_0[3]=0.7246295952243925241;
        wi_0[4]=0.1570673203228566439;
        wi_0[5]=0.00453000990550884564;
        
        xi[0]=-2.3506049736744922228;
        xi[1]=-1.3358490740136969497;
        xi[2]=-0.43607741192761650868;
        xi[3]=0.43607741192761650868;
        xi[4]=1.3358490740136969497;
        xi[5]=2.3506049736744922228;
    }
    else if(no_i==7)
    {
        wi_0[0]=0.00097178124509951915;
        wi_0[1]=0.0545155828191270306;
        wi_0[2]=0.4256072526101278005;
        wi_0[3]=0.81026461755680732676;
        wi_0[4]=0.4256072526101278005;
        wi_0[5]=0.0545155828191270306;
        wi_0[6]=0.00097178124509951915;
        
        xi[0]=-2.6519613568352334924;
        xi[1]=-1.6735516287674714450;
        xi[2]=-0.81628788285896466304;
        xi[3]=0.0;
        xi[4]=0.81628788285896466304;
        xi[5]=1.6735516287674714450;
        xi[6]=2.6519613568352334924;
    }
    else if(no_i==8)
    {
        wi_0[0]=0.000199604072211367619;
        wi_0[1]=0.0170779830074134755;
        wi_0[2]=0.207802325814891880;
        wi_0[3]=0.6611470125582412910;
        wi_0[4]=0.6611470125582412910;
        wi_0[5]=0.207802325814891880;
        wi_0[6]=0.0170779830074134755;
        wi_0[7]=0.000199604072211367619;
        
        xi[0]=-2.9306374202572440192;
        xi[1]=-1.9816567566958429259;
        xi[2]=-1.1571937124467801947;
        xi[3]=-0.38118699020732211685;
        xi[4]=0.38118699020732211685;
        xi[5]=1.1571937124467801947;
        xi[6]=1.9816567566958429259;
        xi[7]=2.9306374202572440192;
    }
    else if(no_i==9)
    {
        wi_0[0]=0.000039606977263264382;
        wi_0[1]=0.0049436242755369472;
        wi_0[2]=0.088474527394376573;
        wi_0[3]=0.432651559002555750;
        wi_0[4]=0.72023521560605095712;
        wi_0[5]=0.432651559002555750;
        wi_0[6]=0.088474527394376573;
        wi_0[7]=0.0049436242755369472;
        wi_0[8]=0.000039606977263264382;
        
        xi[0]=-3.1909932017815276072;
        xi[1]=-2.2665805845318431118;
        xi[2]=-1.4685532892166679317;
        xi[3]=-0.72355101875283757332;
        xi[4]=0.0;
        xi[5]=0.72355101875283757332;
        xi[6]=1.4685532892166679317;
        xi[7]=2.2665805845318431118;
        xi[8]=3.1909932017815276072;
    }
    else if(no_i==10)
    {
        wi_0[0]=7.6404328552326206e-6;
        wi_0[1]=0.0013436457467812327;
        wi_0[2]=0.033874394455481063;
        wi_0[3]=0.240138611082314686;
        wi_0[4]=0.6108626337353257988;
        wi_0[5]=0.6108626337353257988;
        wi_0[6]=0.240138611082314686;
        wi_0[7]=0.033874394455481063;
        wi_0[8]=0.0013436457467812327;
        wi_0[9]=7.6404328552326206e-6;
        
        xi[0]=-3.4361591188377376033;
        xi[1]=-2.5327316742327897964;
        xi[2]=-1.7566836492998817735;
        xi[3]=-1.0366108297895136542;
        xi[4]=-0.34290132722370460879;
        xi[5]=0.34290132722370460879;
        xi[6]=1.0366108297895136542;
        xi[7]=1.7566836492998817735;
        xi[8]=2.5327316742327897964;
        xi[9]=3.4361591188377376033;
    }
    else if(no_i==11)
    {
        wi_0[0]=1.4395603937142582e-6;
        wi_0[1]=0.00034681946632334551;
        wi_0[2]=0.011911395444911532;
        wi_0[3]=0.117227875167708503;
        wi_0[4]=0.429359752356125028;
        wi_0[5]=0.65475928691459177920;
        wi_0[6]=0.429359752356125028;
        wi_0[7]=0.117227875167708503;
        wi_0[8]=0.011911395444911532;
        wi_0[9]=0.00034681946632334551;
        wi_0[10]=1.4395603937142582e-6;
        
        xi[0]=-3.6684708465595825185;
        xi[1]=-2.7832900997816517708;
        xi[2]=-2.0259480158257553352;
        xi[3]=-1.3265570844949328559;
        xi[4]=-0.65680956688209976502;
        xi[5]=0.0;
        xi[6]=0.65680956688209976502;
        xi[7]=1.3265570844949328559;
        xi[8]=2.0259480158257553352;
        xi[9]=2.7832900997816517708;
        xi[10]=3.6684708465595825185;
    }
    else if(no_i==12)
    {
        wi_0[0]=2.6585516843563016e-7;
        wi_0[1]=0.00008573687043587859;
        wi_0[2]=0.003905390584629062;
        wi_0[3]=0.051607985615883930;
        wi_0[4]=0.260492310264161129;
        wi_0[5]=0.5701352362624795783;
        wi_0[6]=0.5701352362624795783;
        wi_0[7]=0.260492310264161129;
        wi_0[8]=0.051607985615883930;
        wi_0[9]=0.003905390584629062;
        wi_0[10]=0.00008573687043587859;
        wi_0[11]=2.6585516843563016e-7;
        
        xi[0]=-3.8897248978697819193;
        xi[1]=-3.0206370251208897717;
        xi[2]=-2.2795070805010599002;
        xi[3]=-1.5976826351526047967;
        xi[4]=-0.94778839124016374370;
        xi[5]=-0.31424037625435911128;
        xi[6]=0.31424037625435911128;
        xi[7]=0.94778839124016374370;
        xi[8]=1.5976826351526047967;
        xi[9]=2.2795070805010599002;
        xi[10]=3.0206370251208897717;
        xi[11]=3.8897248978697819193;
    }
    else if(no_i==13)
    {
        wi_0[0]=4.825731850073131e-8;
        wi_0[1]=0.00002043036040270707;
        wi_0[2]=0.0012074599927193859;
        wi_0[3]=0.020862775296169939;
        wi_0[4]=0.140323320687023438;
        wi_0[5]=0.421616296898543222;
        wi_0[6]=0.60439318792116164234;
        wi_0[7]=0.421616296898543222;
        wi_0[8]=0.140323320687023438;
        wi_0[9]=0.020862775296169939;
        wi_0[10]=0.0012074599927193859;
        wi_0[11]=0.00002043036040270707;
        wi_0[12]=4.825731850073131e-8;
        
        xi[0]=-4.1013375961786396412;
        xi[1]=-3.2466089783724099881;
        xi[2]=-2.5197356856782378834;
        xi[3]=-1.8531076516015121420;
        xi[4]=-1.2200550365907484262;
        xi[5]=-0.60576387917106011308;
        xi[6]=0.0;
        xi[7]=0.60576387917106011308;
        xi[8]=1.2200550365907484262;
        xi[9]=1.8531076516015121420;
        xi[10]=2.5197356856782378834;
        xi[11]=3.2466089783724099881;
        xi[12]=4.1013375961786396412;
    }
    else if(no_i==14)
    {
        wi_0[0]=8.628591168125158e-9;
        wi_0[1]=4.716484355018917e-6;
        wi_0[2]=0.0003550926135519236;
        wi_0[3]=0.007850054726457944;
        wi_0[4]=0.06850553422346521;
        wi_0[5]=0.273105609064246603;
        wi_0[6]=0.5364059097120901498;
        wi_0[7]=0.5364059097120901498;
        wi_0[8]=0.273105609064246603;
        wi_0[9]=0.06850553422346521;
        wi_0[10]=0.007850054726457944;
        wi_0[11]=0.0003550926135519236;
        wi_0[12]=4.716484355018917e-6;
        wi_0[13]=8.628591168125158e-9;
        
        xi[0]=-4.3044485704736318126;
        xi[1]=-3.4626569336022705502;
        xi[2]=-2.7484707249854025686;
        xi[3]=-2.0951832585077168157;
        xi[4]=-1.4766827311411408706;
        xi[5]=-0.87871378732939941611;
        xi[6]=-0.29174551067256207845;
        xi[7]=0.29174551067256207845;
        xi[8]=0.87871378732939941611;
        xi[9]=1.4766827311411408706;
        xi[10]=2.0951832585077168157;
        xi[11]=2.7484707249854025686;
        xi[12]=3.4626569336022705502;
        xi[13]=4.3044485704736318126;
    }
    
    for(int i=0;i<no_i;i++)
    {
        wi_1[i]=wi_0[i]*xi[i];
        wi_2[i]=wi_0[i]*xi[i]*xi[i];
    }
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_hg_1::create_2nd_neigh_lst()
{
    
    
    type0 dx0,dx1,dx2,rsq;
    int iatm,jatm,icomp,jcomp;
    int natms=atoms->natms;
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    type0* x;
    atoms->vectors[x_n].ret(x);
    
    
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
            
            if(rsq<rsq_crd)
            {
                if(neigh_lst_sz[iatm]+1>tmp_neigh_list_size)
                {
                    GROW(tmp_neigh_list,tmp_neigh_list_size,tmp_neigh_list_size+1);
                    tmp_neigh_list_size++;
                }
                tmp_neigh_list[neigh_lst_sz[iatm]]=jatm;
                neigh_lst_sz[iatm]++;
                
            }
            
            CREATE1D(neigh_lst[iatm],neigh_lst_sz[iatm]);
            memcpy(neigh_lst[iatm],tmp_neigh_list,neigh_lst_sz[iatm]*sizeof(int));
            
        }
    }
    
    if(tmp_neigh_list_size)
        delete [] tmp_neigh_list;
    
}
/*--------------------------------------------
 return M_{ij}^{\alpha}
 this fucntion should be symmetric wrt fi &
 fj
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_1::
mat(type0 fi,type0 crdi,type0 fj,type0 crdj
    ,int itype)
{
    if(fi>fj)
        return fi+delta_e[itype];
    else
        return fj+delta_e[itype];
}
/*--------------------------------------------
 dmat(fi,fj)/dfi
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_1::
dmat0(type0 fi,type0 crdi,type0 fj,type0 crdj
      ,int itype)
{
    
    if(fi>fj)
        return 1.0;
    else
        return 0.0;
}

/*--------------------------------------------
 dmat(fi,fj)/dfi
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_1::
dmat1(type0 fi,type0 crdi,type0 fj,type0 crdj
      ,int itype)
{
    
    return 0.0;
}
/*--------------------------------------------
 return M_{ij}^{\alpha}
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_1::calc_ent(type0 x)
{
    if(x<1.0e-8 || 0.99999999<x)
        return 0.0;
    else
        return x*log(x)+(1-x)*log(1-x);
    
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_hg_1::c_d_calc()
{
    
    type0* x;
    atoms->vectors[x_n].ret(x);
    type0* E;
    atoms->vectors[E_n].ret(E);
    type0* dE;
    atoms->vectors[dE_n].ret(dE);
    type0* c;
    atoms->vectors[c_n].ret(c);
    type0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    type0* n;
    atoms->vectors[n_n].ret(n);
    type0* crd;
    atoms->vectors[crd_n].ret(crd);
    
    type0 fi,fj,exp_fi,exp_fj;
    type0 p,p2,p3,p4,tmp0,tmp1,s_ij,crdi,crdj,w_ij,w_ji;
    type0* coef;
    int m;
    int iatm,jatm;
    int icomp,jcomp,istart;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    int natms=atoms->natms;
    for(int i=0;i<natms*no_types;i++)
        c_d[i]=E[i]=0.0;
    
    
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
                crd[iatm+itype]+=c[jcomp+itype];
                if(jatm<natms)
                    crd[jatm+itype]+=c[icomp+itype];
                
                
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    E[icomp+itype]+=c[jcomp+jtype]*rho_phi[istart+type2rho_pair_ji[jtype][itype]];
                    if(jatm<natms)
                        E[jcomp+jtype]+=c[icomp+itype]*rho_phi[istart+type2rho_pair_ij[itype][jtype]];
                }
            }
            
            istart+=stride;
            
        }
        
        for(int itype=0;itype<no_types; itype++)
        {
            p=E[iatm*no_types+itype]*drho_inv;
            m=static_cast<int>(p);
            m=MIN(m,nr-2);
            p-=m;
            p=MIN(p,1.0);
            p2=p*p;
            p3=p2*p;
            p4=p3*p;
            coef=F_arr[itype][m];
            
            
            tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
            tmp1=coef[8]*p3+coef[7]*p2+coef[6]*p+coef[5];
            
            if(E[iatm*no_types+itype]>rho_max)
                tmp0+=tmp1*(E[icomp+itype]-rho_max);
            
            n[icomp+itype]=E[icomp+itype]=tmp0;
            dE[icomp+itype]=tmp1;
        }
    }
    
    atoms->update(dE_n);
    
    
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
                    n[icomp+itype]+=c[jcomp+jtype]*(rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm+jtype]
                                                    +rho_phi[istart+type2phi_pair_ij[itype][jtype]]);
                    
                    if(jatm<natms)
                        n[jcomp+jtype]+=c[icomp+itype]*(rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm+itype]
                                                        +rho_phi[istart+type2phi_pair_ji[jtype][itype]]);
                }
            }
            istart+=stride;
        }
        
        for(int itype=0;itype<no_types;itype++)
        {
            n[icomp+itype]+=c_0[itype];
            //fix this part for change of variables
            n[icomp+itype]-=3.0*kbT*log(x[(3+no_types)*iatm+3+itype]);
        }
        
    }
    
    atoms->update(n_n);
    
    
    atoms->update(crd_n);
    
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        for(int j=0;j<neigh_lst_sz[iatm];j++)
        {
            jatm=neigh_lst[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                crdi=crd[icomp+itype];
                crdj=crd[jcomp+itype];
                fi=n[icomp+itype];
                fj=n[jcomp+itype];
                exp_fi=exp(beta*(fi-mat(fi,crdi,fj,crdj,itype)));
                exp_fj=exp(beta*(fj-mat(fi,crdi,fj,crdj,itype)));
                w_ij=-c[icomp+itype]*(1.0-c[jcomp+itype])*exp_fi;
                w_ji=-c[jcomp+itype]*(1.0-c[icomp+itype])*exp_fj;
                
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
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
type0 ForceField_eam_dmd_hg_1::g_calc(int chk
                                    ,type0 alpha,type0* a,type0* g)
{
    type0* x;
    atoms->vectors[x_n].ret(x);
    type0* E;
    atoms->vectors[E_n].ret(E);
    type0* dE;
    atoms->vectors[dE_n].ret(dE);
    type0* c;
    atoms->vectors[c_n].ret(c);
    type0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    type0* n;
    atoms->vectors[n_n].ret(n);
    type0* s;
    atoms->vectors[s_n].ret(s);
    type0* t;
    atoms->vectors[t_n].ret(t);
    type0* v;
    atoms->vectors[v_n].ret(v);
    type0* crd;
    atoms->vectors[crd_n].ret(crd);
    
    type0 inner,ans;
    type0 fi,fj,exp_fi,exp_fj;
    type0 p,p2,p3,p4,tmp0,tmp1,s_ij;
    type0 w_ij,w_ji,l_ij;
    type0 crdi,crdj;
    type0* coef;
    int m;
    int iatm,jatm;
    int icomp,jcomp,istart;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    int natms=atoms->natms;
    for(int i=0;i<natms*no_types;i++)
        g[i]=E[i]=0.0;
    
    
    
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
                crd[iatm+itype]+=c[jcomp+itype];
                if(jatm<natms)
                    crd[jatm+itype]+=c[icomp+itype];
                
                
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    E[icomp+itype]+=c[jcomp+jtype]*rho_phi[istart+type2rho_pair_ji[jtype][itype]];
                    if(jatm<natms)
                        E[jcomp+jtype]+=c[icomp+itype]*rho_phi[istart+type2rho_pair_ij[itype][jtype]];
                }
            }
            
            istart+=stride;
            
        }
        
        for(int itype=0;itype<no_types; itype++)
        {
            p=E[iatm*no_types+itype]*drho_inv;
            m=static_cast<int>(p);
            m=MIN(m,nr-2);
            p-=m;
            p=MIN(p,1.0);
            p2=p*p;
            p3=p2*p;
            p4=p3*p;
            coef=F_arr[itype][m];
            
            
            tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
            tmp1=coef[8]*p3+coef[7]*p2+coef[6]*p+coef[5];
            
            if(E[iatm*no_types+itype]>rho_max)
                tmp0+=tmp1*(E[icomp+itype]-rho_max);
            
            n[icomp+itype]=E[icomp+itype]=tmp0;
            dE[icomp+itype]=tmp1;
            
            //ddE
            v[2*(icomp+itype)]=(3.0*coef[8]*p2+coef[7]*p+coef[6])*drho_inv;
        }
    }
    
    atoms->update(dE_n);
    
    /*
     claculate f (n)
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
                    n[icomp+itype]+=c[jcomp+jtype]*(rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm+jtype]
                                                    +rho_phi[istart+type2phi_pair_ij[itype][jtype]]);
                    if(jatm<natms)
                        n[jcomp+jtype]+=c[icomp+itype]*(rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm+itype]
                                                        +rho_phi[istart+type2phi_pair_ji[jtype][itype]]);
                }
            }
            istart+=stride;
        }
        
        for(int itype=0;itype<no_types;itype++)
        {
            n[icomp+itype]+=c_0[itype];
            //fix this part for change of variables
            n[icomp+itype]-=3.0*kbT*log(x[(3+no_types)*iatm+3+itype]);
        }
        
    }
    
    atoms->update(n_n);
    
    for(int i=0;i<natms*no_types;i++)
    {
        s[i]=2.0*(c[i]+a[i]);
        c_d[i]=0.0;
    }
    
    /*
     claculate s and c_d
     */
    
    atoms->update(crd_n);
    
    inner=0.0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        for(int j=0;j<neigh_lst_sz[iatm];j++)
        {
            jatm=neigh_lst[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                crdi=crd[icomp+itype];
                crdj=crd[jcomp+itype];
                fi=n[icomp+itype];
                fj=n[jcomp+itype];
                exp_fi=exp(beta*(fi-mat(fi,crdi,fj,crdj,itype)));
                exp_fj=exp(beta*(fj-mat(fi,crdi,fj,crdj,itype)));
                w_ij=-c[icomp+itype]*(1.0-c[jcomp+itype])*exp_fi;
                w_ji=-c[jcomp+itype]*(1.0-c[icomp+itype])*exp_fj;
                
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
        for(int itype=0;itype<no_types;itype++)
        {
            inner+=s[icomp+itype]*s[icomp+itype];
        }
    }
    ans=0.0;
    MPI_Allreduce(&inner,&ans,1,MPI_TYPE0,MPI_SUM,world);
    ans*=0.25;
    
    
    if(chk)
        return ans;
    
    
    atoms->update(s_n);
    
    for(int i=0;i<natms*no_types;i++)
    {
        v[2*i+1]=t[2*i]=t[2*i+1]=0.0;
        g[i]=-s[i];
    }
    
    /*
     claculate t and first part of g
     */
    
    
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        for(int j=0;j<neigh_lst_sz[iatm];j++)
        {
            jatm=neigh_lst[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                
                s_ij=s[icomp+itype]-s[jcomp+itype];
                fi=n[icomp+itype];
                fj=n[jcomp+itype];
                crdi=crd[icomp+itype];
                crdj=crd[jcomp+itype];
                
                exp_fi=exp(beta*(fi-mat(fi,crdi,fj,crdj,itype)));
                exp_fj=exp(beta*(fj-mat(fi,crdi,fj,crdj,itype)));
                
                w_ij=-c[icomp+itype]*(1.0-c[jcomp+itype])*exp_fi;
                w_ji=-c[jcomp+itype]*(1.0-c[icomp+itype])*exp_fj;
                
                l_ij=s_ij*(w_ji-w_ij)*dmat1(fi,crdi,fj,crdj,itype);
                
                t[2*(icomp+itype)]+=s_ij*(w_ij+(w_ji-w_ij)*dmat0(fi,crdi,fj,crdj,itype));
                
                t[2*(icomp+itype)+1]+=l_ij;
                
                g[icomp+itype]-=alpha*((1.0-c[jcomp+itype])*exp_fi+c[jcomp+itype]*exp_fj)*s_ij;
                
                if(jatm<natms)
                {
                    t[2*(jcomp+itype)]-=s_ij*(w_ji+(w_ij-w_ji)*dmat0(fj,crdj,fi,crdi,itype));
                    
                    t[2*(jcomp+itype)+1]+=s_ij*(w_ji-w_ij)*dmat1(fj,crdj,fi,crdi,itype);
                    
                    g[jcomp+itype]+=alpha*((1.0-c[icomp+itype])*exp_fj+c[icomp+itype]*exp_fi)*s_ij;
                    
                }
            }
        }
    }
    
    atoms->update(t_n);
    
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
                g[icomp+itype]+=alpha*beta*t[2*(jcomp+itype)+1];
                if(jatm<natms)
                {
                    g[jcomp+itype]+=alpha*beta*t[2*(icomp+itype)+1];
                }
                
                
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    
                    tmp0=rho_phi[istart+type2phi_pair_ij[itype][jtype]]
                    +rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[icomp+itype]
                    +rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jcomp+jtype];
                    tmp0*=beta*alpha;
                    
                    v[2*(icomp+itype)+1]+=rho_phi[istart+type2rho_pair_ji[jtype][itype]]*t[2*(jcomp+jtype)];
                    g[icomp+itype]+=tmp0*t[2*(jcomp+jtype)];
                    if(jatm<natms)
                    {
                        v[2*(jcomp+jtype)+1]+=rho_phi[istart+type2rho_pair_ij[itype][jtype]]*t[2*(icomp+itype)];
                        g[jcomp+jtype]+=tmp0*t[2*(icomp+itype)];
                    }
                }
            }
            istart+=stride;
        }
    }
    
    atoms->update(v_n);
    
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
                    g[icomp+itype]+=alpha*beta*rho_phi[istart+type2rho_pair_ij[itype][jtype]]
                    *v[2*(jcomp+jtype)]*c[jcomp+jtype]*v[2*(jcomp+jtype)+1];
                    if(jatm<natms)
                        g[jcomp+jtype]+=alpha*beta*rho_phi[istart+type2rho_pair_ji[jtype][itype]]
                        *v[2*(icomp+itype)]*c[icomp+itype]*v[2*(icomp+itype)+1];
                }
            }
            istart+=no_types*no_types;
        }
    }
    
    return ans;
    
}


