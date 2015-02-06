#include "ff_eam_dmd_lg_1.h"
#include "neighbor.h"
#include "atom_types.h"
using namespace MAPP_NS;
#define PI_IN_SQ 0.564189583547756286948079451561
enum{NOT_SET,FUNC_FL,SET_FL,FINNIS_FL};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam_dmd_lg_1::
ForceField_eam_dmd_lg_1(MAPP* mapp) : ForceField(mapp)
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
ForceField_eam_dmd_lg_1::~ForceField_eam_dmd_lg_1()
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
        delete [] wi;
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
void ForceField_eam_dmd_lg_1::
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
    type0 dx0,dx1,dx2,rsq;
    type0 alpha;
    type0 r=0.0,r_inv=0.0;
    int m;
    type0* coef;
    type0** arr_rho_ij;
    type0** arr_rho_ji;
    type0** arr_phi;
    type0 inv_alpha,inv_alpha_sqr,alpha_sqr,tmp0,tmp1;
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
                        
                        if(alpha<alpha_max)
                        {
                            inv_alpha=1.0/alpha;
                            inv_alpha_sqr=inv_alpha*inv_alpha;
                            alpha_sqr=alpha*alpha;
                            
                            phi_it_jt_0=phi_it_jt_1=phi_it_jt_2=0.0;
                            rho_it_jt_0=rho_it_jt_1=rho_it_jt_2=0.0;
                            rho_jt_it_0=rho_jt_it_1=rho_jt_it_2=0.0;
                            
                            
                            arr_rho_ij=rho_r_arr[type2rho[itype][jtype]];
                            arr_rho_ji=rho_r_arr[type2rho[jtype][itype]];
                            arr_phi=phi_r_arr[type2phi[itype][jtype]];
                            for(int i=0;i<no_i;i++)
                            {
                                tmp0=exp(-inv_alpha_sqr*(r-xi[i]*rc)*(r-xi[i]*rc));
                                
                                
                                phi_it_jt_0+=tmp0*arr_phi[0][i];
                                phi_it_jt_1+=tmp0*arr_phi[1][i];
                                phi_it_jt_2+=tmp0*arr_phi[2][i];
                                
                                
                                rho_it_jt_0+=tmp0*arr_rho_ij[0][i];
                                rho_it_jt_1+=tmp0*arr_rho_ij[1][i];
                                rho_it_jt_2+=tmp0*arr_rho_ij[2][i];
                                
                                rho_jt_it_0+=tmp0*arr_rho_ji[0][i];
                                rho_jt_it_1+=tmp0*arr_rho_ji[1][i];
                                rho_jt_it_2+=tmp0*arr_rho_ji[2][i];
                            }
                            
                            
                            rho_it_jt_2*=2.0*r_inv;
                            rho_it_jt_2+=-4.0*rho_it_jt_1;
                            rho_it_jt_2+=(2.0*rsq-alpha_sqr)*r_inv*rho_it_jt_0;
                            rho_it_jt_2*=inv_alpha_sqr*inv_alpha_sqr;
                            rho_it_jt_1*=2.0;
                            rho_it_jt_1-=(2.0*rsq+alpha_sqr)*r_inv*rho_it_jt_0;
                            rho_it_jt_1*=r_inv*inv_alpha*inv_alpha_sqr;
                            rho_it_jt_0*=r_inv*inv_alpha;
                            
                            
                            rho_jt_it_2*=2.0*r_inv;
                            rho_jt_it_2+=-4.0*rho_jt_it_1;
                            rho_jt_it_2+=(2.0*rsq-alpha_sqr)*r_inv*rho_jt_it_0;
                            rho_jt_it_2*=inv_alpha_sqr*inv_alpha_sqr;
                            rho_jt_it_1*=2.0;
                            rho_jt_it_1-=(2.0*rsq+alpha_sqr)*r_inv*rho_jt_it_0;
                            rho_jt_it_1*=r_inv*inv_alpha*inv_alpha_sqr;
                            rho_jt_it_0*=r_inv*inv_alpha;
                            
                            phi_it_jt_2*=2.0*r_inv;
                            phi_it_jt_2+=-4.0*phi_it_jt_1;
                            phi_it_jt_2+=(2.0*rsq-alpha_sqr)*r_inv*phi_it_jt_0;
                            phi_it_jt_2*=inv_alpha_sqr*inv_alpha_sqr;
                            phi_it_jt_1*=2.0;
                            phi_it_jt_1-=(2.0*rsq+alpha_sqr)*r_inv*phi_it_jt_0;
                            phi_it_jt_1*=r_inv*inv_alpha*inv_alpha_sqr;
                            phi_it_jt_0*=r_inv*inv_alpha;
                            
                            
                            
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
type0 ForceField_eam_dmd_lg_1::energy_calc()
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
    type0 dx0,dx1,dx2,rsq;
    type0 alpha;
    type0 r=0.0,r_inv=0.0;
    int m;
    type0* coef;
    type0** arr_rho_ij;
    type0** arr_rho_ji;
    type0** arr_phi;
    type0 inv_alpha,inv_alpha_sqr,alpha_sqr,tmp0,tmp1;
    
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
            r=sqrt(rsq);
            r_inv=1.0/r;
            
            for(int itype=0;itype<no_types; itype++)
            {
                for(int jtype=0;jtype<no_types; jtype++)
                {
                    alpha=sqrt(x[icomp+3+itype]*x[icomp+3+itype]+x[jcomp+3+jtype]*x[jcomp+3+jtype]);
                    
                    if(alpha<alpha_max)
                    {
                        inv_alpha=1.0/alpha;
                        inv_alpha_sqr=inv_alpha*inv_alpha;
                        alpha_sqr=alpha*alpha;
                        
                        phi_it_jt_0=0.0;
                        rho_it_jt_0=0.0;
                        rho_jt_it_0=0.0;
                        
                        
                        arr_rho_ij=rho_r_arr[type2rho[itype][jtype]];
                        arr_rho_ji=rho_r_arr[type2rho[jtype][itype]];
                        arr_phi=phi_r_arr[type2phi[itype][jtype]];
                        for(int i=0;i<no_i;i++)
                        {
                            tmp0=exp(-inv_alpha_sqr*(r-xi[i]*rc)*(r-xi[i]*rc));
                            
                            phi_it_jt_0+=tmp0*arr_phi[0][i];
                            
                            rho_it_jt_0+=tmp0*arr_rho_ij[0][i];
                            
                            rho_jt_it_0+=tmp0*arr_rho_ji[0][i];
                            
                        }
                        
                        
                        rho_it_jt_0*=r_inv*inv_alpha;
                        
                        rho_jt_it_0*=r_inv*inv_alpha;
                        
                        phi_it_jt_0*=r_inv*inv_alpha;
                        
                        
                        
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
void ForceField_eam_dmd_lg_1::init()
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
void ForceField_eam_dmd_lg_1::fin()
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
void ForceField_eam_dmd_lg_1::coef(int narg,char** arg)
{
    type0 kb,T,hbar,r_0;
    if (narg<9)
        error->abort("ff_coef for ff eam_dmd "
                     "should at least have 8 arguments");
    
    
    rsq_crd=atof(arg[1]);
    if(rsq_crd<=0.0)
        error->abort("coordination cutoff for ff eam_dmd "
                     "should be greater than 0.0");
    
    rsq_crd*=rsq_crd;
    
    set_weight_abs(atoi(arg[2]));
    
    alpha_max=atof(arg[3]);
    if(alpha_max<=0.0)
        error->abort("minimum alpha for ff eam_dmd "
                     "should be greater than 0.0");
    
    r_0=atof(arg[4]);
    if(r_0<=0.0)
        error->abort("maximum alpha for ff eam_dmd "
                     "should be greater than 0.0");
    
    
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
        error->abort("wrong coeff command "
                     "for eam Force Field");
    
    rc=(static_cast<type0>(nr)-1.0)*dr;
    set_arrays();
    tmp_clean_up();
    
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
    
    if(alpha_max==0.0)
        error->abort("maximum alpha cannot be zero");
    
    rho_max=(static_cast<type0>(nrho)-1.0)*drho;
    cut_sq_0=rc*rc;
    mod_rc=rc+r_0*alpha_max;
    cut_sq_mod_0=mod_rc*mod_rc;
    
}
/*--------------------------------------------
 read from file and assign DeltaE
 --------------------------------------------*/
void ForceField_eam_dmd_lg_1::setup_delta_e(char* file_name)
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
void ForceField_eam_dmd_lg_1::set_funcfl(int no_files
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
    tmp_allocate();
    
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
            rho_r[type2rho[ityp][0]][i][0]=interpolate(tmp_rho[ityp],nrs[ityp],p,k)*static_cast<type0>(i)*dr;
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
                
                phi_r[type2phi[ityp][jtyp]][i][0]=27.2*0.529*tmp0*tmp1;
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
void ForceField_eam_dmd_lg_1::set_setfl(int no_files
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
    tmp_allocate();
    
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
                rho_r[component][i][0]=tmp[ipos+i]*static_cast<type0>(i)*dr;
            
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
                        phi_r[component][i][0]=tmp[ipos+i];
                    
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
void ForceField_eam_dmd_lg_1::set_fs(int no_files
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
    tmp_allocate();
    
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
                        rho_r[component][i][0]=tmp[ipos+i]*static_cast<type0>(i)*dr;
                    
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
                        phi_r[component][i][0]=tmp[ipos+i];
                    
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
int ForceField_eam_dmd_lg_1::line_read(FILE* file
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
 allocate the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_lg_1::allocate()
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
    
    
    
    CREATE1D(F_arr,no_types);
    for(int ityp=0;ityp<no_types;ityp++)
        CREATE1D(F_arr[ityp],nrho);
    for(int ityp=0;ityp<no_types;ityp++)
        for(int i=0;i<nrho;i++)
            CREATE1D(F_arr[ityp][i],9);
    
    
    
    
    CREATE1D(phi_r_arr,no_types*(no_types+1)/2);
    for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
        CREATE1D(phi_r_arr[ityp],3);
    for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
        for(int i=0;i<3;i++)
            CREATE1D(phi_r_arr[ityp][i],no_i);
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        
        CREATE1D(rho_r_arr,no_types*no_types);
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(rho_r_arr[ityp],3);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<3;i++)
                CREATE1D(rho_r_arr[type2rho[ityp][0]][i],no_i);
        
        for(int ityp=0;ityp<no_types;ityp++)
            for(int jtyp=0;jtyp<no_types;jtyp++)
                rho_r_arr[type2rho[ityp][jtyp]]
                =rho_r_arr[type2rho[ityp][0]];
        
    }
    else if(eam_mode==FINNIS_FL)
    {
        
        CREATE1D(rho_r_arr,no_types*no_types);
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            CREATE1D(rho_r_arr[ityp],3);
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            for(int i=0;i<3;i++)
                CREATE1D(rho_r_arr[ityp][i],no_i);
        
    }
    
    
    allocated=1;
}
/*--------------------------------------------
 clean up the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_lg_1::clean_up()
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
    
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        for(int i=0;i<nrho;i++)
            delete [] F_arr[ityp][i];
        delete [] F_arr[ityp];
    }
    delete [] F_arr;
    
    for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
    {
        for(int i=0;i<3;i++)
            delete [] phi_r_arr[ityp][i];
        
        delete [] phi_r_arr[ityp];
    }
    delete [] phi_r_arr;
    
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        
        for(int ityp=0;ityp<no_types;ityp++)
        {
            for(int i=0;i<3;i++)
                delete [] rho_r_arr[ityp][i];
            delete [] rho_r_arr[ityp];
        }
        delete [] rho_r_arr;
        
    }
    else if(eam_mode==FINNIS_FL)
    {
        
        for(int ityp=0;ityp<no_types*no_types;ityp++)
        {
            for(int i=0;i<3;i++)
                delete [] rho_r_arr[ityp][i];
            delete [] rho_r_arr[ityp];
        }
        delete [] rho_r_arr;
        
    }
    
    allocated=0;
}
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_lg_1::tmp_allocate()
{
    int no_types=atom_types->no_types;
    
    
    CREATE1D(phi_r,no_types*(no_types+1)/2);
    for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
        CREATE1D(phi_r[ityp],nr);
    for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
        for(int i=0;i<nr;i++)
            CREATE1D(phi_r[ityp][i],7);
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        CREATE1D(rho_r,no_types*no_types);
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(rho_r[ityp],nr);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(rho_r[type2rho[ityp][0]][i],7);
        
        
        for(int ityp=0;ityp<no_types;ityp++)
            for(int jtyp=0;jtyp<no_types;jtyp++)
                rho_r[type2rho[ityp][jtyp]]
                =rho_r[type2rho[ityp][0]];
        
    }
    else if(eam_mode==FINNIS_FL)
    {
        CREATE1D(rho_r,no_types*no_types);
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            CREATE1D(rho_r[ityp],nr);
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(rho_r[ityp][i],7);
        
    }
    
}
/*--------------------------------------------
 clean up the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_lg_1::tmp_clean_up()
{
    
    int no_types=atom_types->no_types;
    
    for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
    {
        for(int i=0;i<nr;i++)
            delete [] phi_r[ityp][i];
        delete [] phi_r[ityp];
    }
    delete [] phi_r;
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        for(int ityp=0;ityp<no_types;ityp++)
        {
            for(int i=0;i<nr;i++)
                delete [] rho_r[ityp][i];
            delete [] rho_r[ityp];
        }
        delete [] rho_r;
        
    }
    else if(eam_mode==FINNIS_FL)
    {
        for(int ityp=0;ityp<no_types*no_types;ityp++)
        {
            for(int i=0;i<nr;i++)
                delete [] rho_r[ityp][i];
            delete [] rho_r[ityp];
        }
        delete [] rho_r;
        
        
    }
    
    
}
/*--------------------------------------------
 set the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_lg_1::set_arrays()
{
    for(int i=0;i<no_types;i++)
        interpolate_m(nrho,drho,F_arr[i]);
    
    for(int i=0;i<no_types;i++)
        for(int j=0;j<i+1;j++)
            interpolate(nr,dr,phi_r[type2phi[i][j]]);
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        for(int i=0;i<no_types;i++)
            interpolate(nr,dr,rho_r[i]);
    }
    else if(eam_mode==FINNIS_FL)
    {
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<no_types;j++)
                interpolate(nr,dr,rho_r[type2rho[i][j]]);
    }
    
    type0 rtmp,tmp0,p,p2,p3;
    type0* coef;
    int m;
    
    
    for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
    {
        for(int i=0;i<no_i;i++)
        {
            rtmp=xi[i]*rc;
            p=fabs(rtmp)*dr_inv;
            m=static_cast<int>(p);
            m=MIN(m,nr-2);
            p-=m;
            p=MIN(p,1.0);
            p2=p*p;
            p3=p2*p;
            coef=phi_r[ityp][m];
            tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
            if(rtmp<0.0)
                tmp0*=-1.0;
            tmp0*=wi[i]*rc*PI_IN_SQ;
            phi_r_arr[ityp][0][i]=tmp0;
            tmp0*=rtmp;
            phi_r_arr[ityp][1][i]=tmp0;
            tmp0*=rtmp;
            phi_r_arr[ityp][2][i]=tmp0;
            
        }
    }
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        for(int ityp=0;ityp<no_types;ityp++)
        {
            for(int i=0;i<no_i;i++)
            {
                rtmp=xi[i]*rc;
                p=fabs(rtmp)*dr_inv;
                m=static_cast<int>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);
                p2=p*p;
                p3=p2*p;
                coef=rho_r[ityp][m];
                tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                if(rtmp<0.0)
                    tmp0*=-1.0;
                tmp0*=wi[i]*rc*PI_IN_SQ;
                rho_r_arr[ityp][0][i]=tmp0;
                tmp0*=rtmp;
                rho_r_arr[ityp][1][i]=tmp0;
                tmp0*=rtmp;
                rho_r_arr[ityp][2][i]=tmp0;
                
            }
        }
    }
    else if(eam_mode==FINNIS_FL)
    {
        
        for(int ityp=0;ityp<no_types*no_types;ityp++)
        {
            for(int i=0;i<no_i;i++)
            {
                rtmp=xi[i]*rc;
                p=fabs(rtmp)*dr_inv;
                m=static_cast<int>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);
                p2=p*p;
                p3=p2*p;
                coef=rho_r[ityp][m];
                tmp0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                if(rtmp<0.0)
                    tmp0*=-1.0;
                tmp0*=wi[i]*rc*PI_IN_SQ;
                rho_r_arr[ityp][0][i]=tmp0;
                tmp0*=rtmp;
                rho_r_arr[ityp][1][i]=tmp0;
                tmp0*=rtmp;
                rho_r_arr[ityp][2][i]=tmp0;
                
            }
        }
    }
    
}
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_lg_1::interpolate(int n
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
type0 ForceField_eam_dmd_lg_1::interpolate(type0* arr
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
void ForceField_eam_dmd_lg_1::interpolate_m(int n
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
void ForceField_eam_dmd_lg_1::set_weight_abs(int n)
{
    
    if(n<1 || n>14)
        error->abort("number of gaussian "
                     "points for ff eam_dmd should be "
                     "between 1 & 14");
    
    if(no_i)
    {
        delete [] xi;
        delete [] wi;
    }
    
    no_i=n;
    
    CREATE1D(xi,no_i);
    CREATE1D(wi,no_i);
    
    
    if(no_i==1)
    {
        wi[0]=2.00000000000;
        
        xi[0]=0.00000000000;
    }
    else if(no_i==2)
    {
        wi[0]=1.0000000000;
        wi[1]=1.0000000000;
        
        xi[0]=0.577350269190;
        xi[1]=-0.577350269190;
    }
    else if(no_i==3)
    {
        wi[0]=0.888888888889;
        wi[1]=0.55555555556;
        wi[2]=0.55555555556;
        
        xi[0]=0.00000000000;
        xi[1]=0.774596669241;
        xi[2]=-0.774596669241;
    }
    else if(no_i==4)
    {
        wi[0]=0.65214515486;
        wi[1]=0.65214515486;
        wi[2]=0.34785484514;
        wi[3]=0.34785484514;
        
        xi[0]=0.339981043585;
        xi[1]=-0.339981043585;
        xi[2]=0.861136311594;
        xi[3]=-0.861136311594;
    }
    else if(no_i==5)
    {
        wi[0]=0.568888888889;
        wi[1]=0.47862867050;
        wi[2]=0.47862867050;
        wi[3]=0.23692688506;
        wi[4]=0.23692688506;
        
        xi[0]=0.00000000000;
        xi[1]=0.538469310106;
        xi[2]=-0.538469310106;
        xi[3]=0.906179845939;
        xi[4]=-0.906179845939;
    }
    else if(no_i==6)
    {
        wi[0]=0.17132449238;
        wi[1]=0.17132449238;
        wi[2]=0.467913934573;
        wi[3]=0.467913934573;
        wi[4]=0.36076157305;
        wi[5]=0.36076157305;
        
        xi[0]=0.932469514203;
        xi[1]=-0.932469514203;
        xi[2]=0.238619186083;
        xi[3]=-0.238619186083;
        xi[4]=0.661209386466;
        xi[5]=-0.661209386466;
    }
    else if(no_i==7)
    {
        wi[0]=0.417959183673;
        wi[1]=0.12948496617;
        wi[2]=0.12948496617;
        wi[3]=0.27970539149;
        wi[4]=0.27970539149;
        wi[5]=0.381830050505;
        wi[6]=0.381830050505;
        
        xi[0]=0.00000000000;
        xi[1]=0.949107912343;
        xi[2]=-0.949107912343;
        xi[3]=0.741531185599;
        xi[4]=-0.741531185599;
        xi[5]=0.405845151377;
        xi[6]=-0.405845151377;
    }
    else if(no_i==8)
    {
        wi[0]=0.362683783378;
        wi[1]=0.362683783378;
        wi[2]=0.313706645878;
        wi[3]=0.313706645878;
        wi[4]=0.22238103445;
        wi[5]=0.22238103445;
        wi[6]=0.10122853629;
        wi[7]=0.10122853629;
        
        xi[0]=0.183434642496;
        xi[1]=-0.183434642496;
        xi[2]=0.525532409916;
        xi[3]=-0.525532409916;
        xi[4]=0.796666477414;
        xi[5]=-0.796666477414;
        xi[6]=0.960289856498;
        xi[7]=-0.960289856498;
    }
    else if(no_i==9)
    {
        wi[0]=0.330239355001;
        wi[1]=0.18064816069;
        wi[2]=0.18064816069;
        wi[3]=0.08127438836;
        wi[4]=0.08127438836;
        wi[5]=0.312347077040;
        wi[6]=0.312347077040;
        wi[7]=0.26061069640;
        wi[8]=0.260610696403;
        
        xi[0]=0.00000000000;
        xi[1]=0.836031107327;
        xi[2]=-0.836031107327;
        xi[3]=0.968160239508;
        xi[4]=-0.968160239508;
        xi[5]=0.324253423404;
        xi[6]=-0.324253423404;
        xi[7]=0.613371432701;
        xi[8]=-0.613371432701;
    }
    else if(no_i==10)
    {
        wi[0]=0.295524224715;
        wi[1]=0.295524224715;
        wi[2]=0.269266719310;
        wi[3]=0.269266719310;
        wi[4]=0.219086362516;
        wi[5]=0.21908636252;
        wi[6]=0.14945134915;
        wi[7]=0.14945134915;
        wi[8]=0.06667134431;
        wi[9]=0.06667134431;
        
        xi[0]=-0.148874338982;
        xi[1]=0.148874338982;
        xi[2]=-0.433395394129;
        xi[3]=0.433395394129;
        xi[4]=-0.679409568299;
        xi[5]=0.679409568299;
        xi[6]=-0.865063366689;
        xi[7]=0.865063366689;
        xi[8]=-0.973906528517;
        xi[9]=0.973906528517;
    }
    else if(no_i==11)
    {
        wi[0]=0.272925086778;
        wi[1]=0.262804544510;
        wi[2]=0.262804544510;
        wi[3]=0.233193764592;
        wi[4]=0.233193764592;
        wi[5]=0.186290210928;
        wi[6]=0.18629021093;
        wi[7]=0.12558036946;
        wi[8]=0.12558036946;
        wi[9]=0.05566856712;
        wi[10]=0.05566856712;
        
        xi[0]=0.00000000000;
        xi[1]=-0.269543155952;
        xi[2]=0.269543155952;
        xi[3]=-0.519096129207;
        xi[4]=0.519096129207;
        xi[5]=-0.730152005574;
        xi[6]=0.730152005574;
        xi[7]=-0.887062599768;
        xi[8]=0.887062599768;
        xi[9]=-0.978228658146;
        xi[10]=0.978228658146;
    }
    else if(no_i==12)
    {
        wi[0]=0.249147045813;
        wi[1]=0.249147045813;
        wi[2]=0.233492536538;
        wi[3]=0.233492536538;
        wi[4]=0.203167426723;
        wi[5]=0.203167426723;
        wi[6]=0.160078328543;
        wi[7]=0.16007832854;
        wi[8]=0.10693932600;
        wi[9]=0.10693932600;
        wi[10]=0.04717533639;
        wi[11]=0.04717533639;
        
        xi[0]=-0.125233408511;
        xi[1]=0.125233408511;
        xi[2]=-0.367831498998;
        xi[3]=0.367831498998;
        xi[4]=-0.587317954287;
        xi[5]=0.587317954287;
        xi[6]=-0.769902674194;
        xi[7]=0.769902674194;
        xi[8]=-0.904117256370;
        xi[9]=0.904117256370;
        xi[10]=-0.981560634247;
        xi[11]=0.981560634247;
    }
    else if(no_i==13)
    {
        wi[0]=0.232551553231;
        wi[1]=0.226283180263;
        wi[2]=0.226283180263;
        wi[3]=0.207816047537;
        wi[4]=0.207816047537;
        wi[5]=0.178145980762;
        wi[6]=0.178145980762;
        wi[7]=0.138873510220;
        wi[8]=0.13887351022;
        wi[9]=0.09212149984;
        wi[10]=0.09212149984;
        wi[11]=0.04048400477;
        wi[12]=0.04048400477;
        
        xi[0]=0.00000000000;
        xi[1]=-0.230458315955;
        xi[2]=0.230458315955;
        xi[3]=-0.448492751036;
        xi[4]=0.448492751036;
        xi[5]=-0.642349339440;
        xi[6]=0.642349339440;
        xi[7]=-0.801578090733;
        xi[8]=0.801578090733;
        xi[9]=-0.917598399223;
        xi[10]=0.917598399223;
        xi[11]=-0.984183054719;
        xi[12]=0.984183054719;
    }
    else if(no_i==14)
    {
        wi[0]=0.215263853463;
        wi[1]=0.215263853463;
        wi[2]=0.205198463721;
        wi[3]=0.205198463721;
        wi[4]=0.185538397478;
        wi[5]=0.185538397478;
        wi[6]=0.157203167158;
        wi[7]=0.157203167158;
        wi[8]=0.121518570688;
        wi[9]=0.12151857069;
        wi[10]=0.08015808716;
        wi[11]=0.08015808716;
        wi[12]=0.03511946033;
        wi[13]=0.03511946033;
        
        xi[0]=-0.108054948707;
        xi[1]=0.108054948707;
        xi[2]=-0.319112368928;
        xi[3]=0.319112368928;
        xi[4]=-0.515248636358;
        xi[5]=0.515248636358;
        xi[6]=-0.687292904812;
        xi[7]=0.687292904812;
        xi[8]=-0.827201315070;
        xi[9]=0.827201315070;
        xi[10]=-0.928434883664;
        xi[11]=0.928434883664;
        xi[12]=-0.986283808697;
        xi[13]=0.986283808697;
    }
    
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_lg_1::create_2nd_neigh_lst()
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
type0 ForceField_eam_dmd_lg_1::
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
type0 ForceField_eam_dmd_lg_1::
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
type0 ForceField_eam_dmd_lg_1::
dmat1(type0 fi,type0 crdi,type0 fj,type0 crdj
      ,int itype)
{
    
    return 0.0;
}
/*--------------------------------------------
 return M_{ij}^{\alpha}
 --------------------------------------------*/
type0 ForceField_eam_dmd_lg_1::calc_ent(type0 x)
{
    if(x<1.0e-8 || 0.99999999<x)
        return 0.0;
    else
        return x*log(x)+(1-x)*log(1-x);
    
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_lg_1::c_d_calc()
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
type0 ForceField_eam_dmd_lg_1::g_calc(int chk
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
