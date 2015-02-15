#include "ff_eam_dmd_lg_0.h"
#include "neighbor.h"
#include "atom_types.h"
using namespace MAPP_NS;
#define PI_IN_SQ 0.564189583547756286948079451561
enum{NOT_SET,FUNC_FL,SET_FL,FINNIS_FL};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam_dmd_lg_0::
ForceField_eam_dmd_lg_0(MAPP* mapp) : ForceField(mapp)
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
ForceField_eam_dmd_lg_0::~ForceField_eam_dmd_lg_0()
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
void ForceField_eam_dmd_lg_0::
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
    type0 alpha_sq,alpha_sq_inv,tmp0,tmp1;
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
                        alpha=(x[icomp+3+itype]*x[jcomp+3+jtype])/(x[icomp+3+itype]+x[jcomp+3+jtype]);
                        
                        if(alpha_min<alpha)
                        {
                            alpha_sq=sqrt(alpha);
                            alpha_sq_inv=1.0/alpha_sq;
                            
                            phi_it_jt_0=phi_it_jt_1=phi_it_jt_2=0.0;
                            rho_it_jt_0=rho_it_jt_1=rho_it_jt_2=0.0;
                            rho_jt_it_0=rho_jt_it_1=rho_jt_it_2=0.0;
                            
                            
                            arr_rho_ij=rho_r_arr[type2rho[itype][jtype]];
                            arr_rho_ji=rho_r_arr[type2rho[jtype][itype]];
                            arr_phi=phi_r_arr[type2phi[itype][jtype]];
                            for(int i=0;i<no_i;i++)
                            {
                                tmp0=exp(-alpha*(r-xi[i]*rc)*(r-xi[i]*rc));
                                
                                
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

                            
                            
                            rho_it_jt_2*=-alpha_sq*r_inv;
                            rho_it_jt_2+=2.0*alpha_sq*rho_it_jt_1;
                            rho_it_jt_2+=(0.5*alpha_sq_inv*r_inv-r*alpha_sq)*rho_it_jt_0;
                            rho_it_jt_1*=2.0*alpha*alpha_sq*r_inv;
                            rho_it_jt_1+=(-alpha_sq*r_inv*r_inv-2.0*alpha*alpha_sq)*rho_it_jt_0;
                            rho_it_jt_0*=alpha_sq*r_inv;


                            
                            rho_jt_it_2*=-alpha_sq*r_inv;
                            rho_jt_it_2+=2.0*alpha_sq*rho_jt_it_1;
                            rho_jt_it_2+=(0.5*alpha_sq_inv*r_inv-r*alpha_sq)*rho_jt_it_0;
                            rho_jt_it_1*=2.0*alpha*alpha_sq*r_inv;
                            rho_jt_it_1+=(-alpha_sq*r_inv*r_inv-2.0*alpha*alpha_sq)*rho_jt_it_0;
                            rho_jt_it_0*=alpha_sq*r_inv;
                            
                            phi_it_jt_2*=-alpha_sq*r_inv;
                            phi_it_jt_2+=2.0*alpha_sq*phi_it_jt_1;
                            phi_it_jt_2+=(0.5*alpha_sq_inv*r_inv-r*alpha_sq)*phi_it_jt_0;
                            phi_it_jt_1*=2.0*alpha*alpha_sq*r_inv;
                            phi_it_jt_1+=(-alpha_sq*r_inv*r_inv-2.0*alpha*alpha_sq)*phi_it_jt_0;
                            phi_it_jt_0*=alpha_sq*r_inv;
                            
                            
                            
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
            
            nrgy_strss[0]+=1.5*kbT*c[iatm*no_types+itype]*log(x[icomp+3+itype])
            +c[iatm*no_types+itype]*c_0[itype];
            
            f[icomp+3+itype]-=1.5*kbT*c[iatm*no_types+itype]/x[icomp+3+itype];
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
                    ((x[icomp+3+itype]+x[jcomp+3+jtype])*(x[icomp+3+itype]+x[jcomp+3+jtype]));
                    
                    
                    
                    
                    if(apair!=0.0 || fpair!=0.0)
                    {
                        f[icomp]+=dx0*fpair;
                        f[icomp+1]+=dx1*fpair;
                        f[icomp+2]+=dx2*fpair;
                        f[icomp+3+itype]+=apair*x[jcomp+3+jtype]*x[jcomp+3+jtype];
                        
                        if (jatm<natms)
                        {
                            f[jcomp]-=dx0*fpair;
                            f[jcomp+1]-=dx1*fpair;
                            f[jcomp+2]-=dx2*fpair;
                            f[jcomp+3+jtype]+=apair*x[icomp+3+itype]*x[icomp+3+itype];
                            
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
type0 ForceField_eam_dmd_lg_0::energy_calc()
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
    type0 alpha_sq,alpha_sq_inv,tmp0,tmp1;
    
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
                        alpha=(x[icomp+3+itype]*x[jcomp+3+jtype])/(x[icomp+3+itype]+x[jcomp+3+jtype]);
                        
                        if(alpha_min<alpha)
                        {
                            alpha_sq=sqrt(alpha);
                            alpha_sq_inv=1.0/alpha_sq;
                            
                            phi_it_jt_0=0.0;
                            rho_it_jt_0=0.0;
                            rho_jt_it_0=0.0;
                            
                            
                            arr_rho_ij=rho_r_arr[type2rho[itype][jtype]];
                            arr_rho_ji=rho_r_arr[type2rho[jtype][itype]];
                            arr_phi=phi_r_arr[type2phi[itype][jtype]];
                            for(int i=0;i<no_i;i++)
                            {
                                tmp0=exp(-alpha*(r-xi[i]*rc)*(r-xi[i]*rc));
                                
                                
                                phi_it_jt_0+=tmp0*arr_phi[0][i];
                                
                                rho_it_jt_0+=tmp0*arr_rho_ij[0][i];
                                
                                rho_jt_it_0+=tmp0*arr_rho_ji[0][i];
                            }
                            
                            
                            rho_it_jt_0*=alpha_sq*r_inv;
                            
                            rho_jt_it_0*=alpha_sq*r_inv;
                            
                            phi_it_jt_0*=alpha_sq*r_inv;
                            
                            
                            
                            
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
            
            en+=1.5*kbT*c[iatm*no_types+itype]*log(x[icomp+3+itype])
            +c[iatm*no_types+itype]*c_0[itype];
            
        }
        
    }
    
    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    
    return en_tot;
}
/*--------------------------------------------
 init before running
 --------------------------------------------*/
void ForceField_eam_dmd_lg_0::init()
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
    mu_n=atoms->add<type0>(1,no_types,"mu");
    s_n=atoms->add<type0>(1,no_types,"s");
    t_n=atoms->add<type0>(1,2*no_types,"t");
    v_n=atoms->add<type0>(1,2*no_types,"v");
    crd_n=atoms->add<type0>(1,no_types,"crd");
    
}
/*--------------------------------------------
 fin after running
 --------------------------------------------*/
void ForceField_eam_dmd_lg_0::fin()
{
    atoms->del(crd_n);
    atoms->del(v_n);
    atoms->del(t_n);
    atoms->del(s_n);
    atoms->del(mu_n);
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
void ForceField_eam_dmd_lg_0::coef(int narg,char** arg)
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
    
    alpha_min=atof(arg[3]);
    if(alpha_min<=0.0)
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
        error->abort("unknown file format "
        "for ff eam_dmd: %s",arg[7]);
    
    rc=(static_cast<type0>(nr)-1.0)*dr;
    set_arrays();

    /* test */
    /*
    for(int i=0;i<nr;i++)
    {
        printf("%e ",i*dr);
        for(int j=0;j<no_types;j++)
        {
            for(int k=0;k<no_types;k++)
            {
                printf("%e ",rho_r[type2rho[j][k]][i][0]);
            }
            
        }
        
        for(int j=0;j<no_types;j++)
        {
            for(int k=0;k<no_types;k++)
            {
                printf("%e ",phi_r[type2phi[j][k]][i][0]);
            }
            
        }
        printf("\n");
    }
    */
    /* test */
    
    
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

    
    rho_max=(static_cast<type0>(nrho)-1.0)*drho;
    cut_sq_0=rc*rc;
    mod_rc=rc+r_0/sqrt(alpha_min);
    cut_sq_mod_0=mod_rc*mod_rc;

    
}
/*--------------------------------------------
 read from file and assign DeltaE
 --------------------------------------------*/
void ForceField_eam_dmd_lg_0::setup_delta_e(char* file_name)
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
void ForceField_eam_dmd_lg_0::set_funcfl(int no_files
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
void ForceField_eam_dmd_lg_0::set_setfl(int no_files
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
void ForceField_eam_dmd_lg_0::set_fs(int no_files
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
int ForceField_eam_dmd_lg_0::line_read(FILE* file
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
void ForceField_eam_dmd_lg_0::allocate()
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
void ForceField_eam_dmd_lg_0::clean_up()
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
void ForceField_eam_dmd_lg_0::tmp_allocate()
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
void ForceField_eam_dmd_lg_0::tmp_clean_up()
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
void ForceField_eam_dmd_lg_0::set_arrays()
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
            //tmp0=coef[1]*p+coef[0];
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
                //tmp0=coef[1]*p+coef[0];
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
void ForceField_eam_dmd_lg_0::interpolate(int n
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
type0 ForceField_eam_dmd_lg_0::interpolate(type0* arr
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
void ForceField_eam_dmd_lg_0::interpolate_m(int n
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
void ForceField_eam_dmd_lg_0::set_weight_abs(int n)
{
    
    if(n<1 || n>70)
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
        
        xi[0]=0;
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
        
        xi[0]=0;
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
        
        xi[0]=0;
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
        
        xi[0]=0;
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
        
        xi[0]=0;
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
        
        xi[0]=0;
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
        
        xi[0]=0;
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
    else if(no_i==15)
    {
        wi[0]=0.202578241926;
        wi[1]=0.198431485327;
        wi[2]=0.198431485327;
        wi[3]=0.186161000016;
        wi[4]=0.186161000016;
        wi[5]=0.166269205817;
        wi[6]=0.166269205817;
        wi[7]=0.139570677926;
        wi[8]=0.13957067793;
        wi[9]=0.107159220467;
        wi[10]=0.10715922047;
        wi[11]=0.07036604749;
        wi[12]=0.07036604749;
        wi[13]=0.03075324200;
        wi[14]=0.03075324200;
        
        xi[0]=0;
        xi[1]=-0.201194093997;
        xi[2]=0.201194093997;
        xi[3]=-0.394151347078;
        xi[4]=0.394151347078;
        xi[5]=-0.570972172609;
        xi[6]=0.570972172609;
        xi[7]=-0.724417731360;
        xi[8]=0.724417731360;
        xi[9]=-0.848206583410;
        xi[10]=0.848206583410;
        xi[11]=-0.937273392401;
        xi[12]=0.937273392401;
        xi[13]=-0.987992518020;
        xi[14]=0.987992518020;
    }
    else if(no_i==16)
    {
        wi[0]=0.189450610455;
        wi[1]=0.189450610455;
        wi[2]=0.182603415045;
        wi[3]=0.182603415045;
        wi[4]=0.169156519395;
        wi[5]=0.169156519395;
        wi[6]=0.149595988817;
        wi[7]=0.149595988817;
        wi[8]=0.124628971256;
        wi[9]=0.12462897126;
        wi[10]=0.095158511682;
        wi[11]=0.09515851168;
        wi[12]=0.06225352394;
        wi[13]=0.06225352394;
        wi[14]=0.02715245941;
        wi[15]=0.02715245941;
        
        xi[0]=-0.0950125098376;
        xi[1]=0.0950125098376;
        xi[2]=-0.281603550779;
        xi[3]=0.281603550779;
        xi[4]=-0.458016777657;
        xi[5]=0.458016777657;
        xi[6]=-0.617876244403;
        xi[7]=0.617876244403;
        xi[8]=-0.755404408355;
        xi[9]=0.755404408355;
        xi[10]=-0.865631202388;
        xi[11]=0.865631202388;
        xi[12]=-0.944575023073;
        xi[13]=0.944575023073;
        xi[14]=-0.989400934992;
        xi[15]=0.989400934992;
    }
    else if(no_i==17)
    {
        wi[0]=0.179446470356;
        wi[1]=0.176562705367;
        wi[2]=0.176562705367;
        wi[3]=0.168004102156;
        wi[4]=0.168004102156;
        wi[5]=0.154045761077;
        wi[6]=0.154045761077;
        wi[7]=0.135136368469;
        wi[8]=0.135136368469;
        wi[9]=0.111883847193;
        wi[10]=0.11188384719;
        wi[11]=0.085036148317;
        wi[12]=0.08503614832;
        wi[13]=0.05545952937;
        wi[14]=0.05545952937;
        wi[15]=0.02414830287;
        wi[16]=0.02414830287;
        
        xi[0]=0;
        xi[1]=-0.178484181496;
        xi[2]=0.178484181496;
        xi[3]=-0.351231763454;
        xi[4]=0.351231763454;
        xi[5]=-0.512690537086;
        xi[6]=0.512690537086;
        xi[7]=-0.657671159217;
        xi[8]=0.657671159217;
        xi[9]=-0.781514003897;
        xi[10]=0.781514003897;
        xi[11]=-0.880239153727;
        xi[12]=0.880239153727;
        xi[13]=-0.950675521769;
        xi[14]=0.950675521769;
        xi[15]=-0.990575475314;
        xi[16]=0.990575475314;
    }
    else if(no_i==18)
    {
        wi[0]=0.169142382963;
        wi[1]=0.169142382963;
        wi[2]=0.164276483746;
        wi[3]=0.164276483746;
        wi[4]=0.154684675126;
        wi[5]=0.154684675126;
        wi[6]=0.140642914671;
        wi[7]=0.140642914671;
        wi[8]=0.122555206711;
        wi[9]=0.122555206711;
        wi[10]=0.100942044106;
        wi[11]=0.10094204411;
        wi[12]=0.076425730255;
        wi[13]=0.07642573025;
        wi[14]=0.04971454889;
        wi[15]=0.04971454889;
        wi[16]=0.02161601353;
        wi[17]=0.02161601353;
        
        xi[0]=-0.0847750130417;
        xi[1]=0.0847750130417;
        xi[2]=-0.251886225692;
        xi[3]=0.251886225692;
        xi[4]=-0.411751161463;
        xi[5]=0.411751161463;
        xi[6]=-0.559770831074;
        xi[7]=0.559770831074;
        xi[8]=-0.691687043060;
        xi[9]=0.691687043060;
        xi[10]=-0.803704958973;
        xi[11]=0.803704958973;
        xi[12]=-0.892602466498;
        xi[13]=0.892602466498;
        xi[14]=-0.955823949571;
        xi[15]=0.955823949571;
        xi[16]=-0.991565168421;
        xi[17]=0.991565168421;
    }
    else if(no_i==19)
    {
        wi[0]=0.161054449849;
        wi[1]=0.158968843394;
        wi[2]=0.158968843394;
        wi[3]=0.152766042066;
        wi[4]=0.152766042066;
        wi[5]=0.142606702174;
        wi[6]=0.142606702174;
        wi[7]=0.128753962539;
        wi[8]=0.128753962539;
        wi[9]=0.111566645547;
        wi[10]=0.111566645547;
        wi[11]=0.091490021622;
        wi[12]=0.09149002162;
        wi[13]=0.069044542738;
        wi[14]=0.06904454274;
        wi[15]=0.04481422677;
        wi[16]=0.04481422677;
        wi[17]=0.01946178823;
        wi[18]=0.01946178823;
        
        xi[0]=0;
        xi[1]=-0.160358645640;
        xi[2]=0.160358645640;
        xi[3]=-0.316564099964;
        xi[4]=0.316564099964;
        xi[5]=-0.464570741376;
        xi[6]=0.464570741376;
        xi[7]=-0.600545304662;
        xi[8]=0.600545304662;
        xi[9]=-0.720966177335;
        xi[10]=0.720966177335;
        xi[11]=-0.822714656537;
        xi[12]=0.822714656537;
        xi[13]=-0.903155903615;
        xi[14]=0.903155903615;
        xi[15]=-0.960208152135;
        xi[16]=0.960208152135;
        xi[17]=-0.992406843844;
        xi[18]=0.992406843844;
    }
    else if(no_i==20)
    {
        wi[0]=0.152753387131;
        wi[1]=0.152753387131;
        wi[2]=0.149172986473;
        wi[3]=0.149172986473;
        wi[4]=0.142096109318;
        wi[5]=0.142096109318;
        wi[6]=0.131688638449;
        wi[7]=0.131688638449;
        wi[8]=0.118194531962;
        wi[9]=0.118194531962;
        wi[10]=0.101930119817;
        wi[11]=0.101930119817;
        wi[12]=0.083276741577;
        wi[13]=0.08327674158;
        wi[14]=0.062672048334;
        wi[15]=0.06267204833;
        wi[16]=0.04060142980;
        wi[17]=0.04060142980;
        wi[18]=0.01761400714;
        wi[19]=0.01761400714;
        
        xi[0]=-0.0765265211335;
        xi[1]=0.0765265211335;
        xi[2]=-0.227785851142;
        xi[3]=0.227785851142;
        xi[4]=-0.373706088715;
        xi[5]=0.373706088715;
        xi[6]=-0.510867001951;
        xi[7]=0.510867001951;
        xi[8]=-0.636053680727;
        xi[9]=0.636053680727;
        xi[10]=-0.746331906460;
        xi[11]=0.746331906460;
        xi[12]=-0.839116971822;
        xi[13]=0.839116971822;
        xi[14]=-0.912234428251;
        xi[15]=0.912234428251;
        xi[16]=-0.963971927278;
        xi[17]=0.963971927278;
        xi[18]=-0.993128599185;
        xi[19]=0.993128599185;
    }
    else if(no_i==21)
    {
        wi[0]=0.146081133650;
        wi[1]=0.144524403990;
        wi[2]=0.144524403990;
        wi[3]=0.139887394791;
        wi[4]=0.139887394791;
        wi[5]=0.132268938633;
        wi[6]=0.132268938633;
        wi[7]=0.121831416054;
        wi[8]=0.121831416054;
        wi[9]=0.108797299167;
        wi[10]=0.108797299167;
        wi[11]=0.093444423456;
        wi[12]=0.093444423456;
        wi[13]=0.076100113628;
        wi[14]=0.07610011363;
        wi[15]=0.057134425427;
        wi[16]=0.05713442543;
        wi[17]=0.03695378977;
        wi[18]=0.03695378977;
        wi[19]=0.01601722826;
        wi[20]=0.01601722826;
        
        xi[0]=0;
        xi[1]=-0.145561854161;
        xi[2]=0.145561854161;
        xi[3]=-0.288021316802;
        xi[4]=0.288021316802;
        xi[5]=-0.424342120207;
        xi[6]=0.424342120207;
        xi[7]=-0.551618835887;
        xi[8]=0.551618835887;
        xi[9]=-0.667138804197;
        xi[10]=0.667138804197;
        xi[11]=-0.768439963476;
        xi[12]=0.768439963476;
        xi[13]=-0.853363364583;
        xi[14]=0.853363364583;
        xi[15]=-0.920099334150;
        xi[16]=0.920099334150;
        xi[17]=-0.967226838566;
        xi[18]=0.967226838566;
        xi[19]=-0.993752170620;
        xi[20]=0.993752170620;
    }
    else if(no_i==22)
    {
        wi[0]=0.139251872856;
        wi[1]=0.139251872856;
        wi[2]=0.136541498346;
        wi[3]=0.136541498346;
        wi[4]=0.131173504787;
        wi[5]=0.131173504787;
        wi[6]=0.123252376811;
        wi[7]=0.123252376811;
        wi[8]=0.112932296081;
        wi[9]=0.112932296081;
        wi[10]=0.100414144443;
        wi[11]=0.100414144443;
        wi[12]=0.085941606217;
        wi[13]=0.085941606217;
        wi[14]=0.069796468425;
        wi[15]=0.06979646842;
        wi[16]=0.052293335153;
        wi[17]=0.05229333515;
        wi[18]=0.03377490158;
        wi[19]=0.03377490158;
        wi[20]=0.01462799530;
        wi[21]=0.01462799530;
        
        xi[0]=-0.0697392733197;
        xi[1]=0.0697392733197;
        xi[2]=-0.207860426688;
        xi[3]=0.207860426688;
        xi[4]=-0.341935820892;
        xi[5]=0.341935820892;
        xi[6]=-0.469355837987;
        xi[7]=0.469355837987;
        xi[8]=-0.587640403507;
        xi[9]=0.587640403507;
        xi[10]=-0.694487263187;
        xi[11]=0.694487263187;
        xi[12]=-0.787816805979;
        xi[13]=0.787816805979;
        xi[14]=-0.865812577720;
        xi[15]=0.865812577720;
        xi[16]=-0.926956772187;
        xi[17]=0.926956772187;
        xi[18]=-0.970060497835;
        xi[19]=0.970060497835;
        xi[20]=-0.994294585482;
        xi[21]=0.994294585482;
    }
    else if(no_i==23)
    {
        wi[0]=0.133654572186;
        wi[1]=0.132462039405;
        wi[2]=0.132462039405;
        wi[3]=0.128905722188;
        wi[4]=0.128905722188;
        wi[5]=0.123049084307;
        wi[6]=0.123049084307;
        wi[7]=0.114996640222;
        wi[8]=0.114996640222;
        wi[9]=0.104892091465;
        wi[10]=0.104892091465;
        wi[11]=0.092915766060;
        wi[12]=0.092915766060;
        wi[13]=0.079281411777;
        wi[14]=0.079281411777;
        wi[15]=0.064232421409;
        wi[16]=0.06423242141;
        wi[17]=0.048037671731;
        wi[18]=0.04803767173;
        wi[19]=0.03098800586;
        wi[20]=0.03098800586;
        wi[21]=0.01341185949;
        wi[22]=0.01341185949;
        
        xi[0]=0;
        xi[1]=-0.133256824298;
        xi[2]=0.133256824298;
        xi[3]=-0.264135680970;
        xi[4]=0.264135680970;
        xi[5]=-0.390301038030;
        xi[6]=0.390301038030;
        xi[7]=-0.509501477846;
        xi[8]=0.509501477846;
        xi[9]=-0.619609875764;
        xi[10]=0.619609875764;
        xi[11]=-0.718661363132;
        xi[12]=0.718661363132;
        xi[13]=-0.804888401619;
        xi[14]=0.804888401619;
        xi[15]=-0.876752358270;
        xi[16]=0.876752358270;
        xi[17]=-0.932971086826;
        xi[18]=0.932971086826;
        xi[19]=-0.972542471218;
        xi[20]=0.972542471218;
        xi[21]=-0.994769334998;
        xi[22]=0.994769334998;
    }
    else if(no_i==24)
    {
        wi[0]=0.127938195347;
        wi[1]=0.127938195347;
        wi[2]=0.125837456347;
        wi[3]=0.125837456347;
        wi[4]=0.121670472928;
        wi[5]=0.121670472928;
        wi[6]=0.115505668054;
        wi[7]=0.115505668054;
        wi[8]=0.107444270116;
        wi[9]=0.107444270116;
        wi[10]=0.097618652104;
        wi[11]=0.097618652104;
        wi[12]=0.086190161532;
        wi[13]=0.086190161532;
        wi[14]=0.073346481411;
        wi[15]=0.073346481411;
        wi[16]=0.059298584915;
        wi[17]=0.05929858492;
        wi[18]=0.044277438817;
        wi[19]=0.04427743882;
        wi[20]=0.02853138863;
        wi[21]=0.02853138863;
        wi[22]=0.01234122980;
        wi[23]=0.01234122980;
        
        xi[0]=-0.0640568928626;
        xi[1]=0.0640568928626;
        xi[2]=-0.191118867474;
        xi[3]=0.191118867474;
        xi[4]=-0.315042679696;
        xi[5]=0.315042679696;
        xi[6]=-0.433793507626;
        xi[7]=0.433793507626;
        xi[8]=-0.545421471389;
        xi[9]=0.545421471389;
        xi[10]=-0.648093651937;
        xi[11]=0.648093651937;
        xi[12]=-0.740124191579;
        xi[13]=0.740124191579;
        xi[14]=-0.820001985974;
        xi[15]=0.820001985974;
        xi[16]=-0.886415527004;
        xi[17]=0.886415527004;
        xi[18]=-0.938274552003;
        xi[19]=0.938274552003;
        xi[20]=-0.974728555971;
        xi[21]=0.974728555971;
        xi[22]=-0.995187219997;
        xi[23]=0.995187219997;
    }
    else if(no_i==25)
    {
        wi[0]=0.123176053727;
        wi[1]=0.122242442990;
        wi[2]=0.122242442990;
        wi[3]=0.119455763536;
        wi[4]=0.119455763536;
        wi[5]=0.114858259146;
        wi[6]=0.114858259146;
        wi[7]=0.108519624474;
        wi[8]=0.108519624474;
        wi[9]=0.100535949067;
        wi[10]=0.100535949067;
        wi[11]=0.091028261983;
        wi[12]=0.091028261983;
        wi[13]=0.080140700335;
        wi[14]=0.080140700335;
        wi[15]=0.068038333812;
        wi[16]=0.068038333812;
        wi[17]=0.054904695976;
        wi[18]=0.05490469598;
        wi[19]=0.040939156701;
        wi[20]=0.04093915670;
        wi[21]=0.02635498662;
        wi[22]=0.02635498662;
        wi[23]=0.01139379850;
        wi[24]=0.01139379850;
        
        xi[0]=0;
        xi[1]=-0.122864692611;
        xi[2]=0.122864692611;
        xi[3]=-0.243866883721;
        xi[4]=0.243866883721;
        xi[5]=-0.361172305809;
        xi[6]=0.361172305809;
        xi[7]=-0.473002731446;
        xi[8]=0.473002731446;
        xi[9]=-0.577662930241;
        xi[10]=0.577662930241;
        xi[11]=-0.673566368473;
        xi[12]=0.673566368473;
        xi[13]=-0.759259263037;
        xi[14]=0.759259263037;
        xi[15]=-0.833442628761;
        xi[16]=0.833442628761;
        xi[17]=-0.894991997878;
        xi[18]=0.894991997878;
        xi[19]=-0.942974571229;
        xi[20]=0.942974571229;
        xi[21]=-0.976663921460;
        xi[22]=0.976663921460;
        xi[23]=-0.995556969790;
        xi[24]=0.995556969790;
    }
    else if(no_i==26)
    {
        wi[0]=0.118321415279;
        wi[1]=0.118321415279;
        wi[2]=0.116660443485;
        wi[3]=0.116660443485;
        wi[4]=0.113361816546;
        wi[5]=0.113361816546;
        wi[6]=0.108471840529;
        wi[7]=0.108471840529;
        wi[8]=0.102059161094;
        wi[9]=0.102059161094;
        wi[10]=0.094213800356;
        wi[11]=0.094213800356;
        wi[12]=0.085045894313;
        wi[13]=0.085045894313;
        wi[14]=0.074684149766;
        wi[15]=0.074684149766;
        wi[16]=0.063274046330;
        wi[17]=0.06327404633;
        wi[18]=0.050975825297;
        wi[19]=0.05097582530;
        wi[20]=0.037962383294;
        wi[21]=0.03796238329;
        wi[22]=0.02441785109;
        wi[23]=0.02441785109;
        wi[24]=0.01055137262;
        wi[25]=0.01055137262;
        
        xi[0]=-0.0592300934293;
        xi[1]=0.0592300934293;
        xi[2]=-0.176858820357;
        xi[3]=0.176858820357;
        xi[4]=-0.292004839486;
        xi[5]=0.292004839486;
        xi[6]=-0.403051755123;
        xi[7]=0.403051755123;
        xi[8]=-0.508440714825;
        xi[9]=0.508440714825;
        xi[10]=-0.606692293018;
        xi[11]=0.606692293018;
        xi[12]=-0.696427260420;
        xi[13]=0.696427260420;
        xi[14]=-0.776385948821;
        xi[15]=0.776385948821;
        xi[16]=-0.845445942788;
        xi[17]=0.845445942788;
        xi[18]=-0.902637861984;
        xi[19]=0.902637861984;
        xi[20]=-0.947159066662;
        xi[21]=0.947159066662;
        xi[22]=-0.978385445956;
        xi[23]=0.978385445956;
        xi[24]=-0.995885701146;
        xi[25]=0.995885701146;
    }
    else if(no_i==27)
    {
        wi[0]=0.114220867379;
        wi[1]=0.113476346109;
        wi[2]=0.113476346109;
        wi[3]=0.111252488357;
        wi[4]=0.111252488357;
        wi[5]=0.107578285789;
        wi[6]=0.107578285789;
        wi[7]=0.102501637818;
        wi[8]=0.102501637818;
        wi[9]=0.096088727370;
        wi[10]=0.096088727370;
        wi[11]=0.088423158544;
        wi[12]=0.088423158544;
        wi[13]=0.079604867773;
        wi[14]=0.079604867773;
        wi[15]=0.069748823766;
        wi[16]=0.069748823766;
        wi[17]=0.058983536860;
        wi[18]=0.05898353686;
        wi[19]=0.047449412521;
        wi[20]=0.04744941252;
        wi[21]=0.035297053757;
        wi[22]=0.03529705376;
        wi[23]=0.02268623160;
        wi[24]=0.02268623160;
        wi[25]=0.00979899605;
        wi[26]=0.00979899605;
        
        xi[0]=0;
        xi[1]=-0.113972585610;
        xi[2]=0.113972585610;
        xi[3]=-0.226459365440;
        xi[4]=0.226459365440;
        xi[5]=-0.335993903639;
        xi[6]=0.335993903639;
        xi[7]=-0.441148251750;
        xi[8]=0.441148251750;
        xi[9]=-0.540551564579;
        xi[10]=0.540551564579;
        xi[11]=-0.632907971946;
        xi[12]=0.632907971946;
        xi[13]=-0.717013473739;
        xi[14]=0.717013473739;
        xi[15]=-0.791771639071;
        xi[16]=0.791771639071;
        xi[17]=-0.856207908018;
        xi[18]=0.856207908018;
        xi[19]=-0.909482320677;
        xi[20]=0.909482320677;
        xi[21]=-0.950900557815;
        xi[22]=0.950900557815;
        xi[23]=-0.979923475962;
        xi[24]=0.979923475962;
        xi[25]=-0.996179262889;
        xi[26]=0.996179262889;
    }
    else if(no_i==28)
    {
        wi[0]=0.110047013016;
        wi[1]=0.110047013016;
        wi[2]=0.108711192258;
        wi[3]=0.108711192258;
        wi[4]=0.106055765923;
        wi[5]=0.106055765923;
        wi[6]=0.102112967578;
        wi[7]=0.102112967578;
        wi[8]=0.096930657998;
        wi[9]=0.096930657998;
        wi[10]=0.090571744393;
        wi[11]=0.090571744393;
        wi[12]=0.083113417229;
        wi[13]=0.083113417229;
        wi[14]=0.074646214235;
        wi[15]=0.074646214235;
        wi[16]=0.065272923967;
        wi[17]=0.065272923967;
        wi[18]=0.055107345676;
        wi[19]=0.05510734568;
        wi[20]=0.044272934759;
        wi[21]=0.04427293476;
        wi[22]=0.032901427782;
        wi[23]=0.03290142778;
        wi[24]=0.02113211259;
        wi[25]=0.02113211259;
        wi[26]=0.00912428259;
        wi[27]=0.00912428259;
        
        xi[0]=-0.0550792898840;
        xi[1]=0.0550792898840;
        xi[2]=-0.164569282133;
        xi[3]=0.164569282133;
        xi[4]=-0.272061627635;
        xi[5]=0.272061627635;
        xi[6]=-0.376251516089;
        xi[7]=0.376251516089;
        xi[8]=-0.475874224955;
        xi[9]=0.475874224955;
        xi[10]=-0.569720471811;
        xi[11]=0.569720471811;
        xi[12]=-0.656651094039;
        xi[13]=0.656651094039;
        xi[14]=-0.735610878014;
        xi[15]=0.735610878014;
        xi[16]=-0.805641370917;
        xi[17]=0.805641370917;
        xi[18]=-0.865892522574;
        xi[19]=0.865892522574;
        xi[20]=-0.915633026392;
        xi[21]=0.915633026392;
        xi[22]=-0.954259280629;
        xi[23]=0.954259280629;
        xi[24]=-0.981303165371;
        xi[25]=0.981303165371;
        xi[26]=-0.996442497574;
        xi[27]=0.996442497574;
    }
    else if(no_i==29)
    {
        wi[0]=0.106479381718;
        wi[1]=0.105876155097;
        wi[2]=0.105876155097;
        wi[3]=0.104073310078;
        wi[4]=0.104073310078;
        wi[5]=0.101091273760;
        wi[6]=0.101091273760;
        wi[7]=0.096963834094;
        wi[8]=0.096963834094;
        wi[9]=0.091737757139;
        wi[10]=0.091737757139;
        wi[11]=0.085472257366;
        wi[12]=0.085472257366;
        wi[13]=0.078238327136;
        wi[14]=0.078238327136;
        wi[15]=0.070117933255;
        wi[16]=0.070117933255;
        wi[17]=0.061203090657;
        wi[18]=0.061203090657;
        wi[19]=0.051594826902;
        wi[20]=0.05159482690;
        wi[21]=0.041402062519;
        wi[22]=0.04140206252;
        wi[23]=0.030740492202;
        wi[24]=0.03074049220;
        wi[25]=0.01973208506;
        wi[26]=0.01973208506;
        wi[27]=0.00851690388;
        wi[28]=0.00851690388;
        
        xi[0]=0;
        xi[1]=-0.106278230133;
        xi[2]=0.106278230133;
        xi[3]=-0.211352286166;
        xi[4]=0.211352286166;
        xi[5]=-0.314031637868;
        xi[6]=0.314031637868;
        xi[7]=-0.413152888174;
        xi[8]=0.413152888174;
        xi[9]=-0.507592955124;
        xi[10]=0.507592955124;
        xi[11]=-0.596281797138;
        xi[12]=0.596281797138;
        xi[13]=-0.678214537603;
        xi[14]=0.678214537603;
        xi[15]=-0.752462851734;
        xi[16]=0.752462851734;
        xi[17]=-0.818185487615;
        xi[18]=0.818185487615;
        xi[19]=-0.874637804920;
        xi[20]=0.874637804920;
        xi[21]=-0.921180232953;
        xi[22]=0.921180232953;
        xi[23]=-0.957285595778;
        xi[24]=0.957285595778;
        xi[25]=-0.982545505261;
        xi[26]=0.982545505261;
        xi[27]=-0.996679442261;
        xi[28]=0.996679442261;
    }
    else if(no_i==30)
    {
        wi[0]=0.102852652894;
        wi[1]=0.102852652894;
        wi[2]=0.101762389748;
        wi[3]=0.101762389748;
        wi[4]=0.099593420587;
        wi[5]=0.099593420587;
        wi[6]=0.096368737175;
        wi[7]=0.096368737175;
        wi[8]=0.092122522238;
        wi[9]=0.092122522238;
        wi[10]=0.086899787201;
        wi[11]=0.086899787201;
        wi[12]=0.080755895229;
        wi[13]=0.080755895229;
        wi[14]=0.073755974738;
        wi[15]=0.073755974738;
        wi[16]=0.065974229882;
        wi[17]=0.065974229882;
        wi[18]=0.057493156218;
        wi[19]=0.057493156218;
        wi[20]=0.048402672831;
        wi[21]=0.04840267283;
        wi[22]=0.038799192570;
        wi[23]=0.03879919257;
        wi[24]=0.028784707883;
        wi[25]=0.02878470788;
        wi[26]=0.01846646831;
        wi[27]=0.01846646831;
        wi[28]=0.00796819250;
        wi[29]=0.00796819250;
        
        xi[0]=-0.0514718425553;
        xi[1]=0.0514718425553;
        xi[2]=-0.153869913609;
        xi[3]=0.153869913609;
        xi[4]=-0.254636926168;
        xi[5]=0.254636926168;
        xi[6]=-0.352704725531;
        xi[7]=0.352704725531;
        xi[8]=-0.447033769538;
        xi[9]=0.447033769538;
        xi[10]=-0.536624148142;
        xi[11]=0.536624148142;
        xi[12]=-0.620526182989;
        xi[13]=0.620526182989;
        xi[14]=-0.697850494793;
        xi[15]=0.697850494793;
        xi[16]=-0.767777432105;
        xi[17]=0.767777432105;
        xi[18]=-0.829565762383;
        xi[19]=0.829565762383;
        xi[20]=-0.882560535792;
        xi[21]=0.882560535792;
        xi[22]=-0.926200047429;
        xi[23]=0.926200047429;
        xi[24]=-0.960021864968;
        xi[25]=0.960021864968;
        xi[26]=-0.983668123280;
        xi[27]=0.983668123280;
        xi[28]=-0.996893484075;
        xi[29]=0.996893484075;
    }
    else if(no_i==31)
    {
        wi[0]=0.0997205447934;
        wi[1]=0.099225011227;
        wi[2]=0.099225011227;
        wi[3]=0.097743335386;
        wi[4]=0.097743335386;
        wi[5]=0.095290242912;
        wi[6]=0.095290242912;
        wi[7]=0.091890113894;
        wi[8]=0.091890113894;
        wi[9]=0.087576740608;
        wi[10]=0.087576740608;
        wi[11]=0.082392991762;
        wi[12]=0.082392991762;
        wi[13]=0.076390386599;
        wi[14]=0.076390386599;
        wi[15]=0.069628583235;
        wi[16]=0.069628583235;
        wi[17]=0.062174786561;
        wi[18]=0.062174786561;
        wi[19]=0.054103082425;
        wi[20]=0.054103082425;
        wi[21]=0.045493707527;
        wi[22]=0.04549370753;
        wi[23]=0.036432273912;
        wi[24]=0.03643227391;
        wi[25]=0.027009019185;
        wi[26]=0.02700901918;
        wi[27]=0.01731862079;
        wi[28]=0.01731862079;
        wi[29]=0.00747083158;
        wi[30]=0.00747083158;
        
        xi[0]=0;
        xi[1]=-0.0995553121523;
        xi[2]=0.0995553121523;
        xi[3]=-0.198121199336;
        xi[4]=0.198121199336;
        xi[5]=-0.294718069982;
        xi[6]=0.294718069982;
        xi[7]=-0.388385901608;
        xi[8]=0.388385901608;
        xi[9]=-0.478193782045;
        xi[10]=0.478193782045;
        xi[11]=-0.563249161407;
        xi[12]=0.563249161407;
        xi[13]=-0.642706722924;
        xi[14]=0.642706722924;
        xi[15]=-0.715776784587;
        xi[16]=0.715776784587;
        xi[17]=-0.781733148417;
        xi[18]=0.781733148417;
        xi[19]=-0.839920320146;
        xi[20]=0.839920320146;
        xi[21]=-0.889760029948;
        xi[22]=0.889760029948;
        xi[23]=-0.930756997897;
        xi[24]=0.930756997897;
        xi[25]=-0.962503925093;
        xi[26]=0.962503925093;
        xi[27]=-0.984685909665;
        xi[28]=0.984685909665;
        xi[29]=-0.997087481819;
        xi[30]=0.997087481819;
    }
    else if(no_i==32)
    {
        wi[0]=0.096540088515;
        wi[1]=0.096540088515;
        wi[2]=0.095638720079;
        wi[3]=0.095638720079;
        wi[4]=0.093844399081;
        wi[5]=0.093844399081;
        wi[6]=0.091173878696;
        wi[7]=0.091173878696;
        wi[8]=0.087652093004;
        wi[9]=0.087652093004;
        wi[10]=0.083311924227;
        wi[11]=0.083311924227;
        wi[12]=0.078193895787;
        wi[13]=0.078193895787;
        wi[14]=0.072345794109;
        wi[15]=0.072345794109;
        wi[16]=0.065822222776;
        wi[17]=0.065822222776;
        wi[18]=0.058684093479;
        wi[19]=0.058684093479;
        wi[20]=0.050998059262;
        wi[21]=0.050998059262;
        wi[22]=0.042835898022;
        wi[23]=0.04283589802;
        wi[24]=0.034273862913;
        wi[25]=0.03427386291;
        wi[26]=0.025392065309;
        wi[27]=0.02539206531;
        wi[28]=0.01627439473;
        wi[29]=0.01627439473;
        wi[30]=0.00701861001;
        wi[31]=0.00701861001;
        
        xi[0]=-0.0483076656877;
        xi[1]=0.0483076656877;
        xi[2]=-0.144471961583;
        xi[3]=0.144471961583;
        xi[4]=-0.239287362252;
        xi[5]=0.239287362252;
        xi[6]=-0.331868602282;
        xi[7]=0.331868602282;
        xi[8]=-0.421351276131;
        xi[9]=0.421351276131;
        xi[10]=-0.506899908932;
        xi[11]=0.506899908932;
        xi[12]=-0.587715757241;
        xi[13]=0.587715757241;
        xi[14]=-0.663044266930;
        xi[15]=0.663044266930;
        xi[16]=-0.732182118740;
        xi[17]=0.732182118740;
        xi[18]=-0.794483795968;
        xi[19]=0.794483795968;
        xi[20]=-0.849367613733;
        xi[21]=0.849367613733;
        xi[22]=-0.896321155766;
        xi[23]=0.896321155766;
        xi[24]=-0.934906075938;
        xi[25]=0.934906075938;
        xi[26]=-0.964762255588;
        xi[27]=0.964762255588;
        xi[28]=-0.985611511545;
        xi[29]=0.985611511545;
        xi[30]=-0.997263861849;
        xi[31]=0.997263861849;
    }
    else if(no_i==33)
    {
        wi[0]=0.0937684461602;
        wi[1]=0.093356426066;
        wi[2]=0.093356426066;
        wi[3]=0.092123986643;
        wi[4]=0.092123986643;
        wi[5]=0.090081958661;
        wi[6]=0.090081958661;
        wi[7]=0.087248287619;
        wi[8]=0.087248287619;
        wi[9]=0.083647876067;
        wi[10]=0.083647876067;
        wi[11]=0.079312364795;
        wi[12]=0.079312364795;
        wi[13]=0.074279854844;
        wi[14]=0.074279854844;
        wi[15]=0.068594572819;
        wi[16]=0.068594572819;
        wi[17]=0.062306482530;
        wi[18]=0.062306482530;
        wi[19]=0.055470846632;
        wi[20]=0.055470846632;
        wi[21]=0.048147742819;
        wi[22]=0.048147742819;
        wi[23]=0.040401541332;
        wi[24]=0.04040154133;
        wi[25]=0.032300358632;
        wi[26]=0.03230035863;
        wi[27]=0.023915548102;
        wi[28]=0.02391554810;
        wi[29]=0.01532170151;
        wi[30]=0.01532170151;
        wi[31]=0.00660622785;
        wi[32]=0.00660622785;
        
        xi[0]=0;
        xi[1]=-0.0936310658547;
        xi[2]=0.0936310658547;
        xi[3]=-0.186439298828;
        xi[4]=0.186439298828;
        xi[5]=-0.277609097152;
        xi[6]=0.277609097152;
        xi[7]=-0.366339257748;
        xi[8]=0.366339257748;
        xi[9]=-0.451850017272;
        xi[10]=0.451850017272;
        xi[11]=-0.533389904786;
        xi[12]=0.533389904786;
        xi[13]=-0.610242345836;
        xi[14]=0.610242345836;
        xi[15]=-0.681731959970;
        xi[16]=0.681731959970;
        xi[17]=-0.747230496450;
        xi[18]=0.747230496450;
        xi[19]=-0.806162356274;
        xi[20]=0.806162356274;
        xi[21]=-0.858009652677;
        xi[22]=0.858009652677;
        xi[23]=-0.902316767743;
        xi[24]=0.902316767743;
        xi[25]=-0.938694372611;
        xi[26]=0.938694372611;
        xi[27]=-0.966822909690;
        xi[28]=0.966822909690;
        xi[29]=-0.986455726231;
        xi[30]=0.986455726231;
        xi[31]=-0.997424694246;
        xi[32]=0.997424694246;
    }
    else if(no_i==34)
    {
        wi[0]=0.090956740330;
        wi[1]=0.090956740330;
        wi[2]=0.090203044371;
        wi[3]=0.090203044371;
        wi[4]=0.088701897836;
        wi[5]=0.088701897836;
        wi[6]=0.086465739747;
        wi[7]=0.086465739747;
        wi[8]=0.083513099700;
        wi[9]=0.083513099700;
        wi[10]=0.079868444340;
        wi[11]=0.079868444340;
        wi[12]=0.075561974660;
        wi[13]=0.075561974660;
        wi[14]=0.070629375814;
        wi[15]=0.070629375814;
        wi[16]=0.065111521554;
        wi[17]=0.065111521554;
        wi[18]=0.059054135828;
        wi[19]=0.059054135828;
        wi[20]=0.052507414573;
        wi[21]=0.052507414573;
        wi[22]=0.045525611523;
        wi[23]=0.045525611523;
        wi[24]=0.038166593796;
        wi[25]=0.03816659380;
        wi[26]=0.030491380638;
        wi[27]=0.03049138064;
        wi[28]=0.022563721985;
        wi[29]=0.02256372199;
        wi[30]=0.01445016275;
        wi[31]=0.01445016275;
        wi[32]=0.00622914056;
        wi[33]=0.00622914056;
        
        xi[0]=-0.0455098219531;
        xi[1]=0.0455098219531;
        xi[2]=-0.136152357259;
        xi[3]=0.136152357259;
        xi[4]=-0.225666691616;
        xi[5]=0.225666691616;
        xi[6]=-0.313311081339;
        xi[7]=0.313311081339;
        xi[8]=-0.398359277759;
        xi[9]=0.398359277759;
        xi[10]=-0.480106545190;
        xi[11]=0.480106545190;
        xi[12]=-0.557875500670;
        xi[13]=0.557875500670;
        xi[14]=-0.631021727081;
        xi[15]=0.631021727081;
        xi[16]=-0.698939113216;
        xi[17]=0.698939113216;
        xi[18]=-0.761064876630;
        xi[19]=0.761064876630;
        xi[20]=-0.816884227901;
        xi[21]=0.816884227901;
        xi[22]=-0.865934638335;
        xi[23]=0.865934638335;
        xi[24]=-0.907809677718;
        xi[25]=0.907809677718;
        xi[26]=-0.942162397405;
        xi[27]=0.942162397405;
        xi[28]=-0.968708262533;
        xi[29]=0.968708262533;
        xi[30]=-0.987227816406;
        xi[31]=0.987227816406;
        xi[32]=-0.997571753791;
        xi[33]=0.997571753791;
    }
    else if(no_i==35)
    {
        wi[0]=0.0884867949071;
        wi[1]=0.088140530430;
        wi[2]=0.088140530430;
        wi[3]=0.087104446997;
        wi[4]=0.087104446997;
        wi[5]=0.085386653392;
        wi[6]=0.085386653392;
        wi[7]=0.083000593729;
        wi[8]=0.083000593729;
        wi[9]=0.079964942242;
        wi[10]=0.079964942242;
        wi[11]=0.076303457155;
        wi[12]=0.076303457155;
        wi[13]=0.072044794773;
        wi[14]=0.072044794773;
        wi[15]=0.067222285269;
        wi[16]=0.067222285269;
        wi[17]=0.061873671966;
        wi[18]=0.061873671966;
        wi[19]=0.056040816212;
        wi[20]=0.056040816212;
        wi[21]=0.049769370401;
        wi[22]=0.049769370401;
        wi[23]=0.043108422326;
        wi[24]=0.043108422326;
        wi[25]=0.036110115863;
        wi[26]=0.03611011586;
        wi[27]=0.028829260109;
        wi[28]=0.02882926011;
        wi[29]=0.021322979911;
        wi[30]=0.02132297991;
        wi[31]=0.01365082835;
        wi[32]=0.01365082835;
        wi[33]=0.00588343342;
        wi[34]=0.00588343342;
        
        xi[0]=0;
        xi[1]=-0.0883713432757;
        xi[2]=0.0883713432757;
        xi[3]=-0.176051061166;
        xi[4]=0.176051061166;
        xi[5]=-0.262352941209;
        xi[6]=0.262352941209;
        xi[7]=-0.346601554431;
        xi[8]=0.346601554431;
        xi[9]=-0.428137541518;
        xi[10]=0.428137541518;
        xi[11]=-0.506322773241;
        xi[12]=0.506322773241;
        xi[13]=-0.580545344750;
        xi[14]=0.580545344750;
        xi[15]=-0.650224364666;
        xi[16]=0.650224364666;
        xi[17]=-0.714814501557;
        xi[18]=0.714814501557;
        xi[19]=-0.773810252287;
        xi[20]=0.773810252287;
        xi[21]=-0.826749899092;
        xi[22]=0.826749899092;
        xi[23]=-0.873219125025;
        xi[24]=0.873219125025;
        xi[25]=-0.912854261359;
        xi[26]=0.912854261359;
        xi[27]=-0.945345148208;
        xi[28]=0.945345148208;
        xi[29]=-0.970437616039;
        xi[30]=0.970437616039;
        xi[31]=-0.987935764444;
        xi[32]=0.987935764444;
        xi[33]=-0.997706569100;
        xi[34]=0.997706569100;
    }
    else if(no_i==36)
    {
        wi[0]=0.085983275670;
        wi[1]=0.085983275670;
        wi[2]=0.085346685739;
        wi[3]=0.085346685739;
        wi[4]=0.084078218980;
        wi[5]=0.084078218980;
        wi[6]=0.082187266704;
        wi[7]=0.082187266704;
        wi[8]=0.079687828912;
        wi[9]=0.079687828912;
        wi[10]=0.076598410646;
        wi[11]=0.076598410646;
        wi[12]=0.072941885006;
        wi[13]=0.072941885006;
        wi[14]=0.068745323836;
        wi[15]=0.068745323836;
        wi[16]=0.064039797355;
        wi[17]=0.064039797355;
        wi[18]=0.058860144245;
        wi[19]=0.058860144245;
        wi[20]=0.053244713978;
        wi[21]=0.053244713978;
        wi[22]=0.047235083490;
        wi[23]=0.047235083490;
        wi[24]=0.040875750924;
        wi[25]=0.040875750924;
        wi[26]=0.034213810770;
        wi[27]=0.03421381077;
        wi[28]=0.027298621499;
        wi[29]=0.02729862150;
        wi[30]=0.020181515298;
        wi[31]=0.02018151530;
        wi[32]=0.01291594728;
        wi[33]=0.01291594728;
        wi[34]=0.00556571966;
        wi[35]=0.00556571966;
        
        xi[0]=-0.0430181984737;
        xi[1]=0.0430181984737;
        xi[2]=-0.128736103809;
        xi[3]=0.128736103809;
        xi[4]=-0.213500892317;
        xi[5]=0.213500892317;
        xi[6]=-0.296684995344;
        xi[7]=0.296684995344;
        xi[8]=-0.377672547120;
        xi[9]=0.377672547120;
        xi[10]=-0.455863944433;
        xi[11]=0.455863944433;
        xi[12]=-0.530680285926;
        xi[13]=0.530680285926;
        xi[14]=-0.601567658136;
        xi[15]=0.601567658136;
        xi[16]=-0.668001236586;
        xi[17]=0.668001236586;
        xi[18]=-0.729489171594;
        xi[19]=0.729489171594;
        xi[20]=-0.785576230132;
        xi[21]=0.785576230132;
        xi[22]=-0.835847166992;
        xi[23]=0.835847166992;
        xi[24]=-0.879929800890;
        xi[25]=0.879929800890;
        xi[26]=-0.917497774516;
        xi[27]=0.917497774516;
        xi[28]=-0.948272984400;
        xi[29]=0.948272984400;
        xi[30]=-0.972027691050;
        xi[31]=0.972027691050;
        xi[32]=-0.988586478902;
        xi[33]=0.988586478902;
        xi[34]=-0.997830462484;
        xi[35]=0.997830462484;
    }
    else if(no_i==37)
    {
        wi[0]=0.0837683609931;
        wi[1]=0.083474573626;
        wi[2]=0.083474573626;
        wi[3]=0.082595272236;
        wi[4]=0.082595272236;
        wi[5]=0.081136624508;
        wi[6]=0.081136624508;
        wi[7]=0.079108861838;
        wi[8]=0.079108861838;
        wi[9]=0.076526207571;
        wi[10]=0.076526207571;
        wi[11]=0.073406777248;
        wi[12]=0.073406777248;
        wi[13]=0.069772451556;
        wi[14]=0.069772451556;
        wi[15]=0.065648722873;
        wi[16]=0.065648722873;
        wi[17]=0.061064516523;
        wi[18]=0.061064516523;
        wi[19]=0.056051987998;
        wi[20]=0.056051987998;
        wi[21]=0.050646297655;
        wi[22]=0.050646297655;
        wi[23]=0.044885364662;
        wi[24]=0.044885364662;
        wi[25]=0.038809602502;
        wi[26]=0.038809602502;
        wi[27]=0.032461639848;
        wi[28]=0.03246163985;
        wi[29]=0.025886036991;
        wi[30]=0.02588603699;
        wi[31]=0.019129044489;
        wi[32]=0.01912904449;
        wi[33]=0.01223878010;
        wi[34]=0.01223878010;
        wi[35]=0.00527305728;
        wi[36]=0.00527305728;
        
        xi[0]=0;
        xi[1]=-0.0836704089548;
        xi[2]=0.0836704089548;
        xi[3]=-0.166753930240;
        xi[4]=0.166753930240;
        xi[5]=-0.248667792791;
        xi[6]=0.248667792791;
        xi[7]=-0.328837429884;
        xi[8]=0.328837429884;
        xi[9]=-0.406700509318;
        xi[10]=0.406700509318;
        xi[11]=-0.481710877803;
        xi[12]=0.481710877803;
        xi[13]=-0.553342391862;
        xi[14]=0.553342391862;
        xi[15]=-0.621092608409;
        xi[16]=0.621092608409;
        xi[17]=-0.684486309131;
        xi[18]=0.684486309131;
        xi[19]=-0.743078833982;
        xi[20]=0.743078833982;
        xi[21]=-0.796459200510;
        xi[22]=0.796459200510;
        xi[23]=-0.844252987341;
        xi[24]=0.844252987341;
        xi[25]=-0.886124962155;
        xi[26]=0.886124962155;
        xi[27]=-0.921781437412;
        xi[28]=0.921781437412;
        xi[29]=-0.950972343262;
        xi[30]=0.950972343262;
        xi[31]=-0.973493030056;
        xi[32]=0.973493030056;
        xi[33]=-0.989185963214;
        xi[34]=0.989185963214;
        xi[35]=-0.997944582478;
        xi[36]=0.997944582478;
    }
    else if(no_i==38)
    {
        wi[0]=0.081525029280;
        wi[1]=0.081525029280;
        wi[2]=0.080982493771;
        wi[3]=0.080982493771;
        wi[4]=0.079901033244;
        wi[5]=0.079901033244;
        wi[6]=0.078287844658;
        wi[7]=0.078287844658;
        wi[8]=0.076153663548;
        wi[9]=0.076153663548;
        wi[10]=0.073512692585;
        wi[11]=0.073512692585;
        wi[12]=0.070382507067;
        wi[13]=0.070382507067;
        wi[14]=0.066783937979;
        wi[15]=0.066783937979;
        wi[16]=0.062740933392;
        wi[17]=0.062740933392;
        wi[18]=0.058280399147;
        wi[19]=0.058280399147;
        wi[20]=0.053432019910;
        wi[21]=0.053432019910;
        wi[22]=0.048228061861;
        wi[23]=0.048228061861;
        wi[24]=0.042703158505;
        wi[25]=0.042703158505;
        wi[26]=0.036894081594;
        wi[27]=0.036894081594;
        wi[28]=0.030839500545;
        wi[29]=0.03083950055;
        wi[30]=0.024579739738;
        wi[31]=0.02457973974;
        wi[32]=0.018156577710;
        wi[33]=0.01815657771;
        wi[34]=0.01161344472;
        wi[35]=0.01161344472;
        wi[36]=0.00500288075;
        wi[37]=0.00500288075;
        
        xi[0]=-0.0407851479046;
        xi[1]=0.0407851479046;
        xi[2]=-0.122084025338;
        xi[3]=0.122084025338;
        xi[4]=-0.202570453892;
        xi[5]=0.202570453892;
        xi[6]=-0.281708809790;
        xi[7]=0.281708809790;
        xi[8]=-0.358972440479;
        xi[9]=0.358972440479;
        xi[10]=-0.433847169432;
        xi[11]=0.433847169432;
        xi[12]=-0.505834717928;
        xi[13]=0.505834717928;
        xi[14]=-0.574456021048;
        xi[15]=0.574456021048;
        xi[16]=-0.639254415830;
        xi[17]=0.639254415830;
        xi[18]=-0.699798680379;
        xi[19]=0.699798680379;
        xi[20]=-0.755685903754;
        xi[21]=0.755685903754;
        xi[22]=-0.806544167605;
        xi[23]=0.806544167605;
        xi[24]=-0.852035021932;
        xi[25]=0.852035021932;
        xi[26]=-0.891855739005;
        xi[27]=0.891855739005;
        xi[28]=-0.925741332049;
        xi[29]=0.925741332049;
        xi[30]=-0.953466330934;
        xi[31]=0.953466330934;
        xi[32]=-0.974846328590;
        xi[33]=0.974846328590;
        xi[34]=-0.989739454266;
        xi[35]=0.989739454266;
        xi[36]=-0.998049930536;
        xi[37]=0.998049930536;
    }
    else if(no_i==39)
    {
        wi[0]=0.0795276221394;
        wi[1]=0.079276222568;
        wi[2]=0.079276222568;
        wi[3]=0.078523613287;
        wi[4]=0.078523613287;
        wi[5]=0.077274552545;
        wi[6]=0.077274552545;
        wi[7]=0.075536937323;
        wi[8]=0.075536937323;
        wi[9]=0.073321753414;
        wi[10]=0.073321753414;
        wi[11]=0.070643005971;
        wi[12]=0.070643005971;
        wi[13]=0.067517630966;
        wi[14]=0.067517630966;
        wi[15]=0.063965388139;
        wi[16]=0.063965388139;
        wi[17]=0.060008736089;
        wi[18]=0.060008736089;
        wi[19]=0.055672690341;
        wi[20]=0.055672690341;
        wi[21]=0.050984665292;
        wi[22]=0.050984665292;
        wi[23]=0.045974301109;
        wi[24]=0.045974301109;
        wi[25]=0.040673276848;
        wi[26]=0.040673276848;
        wi[27]=0.035115111498;
        wi[28]=0.035115111498;
        wi[29]=0.029334955984;
        wi[30]=0.02933495598;
        wi[31]=0.023369384832;
        wi[32]=0.02336938483;
        wi[33]=0.017256229094;
        wi[34]=0.01725622909;
        wi[35]=0.01103478894;
        wi[36]=0.01103478894;
        wi[37]=0.00475294469;
        wi[38]=0.00475294469;
        
        xi[0]=0;
        xi[1]=-0.0794438046088;
        xi[2]=0.0794438046088;
        xi[3]=-0.158385339998;
        xi[4]=0.158385339998;
        xi[5]=-0.236325512462;
        xi[6]=0.236325512462;
        xi[7]=-0.312771559248;
        xi[8]=0.312771559248;
        xi[9]=-0.387240163972;
        xi[10]=0.387240163972;
        xi[11]=-0.459260512309;
        xi[12]=0.459260512309;
        xi[13]=-0.528377268660;
        xi[14]=0.528377268660;
        xi[15]=-0.594153454957;
        xi[16]=0.594153454957;
        xi[17]=-0.656173213432;
        xi[18]=0.656173213432;
        xi[19]=-0.714044435895;
        xi[20]=0.714044435895;
        xi[21]=-0.767401242931;
        xi[22]=0.767401242931;
        xi[23]=-0.815906297430;
        xi[24]=0.815906297430;
        xi[25]=-0.859252938000;
        xi[26]=0.859252938000;
        xi[27]=-0.897167119293;
        xi[28]=0.897167119293;
        xi[29]=-0.929409148487;
        xi[30]=0.929409148487;
        xi[31]=-0.955775212325;
        xi[32]=0.955775212325;
        xi[33]=-0.976098709333;
        xi[34]=0.976098709333;
        xi[35]=-0.990251536855;
        xi[36]=0.990251536855;
        xi[37]=-0.998147383066;
        xi[38]=0.998147383066;
    }
    else if(no_i==40)
    {
        wi[0]=0.077505947978;
        wi[1]=0.077505947978;
        wi[2]=0.077039818164;
        wi[3]=0.077039818164;
        wi[4]=0.076110361901;
        wi[5]=0.076110361901;
        wi[6]=0.074723169058;
        wi[7]=0.074723169058;
        wi[8]=0.072886582396;
        wi[9]=0.072886582396;
        wi[10]=0.070611647391;
        wi[11]=0.070611647391;
        wi[12]=0.067912045815;
        wi[13]=0.067912045815;
        wi[14]=0.064804013457;
        wi[15]=0.064804013457;
        wi[16]=0.061306242493;
        wi[17]=0.061306242493;
        wi[18]=0.057439769099;
        wi[19]=0.057439769099;
        wi[20]=0.053227846984;
        wi[21]=0.053227846984;
        wi[22]=0.048695807635;
        wi[23]=0.048695807635;
        wi[24]=0.043870908186;
        wi[25]=0.043870908186;
        wi[26]=0.038782167974;
        wi[27]=0.038782167974;
        wi[28]=0.033460195283;
        wi[29]=0.033460195283;
        wi[30]=0.027937006980;
        wi[31]=0.02793700698;
        wi[32]=0.022245849194;
        wi[33]=0.02224584919;
        wi[34]=0.016421058382;
        wi[35]=0.01642105838;
        wi[36]=0.01049828453;
        wi[37]=0.01049828453;
        wi[38]=0.00452127710;
        wi[39]=0.00452127710;
        
        xi[0]=-0.0387724175061;
        xi[1]=0.0387724175061;
        xi[2]=-0.116084070675;
        xi[3]=0.116084070675;
        xi[4]=-0.192697580701;
        xi[5]=0.192697580701;
        xi[6]=-0.268152185007;
        xi[7]=0.268152185007;
        xi[8]=-0.341994090826;
        xi[9]=0.341994090826;
        xi[10]=-0.413779204372;
        xi[11]=0.413779204372;
        xi[12]=-0.483075801686;
        xi[13]=0.483075801686;
        xi[14]=-0.549467125095;
        xi[15]=0.549467125095;
        xi[16]=-0.612553889668;
        xi[17]=0.612553889668;
        xi[18]=-0.671956684614;
        xi[19]=0.671956684614;
        xi[20]=-0.727318255190;
        xi[21]=0.727318255190;
        xi[22]=-0.778305651427;
        xi[23]=0.778305651427;
        xi[24]=-0.824612230833;
        xi[25]=0.824612230833;
        xi[26]=-0.865959503212;
        xi[27]=0.865959503212;
        xi[28]=-0.902098806969;
        xi[29]=0.902098806969;
        xi[30]=-0.932812808279;
        xi[31]=0.932812808279;
        xi[32]=-0.957916819214;
        xi[33]=0.957916819214;
        xi[34]=-0.977259949984;
        xi[35]=0.977259949984;
        xi[36]=-0.990726238699;
        xi[37]=0.990726238699;
        xi[38]=-0.998237709711;
        xi[39]=0.998237709711;
    }
    else if(no_i==41)
    {
        wi[0]=0.0756955356473;
        wi[1]=0.075478747093;
        wi[2]=0.075478747093;
        wi[3]=0.074829623176;
        wi[4]=0.074829623176;
        wi[5]=0.073751882027;
        wi[6]=0.073751882027;
        wi[7]=0.072251696861;
        wi[8]=0.072251696861;
        wi[9]=0.070337660621;
        wi[10]=0.070337660621;
        wi[11]=0.068020736761;
        wi[12]=0.068020736761;
        wi[13]=0.065314196454;
        wi[14]=0.065314196454;
        wi[15]=0.062233542581;
        wi[16]=0.062233542581;
        wi[17]=0.058796420950;
        wi[18]=0.058796420950;
        wi[19]=0.055022519243;
        wi[20]=0.055022519243;
        wi[21]=0.050933454295;
        wi[22]=0.050933454295;
        wi[23]=0.046552648369;
        wi[24]=0.046552648369;
        wi[25]=0.041905195196;
        wi[26]=0.041905195196;
        wi[27]=0.037017716704;
        wi[28]=0.037017716704;
        wi[29]=0.031918211732;
        wi[30]=0.031918211732;
        wi[31]=0.026635899207;
        wi[32]=0.02663589921;
        wi[33]=0.021201063369;
        wi[34]=0.02120106337;
        wi[35]=0.015644938408;
        wi[36]=0.01564493841;
        wi[37]=0.00999993877;
        wi[38]=0.00999993877;
        wi[39]=0.00430614036;
        wi[40]=0.00430614036;
        
        xi[0]=0;
        xi[1]=-0.0756232589892;
        xi[2]=0.0756232589892;
        xi[3]=-0.150813354864;
        xi[4]=0.150813354864;
        xi[5]=-0.225139605633;
        xi[6]=0.225139605633;
        xi[7]=-0.298176277342;
        xi[8]=0.298176277342;
        xi[9]=-0.369505022640;
        xi[10]=0.369505022640;
        xi[11]=-0.438717277051;
        xi[12]=0.438717277051;
        xi[13]=-0.505416599199;
        xi[14]=0.505416599199;
        xi[15]=-0.569220941610;
        xi[16]=0.569220941610;
        xi[17]=-0.629764839072;
        xi[18]=0.629764839072;
        xi[19]=-0.686701502035;
        xi[20]=0.686701502035;
        xi[21]=-0.739704803070;
        xi[22]=0.739704803070;
        xi[23]=-0.788471145047;
        xi[24]=0.788471145047;
        xi[25]=-0.832721200401;
        xi[26]=0.832721200401;
        xi[27]=-0.872201511692;
        xi[28]=0.872201511692;
        xi[29]=-0.906685944758;
        xi[30]=0.906685944758;
        xi[31]=-0.935976987498;
        xi[32]=0.935976987498;
        xi[33]=-0.959906891730;
        xi[34]=0.959906891730;
        xi[35]=-0.978338673561;
        xi[36]=0.978338673561;
        xi[37]=-0.991167109699;
        xi[38]=0.991167109699;
        xi[39]=-0.998321588575;
        xi[40]=0.998321588575;
    }
    else if(no_i==42)
    {
        wi[0]=0.073864234232;
        wi[1]=0.073864234232;
        wi[2]=0.073460813453;
        wi[3]=0.073460813453;
        wi[4]=0.072656175244;
        wi[5]=0.072656175244;
        wi[6]=0.071454714265;
        wi[7]=0.071454714265;
        wi[8]=0.069862992493;
        wi[9]=0.069862992493;
        wi[10]=0.067889703377;
        wi[11]=0.067889703377;
        wi[12]=0.065545624365;
        wi[13]=0.065545624365;
        wi[14]=0.062843558045;
        wi[15]=0.062843558045;
        wi[16]=0.059798262228;
        wi[17]=0.059798262228;
        wi[18]=0.056426369358;
        wi[19]=0.056426369358;
        wi[20]=0.052746295699;
        wi[21]=0.052746295699;
        wi[22]=0.048778140793;
        wi[23]=0.048778140793;
        wi[24]=0.044543577772;
        wi[25]=0.044543577772;
        wi[26]=0.040065735181;
        wi[27]=0.040065735181;
        wi[28]=0.035369071098;
        wi[29]=0.035369071098;
        wi[30]=0.030479240700;
        wi[31]=0.030479240700;
        wi[32]=0.025422959526;
        wi[33]=0.02542295953;
        wi[34]=0.020227869569;
        wi[35]=0.02022786957;
        wi[36]=0.014922443697;
        wi[37]=0.01492244370;
        wi[38]=0.00953622030;
        wi[39]=0.00953622030;
        wi[40]=0.00410599860;
        wi[41]=0.00410599860;
        
        xi[0]=-0.0369489431654;
        xi[1]=0.0369489431654;
        xi[2]=-0.110645027209;
        xi[3]=0.110645027209;
        xi[4]=-0.183736806565;
        xi[5]=0.183736806565;
        xi[6]=-0.255825079343;
        xi[7]=0.255825079343;
        xi[8]=-0.326516124465;
        xi[9]=0.326516124465;
        xi[10]=-0.395423852043;
        xi[11]=0.395423852043;
        xi[12]=-0.462171912070;
        xi[13]=0.462171912070;
        xi[14]=-0.526395749931;
        xi[15]=0.526395749931;
        xi[16]=-0.587744597485;
        xi[17]=0.587744597485;
        xi[18]=-0.645883388869;
        xi[19]=0.645883388869;
        xi[20]=-0.700494590556;
        xi[21]=0.700494590556;
        xi[22]=-0.751279935689;
        xi[23]=0.751279935689;
        xi[24]=-0.797962053255;
        xi[25]=0.797962053255;
        xi[26]=-0.840285983262;
        xi[27]=0.840285983262;
        xi[28]=-0.878020569812;
        xi[29]=0.878020569812;
        xi[30]=-0.910959724904;
        xi[31]=0.910959724904;
        xi[32]=-0.938923557355;
        xi[33]=0.938923557355;
        xi[34]=-0.961759365338;
        xi[35]=0.961759365338;
        xi[36]=-0.979342508064;
        xi[37]=0.979342508064;
        xi[38]=-0.991577288341;
        xi[39]=0.991577288341;
        xi[40]=-0.998399618990;
        xi[41]=0.998399618990;
    }
    else if(no_i==43)
    {
        wi[0]=0.0722157516938;
        wi[1]=0.072027501971;
        wi[2]=0.072027501971;
        wi[3]=0.071463734253;
        wi[4]=0.071463734253;
        wi[5]=0.070527387765;
        wi[6]=0.070527387765;
        wi[7]=0.069223344194;
        wi[8]=0.069223344194;
        wi[9]=0.067558402229;
        wi[10]=0.067558402229;
        wi[11]=0.065541242126;
        wi[12]=0.065541242126;
        wi[13]=0.063182380449;
        wi[14]=0.063182380449;
        wi[15]=0.060494115250;
        wi[16]=0.060494115250;
        wi[17]=0.057490461957;
        wi[18]=0.057490461957;
        wi[19]=0.054187080319;
        wi[20]=0.054187080319;
        wi[21]=0.050601192784;
        wi[22]=0.050601192784;
        wi[23]=0.046751494754;
        wi[24]=0.046751494754;
        wi[25]=0.042658057198;
        wi[26]=0.042658057198;
        wi[27]=0.038342222194;
        wi[28]=0.038342222194;
        wi[29]=0.033826492087;
        wi[30]=0.033826492087;
        wi[31]=0.029134413261;
        wi[32]=0.029134413261;
        wi[33]=0.024290456614;
        wi[34]=0.02429045661;
        wi[35]=0.019319901424;
        wi[36]=0.01931990142;
        wi[37]=0.014248756432;
        wi[38]=0.01424875643;
        wi[39]=0.00910399664;
        wi[40]=0.00910399664;
        wi[41]=0.00391949025;
        wi[42]=0.00391949025;
        
        xi[0]=0;
        xi[1]=-0.0721529908746;
        xi[2]=0.0721529908746;
        xi[3]=-0.143929809511;
        xi[4]=0.143929809511;
        xi[5]=-0.214956244861;
        xi[6]=0.214956244861;
        xi[7]=-0.284861998033;
        xi[8]=0.284861998033;
        xi[9]=-0.353282612864;
        xi[10]=0.353282612864;
        xi[11]=-0.419861376029;
        xi[12]=0.419861376029;
        xi[13]=-0.484251176786;
        xi[14]=0.484251176786;
        xi[15]=-0.546116316660;
        xi[16]=0.546116316660;
        xi[17]=-0.605134259640;
        xi[18]=0.605134259640;
        xi[19]=-0.660997313751;
        xi[20]=0.660997313751;
        xi[21]=-0.713414235269;
        xi[22]=0.713414235269;
        xi[23]=-0.762111747195;
        xi[24]=0.762111747195;
        xi[25]=-0.806835964137;
        xi[26]=0.806835964137;
        xi[27]=-0.847353716209;
        xi[28]=0.847353716209;
        xi[29]=-0.883453765219;
        xi[30]=0.883453765219;
        xi[31]=-0.914947907206;
        xi[32]=0.914947907206;
        xi[33]=-0.941671956848;
        xi[34]=0.941671956848;
        xi[35]=-0.963486613014;
        xi[36]=0.963486613014;
        xi[37]=-0.980278220980;
        xi[38]=0.980278220980;
        xi[39]=-0.991959557593;
        xi[40]=0.991959557593;
        xi[41]=-0.998472332243;
        xi[42]=0.998472332243;
    }
    else if(no_i==44)
    {
        wi[0]=0.070549157789;
        wi[1]=0.070549157789;
        wi[2]=0.070197685474;
        wi[3]=0.070197685474;
        wi[4]=0.069496491862;
        wi[5]=0.069496491862;
        wi[6]=0.068449070269;
        wi[7]=0.068449070269;
        wi[8]=0.067060638906;
        wi[9]=0.067060638906;
        wi[10]=0.065338114879;
        wi[11]=0.065338114879;
        wi[12]=0.063290079733;
        wi[13]=0.063290079733;
        wi[14]=0.060926736702;
        wi[15]=0.060926736702;
        wi[16]=0.058259859878;
        wi[17]=0.058259859878;
        wi[18]=0.055302735564;
        wi[19]=0.055302735564;
        wi[20]=0.052070096092;
        wi[21]=0.052070096092;
        wi[22]=0.048578046448;
        wi[23]=0.048578046448;
        wi[24]=0.044843984082;
        wi[25]=0.044843984082;
        wi[26]=0.040886512310;
        wi[27]=0.040886512310;
        wi[28]=0.036725347814;
        wi[29]=0.036725347814;
        wi[30]=0.032381222812;
        wi[31]=0.032381222812;
        wi[32]=0.027875782821;
        wi[33]=0.027875782821;
        wi[34]=0.023231481902;
        wi[35]=0.02323148190;
        wi[36]=0.018471481737;
        wi[37]=0.01847148174;
        wi[38]=0.013619586756;
        wi[39]=0.01361958676;
        wi[40]=0.00870048137;
        wi[41]=0.00870048137;
        wi[42]=0.00374540480;
        wi[43]=0.00374540480;
        
        xi[0]=-0.0352892369641;
        xi[1]=0.0352892369641;
        xi[2]=-0.105691901709;
        xi[3]=0.105691901709;
        xi[4]=-0.175568014776;
        xi[5]=0.175568014776;
        xi[6]=-0.244569456928;
        xi[7]=0.244569456928;
        xi[8]=-0.312352466503;
        xi[9]=0.312352466503;
        xi[10]=-0.378579352015;
        xi[11]=0.378579352015;
        xi[12]=-0.442920174525;
        xi[13]=0.442920174525;
        xi[14]=-0.505054391388;
        xi[15]=0.505054391388;
        xi[16]=-0.564672453185;
        xi[17]=0.564672453185;
        xi[18]=-0.621477345904;
        xi[19]=0.621477345904;
        xi[20]=-0.675186070666;
        xi[21]=0.675186070666;
        xi[22]=-0.725531053661;
        xi[23]=0.725531053661;
        xi[24]=-0.772261479249;
        xi[25]=0.772261479249;
        xi[26]=-0.815144539645;
        xi[27]=0.815144539645;
        xi[28]=-0.853966595005;
        xi[29]=0.853966595005;
        xi[30]=-0.888534238286;
        xi[31]=0.888534238286;
        xi[32]=-0.918675259984;
        xi[33]=0.918675259984;
        xi[34]=-0.944239509118;
        xi[35]=0.944239509118;
        xi[36]=-0.965099650422;
        xi[37]=0.965099650422;
        xi[38]=-0.981151833078;
        xi[39]=0.981151833078;
        xi[40]=-0.992316392139;
        xi[41]=0.992316392139;
        xi[42]=-0.998540200637;
        xi[43]=0.998540200637;
    }
    else if(no_i==45)
    {
        wi[0]=0.0690418248292;
        wi[1]=0.068877316978;
        wi[2]=0.068877316978;
        wi[3]=0.068384577379;
        wi[4]=0.068384577379;
        wi[5]=0.067565954164;
        wi[6]=0.067565954164;
        wi[7]=0.066425348450;
        wi[8]=0.066425348450;
        wi[9]=0.064968195751;
        wi[10]=0.064968195751;
        wi[11]=0.063201440074;
        wi[12]=0.063201440074;
        wi[13]=0.061133500831;
        wi[14]=0.061133500831;
        wi[15]=0.058774232719;
        wi[16]=0.058774232719;
        wi[17]=0.056134878760;
        wi[18]=0.056134878760;
        wi[19]=0.053228016731;
        wi[20]=0.053228016731;
        wi[21]=0.050067499238;
        wi[22]=0.050067499238;
        wi[23]=0.046668387718;
        wi[24]=0.046668387718;
        wi[25]=0.043046880709;
        wi[26]=0.043046880709;
        wi[27]=0.039220236729;
        wi[28]=0.039220236729;
        wi[29]=0.035206692202;
        wi[30]=0.035206692202;
        wi[31]=0.031025374935;
        wi[32]=0.031025374935;
        wi[33]=0.026696213968;
        wi[34]=0.026696213968;
        wi[35]=0.022239847551;
        wi[36]=0.02223984755;
        wi[37]=0.017677535258;
        wi[38]=0.01767753526;
        wi[39]=0.013031104992;
        wi[40]=0.01303110499;
        wi[41]=0.00832318930;
        wi[42]=0.00832318930;
        wi[43]=0.00358266316;
        wi[44]=0.00358266316;
        
        xi[0]=0;
        xi[1]=-0.0689869801631;
        xi[2]=0.0689869801631;
        xi[3]=-0.137645205983;
        xi[4]=0.137645205983;
        xi[5]=-0.205647489783;
        xi[6]=0.205647489783;
        xi[7]=-0.272669769752;
        xi[8]=0.272669769752;
        xi[9]=-0.338392654251;
        xi[10]=0.338392654251;
        xi[11]=-0.402502943859;
        xi[12]=0.402502943859;
        xi[13]=-0.464695123920;
        xi[14]=0.464695123920;
        xi[15]=-0.524672820463;
        xi[16]=0.524672820463;
        xi[17]=-0.582150212569;
        xi[18]=0.582150212569;
        xi[19]=-0.636853394453;
        xi[20]=0.636853394453;
        xi[21]=-0.688521680771;
        xi[22]=0.688521680771;
        xi[23]=-0.736908848945;
        xi[24]=0.736908848945;
        xi[25]=-0.781784312594;
        xi[26]=0.781784312594;
        xi[27]=-0.822934220502;
        xi[28]=0.822934220502;
        xi[29]=-0.860162475961;
        xi[30]=0.860162475961;
        xi[31]=-0.893291671753;
        xi[32]=0.893291671753;
        xi[33]=-0.922163936719;
        xi[34]=0.922163936719;
        xi[35]=-0.946641690996;
        xi[36]=0.946641690996;
        xi[37]=-0.966608310397;
        xi[38]=0.966608310397;
        xi[39]=-0.981968715035;
        xi[40]=0.981968715035;
        xi[41]=-0.992649998447;
        xi[42]=0.992649998447;
        xi[43]=-0.998603645182;
        xi[44]=0.998603645182;
    }
    else if(no_i==46)
    {
        wi[0]=0.067518685849;
        wi[1]=0.067518685849;
        wi[2]=0.067210613601;
        wi[3]=0.067210613601;
        wi[4]=0.066595874768;
        wi[5]=0.066595874768;
        wi[6]=0.065677274268;
        wi[7]=0.065677274268;
        wi[8]=0.064459003467;
        wi[9]=0.064459003467;
        wi[10]=0.062946621064;
        wi[11]=0.062946621064;
        wi[12]=0.061147027725;
        wi[13]=0.061147027725;
        wi[14]=0.059068434596;
        wi[15]=0.059068434596;
        wi[16]=0.056720325844;
        wi[17]=0.056720325844;
        wi[18]=0.054113415386;
        wi[19]=0.054113415386;
        wi[20]=0.051259598007;
        wi[21]=0.051259598007;
        wi[22]=0.048171895102;
        wi[23]=0.048171895102;
        wi[24]=0.044864395277;
        wi[25]=0.044864395277;
        wi[26]=0.041352190110;
        wi[27]=0.041352190110;
        wi[28]=0.037651305357;
        wi[29]=0.037651305357;
        wi[30]=0.033778627999;
        wi[31]=0.033778627999;
        wi[32]=0.029751829552;
        wi[33]=0.029751829552;
        wi[34]=0.025589286397;
        wi[35]=0.025589286397;
        wi[36]=0.021309998754;
        wi[37]=0.02130999875;
        wi[38]=0.016933514008;
        wi[39]=0.01693351401;
        wi[40]=0.012479883771;
        wi[41]=0.01247988377;
        wi[42]=0.00796989823;
        wi[43]=0.00796989823;
        wi[44]=0.00343030087;
        wi[45]=0.00343030087;
        
        xi[0]=-0.0337721900161;
        xi[1]=0.0337721900161;
        xi[2]=-0.101162475306;
        xi[3]=0.101162475306;
        xi[4]=-0.168091179467;
        xi[5]=0.168091179467;
        xi[6]=-0.234252922206;
        xi[7]=0.234252922206;
        xi[8]=-0.299345822702;
        xi[9]=0.299345822702;
        xi[10]=-0.363072877021;
        xi[11]=0.363072877021;
        xi[12]=-0.425143313283;
        xi[13]=0.425143313283;
        xi[14]=-0.485273918388;
        xi[15]=0.485273918388;
        xi[16]=-0.543190330262;
        xi[17]=0.543190330262;
        xi[18]=-0.598628289713;
        xi[19]=0.598628289713;
        xi[20]=-0.651334846202;
        xi[21]=0.651334846202;
        xi[22]=-0.701069512020;
        xi[23]=0.701069512020;
        xi[24]=-0.747605359616;
        xi[25]=0.747605359616;
        xi[26]=-0.790730057075;
        xi[27]=0.790730057075;
        xi[28]=-0.830246837066;
        xi[29]=0.830246837066;
        xi[30]=-0.865975394867;
        xi[31]=0.865975394867;
        xi[32]=-0.897752711534;
        xi[33]=0.897752711534;
        xi[34]=-0.925433798807;
        xi[35]=0.925433798807;
        xi[36]=-0.948892363446;
        xi[37]=0.948892363446;
        xi[38]=-0.968021391854;
        xi[39]=0.968021391854;
        xi[40]=-0.982733669804;
        xi[41]=0.982733669804;
        xi[42]=-0.992962348906;
        xi[43]=0.992962348906;
        xi[44]=-0.998663042134;
        xi[45]=0.998663042134;
    }
    else if(no_i==47)
    {
        wi[0]=0.0661351296237;
        wi[1]=0.065990533589;
        wi[2]=0.065990533589;
        wi[3]=0.065557377767;
        wi[4]=0.065557377767;
        wi[5]=0.064837556239;
        wi[6]=0.064837556239;
        wi[7]=0.063834216606;
        wi[8]=0.063834216606;
        wi[9]=0.062551746221;
        wi[10]=0.062551746221;
        wi[11]=0.060995753009;
        wi[12]=0.060995753009;
        wi[13]=0.059173040942;
        wi[14]=0.059173040942;
        wi[15]=0.057091580293;
        wi[16]=0.057091580293;
        wi[17]=0.054760472782;
        wi[18]=0.054760472782;
        wi[19]=0.052189911780;
        wi[20]=0.052189911780;
        wi[21]=0.049391137747;
        wi[22]=0.049391137747;
        wi[23]=0.046376389087;
        wi[24]=0.046376389087;
        wi[25]=0.043158848648;
        wi[26]=0.043158848648;
        wi[27]=0.039752586123;
        wi[28]=0.039752586123;
        wi[29]=0.036172496584;
        wi[30]=0.036172496584;
        wi[31]=0.032434235515;
        wi[32]=0.032434235515;
        wi[33]=0.028554150701;
        wi[34]=0.028554150701;
        wi[35]=0.024549211660;
        wi[36]=0.024549211660;
        wi[37]=0.020436938148;
        wi[38]=0.02043693815;
        wi[39]=0.016235333146;
        wi[40]=0.01623533315;
        wi[41]=0.011962848464;
        wi[42]=0.01196284846;
        wi[43]=0.00763861630;
        wi[44]=0.00763861630;
        wi[45]=0.00328745384;
        wi[46]=0.00328745384;
        
        xi[0]=0;
        xi[1]=-0.0660869239164;
        xi[2]=0.0660869239164;
        xi[3]=-0.131884866555;
        xi[4]=0.131884866555;
        xi[5]=-0.197106110279;
        xi[6]=0.197106110279;
        xi[7]=-0.261465459215;
        xi[8]=0.261465459215;
        xi[9]=-0.324681486338;
        xi[10]=0.324681486338;
        xi[11]=-0.386477764085;
        xi[12]=0.386477764085;
        xi[13]=-0.446584073105;
        xi[14]=0.446584073105;
        xi[15]=-0.504737583864;
        xi[16]=0.504737583864;
        xi[17]=-0.560684005935;
        xi[18]=0.560684005935;
        xi[19]=-0.614178699956;
        xi[20]=0.614178699956;
        xi[21]=-0.664987747390;
        xi[22]=0.664987747390;
        xi[23]=-0.712888973409;
        xi[24]=0.712888973409;
        xi[25]=-0.757672918445;
        xi[26]=0.757672918445;
        xi[27]=-0.799143754168;
        xi[28]=0.799143754168;
        xi[29]=-0.837120139900;
        xi[30]=0.837120139900;
        xi[31]=-0.871436015797;
        xi[32]=0.871436015797;
        xi[33]=-0.901941329439;
        xi[34]=0.901941329439;
        xi[35]=-0.928502693012;
        xi[36]=0.928502693012;
        xi[37]=-0.951003969258;
        xi[38]=0.951003969258;
        xi[39]=-0.969346787327;
        xi[40]=0.969346787327;
        xi[41]=-0.983451003072;
        xi[42]=0.983451003072;
        xi[43]=-0.993255210988;
        xi[44]=0.993255210988;
        xi[45]=-0.998718728584;
        xi[46]=0.998718728584;
    }
    else if(no_i==48)
    {
        wi[0]=0.064737696813;
        wi[1]=0.064737696813;
        wi[2]=0.064466164436;
        wi[3]=0.064466164436;
        wi[4]=0.063924238585;
        wi[5]=0.063924238585;
        wi[6]=0.063114192286;
        wi[7]=0.063114192286;
        wi[8]=0.062039423160;
        wi[9]=0.062039423160;
        wi[10]=0.060704439166;
        wi[11]=0.060704439166;
        wi[12]=0.059114839698;
        wi[13]=0.059114839698;
        wi[14]=0.057277292100;
        wi[15]=0.057277292100;
        wi[16]=0.055199503700;
        wi[17]=0.055199503700;
        wi[18]=0.052890189485;
        wi[19]=0.052890189485;
        wi[20]=0.050359035554;
        wi[21]=0.050359035554;
        wi[22]=0.047616658492;
        wi[23]=0.047616658492;
        wi[24]=0.044674560857;
        wi[25]=0.044674560857;
        wi[26]=0.041545082943;
        wi[27]=0.041545082943;
        wi[28]=0.038241351066;
        wi[29]=0.038241351066;
        wi[30]=0.034777222565;
        wi[31]=0.034777222565;
        wi[32]=0.031167227833;
        wi[33]=0.031167227833;
        wi[34]=0.027426509708;
        wi[35]=0.027426509708;
        wi[36]=0.023570760839;
        wi[37]=0.023570760839;
        wi[38]=0.019616160457;
        wi[39]=0.01961616046;
        wi[40]=0.015579315723;
        wi[41]=0.01557931572;
        wi[42]=0.011477234579;
        wi[43]=0.01147723458;
        wi[44]=0.00732755390;
        wi[45]=0.00732755390;
        wi[46]=0.00315334605;
        wi[47]=0.00315334605;
        
        xi[0]=-0.0323801709629;
        xi[1]=0.0323801709629;
        xi[2]=-0.0970046992095;
        xi[3]=0.0970046992095;
        xi[4]=-0.161222356069;
        xi[5]=0.161222356069;
        xi[6]=-0.224763790395;
        xi[7]=0.224763790395;
        xi[8]=-0.287362487355;
        xi[9]=0.287362487355;
        xi[10]=-0.348755886292;
        xi[11]=0.348755886292;
        xi[12]=-0.408686481991;
        xi[13]=0.408686481991;
        xi[14]=-0.466902904751;
        xi[15]=0.466902904751;
        xi[16]=-0.523160974722;
        xi[17]=0.523160974722;
        xi[18]=-0.577224726084;
        xi[19]=0.577224726084;
        xi[20]=-0.628867396777;
        xi[21]=0.628867396777;
        xi[22]=-0.677872379633;
        xi[23]=0.677872379633;
        xi[24]=-0.724034130924;
        xi[25]=0.724034130924;
        xi[26]=-0.767159032516;
        xi[27]=0.767159032516;
        xi[28]=-0.807066204029;
        xi[29]=0.807066204029;
        xi[30]=-0.843588261624;
        xi[31]=0.843588261624;
        xi[32]=-0.876572020274;
        xi[33]=0.876572020274;
        xi[34]=-0.905879136716;
        xi[35]=0.905879136716;
        xi[36]=-0.931386690707;
        xi[37]=0.931386690707;
        xi[38]=-0.952987703160;
        xi[39]=0.952987703160;
        xi[40]=-0.970591592546;
        xi[41]=0.970591592546;
        xi[42]=-0.984124583723;
        xi[43]=0.984124583723;
        xi[44]=-0.993530172266;
        xi[45]=0.993530172266;
        xi[46]=-0.998771007252;
        xi[47]=0.998771007252;
    }
    else if(no_i==49)
    {
        wi[0]=0.0634632814048;
        wi[1]=0.063335509296;
        wi[2]=0.063335509296;
        wi[3]=0.062952707465;
        wi[4]=0.062952707465;
        wi[5]=0.062316417320;
        wi[6]=0.062316417320;
        wi[7]=0.061429200979;
        wi[8]=0.061429200979;
        wi[9]=0.060294630953;
        wi[10]=0.060294630953;
        wi[11]=0.058917275760;
        wi[12]=0.058917275760;
        wi[13]=0.057302681530;
        wi[14]=0.057302681530;
        wi[15]=0.055457349675;
        wi[16]=0.055457349675;
        wi[17]=0.053388710708;
        wi[18]=0.053388710708;
        wi[19]=0.051105094330;
        wi[20]=0.051105094330;
        wi[21]=0.048615695888;
        wi[22]=0.048615695888;
        wi[23]=0.045930539356;
        wi[24]=0.045930539356;
        wi[25]=0.043060436981;
        wi[26]=0.043060436981;
        wi[27]=0.040016945766;
        wi[28]=0.040016945766;
        wi[29]=0.036812320963;
        wi[30]=0.036812320963;
        wi[31]=0.033459466792;
        wi[32]=0.033459466792;
        wi[33]=0.029971884621;
        wi[34]=0.029971884621;
        wi[35]=0.026363618927;
        wi[36]=0.026363618927;
        wi[37]=0.022649201587;
        wi[38]=0.022649201587;
        wi[39]=0.018843595853;
        wi[40]=0.01884359585;
        wi[41]=0.014962144936;
        wi[42]=0.01496214494;
        wi[43]=0.011020551032;
        wi[44]=0.01102055103;
        wi[45]=0.00703509959;
        wi[46]=0.00703509959;
        wi[47]=0.00302727899;
        wi[48]=0.00302727899;
        
        xi[0]=0;
        xi[1]=-0.0634206849827;
        xi[2]=0.0634206849827;
        xi[3]=-0.126585997270;
        xi[4]=0.126585997270;
        xi[5]=-0.189241592462;
        xi[6]=0.189241592462;
        xi[7]=-0.251135178613;
        xi[8]=0.251135178613;
        xi[9]=-0.312017532120;
        xi[10]=0.312017532120;
        xi[11]=-0.371643501262;
        xi[12]=0.371643501262;
        xi[13]=-0.429772993342;
        xi[14]=0.429772993342;
        xi[15]=-0.486171941452;
        xi[16]=0.486171941452;
        xi[17]=-0.540613246992;
        xi[18]=0.540613246992;
        xi[19]=-0.592877694109;
        xi[20]=0.592877694109;
        xi[21]=-0.642754832419;
        xi[22]=0.642754832419;
        xi[23]=-0.690043824425;
        xi[24]=0.690043824425;
        xi[25]=-0.734554254237;
        xi[26]=0.734554254237;
        xi[27]=-0.776106894345;
        xi[28]=0.776106894345;
        xi[29]=-0.814534427360;
        xi[30]=0.814534427360;
        xi[31]=-0.849682119844;
        xi[32]=0.849682119844;
        xi[33]=-0.881408445573;
        xi[34]=0.881408445573;
        xi[35]=-0.909585655828;
        xi[36]=0.909585655828;
        xi[37]=-0.934100294756;
        xi[38]=0.934100294756;
        xi[39]=-0.954853658674;
        xi[40]=0.954853658674;
        xi[41]=-0.971762200902;
        xi[42]=0.971762200902;
        xi[43]=-0.984757895914;
        xi[44]=0.984757895914;
        xi[45]=-0.993788661944;
        xi[46]=0.993788661944;
        xi[47]=-0.998820150607;
        xi[48]=0.998820150607;
    }
    else if(no_i==50)
    {
        wi[0]=0.062176616655;
        wi[1]=0.062176616655;
        wi[2]=0.061936067421;
        wi[3]=0.061936067421;
        wi[4]=0.061455899590;
        wi[5]=0.061455899590;
        wi[6]=0.060737970842;
        wi[7]=0.060737970842;
        wi[8]=0.059785058704;
        wi[9]=0.059785058704;
        wi[10]=0.058600849813;
        wi[11]=0.058600849813;
        wi[12]=0.057189925648;
        wi[13]=0.057189925648;
        wi[14]=0.055557744806;
        wi[15]=0.055557744806;
        wi[16]=0.053710621889;
        wi[17]=0.053710621889;
        wi[18]=0.051655703070;
        wi[19]=0.051655703070;
        wi[20]=0.049400938449;
        wi[21]=0.049400938449;
        wi[22]=0.046955051304;
        wi[23]=0.046955051304;
        wi[24]=0.044327504339;
        wi[25]=0.044327504339;
        wi[26]=0.041528463090;
        wi[27]=0.041528463090;
        wi[28]=0.038568756613;
        wi[29]=0.038568756613;
        wi[30]=0.035459835615;
        wi[31]=0.035459835615;
        wi[32]=0.032213728224;
        wi[33]=0.032213728224;
        wi[34]=0.028842993581;
        wi[35]=0.028842993581;
        wi[36]=0.025360673570;
        wi[37]=0.025360673570;
        wi[38]=0.021780243170;
        wi[39]=0.021780243170;
        wi[40]=0.018115560713;
        wi[41]=0.01811556071;
        wi[42]=0.014380822761;
        wi[43]=0.01438082276;
        wi[44]=0.010590548384;
        wi[45]=0.01059054838;
        wi[46]=0.00675979920;
        wi[47]=0.00675979920;
        wi[48]=0.00290862255;
        wi[49]=0.00290862255;
        
        xi[0]=-0.0310983383272;
        xi[1]=0.0310983383272;
        xi[2]=-0.0931747015601;
        xi[3]=0.0931747015601;
        xi[4]=-0.154890589998;
        xi[5]=0.154890589998;
        xi[6]=-0.216007236876;
        xi[7]=0.216007236876;
        xi[8]=-0.276288193780;
        xi[9]=0.276288193780;
        xi[10]=-0.335500245419;
        xi[11]=0.335500245419;
        xi[12]=-0.393414311898;
        xi[13]=0.393414311898;
        xi[14]=-0.449806334974;
        xi[15]=0.449806334974;
        xi[16]=-0.504458144907;
        xi[17]=0.504458144907;
        xi[18]=-0.557158304515;
        xi[19]=0.557158304515;
        xi[20]=-0.607702927185;
        xi[21]=0.607702927185;
        xi[22]=-0.655896465685;
        xi[23]=0.655896465685;
        xi[24]=-0.701552468707;
        xi[25]=0.701552468707;
        xi[26]=-0.744494302226;
        xi[27]=0.744494302226;
        xi[28]=-0.784555832900;
        xi[29]=0.784555832900;
        xi[30]=-0.821582070859;
        xi[31]=0.821582070859;
        xi[32]=-0.855429769430;
        xi[33]=0.855429769430;
        xi[34]=-0.885967979524;
        xi[35]=0.885967979524;
        xi[36]=-0.913078556656;
        xi[37]=0.913078556656;
        xi[38]=-0.936656618945;
        xi[39]=0.936656618945;
        xi[40]=-0.956610955243;
        xi[41]=0.956610955243;
        xi[42]=-0.972864385107;
        xi[43]=0.972864385107;
        xi[44]=-0.985354084048;
        xi[45]=0.985354084048;
        xi[46]=-0.994031969432;
        xi[47]=0.994031969432;
        xi[48]=-0.998866404420;
        xi[49]=0.998866404420;
    }
    else if(no_i==51)
    {
        wi[0]=0.0609989248412;
        wi[1]=0.060885464845;
        wi[2]=0.060885464845;
        wi[3]=0.060545506935;
        wi[4]=0.060545506935;
        wi[5]=0.059980315778;
        wi[6]=0.059980315778;
        wi[7]=0.059191993923;
        wi[8]=0.059191993923;
        wi[9]=0.058183473983;
        wi[10]=0.058183473983;
        wi[11]=0.056958507720;
        wi[12]=0.056958507720;
        wi[13]=0.055521652096;
        wi[14]=0.055521652096;
        wi[15]=0.053878252313;
        wi[16]=0.053878252313;
        wi[17]=0.052034421937;
        wi[18]=0.052034421937;
        wi[19]=0.049997020150;
        wi[20]=0.049997020150;
        wi[21]=0.047773626241;
        wi[22]=0.047773626241;
        wi[23]=0.045372511408;
        wi[24]=0.045372511408;
        wi[25]=0.042802607998;
        wi[26]=0.042802607998;
        wi[27]=0.040073476285;
        wi[28]=0.040073476285;
        wi[29]=0.037195268923;
        wi[30]=0.037195268923;
        wi[31]=0.034178693204;
        wi[32]=0.034178693204;
        wi[33]=0.031034971290;
        wi[34]=0.031034971290;
        wi[35]=0.027775798594;
        wi[36]=0.027775798594;
        wi[37]=0.024413300574;
        wi[38]=0.024413300574;
        wi[39]=0.020959988402;
        wi[40]=0.020959988402;
        wi[41]=0.017428714723;
        wi[42]=0.01742871472;
        wi[43]=0.013832634006;
        wi[44]=0.01383263401;
        wi[45]=0.010185191298;
        wi[46]=0.01018519130;
        wi[47]=0.00650033778;
        wi[48]=0.00650033778;
        wi[49]=0.00279680717;
        wi[50]=0.00279680717;
        
        xi[0]=0;
        xi[1]=-0.0609611001506;
        xi[2]=0.0609611001506;
        xi[3]=-0.121695421019;
        xi[4]=0.121695421019;
        xi[5]=-0.181977026957;
        xi[6]=0.181977026957;
        xi[7]=-0.241581666448;
        xi[8]=0.241581666448;
        xi[9]=-0.300287606335;
        xi[10]=0.300287606335;
        xi[11]=-0.357876456688;
        xi[12]=0.357876456688;
        xi[13]=-0.414133983226;
        xi[14]=0.414133983226;
        xi[15]=-0.468850904286;
        xi[16]=0.468850904286;
        xi[17]=-0.521823669366;
        xi[18]=0.521823669366;
        xi[19]=-0.572855216351;
        xi[20]=0.572855216351;
        xi[21]=-0.621755704601;
        xi[22]=0.621755704601;
        xi[23]=-0.668343221175;
        xi[24]=0.668343221175;
        xi[25]=-0.712444457577;
        xi[26]=0.712444457577;
        xi[27]=-0.753895354485;
        xi[28]=0.753895354485;
        xi[29]=-0.792541712099;
        xi[30]=0.792541712099;
        xi[31]=-0.828239763823;
        xi[32]=0.828239763823;
        xi[33]=-0.860856711182;
        xi[34]=0.860856711182;
        xi[35]=-0.890271218030;
        xi[36]=0.890271218030;
        xi[37]=-0.916373862310;
        xi[38]=0.916373862310;
        xi[39]=-0.939067544003;
        xi[40]=0.939067544003;
        xi[41]=-0.958267848614;
        xi[42]=0.958267848614;
        xi[43]=-0.973903368019;
        xi[44]=0.973903368019;
        xi[45]=-0.985915991736;
        xi[46]=0.985915991736;
        xi[47]=-0.994261260437;
        xi[48]=0.994261260437;
        xi[49]=-0.998909990849;
        xi[50]=0.998909990849;
    }
    else if(no_i==52)
    {
        wi[0]=0.059810365745;
        wi[1]=0.059810365745;
        wi[2]=0.059596260171;
        wi[3]=0.059596260171;
        wi[4]=0.059168815466;
        wi[5]=0.059168815466;
        wi[6]=0.058529561772;
        wi[7]=0.058529561772;
        wi[8]=0.057680787453;
        wi[9]=0.057680787453;
        wi[10]=0.056625530902;
        wi[11]=0.056625530902;
        wi[12]=0.055367569669;
        wi[13]=0.055367569669;
        wi[14]=0.053911406933;
        wi[15]=0.053911406933;
        wi[16]=0.052262255384;
        wi[17]=0.052262255384;
        wi[18]=0.050426018566;
        wi[19]=0.050426018566;
        wi[20]=0.048409269744;
        wi[21]=0.048409269744;
        wi[22]=0.046219228373;
        wi[23]=0.046219228373;
        wi[24]=0.043863734259;
        wi[25]=0.043863734259;
        wi[26]=0.041351219501;
        wi[27]=0.041351219501;
        wi[28]=0.038690678310;
        wi[29]=0.038690678310;
        wi[30]=0.035891634835;
        wi[31]=0.035891634835;
        wi[32]=0.032964109090;
        wi[33]=0.032964109090;
        wi[34]=0.029918581147;
        wi[35]=0.029918581147;
        wi[36]=0.026765953747;
        wi[37]=0.026765953747;
        wi[38]=0.023517513554;
        wi[39]=0.023517513554;
        wi[40]=0.020184891508;
        wi[41]=0.020184891508;
        wi[42]=0.016780023396;
        wi[43]=0.01678002340;
        wi[44]=0.013315114982;
        wi[45]=0.01331511498;
        wi[46]=0.009802634579;
        wi[47]=0.00980263458;
        wi[48]=0.00625552396;
        wi[49]=0.00625552396;
        wi[50]=0.00269131695;
        wi[51]=0.00269131695;
        
        xi[0]=-0.0299141097973;
        xi[1]=0.0299141097973;
        xi[2]=-0.0896352446489;
        xi[3]=0.0896352446489;
        xi[4]=-0.149035508607;
        xi[5]=0.149035508607;
        xi[6]=-0.207902264156;
        xi[7]=0.207902264156;
        xi[8]=-0.266024783605;
        xi[9]=0.266024783605;
        xi[10]=-0.323195003435;
        xi[11]=0.323195003435;
        xi[12]=-0.379208269116;
        xi[13]=0.379208269116;
        xi[14]=-0.433864067719;
        xi[15]=0.433864067719;
        xi[16]=-0.486966745698;
        xi[17]=0.486966745698;
        xi[18]=-0.538326209286;
        xi[19]=0.538326209286;
        xi[20]=-0.587758604980;
        xi[21]=0.587758604980;
        xi[22]=-0.635086977695;
        xi[23]=0.635086977695;
        xi[24]=-0.680141904227;
        xi[25]=0.680141904227;
        xi[26]=-0.722762099750;
        xi[27]=0.722762099750;
        xi[28]=-0.762794995194;
        xi[29]=0.762794995194;
        xi[30]=-0.800097283430;
        xi[31]=0.800097283430;
        xi[32]=-0.834535432327;
        xi[33]=0.834535432327;
        xi[34]=-0.865986162846;
        xi[35]=0.865986162846;
        xi[36]=-0.894336890534;
        xi[37]=0.894336890534;
        xi[38]=-0.919486128916;
        xi[39]=0.919486128916;
        xi[40]=-0.941343853641;
        xi[41]=0.941343853641;
        xi[42]=-0.959831826933;
        xi[43]=0.959831826933;
        xi[44]=-0.974883884222;
        xi[45]=0.974883884222;
        xi[46]=-0.986446195652;
        xi[47]=0.986446195652;
        xi[48]=-0.994477590929;
        xi[49]=0.994477590929;
        xi[50]=-0.998951111104;
        xi[51]=0.998951111104;
    }
    else if(no_i==53)
    {
        wi[0]=0.0587187941512;
        wi[1]=0.058617586233;
        wi[2]=0.058617586233;
        wi[3]=0.058314311362;
        wi[4]=0.058314311362;
        wi[5]=0.057810014992;
        wi[6]=0.057810014992;
        wi[7]=0.057106435536;
        wi[8]=0.057106435536;
        wi[9]=0.056205998382;
        wi[10]=0.056205998382;
        wi[11]=0.055111807524;
        wi[12]=0.055111807524;
        wi[13]=0.053827634869;
        wi[14]=0.053827634869;
        wi[15]=0.052357907230;
        wi[16]=0.052357907230;
        wi[17]=0.050707691069;
        wi[18]=0.050707691069;
        wi[19]=0.048882675033;
        wi[20]=0.048882675033;
        wi[21]=0.046889150341;
        wi[22]=0.046889150341;
        wi[23]=0.044733989104;
        wi[24]=0.044733989104;
        wi[25]=0.042424620635;
        wi[26]=0.042424620635;
        wi[27]=0.039969005844;
        wi[28]=0.039969005844;
        wi[29]=0.037375609803;
        wi[30]=0.037375609803;
        wi[31]=0.034653372584;
        wi[32]=0.034653372584;
        wi[33]=0.031811678459;
        wi[34]=0.031811678459;
        wi[35]=0.028860323618;
        wi[36]=0.028860323618;
        wi[37]=0.025809482511;
        wi[38]=0.025809482511;
        wi[39]=0.022669673057;
        wi[40]=0.022669673057;
        wi[41]=0.019451721108;
        wi[42]=0.019451721108;
        wi[43]=0.016166725257;
        wi[44]=0.01616672526;
        wi[45]=0.012826026144;
        wi[46]=0.01282602614;
        wi[47]=0.009441202285;
        wi[48]=0.00944120228;
        wi[49]=0.00602427623;
        wi[50]=0.00602427623;
        wi[51]=0.00259168372;
        wi[52]=0.00259168372;
        
        xi[0]=0;
        xi[1]=-0.0586850543003;
        xi[2]=0.0586850543003;
        xi[3]=-0.117167809072;
        xi[4]=0.117167809072;
        xi[5]=-0.175246662155;
        xi[6]=0.175246662155;
        xi[7]=-0.232721403724;
        xi[8]=0.232721403724;
        xi[9]=-0.289393906452;
        xi[10]=0.289393906452;
        xi[11]=-0.345068808496;
        xi[12]=0.345068808496;
        xi[13]=-0.399554186954;
        xi[14]=0.399554186954;
        xi[15]=-0.452662219462;
        xi[16]=0.452662219462;
        xi[17]=-0.504209831657;
        xi[18]=0.504209831657;
        xi[19]=-0.554019328277;
        xi[20]=0.554019328277;
        xi[21]=-0.601919005714;
        xi[22]=0.601919005714;
        xi[23]=-0.647743743917;
        xi[24]=0.647743743917;
        xi[25]=-0.691335575601;
        xi[26]=0.691335575601;
        xi[27]=-0.732544230808;
        xi[28]=0.732544230808;
        xi[29]=-0.771227654926;
        xi[30]=0.771227654926;
        xi[31]=-0.807252498417;
        xi[32]=0.807252498417;
        xi[33]=-0.840494576546;
        xi[34]=0.840494576546;
        xi[35]=-0.870839297558;
        xi[36]=0.870839297558;
        xi[37]=-0.898182057875;
        xi[38]=0.898182057875;
        xi[39]=-0.922428603043;
        xi[40]=0.922428603043;
        xi[41]=-0.943495353464;
        xi[42]=0.943495353464;
        xi[43]=-0.961309694623;
        xi[44]=0.961309694623;
        xi[45]=-0.975810233715;
        xi[46]=0.975810233715;
        xi[47]=-0.986947035023;
        xi[48]=0.986947035023;
        xi[49]=-0.994681919308;
        xi[50]=0.994681919308;
        xi[51]=-0.998989947776;
        xi[52]=0.998989947776;
    }
    else if(no_i==54)
    {
        wi[0]=0.057617536707;
        wi[1]=0.057617536707;
        wi[2]=0.057426137054;
        wi[3]=0.057426137054;
        wi[4]=0.057043973559;
        wi[5]=0.057043973559;
        wi[6]=0.056472315731;
        wi[7]=0.056472315731;
        wi[8]=0.055713062561;
        wi[9]=0.055713062561;
        wi[10]=0.054768736213;
        wi[11]=0.054768736213;
        wi[12]=0.053642473648;
        wi[13]=0.053642473648;
        wi[14]=0.052338016198;
        wi[15]=0.052338016198;
        wi[16]=0.050859697146;
        wi[17]=0.050859697146;
        wi[18]=0.049212427325;
        wi[19]=0.049212427325;
        wi[20]=0.047401678806;
        wi[21]=0.047401678806;
        wi[22]=0.045433466728;
        wi[23]=0.045433466728;
        wi[24]=0.043314329310;
        wi[25]=0.043314329310;
        wi[26]=0.041051306137;
        wi[27]=0.041051306137;
        wi[28]=0.038651914782;
        wi[29]=0.038651914782;
        wi[30]=0.036124125840;
        wi[31]=0.036124125840;
        wi[32]=0.033476336464;
        wi[33]=0.033476336464;
        wi[34]=0.030717342498;
        wi[35]=0.030717342498;
        wi[36]=0.027856309311;
        wi[37]=0.027856309311;
        wi[38]=0.024902741467;
        wi[39]=0.024902741467;
        wi[40]=0.021866451423;
        wi[41]=0.021866451423;
        wi[42]=0.018757527621;
        wi[43]=0.018757527621;
        wi[44]=0.015586303036;
        wi[45]=0.01558630304;
        wi[46]=0.012363328129;
        wi[47]=0.01236332813;
        wi[48]=0.009099369456;
        wi[49]=0.00909936946;
        wi[50]=0.00580561102;
        wi[51]=0.00580561102;
        wi[52]=0.00249748184;
        wi[53]=0.00249748184;
        
        xi[0]=-0.0288167481993;
        xi[1]=0.0288167481993;
        xi[2]=-0.0863545182632;
        xi[3]=0.0863545182632;
        xi[4]=-0.143605427316;
        xi[5]=0.143605427316;
        xi[6]=-0.200379293606;
        xi[7]=0.200379293606;
        xi[8]=-0.256487520070;
        xi[9]=0.256487520070;
        xi[10]=-0.311743720834;
        xi[11]=0.311743720834;
        xi[12]=-0.365964340372;
        xi[13]=0.365964340372;
        xi[14]=-0.418969263255;
        xi[15]=0.418969263255;
        xi[16]=-0.470582412481;
        xi[17]=0.470582412481;
        xi[18]=-0.520632334386;
        xi[19]=0.520632334386;
        xi[20]=-0.568952768195;
        xi[21]=0.568952768195;
        xi[22]=-0.615383198331;
        xi[23]=0.615383198331;
        xi[24]=-0.659769387632;
        xi[25]=0.659769387632;
        xi[26]=-0.701963889719;
        xi[27]=0.701963889719;
        xi[28]=-0.741826538809;
        xi[29]=0.741826538809;
        xi[30]=-0.779224915346;
        xi[31]=0.779224915346;
        xi[32]=-0.814034785914;
        xi[33]=0.814034785914;
        xi[34]=-0.846140515971;
        xi[35]=0.846140515971;
        xi[36]=-0.875435454066;
        xi[37]=0.875435454066;
        xi[38]=-0.901822286285;
        xi[39]=0.901822286285;
        xi[40]=-0.925213359867;
        xi[41]=0.925213359867;
        xi[42]=-0.945530975165;
        xi[43]=0.945530975165;
        xi[44]=-0.962707645786;
        xi[45]=0.962707645786;
        xi[46]=-0.976686328858;
        xi[47]=0.976686328858;
        xi[48]=-0.987420637397;
        xi[49]=0.987420637397;
        xi[50]=-0.994875117018;
        xi[51]=0.994875117018;
        xi[52]=-0.999026666867;
        xi[53]=0.999026666867;
    }
    else if(no_i==55)
    {
        wi[0]=0.0566029764446;
        wi[1]=0.056512318250;
        wi[2]=0.056512318250;
        wi[3]=0.056240634071;
        wi[4]=0.056240634071;
        wi[5]=0.055788794195;
        wi[6]=0.055788794195;
        wi[7]=0.055158246003;
        wi[8]=0.055158246003;
        wi[9]=0.054351009330;
        wi[10]=0.054351009330;
        wi[11]=0.053369670002;
        wi[12]=0.053369670002;
        wi[13]=0.052217371546;
        wi[14]=0.052217371546;
        wi[15]=0.050897805124;
        wi[16]=0.050897805124;
        wi[17]=0.049415197712;
        wi[18]=0.049415197712;
        wi[19]=0.047774298551;
        wi[20]=0.047774298551;
        wi[21]=0.045980363946;
        wi[22]=0.045980363946;
        wi[23]=0.044039140422;
        wi[24]=0.044039140422;
        wi[25]=0.041956846318;
        wi[26]=0.041956846318;
        wi[27]=0.039740151874;
        wi[28]=0.039740151874;
        wi[29]=0.037396157868;
        wi[30]=0.037396157868;
        wi[31]=0.034932372874;
        wi[32]=0.034932372874;
        wi[33]=0.032356689226;
        wi[34]=0.032356689226;
        wi[35]=0.029677357765;
        wi[36]=0.029677357765;
        wi[37]=0.026902961456;
        wi[38]=0.026902961456;
        wi[39]=0.024042388010;
        wi[40]=0.024042388010;
        wi[41]=0.021104801668;
        wi[42]=0.021104801668;
        wi[43]=0.018099614521;
        wi[44]=0.018099614521;
        wi[45]=0.015036458334;
        wi[46]=0.01503645833;
        wi[47]=0.011925160720;
        wi[48]=0.01192516072;
        wi[49]=0.008775746107;
        wi[50]=0.00877574611;
        wi[51]=0.00559863227;
        wi[52]=0.00559863227;
        wi[53]=0.00240832362;
        wi[54]=0.00240832362;
        
        xi[0]=0;
        xi[1]=-0.0565727538183;
        xi[2]=0.0565727538183;
        xi[3]=-0.112964288059;
        xi[4]=0.112964288059;
        xi[5]=-0.168993963647;
        xi[6]=0.168993963647;
        xi[7]=-0.224482300648;
        xi[8]=0.224482300648;
        xi[9]=-0.279251553201;
        xi[10]=0.279251553201;
        xi[11]=-0.333126278890;
        xi[12]=0.333126278890;
        xi[13]=-0.385933900741;
        xi[14]=0.385933900741;
        xi[15]=-0.437505260037;
        xi[16]=0.437505260037;
        xi[17]=-0.487675158187;
        xi[18]=0.487675158187;
        xi[19]=-0.536282885908;
        xi[20]=0.536282885908;
        xi[21]=-0.583172738026;
        xi[22]=0.583172738026;
        xi[23]=-0.628194512250;
        xi[24]=0.628194512250;
        xi[25]=-0.671203990320;
        xi[26]=0.671203990320;
        xi[27]=-0.712063399987;
        xi[28]=0.712063399987;
        xi[29]=-0.750641856348;
        xi[30]=0.750641856348;
        xi[31]=-0.786815781128;
        xi[32]=0.786815781128;
        xi[33]=-0.820469298559;
        xi[34]=0.820469298559;
        xi[35]=-0.851494606617;
        xi[36]=0.851494606617;
        xi[37]=-0.879792322420;
        xi[38]=0.879792322420;
        xi[39]=-0.905271800744;
        xi[40]=0.905271800744;
        xi[41]=-0.927851424721;
        xi[42]=0.927851424721;
        xi[43]=-0.947458868041;
        xi[44]=0.947458868041;
        xi[45]=-0.964031328593;
        xi[46]=0.964031328593;
        xi[47]=-0.977515735504;
        xi[48]=0.977515735504;
        xi[49]=-0.987868941199;
        xi[50]=0.987868941199;
        xi[51]=-0.995057977847;
        xi[52]=0.995057977847;
        xi[53]=-0.999061419565;
        xi[54]=0.999061419565;
    }
    else if(no_i==56)
    {
        wi[0]=0.055579746307;
        wi[1]=0.055579746307;
        wi[2]=0.055407952503;
        wi[3]=0.055407952503;
        wi[4]=0.055064895902;
        wi[5]=0.055064895902;
        wi[6]=0.054551636871;
        wi[7]=0.054551636871;
        wi[8]=0.053869761866;
        wi[9]=0.053869761866;
        wi[10]=0.053021378524;
        wi[11]=0.053021378524;
        wi[12]=0.052009109152;
        wi[13]=0.052009109152;
        wi[14]=0.050836082618;
        wi[15]=0.050836082618;
        wi[16]=0.049505924683;
        wi[17]=0.049505924683;
        wi[18]=0.048022746794;
        wi[19]=0.048022746794;
        wi[20]=0.046391133373;
        wi[21]=0.046391133373;
        wi[22]=0.044616127653;
        wi[23]=0.044616127653;
        wi[24]=0.042703216085;
        wi[25]=0.042703216085;
        wi[26]=0.040658311385;
        wi[27]=0.040658311385;
        wi[28]=0.038487734259;
        wi[29]=0.038487734259;
        wi[30]=0.036198193872;
        wi[31]=0.036198193872;
        wi[32]=0.033796767116;
        wi[33]=0.033796767116;
        wi[34]=0.031290876747;
        wi[35]=0.031290876747;
        wi[36]=0.028688268474;
        wi[37]=0.028688268474;
        wi[38]=0.025996987058;
        wi[39]=0.025996987058;
        wi[40]=0.023225351563;
        wi[41]=0.023225351563;
        wi[42]=0.020381929882;
        wi[43]=0.020381929882;
        wi[44]=0.017475512911;
        wi[45]=0.017475512911;
        wi[46]=0.014515089278;
        wi[47]=0.01451508928;
        wi[48]=0.011509824340;
        wi[49]=0.01150982434;
        wi[50]=0.008469063163;
        wi[51]=0.00846906316;
        wi[52]=0.00540252225;
        wi[53]=0.00540252225;
        wi[54]=0.00232385538;
        wi[55]=0.00232385538;
        
        xi[0]=-0.0277970352873;
        xi[1]=0.0277970352873;
        xi[2]=-0.0833051868224;
        xi[3]=0.0833051868224;
        xi[4]=-0.138555846810;
        xi[5]=0.138555846810;
        xi[6]=-0.193378238635;
        xi[7]=0.193378238635;
        xi[8]=-0.247602909434;
        xi[9]=0.247602909434;
        xi[10]=-0.301062253867;
        xi[11]=0.301062253867;
        xi[12]=-0.353591032175;
        xi[13]=0.353591032175;
        xi[14]=-0.405026880927;
        xi[15]=0.405026880927;
        xi[16]=-0.455210814878;
        xi[17]=0.455210814878;
        xi[18]=-0.503987718384;
        xi[19]=0.503987718384;
        xi[20]=-0.551206824856;
        xi[21]=0.551206824856;
        xi[22]=-0.596722182771;
        xi[23]=0.596722182771;
        xi[24]=-0.640393106807;
        xi[25]=0.640393106807;
        xi[26]=-0.682084612694;
        xi[27]=0.682084612694;
        xi[28]=-0.721667834450;
        xi[29]=0.721667834450;
        xi[30]=-0.759020422705;
        xi[31]=0.759020422705;
        xi[32]=-0.794026922894;
        xi[33]=0.794026922894;
        xi[34]=-0.826579132143;
        xi[35]=0.826579132143;
        xi[36]=-0.856576433763;
        xi[37]=0.856576433763;
        xi[38]=-0.883926108328;
        xi[39]=0.883926108328;
        xi[40]=-0.908543620421;
        xi[41]=0.908543620421;
        xi[42]=-0.930352880247;
        xi[43]=0.930352880247;
        xi[44]=-0.949286479562;
        xi[45]=0.949286479562;
        xi[46]=-0.965285901905;
        xi[47]=0.965285901905;
        xi[48]=-0.978301709140;
        xi[49]=0.978301709140;
        xi[50]=-0.988293715540;
        xi[51]=0.988293715540;
        xi[52]=-0.995231226081;
        xi[53]=0.995231226081;
        xi[54]=-0.999094343801;
        xi[55]=0.999094343801;
    }
    else if(no_i==57)
    {
        wi[0]=0.0546343287566;
        wi[1]=0.054552803605;
        wi[2]=0.054552803605;
        wi[3]=0.054308471452;
        wi[4]=0.054308471452;
        wi[5]=0.053902061483;
        wi[6]=0.053902061483;
        wi[7]=0.053334786585;
        wi[8]=0.053334786585;
        wi[9]=0.052608339729;
        wi[10]=0.052608339729;
        wi[11]=0.051724888921;
        wi[12]=0.051724888921;
        wi[13]=0.050687070725;
        wi[14]=0.050687070725;
        wi[15]=0.049497982402;
        wi[16]=0.049497982402;
        wi[17]=0.048161172662;
        wi[18]=0.048161172662;
        wi[19]=0.046680631074;
        wi[20]=0.046680631074;
        wi[21]=0.045060776161;
        wi[22]=0.045060776161;
        wi[23]=0.043306442216;
        wi[24]=0.043306442216;
        wi[25]=0.041422864871;
        wi[26]=0.041422864871;
        wi[27]=0.039415665475;
        wi[28]=0.039415665475;
        wi[29]=0.037290834324;
        wi[30]=0.037290834324;
        wi[31]=0.035054712782;
        wi[32]=0.035054712782;
        wi[33]=0.032713974366;
        wi[34]=0.032713974366;
        wi[35]=0.030275604843;
        wi[36]=0.030275604843;
        wi[37]=0.027746881402;
        wi[38]=0.027746881402;
        wi[39]=0.025135350991;
        wi[40]=0.025135350991;
        wi[41]=0.022448807891;
        wi[42]=0.022448807891;
        wi[43]=0.019695270699;
        wi[44]=0.019695270699;
        wi[45]=0.016882959023;
        wi[46]=0.016882959023;
        wi[47]=0.014020270791;
        wi[48]=0.01402027079;
        wi[49]=0.011115763732;
        wi[50]=0.01111576373;
        wi[51]=0.008178160068;
        wi[52]=0.00817816007;
        wi[53]=0.00521653347;
        wi[54]=0.00521653347;
        wi[55]=0.00224375387;
        wi[56]=0.00224375387;
        
        xi[0]=0;
        xi[1]=-0.0546071510016;
        xi[2]=0.0546071510016;
        xi[3]=-0.109051332809;
        xi[4]=0.109051332809;
        xi[5]=-0.163170062591;
        xi[6]=0.163170062591;
        xi[7]=-0.216801828796;
        xi[8]=0.216801828796;
        xi[9]=-0.269786573162;
        xi[10]=0.269786573162;
        xi[11]=-0.321966168395;
        xi[12]=0.321966168395;
        xi[13]=-0.373184890087;
        xi[14]=0.373184890087;
        xi[15]=-0.423289881452;
        xi[16]=0.423289881452;
        xi[17]=-0.472131609518;
        xi[18]=0.472131609518;
        xi[19]=-0.519564311391;
        xi[20]=0.519564311391;
        xi[21]=-0.565446429269;
        xi[22]=0.565446429269;
        xi[23]=-0.609641032909;
        xi[24]=0.609641032909;
        xi[25]=-0.652016228281;
        xi[26]=0.652016228281;
        xi[27]=-0.692445551200;
        xi[28]=0.692445551200;
        xi[29]=-0.730808344745;
        xi[30]=0.730808344745;
        xi[31]=-0.766990119359;
        xi[32]=0.766990119359;
        xi[33]=-0.800882894547;
        xi[34]=0.800882894547;
        xi[35]=-0.832385521150;
        xi[36]=0.832385521150;
        xi[37]=-0.861403983262;
        xi[38]=0.861403983262;
        xi[39]=-0.887851678882;
        xi[40]=0.887851678882;
        xi[41]=-0.911649678521;
        xi[42]=0.911649678521;
        xi[43]=-0.932726961067;
        xi[44]=0.932726961067;
        xi[45]=-0.951020626448;
        xi[46]=0.951020626448;
        xi[47]=-0.966476085172;
        xi[48]=0.966476085172;
        xi[49]=-0.979047226709;
        xi[50]=0.979047226709;
        xi[51]=-0.988696577650;
        xi[52]=0.988696577650;
        xi[53]=-0.995395523678;
        xi[54]=0.995395523678;
        xi[55]=-0.999125565625;
        xi[56]=0.999125565625;
    }
    else if(no_i==58)
    {
        wi[0]=0.053681119863;
        wi[1]=0.053681119863;
        wi[2]=0.053526343304;
        wi[3]=0.053526343304;
        wi[4]=0.053217236447;
        wi[5]=0.053217236447;
        wi[6]=0.052754690526;
        wi[7]=0.052754690526;
        wi[8]=0.052140039184;
        wi[9]=0.052140039184;
        wi[10]=0.051375054618;
        wi[11]=0.051375054618;
        wi[12]=0.050461942480;
        wi[13]=0.050461942480;
        wi[14]=0.049403335509;
        wi[15]=0.049403335509;
        wi[16]=0.048202285945;
        wi[17]=0.048202285945;
        wi[18]=0.046862256729;
        wi[19]=0.046862256729;
        wi[20]=0.045387111515;
        wi[21]=0.045387111515;
        wi[22]=0.043781103534;
        wi[23]=0.043781103534;
        wi[24]=0.042048863330;
        wi[25]=0.042048863330;
        wi[26]=0.040195385410;
        wi[27]=0.040195385410;
        wi[28]=0.038226013846;
        wi[29]=0.038226013846;
        wi[30]=0.036146426867;
        wi[31]=0.036146426867;
        wi[32]=0.033962620493;
        wi[33]=0.033962620493;
        wi[34]=0.031680891254;
        wi[35]=0.031680891254;
        wi[36]=0.029307818044;
        wi[37]=0.029307818044;
        wi[38]=0.026850243182;
        wi[39]=0.026850243182;
        wi[40]=0.024315252725;
        wi[41]=0.024315252725;
        wi[42]=0.021710156140;
        wi[43]=0.021710156140;
        wi[44]=0.019042465462;
        wi[45]=0.019042465462;
        wi[46]=0.016319874235;
        wi[47]=0.016319874235;
        wi[48]=0.013550237113;
        wi[49]=0.01355023711;
        wi[50]=0.010741553533;
        wi[51]=0.01074155353;
        wi[52]=0.007901973850;
        wi[53]=0.00790197385;
        wi[54]=0.00503998161;
        wi[55]=0.00503998161;
        wi[56]=0.00216772325;
        wi[57]=0.00216772325;
        
        xi[0]=-0.0268470123659;
        xi[1]=0.0268470123659;
        xi[2]=-0.0804636302141;
        xi[3]=0.0804636302141;
        xi[4]=-0.133848250595;
        xi[5]=0.133848250595;
        xi[6]=-0.186846951836;
        xi[7]=0.186846951836;
        xi[8]=-0.239306924966;
        xi[9]=0.239306924966;
        xi[10]=-0.291076914311;
        xi[11]=0.291076914311;
        xi[12]=-0.342007653598;
        xi[13]=0.342007653598;
        xi[14]=-0.391952296331;
        xi[15]=0.391952296331;
        xi[16]=-0.440766839187;
        xi[17]=0.440766839187;
        xi[18]=-0.488310537217;
        xi[19]=0.488310537217;
        xi[20]=-0.534446309649;
        xi[21]=0.534446309649;
        xi[22]=-0.579041135130;
        xi[23]=0.579041135130;
        xi[24]=-0.621966435263;
        xi[25]=0.621966435263;
        xi[26]=-0.663098445332;
        xi[27]=0.663098445332;
        xi[28]=-0.702318571154;
        xi[29]=0.702318571154;
        xi[30]=-0.739513731020;
        xi[31]=0.739513731020;
        xi[32]=-0.774576681750;
        xi[33]=0.774576681750;
        xi[34]=-0.807406327913;
        xi[35]=0.807406327913;
        xi[36]=-0.837908013339;
        xi[37]=0.837908013339;
        xi[38]=-0.865993794075;
        xi[39]=0.865993794075;
        xi[40]=-0.891582692022;
        xi[41]=0.891582692022;
        xi[42]=-0.914600928564;
        xi[43]=0.914600928564;
        xi[44]=-0.934982137588;
        xi[45]=0.934982137588;
        xi[46]=-0.952667557519;
        xi[47]=0.952667557519;
        xi[48]=-0.967606202503;
        xi[49]=0.967606202503;
        xi[50]=-0.979755014694;
        xi[51]=0.979755014694;
        xi[52]=-0.989079008248;
        xi[53]=0.989079008248;
        xi[54]=-0.995551476597;
        xi[55]=0.995551476597;
        xi[56]=-0.999155200407;
        xi[57]=0.999155200407;
    }
    else if(no_i==59)
    {
        wi[0]=0.0527980126220;
        wi[1]=0.052724433859;
        wi[2]=0.052724433859;
        wi[3]=0.052503902648;
        wi[4]=0.052503902648;
        wi[5]=0.052137033648;
        wi[6]=0.052137033648;
        wi[7]=0.051624849391;
        wi[8]=0.051624849391;
        wi[9]=0.050968777425;
        wi[10]=0.050968777425;
        wi[11]=0.050170646343;
        wi[12]=0.050170646343;
        wi[13]=0.049232680679;
        wi[14]=0.049232680679;
        wi[15]=0.048157494715;
        wi[16]=0.048157494715;
        wi[17]=0.046948085187;
        wi[18]=0.046948085187;
        wi[19]=0.045607822941;
        wi[20]=0.045607822941;
        wi[21]=0.044140443530;
        wi[22]=0.044140443530;
        wi[23]=0.042550036811;
        wi[24]=0.042550036811;
        wi[25]=0.040841035539;
        wi[26]=0.040841035539;
        wi[27]=0.039018203016;
        wi[28]=0.039018203016;
        wi[29]=0.037086619819;
        wi[30]=0.037086619819;
        wi[31]=0.035051669636;
        wi[32]=0.035051669636;
        wi[33]=0.032919024271;
        wi[34]=0.032919024271;
        wi[35]=0.030694627836;
        wi[36]=0.030694627836;
        wi[37]=0.028384680201;
        wi[38]=0.028384680201;
        wi[39]=0.025995619731;
        wi[40]=0.025995619731;
        wi[41]=0.023534105394;
        wi[42]=0.023534105394;
        wi[43]=0.021006998288;
        wi[44]=0.021006998288;
        wi[45]=0.018421342754;
        wi[46]=0.018421342754;
        wi[47]=0.015784347313;
        wi[48]=0.015784347313;
        wi[49]=0.013103366306;
        wi[50]=0.01310336631;
        wi[51]=0.010385885501;
        wi[52]=0.01038588550;
        wi[53]=0.007639529453;
        wi[54]=0.00763952945;
        wi[55]=0.00487223917;
        wi[56]=0.00487223917;
        wi[57]=0.00209549228;
        wi[58]=0.00209549228;
        
        xi[0]=0;
        xi[1]=-0.0527734840883;
        xi[2]=0.0527734840883;
        xi[3]=-0.105399879016;
        xi[4]=0.105399879016;
        xi[5]=-0.157732505588;
        xi[6]=0.157732505588;
        xi[7]=-0.209625503392;
        xi[8]=0.209625503392;
        xi[9]=-0.260934237343;
        xi[10]=0.260934237343;
        xi[11]=-0.311515700803;
        xi[12]=0.311515700803;
        xi[13]=-0.361228914170;
        xi[14]=0.361228914170;
        xi[15]=-0.409935317810;
        xi[16]=0.409935317810;
        xi[17]=-0.457499158253;
        xi[18]=0.457499158253;
        xi[19]=-0.503787866558;
        xi[20]=0.503787866558;
        xi[21]=-0.548672427808;
        xi[22]=0.548672427808;
        xi[23]=-0.592027740704;
        xi[24]=0.592027740704;
        xi[25]=-0.633732966239;
        xi[26]=0.633732966239;
        xi[27]=-0.673671864505;
        xi[28]=0.673671864505;
        xi[29]=-0.711733118677;
        xi[30]=0.711733118677;
        xi[31]=-0.747810645279;
        xi[32]=0.747810645279;
        xi[33]=-0.781803889862;
        xi[34]=0.781803889862;
        xi[35]=-0.813618107288;
        xi[36]=0.813618107288;
        xi[37]=-0.843164625817;
        xi[38]=0.843164625817;
        xi[39]=-0.870361094293;
        xi[40]=0.870361094293;
        xi[41]=-0.895131711743;
        xi[42]=0.895131711743;
        xi[43]=-0.917407438788;
        xi[44]=0.917407438788;
        xi[45]=-0.937126190353;
        xi[46]=0.937126190353;
        xi[47]=-0.954233009377;
        xi[48]=0.954233009377;
        xi[49]=-0.968680221682;
        xi[50]=0.968680221682;
        xi[51]=-0.980427573957;
        xi[52]=0.980427573957;
        xi[53]=-0.989442365134;
        xi[54]=0.989442365134;
        xi[55]=-0.995699640383;
        xi[56]=0.995699640383;
        xi[57]=-0.999183353909;
        xi[58]=0.999183353909;
    }
    else if(no_i==60)
    {
        wi[0]=0.051907877631;
        wi[1]=0.051907877631;
        wi[2]=0.051767943175;
        wi[3]=0.051767943175;
        wi[4]=0.051488451501;
        wi[5]=0.051488451501;
        wi[6]=0.051070156070;
        wi[7]=0.051070156070;
        wi[8]=0.050514184533;
        wi[9]=0.050514184533;
        wi[10]=0.049822035691;
        wi[11]=0.049822035691;
        wi[12]=0.048995575456;
        wi[13]=0.048995575456;
        wi[14]=0.048037031820;
        wi[15]=0.048037031820;
        wi[16]=0.046948988849;
        wi[17]=0.046948988849;
        wi[18]=0.045734379716;
        wi[19]=0.045734379716;
        wi[20]=0.044396478796;
        wi[21]=0.044396478796;
        wi[22]=0.042938892836;
        wi[23]=0.042938892836;
        wi[24]=0.041365551236;
        wi[25]=0.041365551236;
        wi[26]=0.039680695452;
        wi[27]=0.039680695452;
        wi[28]=0.037888867569;
        wi[29]=0.037888867569;
        wi[30]=0.035994898051;
        wi[31]=0.035994898051;
        wi[32]=0.034003892725;
        wi[33]=0.034003892725;
        wi[34]=0.031921219019;
        wi[35]=0.031921219019;
        wi[36]=0.029752491501;
        wi[37]=0.029752491501;
        wi[38]=0.027503556750;
        wi[39]=0.027503556750;
        wi[40]=0.025180477622;
        wi[41]=0.025180477622;
        wi[42]=0.022789516944;
        wi[43]=0.022789516944;
        wi[44]=0.020337120729;
        wi[45]=0.020337120729;
        wi[46]=0.017829901014;
        wi[47]=0.017829901014;
        wi[48]=0.015274618597;
        wi[49]=0.015274618597;
        wi[50]=0.012678166477;
        wi[51]=0.01267816648;
        wi[52]=0.010047557182;
        wi[53]=0.01004755718;
        wi[54]=0.007389931163;
        wi[55]=0.00738993116;
        wi[56]=0.00471272993;
        wi[57]=0.00471272993;
        wi[58]=0.00202681197;
        wi[59]=0.00202681197;
        
        xi[0]=-0.0259597723012;
        xi[1]=0.0259597723012;
        xi[2]=-0.0778093339495;
        xi[3]=0.0778093339495;
        xi[4]=-0.129449135397;
        xi[5]=0.129449135397;
        xi[6]=-0.180739964873;
        xi[7]=0.180739964873;
        xi[8]=-0.231543551376;
        xi[9]=0.231543551376;
        xi[10]=-0.281722937423;
        xi[11]=0.281722937423;
        xi[12]=-0.331142848268;
        xi[13]=0.331142848268;
        xi[14]=-0.379670056577;
        xi[15]=0.379670056577;
        xi[16]=-0.427173741583;
        xi[17]=0.427173741583;
        xi[18]=-0.473525841762;
        xi[19]=0.473525841762;
        xi[20]=-0.518601400059;
        xi[21]=0.518601400059;
        xi[22]=-0.562278900754;
        xi[23]=0.562278900754;
        xi[24]=-0.604440597049;
        xi[25]=0.604440597049;
        xi[26]=-0.644972828489;
        xi[27]=0.644972828489;
        xi[28]=-0.683766327381;
        xi[29]=0.683766327381;
        xi[30]=-0.720716513356;
        xi[31]=0.720716513356;
        xi[32]=-0.755723775307;
        xi[33]=0.755723775307;
        xi[34]=-0.788693739932;
        xi[35]=0.788693739932;
        xi[36]=-0.819537526162;
        xi[37]=0.819537526162;
        xi[38]=-0.848171984786;
        xi[39]=0.848171984786;
        xi[40]=-0.874519922647;
        xi[41]=0.874519922647;
        xi[42]=-0.898510310810;
        xi[43]=0.898510310810;
        xi[44]=-0.920078476178;
        xi[45]=0.920078476178;
        xi[46]=-0.939166276116;
        xi[47]=0.939166276116;
        xi[48]=-0.955722255840;
        xi[49]=0.955722255840;
        xi[50]=-0.969701788765;
        xi[51]=0.969701788765;
        xi[52]=-0.981067201753;
        xi[53]=0.981067201753;
        xi[54]=-0.989787895222;
        xi[55]=0.989787895222;
        xi[56]=-0.995840525119;
        xi[57]=0.995840525119;
        xi[58]=-0.999210123227;
        xi[59]=0.999210123227;
    }
    else if(no_i==61)
    {
        wi[0]=0.0510811194408;
        wi[1]=0.051014487039;
        wi[2]=0.051014487039;
        wi[3]=0.050814763669;
        wi[4]=0.050814763669;
        wi[5]=0.050482470387;
        wi[6]=0.050482470387;
        wi[7]=0.050018474108;
        wi[8]=0.050018474108;
        wi[9]=0.049423985347;
        wi[10]=0.049423985347;
        wi[11]=0.048700555056;
        wi[12]=0.048700555056;
        wi[13]=0.047850070585;
        wi[14]=0.047850070585;
        wi[15]=0.046874750751;
        wi[16]=0.046874750751;
        wi[17]=0.045777140053;
        wi[18]=0.045777140053;
        wi[19]=0.044560102035;
        wi[20]=0.044560102035;
        wi[21]=0.043226811812;
        wi[22]=0.043226811812;
        wi[23]=0.041780747791;
        wi[24]=0.041780747791;
        wi[25]=0.040225682591;
        wi[26]=0.040225682591;
        wi[27]=0.038565673207;
        wi[28]=0.038565673207;
        wi[29]=0.036805050423;
        wi[30]=0.036805050423;
        wi[31]=0.034948407517;
        wi[32]=0.034948407517;
        wi[33]=0.033000588276;
        wi[34]=0.033000588276;
        wi[35]=0.030966674368;
        wi[36]=0.030966674368;
        wi[37]=0.028851972088;
        wi[38]=0.028851972088;
        wi[39]=0.026661998524;
        wi[40]=0.026661998524;
        wi[41]=0.024402467188;
        wi[42]=0.024402467188;
        wi[43]=0.022079273148;
        wi[44]=0.022079273148;
        wi[45]=0.019698477746;
        wi[46]=0.019698477746;
        wi[47]=0.017266292988;
        wi[48]=0.017266292988;
        wi[49]=0.014789065885;
        wi[50]=0.014789065885;
        wi[51]=0.012273263508;
        wi[52]=0.01227326351;
        wi[53]=0.009725461830;
        wi[54]=0.00972546183;
        wi[55]=0.007152354992;
        wi[56]=0.00715235499;
        wi[57]=0.00456092401;
        wi[58]=0.00456092401;
        wi[59]=0.00196145336;
        wi[60]=0.00196145336;
        
        xi[0]=0;
        xi[1]=-0.0510589067080;
        xi[2]=0.0510589067080;
        xi[3]=-0.101984606562;
        xi[4]=0.101984606562;
        xi[5]=-0.152644240231;
        xi[6]=0.152644240231;
        xi[7]=-0.202905642518;
        xi[8]=0.202905642518;
        xi[9]=-0.252637687169;
        xi[10]=0.252637687169;
        xi[11]=-0.301710628963;
        xi[12]=0.301710628963;
        xi[13]=-0.349996442204;
        xi[14]=0.349996442204;
        xi[15]=-0.397369154726;
        xi[16]=0.397369154726;
        xi[17]=-0.443705176539;
        xi[18]=0.443705176539;
        xi[19]=-0.488883622262;
        xi[20]=0.488883622262;
        xi[21]=-0.532786626503;
        xi[22]=0.532786626503;
        xi[23]=-0.575299651351;
        xi[24]=0.575299651351;
        xi[25]=-0.616311785198;
        xi[26]=0.616311785198;
        xi[27]=-0.655716032095;
        xi[28]=0.655716032095;
        xi[29]=-0.693409590894;
        xi[30]=0.693409590894;
        xi[31]=-0.729294123449;
        xi[32]=0.729294123449;
        xi[33]=-0.763276011172;
        xi[34]=0.763276011172;
        xi[35]=-0.795266599282;
        xi[36]=0.795266599282;
        xi[37]=-0.825182428109;
        xi[38]=0.825182428109;
        xi[39]=-0.852945450848;
        xi[40]=0.852945450848;
        xi[41]=-0.878483237215;
        xi[42]=0.878483237215;
        xi[43]=-0.901729162474;
        xi[44]=0.901729162474;
        xi[45]=-0.922622581383;
        xi[46]=0.922622581383;
        xi[47]=-0.941108986681;
        xi[48]=0.941108986681;
        xi[49]=-0.957140151913;
        xi[50]=0.957140151913;
        xi[51]=-0.970674258833;
        xi[52]=0.970674258833;
        xi[53]=-0.981676011284;
        xi[54]=0.981676011284;
        xi[55]=-0.990116745233;
        xi[56]=0.990116745233;
        xi[57]=-0.995974599815;
        xi[58]=0.995974599815;
        xi[59]=-0.999235597631;
        xi[60]=0.999235597631;
    }
    else if(no_i==62)
    {
        wi[0]=0.050248000375;
        wi[1]=0.050248000375;
        wi[2]=0.050121069569;
        wi[3]=0.050121069569;
        wi[4]=0.049867528595;
        wi[5]=0.049867528595;
        wi[6]=0.049488017920;
        wi[7]=0.049488017920;
        wi[8]=0.048983496221;
        wi[9]=0.048983496221;
        wi[10]=0.048355237963;
        wi[11]=0.048355237963;
        wi[12]=0.047604830184;
        wi[13]=0.047604830184;
        wi[14]=0.046734168478;
        wi[15]=0.046734168478;
        wi[16]=0.045745452215;
        wi[17]=0.045745452215;
        wi[18]=0.044641178977;
        wi[19]=0.044641178977;
        wi[20]=0.043424138258;
        wi[21]=0.043424138258;
        wi[22]=0.042097404410;
        wi[23]=0.042097404410;
        wi[24]=0.040664328882;
        wi[25]=0.040664328882;
        wi[26]=0.039128531752;
        wi[27]=0.039128531752;
        wi[28]=0.037493892582;
        wi[29]=0.037493892582;
        wi[30]=0.035764540623;
        wi[31]=0.035764540623;
        wi[32]=0.033944844379;
        wi[33]=0.033944844379;
        wi[34]=0.032039400582;
        wi[35]=0.032039400582;
        wi[36]=0.030053022574;
        wi[37]=0.030053022574;
        wi[38]=0.027990728163;
        wi[39]=0.027990728163;
        wi[40]=0.025857726954;
        wi[41]=0.025857726954;
        wi[42]=0.023659407209;
        wi[43]=0.023659407209;
        wi[44]=0.021401322278;
        wi[45]=0.021401322278;
        wi[46]=0.019089176659;
        wi[47]=0.019089176659;
        wi[48]=0.016728811790;
        wi[49]=0.016728811790;
        wi[50]=0.014326191824;
        wi[51]=0.014326191824;
        wi[52]=0.011887390117;
        wi[53]=0.01188739012;
        wi[54]=0.009418579428;
        wi[55]=0.00941857943;
        wi[56]=0.006926041902;
        wi[57]=0.00692604190;
        wi[58]=0.00441633346;
        wi[59]=0.00441633346;
        wi[60]=0.00189920568;
        wi[61]=0.00189920568;
        
        xi[0]=-0.0251292914218;
        xi[1]=0.0251292914218;
        xi[2]=-0.0753243954962;
        xi[3]=0.0753243954962;
        xi[4]=-0.125329223616;
        xi[5]=0.125329223616;
        xi[6]=-0.175017459249;
        xi[7]=0.175017459249;
        xi[8]=-0.224263585604;
        xi[9]=0.224263585604;
        xi[10]=-0.272943202697;
        xi[11]=0.272943202697;
        xi[12]=-0.320933341594;
        xi[13]=0.320933341594;
        xi[14]=-0.368112775047;
        xi[15]=0.368112775047;
        xi[16]=-0.414362323717;
        xi[17]=0.414362323717;
        xi[18]=-0.459565157240;
        xi[19]=0.459565157240;
        xi[20]=-0.503607089345;
        xi[21]=0.503607089345;
        xi[22]=-0.546376866300;
        xi[23]=0.546376866300;
        xi[24]=-0.587766447953;
        xi[25]=0.587766447953;
        xi[26]=-0.627671280647;
        xi[27]=0.627671280647;
        xi[28]=-0.665990561335;
        xi[29]=0.665990561335;
        xi[30]=-0.702627492222;
        xi[31]=0.702627492222;
        xi[32]=-0.737489525283;
        xi[33]=0.737489525283;
        xi[34]=-0.770488596055;
        xi[35]=0.770488596055;
        xi[36]=-0.801541346104;
        xi[37]=0.801541346104;
        xi[38]=-0.830569333604;
        xi[39]=0.830569333604;
        xi[40]=-0.857499231512;
        xi[41]=0.857499231512;
        xi[42]=-0.882263012832;
        xi[43]=0.882263012832;
        xi[44]=-0.904798122521;
        xi[45]=0.904798122521;
        xi[46]=-0.925047635636;
        xi[47]=0.925047635636;
        xi[48]=-0.942960401392;
        xi[49]=0.942960401392;
        xi[50]=-0.958491172974;
        xi[51]=0.958491172974;
        xi[52]=-0.971600723372;
        xi[53]=0.971600723372;
        xi[54]=-0.982255949097;
        xi[55]=0.982255949097;
        xi[56]=-0.990429971189;
        xi[57]=0.990429971189;
        xi[58]=-0.996102296316;
        xi[59]=0.996102296316;
        xi[60]=-0.999259859309;
        xi[61]=0.999259859309;
    }
    else if(no_i==63)
    {
        wi[0]=0.0494723666239;
        wi[1]=0.0494118330399;
        wi[2]=0.0494118330399;
        wi[3]=0.0492303804237;
        wi[4]=0.0492303804237;
        wi[5]=0.048928452821;
        wi[6]=0.048928452821;
        wi[7]=0.048506789098;
        wi[8]=0.048506789098;
        wi[9]=0.047966421138;
        wi[10]=0.047966421138;
        wi[11]=0.047308671312;
        wi[12]=0.047308671312;
        wi[13]=0.046535149245;
        wi[14]=0.046535149245;
        wi[15]=0.045647747876;
        wi[16]=0.045647747876;
        wi[17]=0.044648638826;
        wi[18]=0.044648638826;
        wi[19]=0.043540267083;
        wi[20]=0.043540267083;
        wi[21]=0.042325345021;
        wi[22]=0.042325345021;
        wi[23]=0.041006845760;
        wi[24]=0.041006845760;
        wi[25]=0.039587995892;
        wi[26]=0.039587995892;
        wi[27]=0.038072267584;
        wi[28]=0.038072267584;
        wi[29]=0.036463370085;
        wi[30]=0.036463370085;
        wi[31]=0.034765240645;
        wi[32]=0.034765240645;
        wi[33]=0.032982034884;
        wi[34]=0.032982034884;
        wi[35]=0.031118116622;
        wi[36]=0.031118116622;
        wi[37]=0.029178047208;
        wi[38]=0.029178047208;
        wi[39]=0.027166574359;
        wi[40]=0.027166574359;
        wi[41]=0.025088620553;
        wi[42]=0.025088620553;
        wi[43]=0.022949271005;
        wi[44]=0.022949271005;
        wi[45]=0.020753761258;
        wi[46]=0.020753761258;
        wi[47]=0.018507464460;
        wi[48]=0.018507464460;
        wi[49]=0.016215878410;
        wi[50]=0.016215878410;
        wi[51]=0.013884612616;
        wi[52]=0.013884612616;
        wi[53]=0.011519376077;
        wi[54]=0.01151937608;
        wi[55]=0.009125968676;
        wi[56]=0.00912596868;
        wi[57]=0.006710291766;
        wi[58]=0.00671029177;
        wi[59]=0.00427850835;
        wi[60]=0.00427850835;
        wi[61]=0.00183987460;
        wi[62]=0.00183987460;
        
        xi[0]=0;
        xi[1]=-0.0494521871162;
        xi[2]=0.0494521871162;
        xi[3]=-0.0987833564469;
        xi[4]=0.0987833564469;
        xi[5]=-0.147872786358;
        xi[6]=0.147872786358;
        xi[7]=-0.196600346792;
        xi[8]=0.196600346792;
        xi[9]=-0.244846793246;
        xi[10]=0.244846793246;
        xi[11]=-0.292494058586;
        xi[12]=0.292494058586;
        xi[13]=-0.339425541975;
        xi[14]=0.339425541975;
        xi[15]=-0.385526394212;
        xi[16]=0.385526394212;
        xi[17]=-0.430683798795;
        xi[18]=0.430683798795;
        xi[19]=-0.474787247995;
        xi[20]=0.474787247995;
        xi[21]=-0.517728813290;
        xi[22]=0.517728813290;
        xi[23]=-0.559403409486;
        xi[24]=0.559403409486;
        xi[25]=-0.599709051878;
        xi[26]=0.599709051878;
        xi[27]=-0.638547105821;
        xi[28]=0.638547105821;
        xi[29]=-0.675822528115;
        xi[30]=0.675822528115;
        xi[31]=-0.711444099585;
        xi[32]=0.711444099585;
        xi[33]=-0.745324648318;
        xi[34]=0.745324648318;
        xi[35]=-0.777381262990;
        xi[36]=0.777381262990;
        xi[37]=-0.807535495773;
        xi[38]=0.807535495773;
        xi[39]=-0.835713554320;
        xi[40]=0.835713554320;
        xi[41]=-0.861846482364;
        xi[42]=0.861846482364;
        xi[43]=-0.885870328508;
        xi[44]=0.885870328508;
        xi[45]=-0.907726302779;
        xi[46]=0.907726302779;
        xi[47]=-0.927360920622;
        xi[48]=0.927360920622;
        xi[49]=-0.944726134041;
        xi[50]=0.944726134041;
        xi[51]=-0.959779449759;
        xi[52]=0.959779449759;
        xi[53]=-0.972484034698;
        xi[54]=0.972484034698;
        xi[55]=-0.982808810594;
        xi[56]=0.982808810594;
        xi[57]=-0.990728546892;
        xi[58]=0.990728546892;
        xi[59]=-0.996224012778;
        xi[60]=0.996224012778;
        xi[61]=-0.999282984029;
        xi[62]=0.999282984029;
    }
    else if(no_i==64)
    {
        wi[0]=0.0486909570091;
        wi[1]=0.0486909570091;
        wi[2]=0.0485754674415;
        wi[3]=0.0485754674415;
        wi[4]=0.0483447622348;
        wi[5]=0.0483447622348;
        wi[6]=0.0479993885965;
        wi[7]=0.0479993885965;
        wi[8]=0.0475401657148;
        wi[9]=0.0475401657148;
        wi[10]=0.046968182816;
        wi[11]=0.046968182816;
        wi[12]=0.046284796581;
        wi[13]=0.046284796581;
        wi[14]=0.045491627927;
        wi[15]=0.045491627927;
        wi[16]=0.044590558164;
        wi[17]=0.044590558164;
        wi[18]=0.043583724529;
        wi[19]=0.043583724529;
        wi[20]=0.042473515124;
        wi[21]=0.042473515124;
        wi[22]=0.041262563243;
        wi[23]=0.041262563243;
        wi[24]=0.039953741133;
        wi[25]=0.039953741133;
        wi[26]=0.038550153179;
        wi[27]=0.038550153179;
        wi[28]=0.037055128540;
        wi[29]=0.037055128540;
        wi[30]=0.035472213257;
        wi[31]=0.035472213257;
        wi[32]=0.033805161837;
        wi[33]=0.033805161837;
        wi[34]=0.032057928355;
        wi[35]=0.032057928355;
        wi[36]=0.030234657072;
        wi[37]=0.030234657072;
        wi[38]=0.028339672614;
        wi[39]=0.028339672614;
        wi[40]=0.026377469715;
        wi[41]=0.026377469715;
        wi[42]=0.024352702569;
        wi[43]=0.024352702569;
        wi[44]=0.022270173808;
        wi[45]=0.022270173808;
        wi[46]=0.020134823154;
        wi[47]=0.020134823154;
        wi[48]=0.017951715776;
        wi[49]=0.017951715776;
        wi[50]=0.015726030476;
        wi[51]=0.015726030476;
        wi[52]=0.013463047897;
        wi[53]=0.013463047897;
        wi[54]=0.011168139460;
        wi[55]=0.01116813946;
        wi[56]=0.008846759826;
        wi[57]=0.00884675983;
        wi[58]=0.006504457969;
        wi[59]=0.00650445797;
        wi[60]=0.00414703326;
        wi[61]=0.00414703326;
        wi[62]=0.00178328072;
        wi[63]=0.00178328072;
        
        xi[0]=-0.0243502926634;
        xi[1]=0.0243502926634;
        xi[2]=-0.0729931217878;
        xi[3]=0.0729931217878;
        xi[4]=-0.121462819296;
        xi[5]=0.121462819296;
        xi[6]=-0.169644420424;
        xi[7]=0.169644420424;
        xi[8]=-0.217423643740;
        xi[9]=0.217423643740;
        xi[10]=-0.264687162209;
        xi[11]=0.264687162209;
        xi[12]=-0.311322871990;
        xi[13]=0.311322871990;
        xi[14]=-0.357220158338;
        xi[15]=0.357220158338;
        xi[16]=-0.402270157964;
        xi[17]=0.402270157964;
        xi[18]=-0.446366017253;
        xi[19]=0.446366017253;
        xi[20]=-0.489403145707;
        xi[21]=0.489403145707;
        xi[22]=-0.531279464020;
        xi[23]=0.531279464020;
        xi[24]=-0.571895646203;
        xi[25]=0.571895646203;
        xi[26]=-0.611155355172;
        xi[27]=0.611155355172;
        xi[28]=-0.648965471255;
        xi[29]=0.648965471255;
        xi[30]=-0.685236313054;
        xi[31]=0.685236313054;
        xi[32]=-0.719881850172;
        xi[33]=0.719881850172;
        xi[34]=-0.752819907261;
        xi[35]=0.752819907261;
        xi[36]=-0.783972358943;
        xi[37]=0.783972358943;
        xi[38]=-0.813265315123;
        xi[39]=0.813265315123;
        xi[40]=-0.840629296253;
        xi[41]=0.840629296253;
        xi[42]=-0.865999398154;
        xi[43]=0.865999398154;
        xi[44]=-0.889315445995;
        xi[45]=0.889315445995;
        xi[46]=-0.910522137079;
        xi[47]=0.910522137079;
        xi[48]=-0.929569172132;
        xi[49]=0.929569172132;
        xi[50]=-0.946411374858;
        xi[51]=0.946411374858;
        xi[52]=-0.961008799652;
        xi[53]=0.961008799652;
        xi[54]=-0.973326827790;
        xi[55]=0.973326827790;
        xi[56]=-0.983336253885;
        xi[57]=0.983336253885;
        xi[58]=-0.991013371477;
        xi[59]=0.991013371477;
        xi[60]=-0.996340116772;
        xi[61]=0.996340116772;
        xi[62]=-0.999305041736;
        xi[63]=0.999305041736;
    }
    else if(no_i==65)
    {
        wi[0]=0.0479618493945;
        wi[1]=0.0479066925005;
        wi[2]=0.0479066925005;
        wi[3]=0.0477413486812;
        wi[4]=0.0477413486812;
        wi[5]=0.0474661982329;
        wi[6]=0.0474661982329;
        wi[7]=0.0470818740105;
        wi[8]=0.0470818740105;
        wi[9]=0.0465892599722;
        wi[10]=0.0465892599722;
        wi[11]=0.045989489147;
        wi[12]=0.045989489147;
        wi[13]=0.045283941026;
        wi[14]=0.045283941026;
        wi[15]=0.044474238395;
        wi[16]=0.044474238395;
        wi[17]=0.043562243596;
        wi[18]=0.043562243596;
        wi[19]=0.042550054247;
        wi[20]=0.042550054247;
        wi[21]=0.041439998417;
        wi[22]=0.041439998417;
        wi[23]=0.040234629273;
        wi[24]=0.040234629273;
        wi[25]=0.038936719204;
        wi[26]=0.038936719204;
        wi[27]=0.037549253448;
        wi[28]=0.037549253448;
        wi[29]=0.036075423226;
        wi[30]=0.036075423226;
        wi[31]=0.034518618399;
        wi[32]=0.034518618399;
        wi[33]=0.032882419676;
        wi[34]=0.032882419676;
        wi[35]=0.031170590380;
        wi[36]=0.031170590380;
        wi[37]=0.029387067789;
        wi[38]=0.029387067789;
        wi[39]=0.027535954088;
        wi[40]=0.027535954088;
        wi[41]=0.025621506938;
        wi[42]=0.025621506938;
        wi[43]=0.023648129691;
        wi[44]=0.023648129691;
        wi[45]=0.021620361285;
        wi[46]=0.021620361285;
        wi[47]=0.019542865837;
        wi[48]=0.019542865837;
        wi[49]=0.017420421998;
        wi[50]=0.017420421998;
        wi[51]=0.015257912146;
        wi[52]=0.015257912146;
        wi[53]=0.013060311640;
        wi[54]=0.013060311640;
        wi[55]=0.010832678790;
        wi[56]=0.01083267879;
        wi[57]=0.008580148267;
        wi[58]=0.00858014827;
        wi[59]=0.006307942579;
        wi[60]=0.00630794258;
        wi[61]=0.00402152417;
        wi[62]=0.00402152417;
        wi[63]=0.00172925825;
        wi[64]=0.00172925825;
        
        xi[0]=0;
        xi[1]=-0.0479434623532;
        xi[2]=0.0479434623532;
        xi[3]=-0.0957766532092;
        xi[4]=0.0957766532092;
        xi[5]=-0.143389554699;
        xi[6]=0.143389554699;
        xi[7]=-0.190672655626;
        xi[8]=0.190672655626;
        xi[9]=-0.237517203346;
        xi[10]=0.237517203346;
        xi[11]=-0.283815453902;
        xi[12]=0.283815453902;
        xi[13]=-0.329460919837;
        xi[14]=0.329460919837;
        xi[15]=-0.374348615122;
        xi[16]=0.374348615122;
        xi[17]=-0.418375296623;
        xi[18]=0.418375296623;
        xi[19]=-0.461439701569;
        xi[20]=0.461439701569;
        xi[21]=-0.503442780455;
        xi[22]=0.503442780455;
        xi[23]=-0.544287924862;
        xi[24]=0.544287924862;
        xi[25]=-0.583881189660;
        xi[26]=0.583881189660;
        xi[27]=-0.622131509085;
        xi[28]=0.622131509085;
        xi[29]=-0.658950906194;
        xi[30]=0.658950906194;
        xi[31]=-0.694254695214;
        xi[32]=0.694254695214;
        xi[33]=-0.727961676329;
        xi[34]=0.727961676329;
        xi[35]=-0.759994322442;
        xi[36]=0.759994322442;
        xi[37]=-0.790278957492;
        xi[38]=0.790278957492;
        xi[39]=-0.818745925923;
        xi[40]=0.818745925923;
        xi[41]=-0.845329752900;
        xi[42]=0.845329752900;
        xi[43]=-0.869969294926;
        xi[44]=0.869969294926;
        xi[45]=-0.892607880505;
        xi[46]=0.892607880505;
        xi[47]=-0.913193440543;
        xi[48]=0.913193440543;
        xi[49]=-0.931678628229;
        xi[50]=0.931678628229;
        xi[51]=-0.948020928168;
        xi[52]=0.948020928168;
        xi[53]=-0.962182754718;
        xi[54]=0.962182754718;
        xi[55]=-0.974131539834;
        xi[56]=0.974131539834;
        xi[57]=-0.983839812187;
        xi[58]=0.983839812187;
        xi[59]=-0.991285276177;
        xi[60]=0.991285276177;
        xi[61]=-0.996450948062;
        xi[62]=0.996450948062;
        xi[63]=-0.999326097075;
        xi[64]=0.999326097075;
    }
    else if(no_i==66)
    {
        wi[0]=0.0472274812630;
        wi[1]=0.0472274812630;
        wi[2]=0.0471220982876;
        wi[3]=0.0471220982876;
        wi[4]=0.0469115674876;
        wi[5]=0.0469115674876;
        wi[6]=0.0465963586396;
        wi[7]=0.0465963586396;
        wi[8]=0.0461771750979;
        wi[9]=0.0461771750979;
        wi[10]=0.0456549522253;
        wi[11]=0.0456549522253;
        wi[12]=0.0450308553054;
        wi[13]=0.0450308553054;
        wi[14]=0.044306276943;
        wi[15]=0.044306276943;
        wi[16]=0.043482833957;
        wi[17]=0.043482833957;
        wi[18]=0.042562363770;
        wi[19]=0.042562363770;
        wi[20]=0.041546920313;
        wi[21]=0.041546920313;
        wi[22]=0.040438769439;
        wi[23]=0.040438769439;
        wi[24]=0.039240383867;
        wi[25]=0.039240383867;
        wi[26]=0.037954437666;
        wi[27]=0.037954437666;
        wi[28]=0.036583800288;
        wi[29]=0.036583800288;
        wi[30]=0.035131530165;
        wi[31]=0.035131530165;
        wi[32]=0.033600867886;
        wi[33]=0.033600867886;
        wi[34]=0.031995228964;
        wi[35]=0.031995228964;
        wi[36]=0.030318196219;
        wi[37]=0.030318196219;
        wi[38]=0.028573511783;
        wi[39]=0.028573511783;
        wi[40]=0.026765068754;
        wi[41]=0.026765068754;
        wi[42]=0.024896902515;
        wi[43]=0.024896902515;
        wi[44]=0.022973181735;
        wi[45]=0.022973181735;
        wi[46]=0.020998199092;
        wi[47]=0.020998199092;
        wi[48]=0.018976361723;
        wi[49]=0.018976361723;
        wi[50]=0.016912181472;
        wi[51]=0.016912181472;
        wi[52]=0.014810265003;
        wi[53]=0.014810265003;
        wi[54]=0.012675303981;
        wi[55]=0.012675303981;
        wi[56]=0.010512065988;
        wi[57]=0.01051206599;
        wi[58]=0.008325388766;
        wi[59]=0.00832538877;
        wi[60]=0.006120192018;
        wi[61]=0.00612019202;
        wi[62]=0.00390162564;
        wi[63]=0.00390162564;
        wi[64]=0.00167765374;
        wi[65]=0.00167765374;
        
        xi[0]=-0.0236181333859;
        xi[1]=0.0236181333859;
        xi[2]=-0.0708016988681;
        xi[3]=0.0708016988681;
        xi[4]=-0.117827278079;
        xi[5]=0.117827278079;
        xi[6]=-0.164589938565;
        xi[7]=0.164589938565;
        xi[8]=-0.210985334548;
        xi[9]=0.210985334548;
        xi[10]=-0.256909939762;
        xi[11]=0.256909939762;
        xi[12]=-0.302261278456;
        xi[13]=0.302261278456;
        xi[14]=-0.346938154066;
        xi[15]=0.346938154066;
        xi[16]=-0.390840875016;
        xi[17]=0.390840875016;
        xi[18]=-0.433871477173;
        xi[19]=0.433871477173;
        xi[20]=-0.475933942444;
        xi[21]=0.475933942444;
        xi[22]=-0.516934413028;
        xi[23]=0.516934413028;
        xi[24]=-0.556781400852;
        xi[25]=0.556781400852;
        xi[26]=-0.595385991714;
        xi[27]=0.595385991714;
        xi[28]=-0.632662043691;
        xi[29]=0.632662043691;
        xi[30]=-0.668526379351;
        xi[31]=0.668526379351;
        xi[32]=-0.702898971359;
        xi[33]=0.702898971359;
        xi[34]=-0.735703121050;
        xi[35]=0.735703121050;
        xi[36]=-0.766865629575;
        xi[37]=0.766865629575;
        xi[38]=-0.796316961246;
        xi[39]=0.796316961246;
        xi[40]=-0.823991398696;
        xi[41]=0.823991398696;
        xi[42]=-0.849827189536;
        xi[43]=0.849827189536;
        xi[44]=-0.873766684167;
        xi[45]=0.873766684167;
        xi[46]=-0.895756464447;
        xi[47]=0.895756464447;
        xi[48]=-0.915747462945;
        xi[49]=0.915747462945;
        xi[50]=-0.933695072545;
        xi[51]=0.933695072545;
        xi[52]=-0.949559246206;
        xi[53]=0.949559246206;
        xi[54]=-0.963304586859;
        xi[55]=0.963304586859;
        xi[56]=-0.974900427746;
        xi[57]=0.974900427746;
        xi[58]=-0.984320904935;
        xi[59]=0.984320904935;
        xi[60]=-0.991545030383;
        xi[61]=0.991545030383;
        xi[62]=-0.996556821090;
        xi[63]=0.996556821090;
        xi[64]=-0.999346209872;
        xi[65]=0.999346209872;
    }
    else if(no_i==67)
    {
        wi[0]=0.0465408367036;
        wi[1]=0.0464904381603;
        wi[2]=0.0464904381603;
        wi[3]=0.0463393516824;
        wi[4]=0.0463393516824;
        wi[5]=0.0460879044898;
        wi[6]=0.0460879044898;
        wi[7]=0.0457366411611;
        wi[8]=0.0457366411611;
        wi[9]=0.0452863224547;
        wi[10]=0.0452863224547;
        wi[11]=0.0447379236609;
        wi[12]=0.0447379236609;
        wi[13]=0.0440926324898;
        wi[14]=0.0440926324898;
        wi[15]=0.0433518464987;
        wi[16]=0.0433518464987;
        wi[17]=0.042517170066;
        wi[18]=0.042517170066;
        wi[19]=0.041590410915;
        wi[20]=0.041590410915;
        wi[21]=0.040573576202;
        wi[22]=0.040573576202;
        wi[23]=0.039468868164;
        wi[24]=0.039468868164;
        wi[25]=0.038278679356;
        wi[26]=0.038278679356;
        wi[27]=0.037005587463;
        wi[28]=0.037005587463;
        wi[29]=0.035652349723;
        wi[30]=0.035652349723;
        wi[31]=0.034221896950;
        wi[32]=0.034221896950;
        wi[33]=0.032717327192;
        wi[34]=0.032717327192;
        wi[35]=0.031141899019;
        wi[36]=0.031141899019;
        wi[37]=0.029499024471;
        wi[38]=0.029499024471;
        wi[39]=0.027792261662;
        wi[40]=0.027792261662;
        wi[41]=0.026025307086;
        wi[42]=0.026025307086;
        wi[43]=0.024201987610;
        wi[44]=0.024201987610;
        wi[45]=0.022326252196;
        wi[46]=0.022326252196;
        wi[47]=0.020402163371;
        wi[48]=0.020402163371;
        wi[49]=0.018433888457;
        wi[50]=0.018433888457;
        wi[51]=0.016425690623;
        wi[52]=0.016425690623;
        wi[53]=0.014381919827;
        wi[54]=0.014381919827;
        wi[55]=0.012307003849;
        wi[56]=0.012307003849;
        wi[57]=0.010205440034;
        wi[58]=0.01020544003;
        wi[59]=0.008081790299;
        wi[60]=0.00808179030;
        wi[61]=0.005940693178;
        wi[62]=0.00594069318;
        wi[63]=0.00378700830;
        wi[64]=0.00378700830;
        wi[65]=0.00162832504;
        wi[66]=0.00162832504;
        
        xi[0]=0;
        xi[1]=-0.0465240359760;
        xi[2]=0.0465240359760;
        xi[3]=-0.0929473112521;
        xi[4]=0.0929473112521;
        xi[5]=-0.139169283354;
        xi[6]=0.139169283354;
        xi[7]=-0.185089845784;
        xi[8]=0.185089845784;
        xi[9]=-0.230609544832;
        xi[10]=0.230609544832;
        xi[11]=-0.275629794970;
        xi[12]=0.275629794970;
        xi[13]=-0.320053092363;
        xi[14]=0.320053092363;
        xi[15]=-0.363783226044;
        xi[16]=0.363783226044;
        xi[17]=-0.406725486286;
        xi[18]=0.406725486286;
        xi[19]=-0.448786869720;
        xi[20]=0.448786869720;
        xi[21]=-0.489876280761;
        xi[22]=0.489876280761;
        xi[23]=-0.529904728900;
        xi[24]=0.529904728900;
        xi[25]=-0.568785521441;
        xi[26]=0.568785521441;
        xi[27]=-0.606434451256;
        xi[28]=0.606434451256;
        xi[29]=-0.642769979157;
        xi[30]=0.642769979157;
        xi[31]=-0.677713410498;
        xi[32]=0.677713410498;
        xi[33]=-0.711189065610;
        xi[34]=0.711189065610;
        xi[35]=-0.743124443704;
        xi[36]=0.743124443704;
        xi[37]=-0.773450379900;
        xi[38]=0.773450379900;
        xi[39]=-0.802101195022;
        xi[40]=0.802101195022;
        xi[41]=-0.829014837853;
        xi[42]=0.829014837853;
        xi[43]=-0.854133019537;
        xi[44]=0.854133019537;
        xi[45]=-0.877401339832;
        xi[46]=0.877401339832;
        xi[47]=-0.898769404961;
        xi[48]=0.898769404961;
        xi[49]=-0.918190936813;
        xi[50]=0.918190936813;
        xi[51]=-0.935623873271;
        xi[52]=0.935623873271;
        xi[53]=-0.951030459540;
        xi[54]=0.951030459540;
        xi[55]=-0.964377330428;
        xi[56]=0.964377330428;
        xi[57]=-0.975635583923;
        xi[58]=0.975635583923;
        xi[59]=-0.984780847756;
        xi[60]=0.984780847756;
        xi[61]=-0.991793347082;
        xi[62]=0.991793347082;
        xi[63]=-0.996658027206;
        xi[64]=0.996658027206;
        xi[65]=-0.999365435553;
        xi[66]=0.999365435553;
    }
    else if(no_i==68)
    {
        wi[0]=0.0458493873873;
        wi[1]=0.0458493873873;
        wi[2]=0.0457529654161;
        wi[3]=0.0457529654161;
        wi[4]=0.0455603242506;
        wi[5]=0.0455603242506;
        wi[6]=0.0452718690184;
        wi[7]=0.0452718690184;
        wi[8]=0.0448882063454;
        wi[9]=0.0448882063454;
        wi[10]=0.0444101430804;
        wi[11]=0.0444101430804;
        wi[12]=0.0438386845980;
        wi[13]=0.0438386845980;
        wi[14]=0.0431750326846;
        wi[15]=0.0431750326846;
        wi[16]=0.0424205830111;
        wi[17]=0.0424205830111;
        wi[18]=0.041576922197;
        wi[19]=0.041576922197;
        wi[20]=0.040645824476;
        wi[21]=0.040645824476;
        wi[22]=0.039629247961;
        wi[23]=0.039629247961;
        wi[24]=0.038529330529;
        wi[25]=0.038529330529;
        wi[26]=0.037348385326;
        wi[27]=0.037348385326;
        wi[28]=0.036088895900;
        wi[29]=0.036088895900;
        wi[30]=0.034753510980;
        wi[31]=0.034753510980;
        wi[32]=0.033345038904;
        wi[33]=0.033345038904;
        wi[34]=0.031866441717;
        wi[35]=0.031866441717;
        wi[36]=0.030320828939;
        wi[37]=0.030320828939;
        wi[38]=0.028711451027;
        wi[39]=0.028711451027;
        wi[40]=0.027041692546;
        wi[41]=0.027041692546;
        wi[42]=0.025315065045;
        wi[43]=0.025315065045;
        wi[44]=0.023535199686;
        wi[45]=0.023535199686;
        wi[46]=0.021705839610;
        wi[47]=0.021705839610;
        wi[48]=0.019830832088;
        wi[49]=0.019830832088;
        wi[50]=0.017914120458;
        wi[51]=0.017914120458;
        wi[52]=0.015959735910;
        wi[53]=0.015959735910;
        wi[54]=0.013971789174;
        wi[55]=0.013971789174;
        wi[56]=0.011954462320;
        wi[57]=0.011954462320;
        wi[58]=0.009912001252;
        wi[59]=0.00991200125;
        wi[60]=0.007848711393;
        wi[61]=0.00784871139;
        wi[62]=0.005768969919;
        wi[63]=0.00576896992;
        wi[64]=0.00367736660;
        wi[65]=0.00367736660;
        wi[66]=0.00158114026;
        wi[67]=0.00158114026;
        
        xi[0]=-0.0229287128254;
        xi[1]=0.0229287128254;
        xi[2]=-0.0687379190384;
        xi[3]=0.0687379190384;
        xi[4]=-0.114402568344;
        xi[5]=0.114402568344;
        xi[6]=-0.159826627278;
        xi[7]=0.159826627278;
        xi[8]=-0.204914568342;
        xi[9]=0.204914568342;
        xi[10]=-0.249571570898;
        xi[11]=0.249571570898;
        xi[12]=-0.293703720581;
        xi[13]=0.293703720581;
        xi[14]=-0.337218206796;
        xi[15]=0.337218206796;
        xi[16]=-0.380023517908;
        xi[17]=0.380023517908;
        xi[18]=-0.422029633686;
        xi[19]=0.422029633686;
        xi[20]=-0.463148214623;
        xi[21]=0.463148214623;
        xi[22]=-0.503292787709;
        xi[23]=0.503292787709;
        xi[24]=-0.542378928291;
        xi[25]=0.542378928291;
        xi[26]=-0.580324437616;
        xi[27]=0.580324437616;
        xi[28]=-0.617049515700;
        xi[29]=0.617049515700;
        xi[30]=-0.652476929145;
        xi[31]=0.652476929145;
        xi[32]=-0.686532173566;
        xi[33]=0.686532173566;
        xi[34]=-0.719143630276;
        xi[35]=0.719143630276;
        xi[36]=-0.750242716903;
        xi[37]=0.750242716903;
        xi[38]=-0.779764031619;
        xi[39]=0.779764031619;
        xi[40]=-0.807645490690;
        xi[41]=0.807645490690;
        xi[42]=-0.833828459042;
        xi[43]=0.833828459042;
        xi[44]=-0.858257873583;
        xi[45]=0.858257873583;
        xi[46]=-0.880882359017;
        xi[47]=0.880882359017;
        xi[48]=-0.901654335920;
        xi[49]=0.901654335920;
        xi[50]=-0.920530120848;
        xi[51]=0.920530120848;
        xi[52]=-0.937470018313;
        xi[53]=0.937470018313;
        xi[54]=-0.952438404480;
        xi[55]=0.952438404480;
        xi[56]=-0.965403802582;
        xi[57]=0.965403802582;
        xi[58]=-0.976338950400;
        xi[59]=0.976338950400;
        xi[60]=-0.985220861440;
        xi[61]=0.985220861440;
        xi[62]=-0.992030887742;
        xi[63]=0.992030887742;
        xi[64]=-0.996754836677;
        xi[65]=0.996754836677;
        xi[66]=-0.999383825530;
        xi[67]=0.999383825530;
    }
    else if(no_i==69)
    {
        wi[0]=0.0452016023771;
        wi[1]=0.0451554302361;
        wi[2]=0.0451554302361;
        wi[3]=0.0450170081404;
        wi[4]=0.0450170081404;
        wi[5]=0.0447866188783;
        wi[6]=0.0447866188783;
        wi[7]=0.0444647331220;
        wi[8]=0.0444647331220;
        wi[9]=0.0440520084659;
        wi[10]=0.0440520084659;
        wi[11]=0.0435492880829;
        wi[12]=0.0435492880829;
        wi[13]=0.0429575990023;
        wi[14]=0.0429575990023;
        wi[15]=0.0422781500113;
        wi[16]=0.0422781500113;
        wi[17]=0.0415123291857;
        wi[18]=0.0415123291857;
        wi[19]=0.040661701054;
        wi[20]=0.040661701054;
        wi[21]=0.039728003402;
        wi[22]=0.039728003402;
        wi[23]=0.038713143720;
        wi[24]=0.038713143720;
        wi[25]=0.037619195312;
        wi[26]=0.037619195312;
        wi[27]=0.036448393051;
        wi[28]=0.036448393051;
        wi[29]=0.035203128822;
        wi[30]=0.035203128822;
        wi[31]=0.033885946631;
        wi[32]=0.033885946631;
        wi[33]=0.032499537410;
        wi[34]=0.032499537410;
        wi[35]=0.031046733518;
        wi[36]=0.031046733518;
        wi[37]=0.029530502958;
        wi[38]=0.029530502958;
        wi[39]=0.027953943312;
        wi[40]=0.027953943312;
        wi[41]=0.026320275417;
        wi[42]=0.026320275417;
        wi[43]=0.024632836785;
        wi[44]=0.024632836785;
        wi[45]=0.022895074791;
        wi[46]=0.022895074791;
        wi[47]=0.021110539640;
        wi[48]=0.021110539640;
        wi[49]=0.019282877129;
        wi[50]=0.019282877129;
        wi[51]=0.017415821232;
        wi[52]=0.017415821232;
        wi[53]=0.015513186543;
        wi[54]=0.015513186543;
        wi[55]=0.013578860649;
        wi[56]=0.013578860649;
        wi[57]=0.011616796611;
        wi[58]=0.011616796611;
        wi[59]=0.009631006150;
        wi[60]=0.00963100615;
        wi[61]=0.007625555931;
        wi[62]=0.00762555593;
        wi[63]=0.005604579928;
        wi[64]=0.00560457993;
        wi[65]=0.00357241674;
        wi[66]=0.00357241674;
        wi[67]=0.00153597695;
        wi[68]=0.00153597695;
        
        xi[0]=0;
        xi[1]=-0.0451862106151;
        xi[2]=0.0451862106151;
        xi[3]=-0.0902801083928;
        xi[4]=0.0902801083928;
        xi[5]=-0.135189569086;
        xi[6]=0.135189569086;
        xi[7]=-0.179822845241;
        xi[8]=0.179822845241;
        xi[9]=-0.224088753634;
        xi[10]=0.224088753634;
        xi[11]=-0.267896861554;
        xi[12]=0.267896861554;
        xi[13]=-0.311157671547;
        xi[14]=0.311157671547;
        xi[15]=-0.353782804262;
        xi[16]=0.353782804262;
        xi[17]=-0.395685178997;
        xi[18]=0.395685178997;
        xi[19]=-0.436779191606;
        xi[20]=0.436779191606;
        xi[21]=-0.476980889379;
        xi[22]=0.476980889379;
        xi[23]=-0.516208142557;
        xi[24]=0.516208142557;
        xi[25]=-0.554380812115;
        xi[26]=0.554380812115;
        xi[27]=-0.591420913485;
        xi[28]=0.591420913485;
        xi[29]=-0.627252775870;
        xi[30]=0.627252775870;
        xi[31]=-0.661803196841;
        xi[32]=0.661803196841;
        xi[33]=-0.695001591883;
        xi[34]=0.695001591883;
        xi[35]=-0.726780138596;
        xi[36]=0.726780138596;
        xi[37]=-0.757073915256;
        xi[38]=0.757073915256;
        xi[39]=-0.785821033447;
        xi[40]=0.785821033447;
        xi[41]=-0.812962764500;
        xi[42]=0.812962764500;
        xi[43]=-0.838443659481;
        xi[44]=0.838443659481;
        xi[45]=-0.862211662473;
        xi[46]=0.862211662473;
        xi[47]=-0.884218216947;
        xi[48]=0.884218216947;
        xi[49]=-0.904418364986;
        xi[50]=0.904418364986;
        xi[51]=-0.922770839180;
        xi[52]=0.922770839180;
        xi[53]=-0.939238147042;
        xi[54]=0.939238147042;
        xi[55]=-0.953786647805;
        xi[56]=0.953786647805;
        xi[57]=-0.966386621643;
        xi[58]=0.966386621643;
        xi[59]=-0.977012331622;
        xi[60]=0.977012331622;
        xi[61]=-0.985642080023;
        xi[62]=0.985642080023;
        xi[63]=-0.992258266717;
        xi[64]=0.992258266717;
        xi[65]=-0.996847500485;
        xi[66]=0.996847500485;
        xi[67]=-0.999401427543;
        xi[68]=0.999401427543;
    }
    else if(no_i==70)
    {
        wi[0]=0.0445494171598;
        wi[1]=0.0445494171598;
        wi[2]=0.0444609684172;
        wi[3]=0.0444609684172;
        wi[4]=0.0442842465391;
        wi[5]=0.0442842465391;
        wi[6]=0.0440196023902;
        wi[7]=0.0440196023902;
        wi[8]=0.0436675613972;
        wi[9]=0.0436675613972;
        wi[10]=0.0432288225051;
        wi[11]=0.0432288225051;
        wi[12]=0.0427042567894;
        wi[13]=0.0427042567894;
        wi[14]=0.0420949057273;
        wi[15]=0.0420949057273;
        wi[16]=0.0414019791290;
        wi[17]=0.0414019791290;
        wi[18]=0.0406268527368;
        wi[19]=0.0406268527368;
        wi[20]=0.0397710654928;
        wi[21]=0.0397710654928;
        wi[22]=0.038836316484;
        wi[23]=0.038836316484;
        wi[24]=0.037824461569;
        wi[25]=0.037824461569;
        wi[26]=0.036737509694;
        wi[27]=0.036737509694;
        wi[28]=0.035577618901;
        wi[29]=0.035577618901;
        wi[30]=0.034347092050;
        wi[31]=0.034347092050;
        wi[32]=0.033048372239;
        wi[33]=0.033048372239;
        wi[34]=0.031684037961;
        wi[35]=0.031684037961;
        wi[36]=0.030256797980;
        wi[37]=0.030256797980;
        wi[38]=0.028769485956;
        wi[39]=0.028769485956;
        wi[40]=0.027225054819;
        wi[41]=0.027225054819;
        wi[42]=0.025626570908;
        wi[43]=0.025626570908;
        wi[44]=0.023977207889;
        wi[45]=0.023977207889;
        wi[46]=0.022280240452;
        wi[47]=0.022280240452;
        wi[48]=0.020539037824;
        wi[49]=0.020539037824;
        wi[50]=0.018757057093;
        wi[51]=0.018757057093;
        wi[52]=0.016937836376;
        wi[53]=0.016937836376;
        wi[54]=0.015084987865;
        wi[55]=0.015084987865;
        wi[56]=0.013202190815;
        wi[57]=0.013202190815;
        wi[58]=0.011293184650;
        wi[59]=0.011293184650;
        wi[60]=0.009361762770;
        wi[61]=0.00936176277;
        wi[62]=0.007411769363;
        wi[63]=0.00741176936;
        wi[64]=0.005447111874;
        wi[65]=0.00544711187;
        wi[66]=0.00347189489;
        wi[67]=0.00347189489;
        wi[68]=0.00149272129;
        wi[69]=0.00149272129;
        
        xi[0]=-0.0222783952861;
        xi[1]=0.0222783952861;
        xi[2]=-0.0667909541676;
        xi[3]=0.0667909541676;
        xi[4]=-0.111170905794;
        xi[5]=0.111170905794;
        xi[6]=-0.155330137882;
        xi[7]=0.155330137882;
        xi[8]=-0.199180976365;
        xi[9]=0.199180976365;
        xi[10]=-0.242636359464;
        xi[11]=0.242636359464;
        xi[12]=-0.285610010540;
        xi[13]=0.285610010540;
        xi[14]=-0.328016609390;
        xi[15]=0.328016609390;
        xi[16]=-0.369771961638;
        xi[17]=0.369771961638;
        xi[18]=-0.410793165903;
        xi[19]=0.410793165903;
        xi[20]=-0.450998778382;
        xi[21]=0.450998778382;
        xi[22]=-0.490308974558;
        xi[23]=0.490308974558;
        xi[24]=-0.528645707680;
        xi[25]=0.528645707680;
        xi[26]=-0.565932863719;
        xi[27]=0.565932863719;
        xi[28]=-0.602096412485;
        xi[29]=0.602096412485;
        xi[30]=-0.637064554610;
        xi[31]=0.637064554610;
        xi[32]=-0.670767864094;
        xi[33]=0.670767864094;
        xi[34]=-0.703139426152;
        xi[35]=0.703139426152;
        xi[36]=-0.734114970061;
        xi[37]=0.734114970061;
        xi[38]=-0.763632996772;
        xi[39]=0.763632996772;
        xi[40]=-0.791634901008;
        xi[41]=0.791634901008;
        xi[42]=-0.818065087625;
        xi[43]=0.818065087625;
        xi[44]=-0.842871081999;
        xi[45]=0.842871081999;
        xi[46]=-0.866003634214;
        xi[47]=0.866003634214;
        xi[48]=-0.887416816863;
        xi[49]=0.887416816863;
        xi[50]=-0.907068116261;
        xi[51]=0.907068116261;
        xi[52]=-0.924918516898;
        xi[53]=0.924918516898;
        xi[54]=-0.940932579004;
        xi[55]=0.940932579004;
        xi[56]=-0.955078509114;
        xi[57]=0.955078509114;
        xi[58]=-0.967328223665;
        xi[59]=0.967328223665;
        xi[60]=-0.977657405958;
        xi[61]=0.977657405958;
        xi[62]=-0.986045558070;
        xi[63]=0.986045558070;
        xi[64]=-0.992476055212;
        xi[65]=0.992476055212;
        xi[66]=-0.996936251962;
        xi[67]=0.996936251962;
        xi[68]=-0.999418285974;
        xi[69]=0.999418285974;
    }

    
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_lg_0::create_2nd_neigh_lst()
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
type0 ForceField_eam_dmd_lg_0::
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
type0 ForceField_eam_dmd_lg_0::
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
type0 ForceField_eam_dmd_lg_0::
dmat1(type0 fi,type0 crdi,type0 fj,type0 crdj
      ,int itype)
{
    
    return 0.0;
}
/*--------------------------------------------
 return M_{ij}^{\alpha}
 --------------------------------------------*/
type0 ForceField_eam_dmd_lg_0::calc_ent(type0 x)
{
    if(x<1.0e-8 || 0.99999999<x)
        return 0.0;
    else
        return x*log(x)+(1-x)*log(1-x);
    
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_lg_0::c_d_calc()
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
    type0* mu;
    atoms->vectors[mu_n].ret(mu);
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
            
            mu[icomp+itype]=E[icomp+itype]=tmp0;
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
                    mu[icomp+itype]+=c[jcomp+jtype]*(rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm+jtype]
                                                    +rho_phi[istart+type2phi_pair_ij[itype][jtype]]);
                    
                    if(jatm<natms)
                        mu[jcomp+jtype]+=c[icomp+itype]*(rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm+itype]
                                                        +rho_phi[istart+type2phi_pair_ji[jtype][itype]]);
                }
            }
            istart+=stride;
        }
        
        for(int itype=0;itype<no_types;itype++)
        {
            mu[icomp+itype]+=c_0[itype];
            //fix this part for change of variables
            mu[icomp+itype]+=1.5*kbT*log(x[(3+no_types)*iatm+3+itype]);
        }
        
    }
    
    atoms->update(mu_n);
    
    
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
                fi=mu[icomp+itype];
                fj=mu[jcomp+itype];
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
type0 ForceField_eam_dmd_lg_0::g_calc(int chk
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
    type0* mu;
    atoms->vectors[mu_n].ret(mu);
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
            
            mu[icomp+itype]=E[icomp+itype]=tmp0;
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
                    mu[icomp+itype]+=c[jcomp+jtype]*(rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jatm+jtype]
                                                    +rho_phi[istart+type2phi_pair_ij[itype][jtype]]);
                    if(jatm<natms)
                        mu[jcomp+jtype]+=c[icomp+itype]*(rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[iatm+itype]
                                                        +rho_phi[istart+type2phi_pair_ji[jtype][itype]]);
                }
            }
            istart+=stride;
        }
        
        for(int itype=0;itype<no_types;itype++)
        {
            mu[icomp+itype]+=c_0[itype];
            //fix this part for change of variables
            mu[icomp+itype]+=1.5*kbT*log(x[(3+no_types)*iatm+3+itype]);
        }
        
    }
    
    atoms->update(mu_n);
    
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
                fi=mu[icomp+itype];
                fj=mu[jcomp+itype];
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
                fi=mu[icomp+itype];
                fj=mu[jcomp+itype];
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
