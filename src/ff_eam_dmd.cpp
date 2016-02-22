#include <stdlib.h>
#include "ff_eam_dmd.h"
#include "neighbor.h"
#include "error.h"
#include "memory.h"
#include "atom_types.h"
#include "xmath.h"
#include <limits>
#define PI_IN_SQ 0.564189583547756286948079451561
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam_dmd::
ForceField_eam_dmd(MAPP* mapp): ForceFieldDMD(mapp)
{
    if(mapp->mode!=DMD_mode)
        error->abort("ff eam_dmd works only "
        "for dmd mode");

    eam_reader=NULL;
    no_types=0;
    no_i=0;
    max_pairs=0;
    phi_psi_sz_sz=0;
    n_phi_psi=0;
    M_N_sz_sz=0;
    min_log=log(std::numeric_limits<type0>::min());
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
    
    type0* x=mapp->x->begin();
    dmd_type* type=mapp->ctype->begin();
    type0* fvec=f->begin();
    type0* E=E_ptr->begin();
    type0* dE=dE_ptr->begin();
    type0* mu=mu_ptr->begin();;
    type0* c=mapp->c->begin();
    
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
    
    int natms=atoms->natms;
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
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            if(rsq<cut_sq_mod_0)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                
                for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && c[ic_dim]>=0.0;ic_dim++)
                {
                    itype=type[ic_dim];
                    for(int jc_dim=jatm*c_dim,jtype;jc_dim<(jatm+1)*c_dim && c[jc_dim]>=0.0;jc_dim++)
                    {
                        jtype=type[jc_dim];
                        
                        if(rsq>cut_sq[COMP(itype,jtype)])
                            continue;
                        
                        type0 alpha_ij=sqrt(x[3*iatm+ic_dim+3]*x[3*iatm+ic_dim+3]+x[3*jatm+jc_dim+3]*x[3*jatm+jc_dim+3]);
                        
                        if(alpha_ij>=alpha_max)
                            continue;
                        
                        type0 upper=(r+cut[COMP(itype,jtype)])/alpha_ij;
                        type0 lower=(r-cut[COMP(itype,jtype)])/alpha_ij;
                        
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
                        
                        
                        E[ic_dim]+=c[jc_dim]*rho_phi[rho_ji];
                        
                        if(jatm<natms)
                        {
                            E[jc_dim]+=c[ic_dim]*rho_phi[rho_ij];
                            nrgy_strss_lcl[0]+=c[ic_dim]*c[jc_dim]*rho_phi[phi_ij];
                        }
                        else
                            nrgy_strss_lcl[0]+=0.5*c[ic_dim]*c[jc_dim]*rho_phi[phi_ij];
                    }
                }
            }
            
            
            istart+=stride;
        }
        
        c_iv=1.0;
        for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && c[ic_dim]>=0.0;ic_dim++)
        {
            
            itype=type[ic_dim];
            c_iv-=c[ic_dim];
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
            
            nrgy_strss_lcl[0]+=c[ic_dim]*(tmp0+c_0[itype]-3.0*kbT*log(x[3*iatm+ic_dim+3]));
            mu[ic_dim]=tmp0-3.0*kbT*log(x[3*iatm+ic_dim+3]);
            nrgy_strss_lcl[0]+=kbT*calc_ent(c[ic_dim]);
        }
        nrgy_strss_lcl[0]+=kbT*calc_ent(c_iv);
    }
    
    atoms->update(dE_ptr);

    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+c_dim)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+c_dim)*jatm;
            
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            r_inv=1.0/sqrt(dx0*dx0+dx1*dx1+dx2*dx2);
            
            for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && c[ic_dim]>=0.0;ic_dim++)
            {
                itype=type[ic_dim];
                for(int jc_dim=jatm*c_dim,jtype;jc_dim<(jatm+1)*c_dim && c[jc_dim]>=0.0;jc_dim++)
                {
                    jtype=type[jc_dim];
                    
                    fpair=-(drho_phi_dr[istart+type2rho_pair_ji[jtype][itype]]*dE[ic_dim]
                    +drho_phi_dr[istart+type2rho_pair_ij[itype][jtype]]*dE[jc_dim]
                    +drho_phi_dr[istart+type2phi_pair_ij[itype][jtype]])*c[ic_dim]*c[jc_dim];
                    
                    apair=-(drho_phi_dalpha[istart+type2rho_pair_ji[jtype][itype]]*dE[ic_dim]
                    +drho_phi_dalpha[istart+type2rho_pair_ij[itype][jtype]]*dE[jc_dim]
                    +drho_phi_dalpha[istart+type2phi_pair_ij[itype][jtype]])*c[ic_dim]*c[jc_dim];
                    
                    mu[ic_dim]+=c[jc_dim]*(rho_phi[istart+type2phi_pair_ij[itype][jtype]]+
                    rho_phi[istart+type2rho_pair_ij[itype][jtype]]*dE[jc_dim]);
                    
                    if(jatm<natms)
                        mu[jc_dim]+=c[ic_dim]*(rho_phi[istart+type2phi_pair_ji[jtype][itype]]+
                        rho_phi[istart+type2rho_pair_ji[jtype][itype]]*dE[ic_dim]);
                    
                    if(apair==0.0 && fpair==0.0)
                        continue;
                
                    
                    fvec[icomp]+=dx0*fpair;
                    fvec[icomp+1]+=dx1*fpair;
                    fvec[icomp+2]+=dx2*fpair;
                    fvec[3*iatm+ic_dim+3]+=apair*x[3*iatm+ic_dim+3];
                    if(jatm<natms)
                    {
                        fvec[jcomp]-=dx0*fpair;
                        fvec[jcomp+1]-=dx1*fpair;
                        fvec[jcomp+2]-=dx2*fpair;
                        fvec[3*jatm+jc_dim+3]+=apair*x[3*jatm+jc_dim+3];
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
    
    type0* alpha=x+dim;
    type0* fvec_alpha=fvec+dim;
    for(iatm=0;iatm<natms;iatm++)
    {
        for(int i=0;i<c_dim;i++)
            if(c[i]>=0.0)
                fvec_alpha[i]+=3.0*kbT*c[i]/alpha[i];
        
        fvec_alpha+=x_dim;
        alpha+=x_dim;
        c+=c_dim;
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
    type0* x=mapp->x->begin();
    dmd_type* type=mapp->ctype->begin();
    type0* E=E_ptr->begin();
    type0* c=mapp->c->begin();
    
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
    
    int natms=atoms->natms;
    for(int i=0;i<natms*c_dim;i++) E[i]=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+c_dim)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+c_dim)*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            if(rsq<cut_sq_mod_0)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                
                for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && c[ic_dim]>=0.0;ic_dim++)
                {
                    itype=type[ic_dim];
                    for(int jc_dim=jatm*c_dim,jtype;jc_dim<(jatm+1)*c_dim && c[jc_dim]>=0.0;jc_dim++)
                    {
                        jtype=type[jc_dim];
                        
                        if(rsq>cut_sq[COMP(itype,jtype)])
                            continue;
                        
                        alpha=sqrt(x[3*iatm+ic_dim+3]*x[3*iatm+ic_dim+3]+x[3*jatm+jc_dim+3]*x[3*jatm+jc_dim+3]);
                        
                        if(alpha>=alpha_max)
                            continue;
                        
                        upper=(r+cut[COMP(itype,jtype)])/alpha;
                        lower=(r-cut[COMP(itype,jtype)])/alpha;
                        
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
                        
                        E[ic_dim]+=c[jc_dim]*rho_jt_it_0;
                        
                        if(jatm<natms)
                        {
                            E[jc_dim]+=c[ic_dim]*rho_it_jt_0;
                            en+=c[ic_dim]*c[jc_dim]*phi_it_jt_0;
                        }
                        else
                            en+=0.5*c[ic_dim]*c[jc_dim]*phi_it_jt_0;
                    }
                }
            }
        }
        
        c_iv=1.0;
        for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && c[ic_dim]>=0.0;ic_dim++)
        {
            itype=type[ic_dim];
            c_iv-=c[ic_dim];
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

            en+=c[ic_dim]*(tmp0+c_0[itype])
            -3.0*kbT*c[ic_dim]*log(x[3*iatm+ic_dim+3]);
            
            en+=kbT*calc_ent(c[ic_dim]);
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
    neighbor->pair_wise=true;
    
    if(kbT==-1.0)
        error->abort("temperature for the siumulation has not been set");
    
    
    cv_ptr=new Vec<type0>(atoms,1);
    E_ptr=new Vec<type0>(atoms,c_dim);
    dE_ptr=new Vec<type0>(atoms,c_dim);
    mu_ptr=new Vec<type0>(atoms,no_types,"mu");
    mu_ptr->change_dimension(0.0,no_types,no_types-c_dim);

    crd_ptr=new Vec<type0>(atoms,c_dim);
    s_ptr=new Vec<type0>(atoms,c_dim);
    
    M_ptr=new Vec<type0>(atoms,c_dim);
    N_ptr=new Vec<type0>(atoms,c_dim);

    t_ptr=new Vec<type0>(atoms,c_dim);
    
    type0* mu=mu_ptr->begin();
    for(int i=0;i<c_dim*atoms->natms;i++)
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
        delete [] N_IJ;
        delete [] N_JI;
        delete [] M_IJ;
        delete [] M_JI;
        M_N_sz_sz=0;
    }
    
    delete t_ptr;
    delete s_ptr;
    delete crd_ptr;
    delete mu_ptr;
    delete dE_ptr;
    delete E_ptr;
    delete cv_ptr;
    delete M_ptr;
    delete N_ptr;
}
/*--------------------------------------------
 ff_coef
 --------------------------------------------*/
void ForceField_eam_dmd::coef(int nargs,char** args)
{
    if(nargs!=2)
        error->abort("wrong coeff command "
        "for ff dmd");
    
    cut_off_alloc();
    allocate();

    no_types=atom_types->no_types;
    c_dim=mapp->c->dim;
    x_dim=mapp->x->dim;
    dim=atoms->dimension;
    
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
    CREATE1D(f_t,no_types);
    CREATE1D(Q_nrm,no_types);
}
/*-------------------------------------------
 
 --------------------------------------------*/
void ForceField_eam_dmd::deallocate()
{
    if(!no_types)
        return;
    
    delete [] c_0;
    delete [] f_t;
    delete [] Q_nrm;
    no_types=0;

}
/*-------------------------------------------

 --------------------------------------------*/
void ForceField_eam_dmd::calc_cv(type0* c,type0* cv,int n)
{
    for(int i=0;i<n;i++)
    {
        cv[i]=1.0;
        for(int j=0;j<c_dim;j++)
            if(c[i*c_dim+j]>=0.0)
                cv[i]-=c[i*c_dim+j];
    }
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_eam_dmd::read_file(char* file_name)
{
    int nargs;
    int tot_no_types;
    int icmp;
    type0 tmp0,alpha_coef=0.0;
    int tmp1;
    
    FILE* fp=NULL;
    char** args=NULL;
    int args_cpcty=0;
    int* type_ref;
    char* line;
    char* strtmp;
    
    delete eam_reader;
    eam_reader=new EAMFileReader(mapp);

    mapp->open_file(fp,file_name,"r");
    
    CREATE1D(line,MAXCHAR);
    CREATE1D(strtmp,MAXCHAR);
    nargs=0;
    while(nargs==0 && mapp->read_line(fp,line) !=-1)
        nargs=mapp->parse_line(line,args,args_cpcty);
    
    if(nargs==0)
        error->abort("%s file ended immaturely",file_name);

    if(nargs<no_types)
        error->abort("the number of atoms in %s file"
        " is less than the number of atom types present in the system",file_name);
    
    tot_no_types=nargs;
    
    CREATE1D(type_ref,tot_no_types);
    for(int ityp=0;ityp<tot_no_types;ityp++)
        type_ref[ityp]=atom_types->find_type_exist(args[ityp]);
    
    if(args_cpcty)
        delete [] args;
    
    for(int i=0;i<no_types;i++)
        rsq_crd[i]=f_t[i]=Q_nrm[i]=-1.0;
    
    alpha_min=alpha_max=0.0;
    kbT=beta=-1.0;

    while(mapp->read_line(fp,line)!=-1)
    {
        if(mapp->hash_remover(line)==0)
            continue;
        
        if(sscanf(line,"degree = %d",&tmp1)==1)
        {
            alpha_coef=set_weight_abs(tmp1);
        }
        else if(sscanf(line,"alpha_max = %lf",&tmp0)==1)
        {
            if(tmp0<=0.0)
                error->abort("maximum alpha for ff eam_dmd "
                "should be greater than 0.0");
            alpha_max=tmp0;
        }
        else if(sscanf(line,"T = %lf",&tmp0)==1)
        {
            if(tmp0<=0.0)
                error->abort("temperature for ff eam_dmd "
                "should be greater than 0.0");
            set_temp(tmp0);
        }
        else if(sscanf(line,"r_crd(%d) = %lf",&icmp,&tmp0)==2)
        {
            if(icmp<0 || icmp>tot_no_types-1)
                error->abort("wrong component in %s "
                "file for r_crd(%d)",file_name,icmp);
            
            if(tmp0<=0.0)
                error->abort("r_crd(%d) for ff eam_dmd "
                "in %s file should be greater than 0.0",icmp,file_name);
            
            if(type_ref[icmp]!=-1)
            {
                rsq_crd[type_ref[icmp]]=tmp0*tmp0;
            }
        }
        else if(sscanf(line,"f_t(%d) = %lf",&icmp,&tmp0)==2)
        {
            if(icmp<0 || icmp>tot_no_types-1)
                error->abort("wrong component in %s "
                "file for f_t(%d)",file_name,icmp);
            
            if(tmp0<=0.0)
                error->abort("f_t(%d) for ff eam_dmd "
                "in %s file should be greater than 0.0",icmp,file_name);
            
            if(type_ref[icmp]!=-1)
            {
                f_t[type_ref[icmp]]=tmp0;
            }
        }
        else if(sscanf(line,"Q(%d) = %lf",&icmp,&tmp0)==2)
        {
            if(icmp<0 || icmp>tot_no_types-1)
                error->abort("wrong component in %s "
                "file for Q(%d)",file_name,icmp);
            
            if(tmp0<=0.0)
                error->abort("Q(%d) for ff eam_dmd "
                "in %s file should be greater than 0.0",icmp,file_name);
            
            if(type_ref[icmp]!=-1)
            {
                Q_nrm[type_ref[icmp]]=tmp0;
            }
        }
        else if(sscanf(line,"file_format = %s",strtmp)==1)
        {
            eam_reader->file_format(strtmp);
        }
        else if(sscanf(line,"file(%d) = %s",&icmp,strtmp)==2)
        {
            if(icmp<0 || icmp>tot_no_types-1)
                error->abort("wrong component in %s "
                "file for file(%d)",file_name,icmp);
            if(type_ref[icmp]!=-1)
                eam_reader->add_file(strtmp,type_ref[icmp]);
        }
        else
        {
            error->abort("invalid line in %s file: %s",file_name,line);
        }
    }

    delete [] line;
    delete [] strtmp;
    
    if(tot_no_types)
        delete [] type_ref;
    
    
    if(no_i==0)
        error->abort("degree was not set by %s file ",file_name);
   
    if(alpha_max==0.0)
        error->abort("miximum alpha was not set by %s file ",file_name);
    
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
    
    for(int i=0;i<no_types;i++)
        Q_nrm[i]/=f_t[i]*sqrt(f_t[i]);
    
    if(alpha_max<=alpha_min)
        error->abort("maximum alpha should be less "
        "than minimum alpha for ff eam_dmd");
    
    type0* x=mapp->x->begin();
    int* id=mapp->id->begin();
    dmd_type* type=mapp->ctype->begin();
    type0 alpha_bound=alpha_max/sqrt(2.0);
    for(int iatm=0;iatm<atoms->natms;iatm++)
    {
        for(int ic_dim=iatm*c_dim;ic_dim<(iatm+1)*c_dim;ic_dim++)
            if(x[3*iatm+ic_dim+3]>=alpha_bound)
                error->abort("alpha_%s of atom %d is greater than "
                "the maximum allowed %lf (alpha_max/sqrt(2.0)) set "
                "by %s file",atom_types->atom_names[type[ic_dim]]
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
        for(int jtype=itype;jtype<no_types;jtype++)
        {
            type0 tmp=sqrt(eam_reader->cut_sq[COMP(itype,jtype)]);
            cut_max=MAX(cut_max,tmp);
            cut[COMP(itype,jtype)]=tmp;
            cut_sq[COMP(itype,jtype)]=(tmp+xi[no_i-1]*alpha_max)*(tmp+xi[no_i-1]*alpha_max);
        }
    
    type0 mod_rc=cut_max+xi[no_i-1]*alpha_max;
    cut_sq_mod_0=mod_rc*mod_rc;
    
    for(int itype=0;itype<no_types;itype++)
        if(rsq_crd[itype]>cut_sq[COMP(itype,itype)])
            error->abort("r_crd(%s) set by file %s should be less than %lf"
            ,file_name,atom_types->atom_names[itype],sqrt(cut_sq[COMP(itype,itype)]));
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
    
    XMath* xmath=new XMath();
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
 return M_{ij}^{\alpha}
 this fucntion should be symmetric wrt fi &
 fj
 --------------------------------------------*/
type0 ForceField_eam_dmd::
mat(type0 fi,type0 fj,int itype)
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
type0 ForceField_eam_dmd::
dmat(type0 fi,type0 fj,int itype)
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
 return M_{ij}^{\alpha}
 --------------------------------------------*/
inline type0 ForceField_eam_dmd::mod_log(type0 x)
{
    if(x==0.0)
        return min_log;
    return log(x);
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd::dc()
{
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    dmd_type* type=mapp->ctype->begin();
    type0* cv=cv_ptr->begin();
    type0* mu=mu_ptr->begin();
    type0 fi,fj,exp_fi,exp_fj;
    type0 tmp0,s_ij,w_ij,w_ji;
    calc_mu();
    
    int natms=atoms->natms;
    for(int i=0;i<natms*c_dim;i++)
        c_d[i]=0.0;
    
    int** neighbor_list_2nd=neighbor->neighbor_list_2nd;
    int* neighbor_list_size_2nd=neighbor->neighbor_list_size_2nd;
    
    for(int ic_dim=0,itype;ic_dim<natms*c_dim;ic_dim++)
    {
        itype=type[ic_dim];
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            jc_dim=neighbor_list_2nd[ic_dim][j];
            
            fi=mu[ic_dim];
            fj=mu[jc_dim];
            
            tmp0=mat(fi,fj,itype);
            exp_fi=exp(beta*(fi-tmp0));
            exp_fj=exp(beta*(fj-tmp0));
            w_ij=-c[ic_dim]*cv[jc_dim/c_dim]*exp_fi;
            w_ji=-c[jc_dim]*cv[ic_dim/c_dim]*exp_fj;
            
            s_ij=w_ij-w_ji;
            
            c_d[ic_dim]+=s_ij;
            if(jc_dim<natms*c_dim)
            {
                c_d[jc_dim]-=s_ij;
            }
        }
    }
}
/*--------------------------------------------
 calculate norm of d^2c/dt^2
 --------------------------------------------*/
type0 ForceField_eam_dmd::ddc_norm()
{
    type0* c=mapp->c->begin();
    dmd_type* type=mapp->ctype->begin();
    type0* dE=dE_ptr->begin();
    type0* ddE=E_ptr->begin();;
    type0* cv=cv_ptr->begin();
    type0* s=s_ptr->begin();
    type0* c_d=mapp->c_d->begin();
    type0* mu=mu_ptr->begin();
    type0* t=t_ptr->begin();
    
    type0 fi,fj,exp_fi,exp_fj,norm_lcl;
    type0 tmp0,s_ij,w_ij,w_ji;
    int istart;
    
    
    int natms=atoms->natms;
    
    calc_mu();
    
    for(int i=0;i<natms*c_dim;i++)
        s[i]=c_d[i]=0.0;
    
    int** neighbor_list_2nd=neighbor->neighbor_list_2nd;
    int* neighbor_list_size_2nd=neighbor->neighbor_list_size_2nd;

    for(int ic_dim=0,itype;ic_dim<natms*c_dim;ic_dim++)
    {
        itype=type[ic_dim];
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            jc_dim=neighbor_list_2nd[ic_dim][j];
            
            fi=mu[ic_dim];
            fj=mu[jc_dim];
            
            tmp0=mat(fi,fj,itype);
            exp_fi=exp(beta*(fi-tmp0));
            exp_fj=exp(beta*(fj-tmp0));
            w_ij=-c[ic_dim]*cv[jc_dim/c_dim]*exp_fi;
            w_ji=-c[jc_dim]*cv[ic_dim/c_dim]*exp_fj;
            
            s_ij=w_ij-w_ji;
            
            s[ic_dim]+=s_ij;
            c_d[ic_dim]+=s_ij;
            if(jc_dim<natms*c_dim)
            {
                s[jc_dim]-=s_ij;
                c_d[jc_dim]-=s_ij;
            }
        }
    }

    atoms->update(s_ptr);
    
    for(int i=0;i<natms*c_dim;i++)
        t[i]=0.0;
    
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            t[ic_dim]+=c[ic_dim]*ddE[ic_dim]*psi_JI[istart]*s[jc_dim];
            if(jc_dim<phi_psi_sz_sz)
                t[jc_dim]+=c[jc_dim]*ddE[jc_dim]*psi_IJ[istart]*s[ic_dim];
            istart++;
        }

    
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        t[ic_dim]*=c[ic_dim]*ddE[ic_dim];

    
    atoms->update(t_ptr);
    
    
    for(int i=0;i<natms*c_dim;i++)
        ddE[i]=0.0;
    
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            tmp0=psi_IJ[istart]*dE[ic_dim]+psi_JI[istart]*dE[jc_dim]+phi_IJ[istart];
            
            ddE[ic_dim]+=psi_IJ[istart]*t[jc_dim]+tmp0*s[jc_dim];
            if(jc_dim<phi_psi_sz_sz)
                ddE[jc_dim]+=psi_JI[istart]*t[ic_dim]+tmp0*s[ic_dim];
            
            istart++;
        }
    
    atoms->update(E_ptr);

    for(int i=0;i<natms*c_dim;i++)
        dE[i]=0.0;
    
    norm_lcl=0.0;

    for(int ic_dim=0,itype;ic_dim<natms*c_dim;ic_dim++)
    {
        itype=type[ic_dim];
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            jc_dim=neighbor_list_2nd[ic_dim][j];
            
            fi=mu[ic_dim];
            fj=mu[jc_dim];
            tmp0=mat(fi,fj,itype);
            exp_fi=exp(beta*(fi-tmp0));
            exp_fj=exp(beta*(fj-tmp0));
            w_ij=-c[ic_dim]*cv[jc_dim/c_dim]*exp_fi;
            w_ji=-c[jc_dim]*cv[ic_dim/c_dim]*exp_fj;
            s_ij=w_ij-w_ji;
            
            tmp0=(s[jc_dim]*c[ic_dim]-s[ic_dim]*cv[jc_dim/c_dim])*exp_fi
            -(s[ic_dim]*c[jc_dim]-s[jc_dim]*cv[ic_dim/c_dim])*exp_fj;
            
            tmp0+=beta*(ddE[ic_dim]*w_ij-ddE[jc_dim]*w_ji);
            
            tmp0+=-beta*(w_ij-w_ji)*(
            ddE[ic_dim]*dmat(fi,fj,itype)
            +ddE[jc_dim]*dmat(fj,fi,itype));
            
            dE[ic_dim]+=tmp0;
            if(jc_dim<natms*c_dim)
                dE[jc_dim]-=tmp0;
        }
        norm_lcl+=dE[ic_dim]*dE[ic_dim];
    }
    
    MPI_Allreduce(&norm_lcl,&tmp0,1,MPI_TYPE0,MPI_SUM,world);
    return sqrt(tmp0);
}
/*--------------------------------------------
 calculate norm of d^2c/dt^2
 --------------------------------------------*/
void ForceField_eam_dmd::ddc(type0* ddc_)
{
    type0* c=mapp->c->begin();
    dmd_type* type=mapp->ctype->begin();
    type0* dE=dE_ptr->begin();
    type0* ddE=E_ptr->begin();;
    type0* cv=cv_ptr->begin();
    type0* s=s_ptr->begin();
    type0* mu=mu_ptr->begin();
    type0* t=t_ptr->begin();
    
    type0 fi,fj,exp_fi,exp_fj;
    type0 tmp0,s_ij,w_ij,w_ji;
    int istart;
    
    
    int natms=atoms->natms;
    
    calc_mu();
    
    for(int i=0;i<natms*c_dim;i++)
        s[i]=0.0;
    
    int** neighbor_list_2nd=neighbor->neighbor_list_2nd;
    int* neighbor_list_size_2nd=neighbor->neighbor_list_size_2nd;
    
    for(int ic_dim=0,itype;ic_dim<natms*c_dim;ic_dim++)
    {
        itype=type[ic_dim];
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            jc_dim=neighbor_list_2nd[ic_dim][j];
            
            fi=mu[ic_dim];
            fj=mu[jc_dim];
            
            tmp0=mat(fi,fj,itype);
            exp_fi=exp(beta*(fi-tmp0));
            exp_fj=exp(beta*(fj-tmp0));
            w_ij=-c[ic_dim]*cv[jc_dim/c_dim]*exp_fi;
            w_ji=-c[jc_dim]*cv[ic_dim/c_dim]*exp_fj;
            
            s_ij=w_ij-w_ji;
            
            s[ic_dim]+=s_ij;
            if(jc_dim<natms*c_dim)
            {
                s[jc_dim]-=s_ij;
            }
        }
    }
    
    atoms->update(s_ptr);
    
    for(int i=0;i<natms*c_dim;i++)
        t[i]=0.0;
    
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            t[ic_dim]+=c[ic_dim]*ddE[ic_dim]*psi_JI[istart]*s[jc_dim];
            if(jc_dim<phi_psi_sz_sz)
                t[jc_dim]+=c[jc_dim]*ddE[jc_dim]*psi_IJ[istart]*s[ic_dim];
            istart++;
        }
    
    
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        t[ic_dim]*=c[ic_dim]*ddE[ic_dim];
    
    
    atoms->update(t_ptr);
    
    
    for(int i=0;i<natms*c_dim;i++)
        ddE[i]=0.0;
    
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            tmp0=psi_IJ[istart]*dE[ic_dim]+psi_JI[istart]*dE[jc_dim]+phi_IJ[istart];
            
            ddE[ic_dim]+=psi_IJ[istart]*t[jc_dim]+tmp0*s[jc_dim];
            if(jc_dim<phi_psi_sz_sz)
                ddE[jc_dim]+=psi_JI[istart]*t[ic_dim]+tmp0*s[ic_dim];
            
            istart++;
        }
    
    atoms->update(E_ptr);
    
    for(int i=0;i<natms*c_dim;i++)
        ddc_[i]=0.0;
    
    
    for(int ic_dim=0,itype;ic_dim<natms*c_dim;ic_dim++)
    {
        itype=type[ic_dim];
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            jc_dim=neighbor_list_2nd[ic_dim][j];
            
            fi=mu[ic_dim];
            fj=mu[jc_dim];
            tmp0=mat(fi,fj,itype);
            exp_fi=exp(beta*(fi-tmp0));
            exp_fj=exp(beta*(fj-tmp0));
            w_ij=-c[ic_dim]*cv[jc_dim/c_dim]*exp_fi;
            w_ji=-c[jc_dim]*cv[ic_dim/c_dim]*exp_fj;
            s_ij=w_ij-w_ji;
            
            tmp0=(s[jc_dim]*c[ic_dim]-s[ic_dim]*cv[jc_dim/c_dim])*exp_fi
            -(s[ic_dim]*c[jc_dim]-s[jc_dim]*cv[ic_dim/c_dim])*exp_fj;
            
            tmp0+=beta*(ddE[ic_dim]*w_ij-ddE[jc_dim]*w_ji);
            
            tmp0+=-beta*(w_ij-w_ji)*(
            ddE[ic_dim]*dmat(fi,fj,itype)
            +ddE[jc_dim]*dmat(fj,fi,itype));
            
            ddc_[ic_dim]+=tmp0;
            if(jc_dim<natms*c_dim)
                ddc_[jc_dim]-=tmp0;
        }
    }
}
/*--------------------------------------------
 calculate mu crd dE ddE and local energy
 --------------------------------------------*/
void ForceField_eam_dmd::calc_mu()
{
    int natms=atoms->natms;
    type0* c=mapp->c->begin();
    dmd_type* type=mapp->ctype->begin();
    type0* x=mapp->x->begin();
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
    
    for(int i=0;i<natms+atoms->natms_ph;i++)
    {
        cv[i]=1.0;
        for(int j=0;j<c_dim;j++)
            if(c[i*c_dim+j]>=0.0)
                cv[i]-=c[i*c_dim+j];
    }
    
    int istart=0;
    for(int ic_dim=0;ic_dim<c_dim*natms;ic_dim++)
    {
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            E[ic_dim]+=c[jc_dim]*psi_JI[istart];
            mu[ic_dim]+=c[jc_dim]*phi_IJ[istart];
            
            
            if(jc_dim<phi_psi_sz_sz)
            {
                E[jc_dim]+=c[ic_dim]*psi_IJ[istart];
                mu[jc_dim]+=c[ic_dim]*phi_IJ[istart];
                nrgy_strss_lcl[0]+=c[ic_dim]*c[jc_dim]*phi_IJ[istart];
            }
            else
                nrgy_strss_lcl[0]+=0.5*c[ic_dim]*c[jc_dim]*phi_IJ[istart];
                
            istart++;
        }
        if(c[ic_dim]>=0.0)
        {
            itype=type[ic_dim];
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
            tmp2=tmp0+c_0[itype]-3.0*kbT*log(x[3*(ic_dim/c_dim +1)+ic_dim]);
            nrgy_strss_lcl[0]+=c[ic_dim]*tmp2;
            mu[ic_dim]+=tmp2;
            nrgy_strss_lcl[0]+=kbT*calc_ent(c[ic_dim]);
        }
        if(ic_dim%c_dim==c_dim-1)
            nrgy_strss_lcl[0]+=kbT*calc_ent(cv[ic_dim/c_dim]);
    }
    
    atoms->update(dE_ptr);
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
    {
        for(int j=0,jc_dim;j<phi_psi_sz[ic_dim];j++)
        {
            jc_dim=phi_psi_cmp[istart];
            mu[ic_dim]+=psi_IJ[istart]*c[jc_dim]*dE[jc_dim];
            
            if(jc_dim<phi_psi_sz_sz)
            {
                mu[jc_dim]+=psi_JI[istart]*c[ic_dim]*dE[ic_dim];
            }
            istart++;
        }
    }
    
    atoms->update(mu_ptr);
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
        delete [] N_IJ;
        delete [] N_JI;
        delete [] M_IJ;
        delete [] M_JI;
    }
    
    
    int no_pairs=neighbor->no_pairs;
    int natms=atoms->natms;
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
    
    
    type0* c=mapp->c->begin();
    dmd_type* type=mapp->ctype->begin();
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
        
        for(int ic_dim=iatm*c_dim,itype;ic_dim<(iatm+1)*c_dim && c[ic_dim]>=0.0;ic_dim++)
        {
            itype=type[ic_dim];
            
            iistart=istart;
            for(int j=0,jatm;j<neighbor_list_size[iatm];j++)
            {
                jatm=neighbor_list[iatm][j];

                for(int jc_dim=jatm*c_dim,jtype;jc_dim<(jatm+1)*c_dim && c[jc_dim]>=0.0;jc_dim++)
                {
                    jtype=type[jc_dim];
                    
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
    
    CREATE1D(M_IJ,M_N_sz_sz);
    CREATE1D(M_JI,M_N_sz_sz);
    CREATE1D(N_IJ,M_N_sz_sz);
    CREATE1D(N_JI,M_N_sz_sz);
    
    
}
/*--------------------------------------------
 create the sparse matrices
 --------------------------------------------*/
type0 ForceField_eam_dmd::update_J(type0 alpha_inv,type0* a,type0* g)
{
    calc_mu();
    
    type0* c=mapp->c->begin();
    dmd_type* type=mapp->ctype->begin();
    type0* mu=mu_ptr->begin();
    type0* cv=cv_ptr->begin();

    int istart;
    type0 tmp0;
    

    type0* M=M_ptr->begin();
    type0* N=N_ptr->begin();
    type0* c_d=mapp->c_d->begin();
    for(int i=0;i<phi_psi_sz_sz;i++)
    {
        c_d[i]=M[i]=0.0;
        N[i]=alpha_inv;
        if(c[i]>=0.0)
            g[i]=-alpha_inv*c[i]-a[i];
        else
            g[i]=0.0;
    }
    
    type0 fi,fj,exp_fi,exp_fj,w_ij,w_ji,ans_lcl=0.0,dc_ij;
    int** neighbor_list_2nd=neighbor->neighbor_list_2nd;
    int* neighbor_list_size_2nd=neighbor->neighbor_list_size_2nd;

    istart=0;
    for(int ic_dim=0,itype;ic_dim<phi_psi_sz_sz;ic_dim++)
    {
        itype=type[ic_dim];
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            jc_dim=neighbor_list_2nd[ic_dim][j];
            
            fi=mu[ic_dim];
            fj=mu[jc_dim];
            tmp0=mat(fi,fj,itype);
            exp_fi=exp(beta*(fi-tmp0));
            exp_fj=exp(beta*(fj-tmp0));
            w_ij=-c[ic_dim]*cv[jc_dim/c_dim]*exp_fi;
            w_ji=-c[jc_dim]*cv[ic_dim/c_dim]*exp_fj;
            dc_ij=w_ij-w_ji;
            
            
            M_IJ[istart]=beta*(w_ji+dc_ij*dmat(fj,fi,itype));
            N_IJ[istart]=-(c[ic_dim]*exp_fi+cv[ic_dim/c_dim]*exp_fj);

            M_JI[istart]=beta*(w_ij-dc_ij*dmat(fi,fj,itype));
            N_JI[istart]=-(c[jc_dim]*exp_fj+cv[jc_dim/c_dim]*exp_fi);

            M[ic_dim]-=M_JI[istart];
            N[ic_dim]-=N_JI[istart];
            g[ic_dim]+=dc_ij;
            c_d[ic_dim]+=dc_ij;
            
            if(jc_dim<phi_psi_sz_sz)
            {
                M[jc_dim]-=M_IJ[istart];
                N[jc_dim]-=N_IJ[istart];
                g[jc_dim]-=dc_ij;
                c_d[jc_dim]-=dc_ij;
            }
            
            istart++;
        }
        if(c[ic_dim]>=0.0)
            ans_lcl+=g[ic_dim]*g[ic_dim];
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
    type0 tmp0;
    type0* c=mapp->c->begin();
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
            Ax[ic_dim]+=c[ic_dim]*ddE[ic_dim]*psi_JI[istart]*x[jc_dim];
            if(jc_dim<phi_psi_sz_sz)
                Ax[jc_dim]+=c[jc_dim]*ddE[jc_dim]*psi_IJ[istart]*x[ic_dim];
            istart++;
        }
    
    atoms->update(Ax_ptr);
    
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
    
    atoms->update(crd_ptr);
    
    type0* M=M_ptr->begin();
    type0* N=N_ptr->begin();
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        Ax[ic_dim]=M[ic_dim]*b0[ic_dim]+N[ic_dim]*x[ic_dim];
    
    int** neighbor_list_2nd=neighbor->neighbor_list_2nd;
    int* neighbor_list_size_2nd=neighbor->neighbor_list_size_2nd;
    
    istart=0;
    for(int ic_dim=0;ic_dim<phi_psi_sz_sz;ic_dim++)
        for(int j=0,jc_dim;j<neighbor_list_size_2nd[ic_dim];j++)
        {
            jc_dim=neighbor_list_2nd[ic_dim][j];
            Ax[ic_dim]+=N_IJ[istart]*x[jc_dim]+M_IJ[istart]*b0[jc_dim];
            
            if(jc_dim<phi_psi_sz_sz)
                Ax[jc_dim]+=N_JI[istart]*x[ic_dim]+M_JI[istart]*b0[ic_dim];
            
            istart++;
        }
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceField_eam_dmd::
force_calc_static(bool st_clc)
{
    type0* x=mapp->x->begin();
    type0* fvec=f->begin();
    type0* dE=dE_ptr->begin();
    type0* c=mapp->c->begin();
    
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
            fpair=-(psi_r_JI[istart]*dE[ic_dim]+psi_r_IJ[istart]*dE[jc_dim]+phi_r_IJ[istart])*c[ic_dim]*c[jc_dim];
            apair=-(psi_alpha_JI[istart]*dE[ic_dim]+psi_alpha_IJ[istart]*dE[jc_dim]+phi_alpha_IJ[istart])*c[ic_dim]*c[jc_dim];
            
            j_comp=(3+c_dim)*(jc_dim/c_dim);
            J_comp=3*(jc_dim/c_dim)+jc_dim;
            
            dx0=x[i_comp]-x[j_comp];
            dx1=x[i_comp+1]-x[j_comp+1];
            dx2=x[i_comp+2]-x[j_comp+2];
            
            fvec[i_comp]+=dx0*fpair;
            fvec[i_comp+1]+=dx1*fpair;
            fvec[i_comp+2]+=dx2*fpair;
            fvec[I_comp]+=apair*x[I_comp];
            
            if(jc_dim<phi_psi_sz_sz)
            {
                fvec[j_comp]-=dx0*fpair;
                fvec[j_comp+1]-=dx1*fpair;
                fvec[j_comp+2]-=dx2*fpair;
                fvec[J_comp]+=apair*x[J_comp];
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
        
        if(c[ic_dim]>=0.0)
            fvec[I_comp]+=3.0*kbT*c[ic_dim]/x[I_comp];
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


