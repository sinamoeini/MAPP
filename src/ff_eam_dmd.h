#ifdef FF_Style
    FFStyle(ForceField_eam_dmd,eam_dmd)
#else
#ifndef __MAPP__ff_eam_dmd__
#define __MAPP__ff_eam_dmd__
#include "ff.h"
#include "eam_file_reader.h"
namespace MAPP_NS
{    
    class ForceField_eam_dmd : public ForceFieldDMD
    {
    private:
        int no_types;
        type0 min_log;
    protected:
        int stride;
        int nr,nrho;
        type0 dr,drho,dr_inv,drho_inv,rho_max;
        type0*** F_arr;
        type0*** phi_r_arr;
        type0*** rho_r_arr;
        int** type2rho;
        int** type2phi;
        int** type2rho_pair_ij;
        int** type2rho_pair_ji;
        int** type2phi_pair_ij;
        int** type2phi_pair_ji;
        
        Vec<type0>* cv_ptr;
        Vec<type0>* E_ptr;
        Vec<type0>* dE_ptr;
        Vec<type0>* mu_ptr;
        Vec<type0>* crd_ptr;        
        Vec<type0>* s_ptr;
        Vec<type0>* t_ptr;

        EAMFileReader* eam_reader;

        void read_file(char*);
        
        /*--------------------------------------------*/
        
        /*--------------------------------------------*/

        
        int c_dim;
        int dim;
        int x_dim;
        type0 cut_sq_mod_0;
        type0 kbT,beta;
        type0* c_0;
        type0* f_t;
        type0* Q_nrm;
        
        /*--------------------------------------------*/
        
        type0* rho_phi;
        type0* drho_phi_dr;
        type0* drho_phi_dalpha;
        int max_pairs;
        
        /*--------------------------------------------*/

        type0* xi;
        type0* wi_0;
        type0* wi_1;
        type0* wi_2;
        int no_i;
        type0 set_weight_abs(int);
        
        /*--------------------------------------------*/
        
        type0 mat(type0,type0,type0,type0,int);
        type0 dmat0(type0,type0,type0,type0,int);
        type0 dmat1(type0,type0,type0,type0,int);
        type0 calc_ent(type0);
        type0 mod_log(type0);
        void calc_cv(type0*,type0*,int);
        
        void dc_crd();
        void dc_ncrd();
        type0 dc_en_proj_crd(bool,type0*,type0&);
        type0 dc_en_proj_ncrd(bool,type0*,type0&);
        type0 imp_cost_grad_crd(bool,type0,type0*,type0*);
        type0 imp_cost_grad_ncrd(bool,type0,type0*,type0*);
        type0 ddc_norm_crd();
        type0 ddc_norm_ncrd();
        type0 dc_norm_grad_crd(bool,type0*,type0*);
        type0 dc_norm_grad_ncrd(bool,type0*,type0*);
        
        void calc_mu_crd();
        void calc_mu_ncrd();
        
        void force_calc(bool,type0*);

        type0 energy_calc();
        void dc();
        type0 dc_en_proj(bool,type0*,type0&);
        type0 ddc_norm();
        type0 imp_cost_grad(bool,type0,type0*,type0*);
        type0 dc_norm_grad(bool,type0*,type0*);
        type0 en_grad(bool,type0*,type0*);
        void enst_calc(bool,type0*);
    public:
        ForceField_eam_dmd(MAPP *);
        ~ForceField_eam_dmd();
        void init();
        void fin();
        void coef(int,char**);
        


    };
    
    
    
}

#endif
#endif 




/*auto f0=[this]
 (type0& r_ij,type0& r_inv_ij,type0 alpha_ij,
 int rho_ij,int rho_ji,int phi_ij,
 type0& ci, type0& cj,
 type0& Ei, type0& Ej,
 int itype,int jtype,
 bool jown)->void
 {
 if(alpha_ij>=alpha_max)
 return;
 
 type0 upper=(r_ij+rc)/alpha_ij;
 type0 lower=(r_ij-rc)/alpha_ij;
 
 if(lower>=xi[no_i-1])
 return;
 
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
 rtmp=r_ij-xi[i]*alpha_ij;
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
 
 tmp0=PI_IN_SQ*r_inv_ij;
 
 rho_phi[rho_ij]*=tmp0;
 drho_phi_dr[rho_ij]*=tmp0;
 drho_phi_dr[rho_ij]-=rho_phi[rho_ij]*r_inv_ij;
 drho_phi_dalpha[rho_ij]*=-tmp0;
 
 if(itype!=jtype)
 {
 rho_phi[rho_ji]*=tmp0;
 drho_phi_dr[rho_ji]*=tmp0;
 drho_phi_dr[rho_ji]-=rho_phi[rho_ji]*r_inv_ij;
 drho_phi_dalpha[rho_ji]*=-tmp0;
 }
 
 
 rho_phi[phi_ij]*=tmp0;
 drho_phi_dr[phi_ij]*=tmp0;
 drho_phi_dr[phi_ij]-=rho_phi[phi_ij]*r_inv_ij;
 drho_phi_dalpha[phi_ij]*=-tmp0;
 
 Ei+=cj*rho_phi[rho_ji];
 
 if(jown)
 {
 Ej+=ci*rho_phi[rho_ij];
 nrgy_strss[0]+=ci*cj*rho_phi[phi_ij];
 }
 else
 nrgy_strss[0]+=0.5*ci*cj*rho_phi[phi_ij];
 };*/

/*
 auto f0=[this]
 (type0& r_ij,type0& r_inv_ij,type0 alpha_ij,
 int rho_ij,int rho_ji,int phi_ij,
 type0& ci, type0& cj,
 type0& Ei, type0& Ej,
 int itype,int jtype,
 bool jown)->void
 {
 if(alpha_ij>=alpha_max)
 return;
 
 type0 upper=(r_ij+rc)/alpha_ij;
 type0 lower=(r_ij-rc)/alpha_ij;
 
 if(lower>=xi[no_i-1])
 return;
 
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
 type0 p,tmp0,tmp1,tmp2;
 type0* coef;
 for(int i=0;i<no_i;i++)
 {
 if(xi[i]<=lower && xi[i]>=upper)
 continue;
 rtmp=r_ij-xi[i]*alpha_ij;
 p=fabs(rtmp)*dr_inv;
 m=static_cast<int>(p);
 m=MIN(m,nr-2);
 p-=m;
 p=MIN(p,1.0);
 
 coef=rho_r_arr[type2rho[itype][jtype]][m];
 tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
 if(rtmp<0.0)
 tmp0*=-1.0;
 
 rho_phi[rho_ij]+=wi_0[i]*tmp0;
 drho_phi_dr[rho_ij]+=wi_1[i]*tmp0;
 drho_phi_dalpha[rho_ij]+=wi_2[i]*tmp0;
 
 if(itype!=jtype)
 {
 coef=rho_r_arr[type2rho[jtype][itype]][m];
 tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
 if(rtmp<0.0)
 tmp0*=-1.0;
 rho_phi[rho_ji]+=wi_0[i]*tmp0;
 drho_phi_dr[rho_ji]+=wi_1[i]*tmp0;
 drho_phi_dalpha[rho_ji]+=wi_2[i]*tmp0;
 }
 
 
 
 coef=phi_r_arr[type2phi[itype][jtype]][m];
 tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
 if(rtmp<0.0)
 tmp0*=-1.0;
 
 rho_phi[phi_ij]+=wi_0[i]*tmp0;
 drho_phi_dr[phi_ij]+=wi_1[i]*tmp0;
 drho_phi_dalpha[phi_ij]+=wi_2[i]*tmp0;
 }
 
 tmp0=PI_IN_SQ*r_inv_ij;
 tmp1=tmp0*r_inv_ij;
 tmp2=tmp0/alpha_ij;
 
 drho_phi_dalpha[rho_ij]*=2.0*tmp2;
 drho_phi_dalpha[rho_ij]-=tmp2*rho_phi[rho_ij];
 drho_phi_dr[rho_ij]*=-2.0*tmp2;
 drho_phi_dr[rho_ij]-=tmp1*rho_phi[rho_ij];
 rho_phi[rho_ij]*=tmp0;
 
 
 
 if(itype!=jtype)
 {
 drho_phi_dalpha[rho_ji]*=2.0*tmp2;;
 drho_phi_dalpha[rho_ji]-=tmp2*rho_phi[rho_ji];
 drho_phi_dr[rho_ji]*=-2.0*tmp2;
 drho_phi_dr[rho_ji]-=tmp1*rho_phi[rho_ji];
 rho_phi[rho_ji]*=tmp0;
 }
 
 
 drho_phi_dalpha[phi_ij]*=2.0*tmp2;;
 drho_phi_dalpha[phi_ij]-=tmp2*rho_phi[phi_ij];
 drho_phi_dr[phi_ij]*=-2.0*tmp2;
 drho_phi_dr[phi_ij]-=tmp1*rho_phi[phi_ij];
 rho_phi[phi_ij]*=tmp0;
 
 Ei+=cj*rho_phi[rho_ji];
 
 if(jown)
 {
 Ej+=ci*rho_phi[rho_ij];
 nrgy_strss[0]+=ci*cj*rho_phi[phi_ij];
 }
 else
 nrgy_strss[0]+=0.5*ci*cj*rho_phi[phi_ij];
 };
 */
/*
 f0(r,r_inv,sqrt(x[3*iatm+ic_dim+3]*x[3*iatm+ic_dim+3]+x[3*jatm+jc_dim+3]*x[3*jatm+jc_dim+3]),
 istart+type2rho_pair_ij[itype][jtype],
 istart+type2rho_pair_ji[jtype][itype],
 istart+type2phi_pair_ij[itype][jtype],
 c[ic_dim], c[jc_dim],
 E[ic_dim], E[jc_dim],
 type[ic_dim],type[jc_dim],
 jatm<natms);*/