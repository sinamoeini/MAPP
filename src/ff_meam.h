
#ifdef FF_Style
    FFStyle(ForceField_meam,meam)
#else
#ifndef __MAPP__ff_meam__
#define __MAPP__ff_meam__

#include <iostream>
#include "ff.h"
#include "atoms.h"
namespace MAPP_NS
{
    class ForceField_meam : public ForceField
    {
    private:
        
        type0* delij;
        type0* delik;
        type0* deljk;
        type0* s;
        type0* si;
        type0* sj;
        type0* drho1drm1;
        type0* drho1drm2;
        type0* drho2drm1;
        type0* drho2drm2;
        type0* drho3drm1;
        type0* drho3drm2;
        type0* drhodrm1;
        type0* drhodrm2;
        type0* dUdrijm;
        type0* fi;
        type0* fj;
        type0* v;
        
    protected:
        /*---------------------------------*/
        int rho_dim;
        int rho_n,rho_vec_n,type_n,x_n,f_n;
        /*---------------------------------*/
        
        
        
        /*---------------------------------*/
        int* v2d;
        int* v3d;
        int*** vind3d;
        int** vind2d;
        /*---------------------------------*/
        
        /*
         rho dimension is 38!!!!
         **** rho functions (4)
         00: rho0
         01: rho1
         02: rho2
         03: rho3
         
         **** rho_1 terms (3)
         04: Arho1_1
         05: Arho1_2
         06: Arho1_3
         
         **** rho_2 terms (7)
         ******** 1st terms (6)
         07: Arho2_1
         08: Arho2_2
         09: Arho2_3
         10: Arho2_4
         11: Arho2_5
         12: Arho2_6
         ******** 2nd term (1)
         13: Arho2b
         
         **** rho_3 terms
         ******** 1st terms (10)
         14: Arho3_1
         15: Arho3_2
         16: Arho3_3
         17: Arho3_4
         18: Arho3_5
         19: Arho3_6
         20: Arho3_7
         21: Arho3_8
         22: Arho3_9
         23: Arho3_10
         ******** 2nd terms (3)
         24: Arho3b_1
         25: Arho3b_2
         26: Arho3b_3
         
         **** t_ave terms (3)
         27: t_ave_1
         28: t_ave_2
         29: t_ave_3
         
         **** tsq_ave terms (3)
         30: tsq_ave_1
         31: tsq_ave_2
         32: tsq_ave_3
         
         **** Gamma function (1)
         33: gamma
         
         **** dGamma functions (3)
         34: dgamma1
         35: dgamma2
         36: dgamma3
         
         **** fhop term  (1)
         37: fhop
         */
        
        /*---------------------------------*/
        int emb_lin_neg;//checked
        int ialloy;//checked
        int mix_ref_t;//checked
        int bkgd_dyn;//checked
        int augt1;//checked
        int erose_form;//checked
        int nr;
        type0 delr_meam;//checked
        type0 delr_meam_inv;//checked
        type0 rc_meam;//checked
        type0 gsmooth_factor;//checked
        type0 dr;
        type0 dr_inv;
        
        type0*** c_min;//checked
        type0*** c_max;//checked
        type0*** phirar;
        
        type0** re_meam;//checked
        type0** ebound_meam;//checked
        type0** Ec_meam;//checked
        type0** alpha_meam;//checked
        type0** delta_meam;//checked
        type0** attrac_meam;//checked
        type0** repuls_meam;//checked
        int** lattice;//checked
        int** nn2_meam;//checked
        int** zbl_meam;//checked
        
        type0* rho0_meam;//checked
        type0* beta0_meam;//checked
        type0* beta1_meam;//checked
        type0* beta2_meam;//checked
        type0* beta3_meam;//checked
        type0* t0_meam;//checked
        type0* t1_meam;//checked
        type0* t2_meam;//checked
        type0* t3_meam;//checked
        type0* Z_meam;//checked
        type0* rho_ref_meam;//checked
        type0* A_meam;//checked
        int* ibar_meam;//checked
        int* ielt_meam;//checked
        int* type_ref;//not set yet
        
        void get_dens_ref(type0,int,int,type0&,type0&,type0&,type0&,type0&,type0&,type0&,type0&);
        void get_sijk(type0,int,int,int,type0&);
        void get_tavref(type0&,type0&,type0&,type0&,type0&,type0&,type0,type0,type0,type0,type0,type0,type0,int,int,int);
        type0 phi_meam(type0,int,int);
        /*---------------------------------*/
        
        
        /*---------------------------------*/
        int max_pairs;
        type0* scrfcn;//not set yet
        type0* dscrfcn;//not set yet
        type0* fcpair;//not set yet
        /*---------------------------------*/
                

        void reset();

        void fcut(type0,type0&);
        void dfcut(type0,type0&,type0&);
        void dCfunc(type0,type0,type0,type0&);
        void dCfunc2(type0,type0,type0,type0&,type0&);
        
        void G_gam(type0,int,type0,type0&);
        void dG_gam(type0,int,type0,type0&,type0&);
        void get_shpfcn(type0*,int);
        
        void read_global(char*);
        void read_local(char*);
        
        void get_Zij(type0&,int);
        void get_Zij2(type0&,type0&,type0&,int,type0,type0);
        type0 zbl(type0,type0,type0);
        type0 erose(type0,type0,type0,type0,type0,type0,int);
        void alloy_params();
        void compute_reference_density();
        void compute_pair_meam();
        void compute_phi(type0,int,int,type0&);
        void compute_phi_dphi(type0,int,int,type0&,type0&);
        
        type0 third,sixth;
        
    public:
        ForceField_meam(MAPP *);
        ~ForceField_meam();
        void force_calc1(int,type0*);
        void force_calc(int,type0*);
        type0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);
        
        void create_2nd_neigh_lst(){};
        type0 g_calc(int,type0,type0*,type0*){return 0.0;};
        void c_d_calc(){};
    };
}


#endif
#endif
