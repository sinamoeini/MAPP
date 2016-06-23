#ifdef FF_Style
    FFStyle(ForceField_eam_dmd_II,eam_dmd_II)
#else
#ifndef __MAPP__ff_eam_dmd_II__
#define __MAPP__ff_eam_dmd_II__
#include "ff.h"
#include "eam_file_reader.h"
namespace MAPP_NS
{    
    class ForceField_eam_dmd_II : public ForceFieldDMD
    {
    private:
        int no_types;
        void allocate();
        void deallocate();
        
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
        int c_dim;
        int dim;
        int x_dim;
        type0 cut_sq_mod_0;
        type0 kbT,beta;
        type0* c_0;
        type0* c_1;
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
        
        type0 calc_ent(type0);
        void calc_mu();

        
        type0* psi_IJ;
        type0* psi_JI;
        type0* phi_IJ;
        
        type0* psi_r_IJ;
        type0* psi_r_JI;
        type0* phi_r_IJ;
        type0* psi_alpha_IJ;
        type0* psi_alpha_JI;
        type0* phi_alpha_IJ;
        
        
        int* phi_psi_cmp;
        int* phi_psi_sz;
        int phi_psi_sz_sz;
        int n_phi_psi;
        type0* N_IJ;
        type0* M_IJ;
        int M_N_sz_sz;
        void force_calc_static(bool);
        /*--------------------------------------------*/
        void calc_Q_all(type0&,type0&,type0&,
                        type0&,type0&,type0&,
                        int&,type0&,type0&,
                        type0&,
                        type0&,type0&,type0&,type0&);
        
        void calc_Q_all(type0&,type0&,type0&,
                        type0&,type0&,type0&,
                        int&,type0&,type0&,
                        type0&);
        Vec<type0>* x_tmp_ptr;
        type0 alpha_tmp,alpha_inv_tmp;
        /*--------------------------------------------*/
        
    protected:
        void force_calc(bool);
        type0 energy_calc();
        void dc();
        type0 ddc_norm();
        void ddc(type0*);
    public:
        ForceField_eam_dmd_II(MAPP *);
        ~ForceField_eam_dmd_II();
        void init();
        void fin();
        void coef(int,char**);
        void set_temp(type0);
        
        void operator()(Vec<type0>*,Vec<type0>*);
        void init_static();
        type0 update_J(type0,type0*,type0*);
    };
}
#endif
#endif 


