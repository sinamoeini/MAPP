#ifdef FF_Style
    FFStyle(ForceField_eam_dmd_lg_0,eam_dmd_lg_0)
#else
#ifndef __MAPP__ff_eam_dmd_lg_0__
#define __MAPP__ff_eam_dmd_lg_0__

#include <stdio.h>
#include <iostream>
#include "ff.h"
#include "atoms.h"
namespace MAPP_NS
{
    class ForceField_eam_dmd_lg_0 : public ForceField
    {
    private:
        int no_types;
    protected:
        type0*** F_arr;
        type0*** phi_r_arr;
        type0*** rho_r_arr;
        type0*** phi_r;
        type0*** rho_r;
        
        
        
        int eam_mode;
        int nr,nrho;
        type0 dr,drho,dr_inv,drho_inv,rho_max;
        
        void set_setfl(int,char**);
        void set_funcfl(int,char**);
        void set_fs(int,char**);
        int line_read(FILE*,char*&);
        void setup_delta_e(char*);
        
        void allocate();
        int allocated;
        void clean_up();
        void tmp_allocate();
        void tmp_clean_up();
        
        
        void set_arrays();
        void interpolate(int,type0,type0**);
        type0 interpolate(type0*,int,type0,int);
        void interpolate_m(int,type0,type0**);
        
        int** type2rho;
        int** type2phi;
        int** type2rho_pair_ij;
        int** type2rho_pair_ji;
        int** type2phi_pair_ij;
        int** type2phi_pair_ji;
        int stride;
        
        int f_n,type_n,x_n,c_n,E_n,dE_n,ddE_n,
        n_n,s_n,crd_n,t_n,v_n,c_d_n,m_n;
        type0 cut_sq_0;
        type0 cut_sq_mod_0;
        type0 rc,mod_rc,kbT,beta;
        type0* c_0;
        type0* delta_e;
        
        /*--------------------------------------------*/
        type0* rho_phi;
        type0* drho_phi_dr;
        type0* drho_phi_dalpha;
        
        int max_pairs;
        /*--------------------------------------------*/
        
        
        type0* xi;
        type0* wi;
        
        type0 alpha_min;
        int no_i;
        void set_weight_abs(int);
        void rho_calc(type0,type0,int,int);
        void phi_calc(type0,type0,int,int);
        /*--------------------------------------------*/
        type0 r_crd,rsq_crd;
        int** neigh_lst;
        int* neigh_lst_sz;
        int neigh_lst_sz_sz;
        /*--------------------------------------------*/
        
        type0 mat(type0,type0,type0,type0,int);
        type0 dmat0(type0,type0,type0,type0,int);
        type0 dmat1(type0,type0,type0,type0,int);
        type0 calc_ent(type0);
        
    public:
        ForceField_eam_dmd_lg_0(MAPP *);
        ~ForceField_eam_dmd_lg_0();
        void force_calc(int,type0*);
        type0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);
        
        void create_2nd_neigh_lst();
        type0 g_calc(int,type0,type0*,type0*);
        void c_d_calc();
    };
    
    
}


#endif
#endif