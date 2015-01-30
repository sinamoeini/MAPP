/*--------------------------------------------
 Created by Sina on 06/20/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifdef MD_Style
    MDStyle(MD_nh,nh)
#else
#ifndef __MAPP__md_nh__
#define __MAPP__md_nh__
#include "init.h"
#include "atoms.h"
#include "thermo_dynamics.h"
#include "md.h"
#include <iostream>
namespace MAPP_NS {
    class MD_nh : public MD{
    private:
        int x_n,x_d_n,f_n,type_n,id_n,dof_n;
        TYPE0 MTK_1,MTK_2;
        TYPE0 ke_cur,ke_tar,t_tar,t_cur;
        TYPE0* ke_curr;
        TYPE0 t_freq,*tau_freq,tau_freq_m;
        TYPE0 dt2,dt4,dt8;
        TYPE0 no_dof;
        int chk_stress;
        int chk_create_vel,seed;
        TYPE0** M1;
        TYPE0** M2;
        TYPE0** H0;
        TYPE0** H0_inv;
        
        TYPE0* tmp_fac;
        TYPE0* tmp_ke_curr;
        TYPE0* tmp_x;
        
        int* chk_tau;
        TYPE0* tau_tar;
        TYPE0* v_per_atm;
        
        int no_it_eta,no_ch_eta;
        TYPE0* eta_m;
        TYPE0* eta_d;
        TYPE0* eta_dd;

        int no_it_peta,no_ch_peta;
        TYPE0* peta_m;
        TYPE0* peta_d;
        TYPE0* peta_dd;
        
        TYPE0* omega_m;
        TYPE0* omega_d;
                
        void update_H(TYPE0);
        void update_x(TYPE0);
        void update_x_d(TYPE0);
        void update_NH_T(TYPE0);
        void update_NH_tau(TYPE0);
        void update_omega_d(TYPE0);
        void update_x_d_xpnd(TYPE0);
        
        void zero_f();
        void couple();
        
        void create_vel(int,TYPE0);
        void init_vel(TYPE0);
        
        TYPE0* x_ave;
        TYPE0* x_ave_tot;
        
        int omega_denom;
        int x_dim,x_d_dim,f_dim,dof_dim;
    protected:
    public:
        MD_nh(MAPP *,int,char**);
        ~MD_nh();
        void init();
        void fin();
        void run(int);
    };
}
#endif
#endif
