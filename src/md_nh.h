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
#include "md.h"
#include <iostream>
namespace MAPP_NS {
    class MD_nh : public MD{
    private:
        int x_d_n,f_n,type_n,id_n;
        bool dof_xst;
        type0 MTK_1,MTK_2;
        type0 ke_cur,ke_tar,t_tar,t_cur;
        type0* ke_curr;
        type0 t_freq,*tau_freq,tau_freq_m;
        type0 dt2,dt4,dt8;
        type0 no_dof;
        int chk_stress;
        int chk_create_vel,seed;
        type0** M1;
        type0** M2;
        
        type0* tmp_fac;
        type0* tmp_ke_curr;
        type0* tmp_x;
        
        int* chk_tau;
        type0* tau_tar;
        type0* v_per_atm;
        
        int no_it_eta,no_ch_eta;
        type0* eta_m;
        type0* eta_d;
        type0* eta_dd;

        int no_it_peta,no_ch_peta;
        type0* peta_m;
        type0* peta_d;
        type0* peta_dd;
        
        type0* omega_m;
        type0* omega_d;
                
        void update_H(type0);
        void update_x(type0);
        void update_x_d(type0);
        void update_NH_T(type0);
        void update_NH_tau(type0);
        void update_omega_d(type0);
        void update_x_d_xpnd(type0);
        
        void zero_f();
        void couple();
        
        void create_vel(int,type0);
        void init_vel(type0);
        
        type0* x_ave;
        type0* x_ave_tot;
        
        int omega_denom;
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
