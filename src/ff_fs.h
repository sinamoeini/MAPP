/*--------------------------------------------
 Created by Sina on 01/29/14.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifdef FF_Style
    FFStyle(ForceField_fs,fs)
#else
#ifndef __MAPP__ff_fs__
#define __MAPP__ff_fs__
#include "ff.h"
namespace MAPP_NS {
    class ForceField_fs : public ForceField{
    private:
        int x_n,f_n,type_n,rho_n;
        int arr_size;
        
        type0** mat_t_1;
        type0** mat_t_2;
        type0* mat_A;
        
        type0* cut_phi;
        type0* cut_rho;
        type0* mat_k_1;
        type0* mat_k_2;
        type0* mat_k_3;
        
        /*--------------------------------------------*/
        type0* drhoi_dr;
        type0* drhoj_dr;
        int max_pairs;
        /*--------------------------------------------*/
        
        void read_file(char*);
        int read_line(FILE*,char*&);
        
        
    protected:
    public:
        ForceField_fs(MAPP *);
        ~ForceField_fs();
        void force_calc(int,type0*);
        void force_calc0(int,type0*);
        type0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);
        int shift;

        void create_2nd_neigh_lst(){};
        type0 g_calc(int,type0,type0*,type0*,type0*){return 0.0;};
        void c_d_calc(int,type0*){};
    };
}
#endif
#endif
