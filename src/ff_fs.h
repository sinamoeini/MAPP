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
namespace MAPP_NS
{
    class ForceField_fs: public ForceFieldMD
    {
    private:
        int no_types;
        
        type0** mat_t_1;
        type0** mat_t_2;
        type0* mat_A;
        
        type0** cut_phi;
        type0** cut_rho;
        type0** mat_k_1;
        type0** mat_k_2;
        type0** mat_k_3;
        
        /*--------------------------------------------*/
        type0* drhoi_dr;
        type0* drhoj_dr;
        int max_pairs;
        /*--------------------------------------------*/
        
        void read_file(char*);        
        Vec<type0>* rho_ptr;

        void allocate();
        void deallocate();
    protected:
        void force_calc(bool);
        type0 energy_calc();
        void xchng_energy(GCMC*);
    public:
        ForceField_fs(MAPP *);
        ~ForceField_fs();
        void init();
        void fin();
        void init_xchng();
        void fin_xchng();
        void coef(int,char**);
    };
}
#endif
#endif
