/*--------------------------------------------
 Created by Sina on 01/29/14.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifdef FF_Style
    FFStyle(ForceField_fsm,fsm)
#else
#ifndef __MAPP__ff_fsm__
#define __MAPP__ff_fsm__
#include "ff.h"
namespace MAPP_NS
{
    class ForceField_fsm: public ForceFieldMD
    {
    private:
        int arr_size;
        int no_types;
        
        type0** mat_t_1;
        type0** mat_t_2;
        type0* mat_A;
        
        type0* cut_sq_phi;
        type0* cut_sq_rho;
        type0* mat_k_1;
        type0* mat_k_2;
        type0* mat_k_3;
        
        void read_file(char*);
        
        Vec<type0>* rho_ptr;
    protected:
        void force_calc(bool,type0*);
        type0 energy_calc();
    public:
        ForceField_fsm(MAPP *);
        ~ForceField_fsm();
        void init();
        void fin();
        void coef(int,char**);

    };
}
#endif
#endif
