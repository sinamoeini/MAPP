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
namespace MAPP_NS {
    class ForceField_fsm : public ForceField{
    private:
        int x_n,f_n,type_n,rho_n;
        int arr_size;
        
        type0** mat_t_1;
        type0** mat_t_2;
        type0* mat_A;
        
        type0* cut_sq_phi;
        type0* cut_sq_rho;
        type0* mat_k_1;
        type0* mat_k_2;
        type0* mat_k_3;
        
        void read_file(char*);
        int read_line(FILE*,char*&);
        
    protected:
    public:
        ForceField_fsm(MAPP *);
        ~ForceField_fsm();
        void force_calc(int,type0*);
        type0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);
        int shift;

        void create_2nd_neigh_lst(){};
        type0 g_calc(int,type0,type0*,type0*){return 0.0;};
        void c_d_calc(){};
    };
}
#endif
#endif
