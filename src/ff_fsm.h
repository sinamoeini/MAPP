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
        int no_types;
        
        type0** mat_t_1;
        type0** mat_t_2;
        type0** cut_sq_phi;
        type0** cut_sq_rho;
        type0** mat_k_1;
        type0** mat_k_2;
        type0** mat_k_3;
        type0* mat_A;
        
        void read_file(char*);
        void allocate();
        void deallocate();
        
        Vec<type0>* rho_ptr;
    protected:
        void force_calc(bool);
        type0 energy_calc();
        void pre_xchng_energy(GCMC*);
        type0 xchng_energy(GCMC*);
        void post_xchng_energy(GCMC*);
    public:
        ForceField_fsm();
        ~ForceField_fsm();
        void init();
        void fin();
        void init_xchng();
        void fin_xchng();
        void coef(int,char**);


    };
}
#endif
#endif
