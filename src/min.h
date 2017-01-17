#ifndef __MAPP__min__
#define __MAPP__min__
#include "type_def.h"
#include "vec.h"
namespace MAPP_NS
{
    class Min
    {
    private:
        class ForceFieldDMD* forcefield_dmd;
        type0 max_dx;
    protected:
        bool chng_box;
        int err;
        type0 f_h;
        type0 curr_energy;
        type0*& nrgy_strss;
        
        int pe_idx;
        int stress_idx;
        bool sts_flag;
        
        void prepare_affine_h();
        type0 calc_ndofs();
        type0 ndofs;
        
        class ThermoDynamics* thermo;
        
        int max_iter;
        bool affine;
        type0 energy_tolerance;
        
        bool** H_dof;
        
        int x_dim;
        
                
        VecTens<type0> h;
        VecTens<type0> x;
        VecTens<type0> x0;
        VecTens<type0> f;
        VecTens<type0> f0;
        
        class Write*& __write__;
        MPI_Comm& world;
    public:
        Min();
        virtual ~Min();
        void print_error();
        virtual void run()=0;
        virtual void init();
        virtual void fin();
        
        void force_calc();
        type0 calc_ave_f_norm();
        bool output_flag;

        type0 F(type0);
        type0 dF(type0,type0&);
        void ls_prep(type0&,type0&,type0&);
        void F_reset();

        static LineSearch<Min>* ls;
    };
}

#endif
