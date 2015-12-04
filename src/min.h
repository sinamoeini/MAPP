/*--------------------------------------------
 Created by Sina on 2/1/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__min__
#define __MAPP__min__
#include "init.h"
#include "ls.h"
#include "ff.h"
#include "atoms.h"
#include "ls.h"
#include "vec.h"
#include "xmath.h"
#include "thermo_dynamics.h"
namespace MAPP_NS
{
    class Min : protected InitPtrs
    {
    private:
        XMath* xmath;
        ForceFieldDMD* forcefield_dmd;
        type0 max_dx;
    protected:
        VecLst* vecs_comm;
        bool chng_box;
        int dim;
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
        
        ThermoDynamics* thermo;
        
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
        
    public:
        Min(MAPP *);
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

    };


}

#endif