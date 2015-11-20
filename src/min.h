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
#include "xmath.h"
#include "thermo_dynamics.h"
namespace MAPP_NS
{
    class Min : protected InitPtrs
    {
    private:
        int ns_alloc;
        XMath* xmath;
        ForceFieldDMD* forcefield_dmd;
        type0 max_dx;
    protected:
        VecLst* vecs_comm;
        int chng_box;
        int dim;
        int err;
        type0 curr_energy;
        type0* nrgy_strss;
        
        int pe_idx;
        int stress_idx;
        int sts_flag;
        
        void force_calc();
        void prepare_affine_h(type0* x,type0* h);
        void zero_f();
        
        ThermoDynamics* thermo;
    public:
        Min(MAPP *);
        virtual ~Min();
        void print_error();
        virtual void run()=0;
        virtual void init();
        virtual void fin();
        void rectify(type0*);
        int max_iter;
        int affine;
        type0 energy_tolerance;

        int** H_dof;
        type0** N;
        
        int x_dim;
        
        Vec<type0>* h_ptr;
        type0** h_H;

        type0** H_prev;
        Vec<type0>* x_prev_ptr;
        
        type0** f_H_prev;
        Vec<type0>* f_prev_ptr;
        
        //Vec<type0>* f_ptr;
        type0** f_H;
        
        type0 F(type0);
        type0 dF(type0,type0&);
        void ls_prep(type0&,type0&,type0&);
        void F_reset();

    };


}

#endif