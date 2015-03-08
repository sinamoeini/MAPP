/*--------------------------------------------
 Created by Sina on 2/1/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__min__
#define __MAPP__min__

#include <iostream>
#include "init.h"
#include "atoms.h"
#include "line_search.h"
#include "thermo_dynamics.h"

namespace MAPP_NS {

    class Min : protected InitPtrs
    {
    private:
        int ns_alloc;
    protected:
        VecLst* vecs_comm;
        LineSearch* line_search;
        int chng_box;
        int dim;
        int err;
        type0 curr_energy;
        ThermoDynamics* thermo;
        
        type0* nrgy_strss;
        
        int pe_idx;
        int stress_idx;
        char* dof;
        int dof_n,cdof_n;
        int id_n;
        int c_type_n;
        
        void reg_h_H(type0**,type0**);
    public:
        Min(MAPP *);
        virtual ~Min();
        void print_error();
        virtual void run()=0;
        virtual void init()=0;
        virtual void fin()=0;
        void rectify_f(type0*);
        int max_iter;
        type0 energy_tolerance;
        int** H_dof;

    };


}
#endif /* defined(__MAPP__min__) */
