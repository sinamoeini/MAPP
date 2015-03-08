/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__ff__
#define __MAPP__ff__
#include <iostream>
#include "init.h"
#include "vecmath.h"
#include "atoms.h"
namespace MAPP_NS {
    class ForceField : protected InitPtrs{
    private:
        int cut_alloc;
        int ns_alloc;

        
    protected:
        type0* nrgy_strss;
    public:
        ForceField(MAPP *);
        virtual ~ForceField();

        
        virtual void init(){};
        virtual void fin(){};
        virtual void coef(int,char**)=0;
        type0* cut_sq;
        type0* cut_sk_sq;


        virtual void force_calc(int,type0*)=0;
        virtual type0 energy_calc()=0;
        virtual void create_2nd_neigh_lst()=0;
        virtual type0 g_calc(int,type0,type0*,type0*)=0;
        virtual void c_d_calc()=0;
        
        void force_calc_timer(int,type0*);
        type0 energy_calc_timer();
        type0 g_calc_timer(int,type0,type0*,type0*);
        void c_d_calc_timer();
        void create_2nd_neigh_lst_timer();
        
    };
}
#endif

