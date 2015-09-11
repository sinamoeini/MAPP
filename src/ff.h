/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__ff__
#define __MAPP__ff__
#include "init.h"
#include "atoms.h"
namespace MAPP_NS {
    class ForceField : protected InitPtrs{
    private:
        int cut_alloc;
        int ns_alloc;

        
    protected:
        type0* nrgy_strss;
        type0* image;
        int image_size;
        void image_calc();
        
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
        virtual type0 g_calc(int,type0,type0*,type0*,type0*)=0;
        virtual void c_d_calc(int,type0*)=0;
        virtual type0 c_dd_norm()=0;
        
        void force_calc_timer(int,type0*);
        type0 energy_calc_timer();
        type0 g_calc_timer(int,type0,type0*,type0*,type0*);
        void c_d_calc_timer(int,type0*);
        void create_2nd_neigh_lst_timer();
        
        type0 alpha_min,alpha_max;
        
    };
}
#endif

