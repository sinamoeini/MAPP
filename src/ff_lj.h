/*--------------------------------------------
 Created by Sina on 07/15/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifdef FF_Style
    FFStyle(ForceField_lj,lj)
#else
#ifndef __MAPP__ff_lj__
#define __MAPP__ff_lj__
#include "ff.h"
namespace MAPP_NS
{
    class ForceField_lj: public ForceFieldMD
    {
    private:
        void allocate();
        void deallocate();
        int no_types;
        type0** sigma;
        type0** epsilon;
        type0** offset;
        void read_file(char*);
        int shift;
    protected:
        void force_calc(bool);
        type0 energy_calc();
        void pre_xchng_energy(GCMC*);
        type0 xchng_energy(GCMC*);
        void post_xchng_energy(GCMC*);
    public:
        ForceField_lj(MAPP *);
        ~ForceField_lj();

        void init();
        void fin();
        void init_xchng();
        void fin_xchng();
        void coef(int,char**);
    };
}
#endif
#endif
