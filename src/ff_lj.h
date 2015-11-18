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
        int x_n,f_n,type_n;
        int arr_size;
        type0* sigma;
        type0* epsilon;
        type0* offset;
        void read_file(char*);
        int shift;
    protected:
        void force_calc(int,type0*);
        type0 energy_calc();
    public:
        ForceField_lj(MAPP *);
        ~ForceField_lj();

        void init();
        void fin();
        void coef(int,char**);


    };
}
#endif
#endif
