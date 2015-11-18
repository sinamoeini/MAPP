#ifndef __MAPP__eam_file_reader__
#define __MAPP__eam_file_reader__
#include "init.h"
#include <stdio.h>

namespace MAPP_NS
{
    class EAMFileReader:protected InitPtrs
    {
    private:
        int no_types;
        int allocated;
        void allocate();
        void deallocate();
        void set_arrays();
        type0 interpolate(type0*,int,type0,int);
        void interpolate_4(int,type0,type0**);
        void interpolate_5(int,type0,type0**);
        void interpolate_7(int,type0,type0**);
        void set_setfl();
        void set_funcfl();
        void set_fs();
        
        int eam_mode;
        
        int nfiles;
        char** files;
        int* types;
    protected:
    public:
        EAMFileReader(MAPP*);
        ~EAMFileReader();
        void setup();
        void set_cut_sq();
        int nr,nrho;
        type0 dr,drho,dr_inv,drho_inv,rho_max,rc;
        type0*** F_arr;
        type0*** phi_r_arr;
        type0*** rho_r_arr;
        type0*** rho_arr;
        type0* cut_sq;

        
        int** type2rho;
        int** type2phi;
        int** type2rho_pair_ij;
        int** type2rho_pair_ji;
        int** type2phi_pair_ij;
        int** type2phi_pair_ji;
        int stride;
        
        
        void file_format(const char*);
        void add_file(const char*,int);

    };
}

#endif
