#ifdef FF_Style
    FFStyle(ForceField_eam,eam)
#else
#ifndef __MAPP__ff_eam__
#define __MAPP__ff_eam__
#include "ff.h"
namespace MAPP_NS
{
    class ForceField_eam : public ForceField
    {
    private:
        int no_types;
    protected:
        type0*** F_arr;
        type0*** phi_r_arr;
        type0*** rho_arr;
        

        
        int eam_mode;
        int nr,nrho;
        type0 dr,drho,dr_inv,drho_inv,rho_max;
        
        
        void set_setfl(int,char**);
        void set_funcfl(int,char**);
        void set_fs(int,char**);
        int line_read(FILE*,char*&);
        
        void allocate();
        int allocated;
        void clean_up();
        
        void set_arrays();
        void interpolate(int,type0,type0**);
        type0 interpolate(type0*,int,type0,int);
        
        
        int** type2rho;
        int** type2phi;
        
        int rho_n,f_n,type_n,x_n,dF_n;
        type0 cut_sq_0;
        
        /*--------------------------------------------*/
        type0* drhoi_dr;
        type0* drhoj_dr;
        int max_pairs;
        /*--------------------------------------------*/
        
    public:
        ForceField_eam(MAPP *);
        ~ForceField_eam();
        void force_calc(int,type0*);
        type0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);

        void create_2nd_neigh_lst(){};
        type0 g_calc(int,type0,type0*,type0*,type0*){return 0.0;};
        type0 c_dd_norm(){return 0.0;};
        void c_d_calc(int,type0*){};
    };
    

    
}
#endif
#endif
