#ifdef DMD_Style
    DMDStyle(DMD_bdf_f,bdf_f)
#else
#ifndef __MAPP__dmd_bdf_f__
#define __MAPP__dmd_bdf_f__
#include "dmd.h"
namespace MAPP_NS {
    class DMD_bdf_f :public DMDImplicit
    {
    private:
        type0* l;
        void update_l(int);
        void update_A_bar();

        type0** A;
        type0** A_bar;
        
        type0* dy;
        type0* z;
        
        
        type0 lo_err_fac[2];
        type0 hi_err_fac[2];
        
        void calc_err_coef_flc();
        void calc_err_coef_vc();
        
        bool lo_avail,hi_avail;
        type0 eta[3];
        
        void update_z_vc();
        
        void update_z_flc();
        
        void eta_calc();
        
        Vec<type0>** vecs_1;
        
        void beta_calc_flc();
        void beta_calc_vc();
    protected:
        // stuff for book keeping
        
        type0* t;
        
        void allocate();
        void deallocate();
        
        //functions to be accessed by run()
        void restart();
        void update_for_next();
        bool interpolate();
        void ord_dt(type0&);
        void err_calc();
        void start();
        void err_fac_calc();
    public:
        DMD_bdf_f(MAPP*,int,char**);
        ~DMD_bdf_f();

    };
}

#endif 
#endif