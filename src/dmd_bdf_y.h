#ifdef DMD_Style
    DMDStyle(DMD_bdf_y,bdf_y)
#else
#ifndef __MAPP__dmd_bdf_y__
#define __MAPP__dmd_bdf_y__
#include "dmd.h"
namespace MAPP_NS {
    class DMD_bdf_y :public DMDImplicit
    {
    private:
        type0 err_est(int);
    protected:
        // stuff for book keeping
        type0* t;
        type0** y;
        type0* dy;
        Vec<type0>** vecs_1;

        // stuff for coefficients
        type0* alpha_y;
        type0* dalpha_y;
        
        //stuff for error calculation
        type0* alph_err;
        
        void allocate();
        void deallocate();
        
        //functions to be accessed by run()
        void restart();
        void update_for_next();        
        bool interpolate();
        void interpolate_fail();
        void ord_dt(type0&);
        void err_calc();
        void err_fac_calc();        
        void start();        
    public:
        DMD_bdf_y(MAPP*,int,char**);
        ~DMD_bdf_y();

    };
}

#endif 
#endif