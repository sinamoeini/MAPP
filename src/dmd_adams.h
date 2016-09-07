#ifdef DMD_Style
    DMDStyle(DMD_adams,adams)
#else
#ifndef __MAPP__dmd_adams__
#define __MAPP__dmd_adams__
#include "dmd.h"
namespace MAPP_NS
{
    class DMD_adams :public DMDImplicit
    {
    private:
        type0** xi; type0** wi;
        void ratio_calc(type0&,type0&);
    protected:
        // stuff for book keeping
        type0* t;
        type0** dy;
        type0* y;
        type0* e_n;
        Vec<type0>** vecs_1;
        
        // stuff for coefficients
        type0 alpha_y;
        type0* alpha_dy;
        type0 dalpha_y;
        type0* dalpha_dy;
        
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
        DMD_adams(int,char**);
        ~DMD_adams();
    };
}

#endif
#endif
