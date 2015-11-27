#ifdef DMD_Style
    DMDStyle(DMD_adams,adams)
#else
#ifndef __MAPP__dmd_adams__
#define __MAPP__dmd_adams__
#include "dmd.h"
namespace MAPP_NS {
    class DMD_adams :public DMDImplicit
    {
    private:
    protected:
        // stuff for book keeping
        type0* t;
        type0** dy;
        type0* y;
        Vec<type0>** vecs_1;
        void reset();
        
        // stuff for coefficients
        type0 alpha_y;
        type0* alpha_dy;
        type0 dalpha_y;
        type0* dalpha_dy;
        
        //stuff for error calculation
        type0* e_n;
        
        void interpolate(type0&,int&);
        void ord_dt(type0&,int&,type0);
        void ratio_calc(int,type0,type0&,type0&);
        void allocate();
        void deallocate();
        inline void fail_stp_adj(type0,type0&,int&);
        inline void init_stp_adj(type0&);
        
        type0** xi;
        type0** wi;
        type0 est_dt();        
    public:
        DMD_adams(MAPP *,int,char**);
        ~DMD_adams();
        void run();
        void init();
        void fin();
    };
}

#endif
#endif
