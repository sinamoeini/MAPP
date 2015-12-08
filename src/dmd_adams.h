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
        void ratio_calc(int,type0,type0&,type0&);
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
        void reset();
        type0 est_dt();
        void interpolate(type0&,int&);
        void ord_dt(type0,type0&,int&);
        void fail_stp_adj(type0,type0,type0&,int&);        
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
