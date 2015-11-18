#ifdef Clock_Style
    ClockStyle(Clock_adams,adams)
#else
#ifndef __MAPP__clock_adams__
#define __MAPP__clock_adams__
#include "clock.h"
namespace MAPP_NS {
    class Clock_adams :public ClockImplicit
    {
    private:
    protected:
        // stuff for book keeping
        type0* t;
        type0** dy;
        type0* y;
        
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
    public:
        Clock_adams(MAPP *,int,char**);
        ~Clock_adams();
        void run();
        void init();
        void fin();
    };
}

#endif
#endif
