#ifdef Clock_Style
    ClockStyle(Clock_mbdf,mbdf)
#else
#ifndef __MAPP__clock_mbdf__
#define __MAPP__clock_mbdf__
#include "clock.h"
namespace MAPP_NS {
    class Clock_mbdf :public ClockImplicit
    {
    private:
    protected:
        // stuff for book keeping
        type0* t;
        type0** y;
        type0* dy;

        // stuff for coefficients
        type0* alpha_y;
        type0* dalpha_y;
        
        //stuff for error calculation
        type0* alph_err;
        
        void interpolate(int&,int&,type0&,int&);
        void ord_dt(int&,int&,type0&,int&);
        type0 err_est(int,type0);
        
        void allocate();
        void deallocate();
        inline void fail_stp_adj(int&,int&,type0,type0&,int&);
        inline void init_stp_adj(type0&);
        inline type0 precond_rat_adj(type0);
    public:
        Clock_mbdf(MAPP *,int,char**);
        ~Clock_mbdf();
        void run();
        void init();
        void fin();
    };
}

#endif 
#endif