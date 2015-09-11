#ifdef Clock_Style
    ClockStyle(Clock_bdf,bdf)
#else
#ifndef __MAPP__clock_bdf__
#define __MAPP__clock_bdf__
#include "clock.h"
namespace MAPP_NS {
    class Clock_bdf :public Clock
    {
    private:
    protected:
        // inputs by user
        int max_step,max_order;
        type0 min_del_t,initial_del_t,max_t;

        
        // stuff for book keeping
        type0* t;
        type0** y;
        type0* dy;
        
        // stuff for coefficients
        type0 alpha_dy;
        type0* alpha_y;
        type0 dalpha_dy;
        type0* dalpha_y;
        
        
        //stuff for error calculation
        type0* e_n;
        
        
        void interpolate(type0&,int&);
        void ord_dt(type0&,int&,type0);
        void ratio_calc(int,type0,type0&,type0&);
        
        int old_comm_mode;
        type0 old_skin;
        inline void fail_stp_adj(type0,type0&,int&);
        inline void init_stp_adj(type0&);        
    public:
        Clock_bdf(MAPP *,int,char**);
        ~Clock_bdf();
        void run();
        void init();
        void fin();
    };
}

#endif 
#endif 