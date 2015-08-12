#ifdef Clock_Style
    ClockStyle(Clock_adams,adams)
#else
#ifndef __MAPP__clock_adams__
#define __MAPP__clock_adams__
#include "clock.h"
namespace MAPP_NS {
    class Clock_adams :public Clock
    {
    private:
    protected:
        // inputs by user
        int max_step,max_order;
        type0 min_del_t,max_del_t,initial_del_t,max_t;
        
        
        
        
        
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
        type0 err_coef_dy,err_coef_y;
        type0* err_coef;
        type0* e_n;
        
        
        type0 solve(type0,int,type0&);
        void interpolate(type0&,int&);
        void ord_dt(type0&,int&,int,type0);
        
        int old_comm_mode;
        type0 old_skin;
        
        
        
        type0** xi;
        type0** wi;
        void quad(int,type0*,type0*);
        
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
