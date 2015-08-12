#ifdef Clock_Style
    ClockStyle(Clock_cn,cn)
#else
#ifndef __MAPP__clock_cn__
#define __MAPP__clock_cn__
#include "clock.h"
namespace MAPP_NS {
    class Clock_cn :public Clock
    {
    private:
    protected:
        // inputs by user
        int max_step,max_order;
        type0 min_del_t,max_del_t,initial_del_t,max_t;

        // stuff for book keeping
        type0* t;
        type0* y;
        type0** dy;
        type0* dy0;
        
        void interpolate(type0&);
        void ord_dt(type0&,type0);
        
        int old_comm_mode;
        type0 old_skin;
        
    public:
        Clock_cn(MAPP *,int,char**);
        ~Clock_cn();
        void run();
        void init();
        void fin();
        
    };
}
#endif
#endif 
