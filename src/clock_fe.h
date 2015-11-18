#ifdef Clock_Style
    ClockStyle(Clock_fe,fe)
#else
#ifndef __MAPP__clock_fe__
#define __MAPP__clock_fe__
#include "clock.h"
namespace MAPP_NS
{
    class Clock_fe :public ClockExplicit
    {
    private:
    protected:
        type0* y;
        type0* dy;
        void ord_dt(type0&,type0);
        void interpolate_n_err(type0&,type0&);
        
        void allocate();
        void deallocate();
        inline void fail_stp_adj(type0,type0&);
        inline void init_stp_adj(type0&);
    public:
        Clock_fe(MAPP *,int,char**);
        ~Clock_fe();
        void run();
        void init();
        void fin();
        
    };
}
#endif
#endif 
