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
        int order;
        type0* beta;
        type0 ave_err;
        
        type0 delta_t;
        VecLst* vecs_comm;
        
        
        type0** y;
        type0* a;
        type0* g;
        type0* c0;
        type0* g0;
        type0* h;
        int c_n,c_d_n;
        int tot_dim;
        type0 solve(type0);
        void line_search();
        type0 gamma_red,slope;
        int no_steps;
        
        type0* t;
        type0* beta_mod;
        type0 tot_time;
        void interpolate(type0);
        
        
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
