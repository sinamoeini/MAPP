#ifdef Clock_Style
    ClockStyle(Clock_fe,fe)
#else
#ifndef __MAPP__clock_fe__
#define __MAPP__clock_fe__

#include <stdio.h>
#include "clock.h"
namespace MAPP_NS {
    class Clock_fe :public Clock
    {
    private:
    protected:
        type0 delta_t;
        VecLst* vecs_comm;
        int dof_tot,dof_lcl;
        int c_n,c_d_n;
        int no_steps;
        
        type0 a_tol,e_tol;
        type0 min_del_t;
        type0 max_del_t;
        type0 eq_ratio;
        type0* y;
        type0* dy;
        void ord_dt(type0&);
        void solve(type0&);
        type0 err;
        
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
