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
        TYPE0 delta_t;
        VecLst* vecs_comm;
        int dof_tot,dof_lcl;
        int c_n,c_d_n;
        int no_steps;
        
        TYPE0 a_tol,e_tol;
        TYPE0 min_del_t;
        TYPE0 max_del_t;
        TYPE0 eq_ratio;
        TYPE0* y;
        TYPE0* dy;
        void ord_dt(TYPE0&);
        void solve(TYPE0&);
        TYPE0 err;
        
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
