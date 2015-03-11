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
        int c_n,c_d_n,dof_tot,dof_lcl;
        int max_iter,no_steps;
        type0 min_gamma,gamma_red,slope;
        type0 m_tol,a_tol,e_tol;
        type0 min_del_t,max_del_t,initial_del_t;
        type0 eq_ratio;

        type0 beta;
        type0 err,err_prefac;
        
        
        type0* t;
        
        type0** dy;
        
        type0* y_0;
        type0* y;
        type0* dy0;
        type0* a;
        type0* g;
        type0* c0;
        type0* g0;
        type0* h;

        type0 solve(type0);
        int interpolate(type0);
        void ord_dt(type0&);
        
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
