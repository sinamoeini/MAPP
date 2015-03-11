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
        int max_iter,no_steps,max_order;
        type0 min_gamma,gamma_red,slope;
        type0 m_tol,a_tol,e_tol;
        type0 min_del_t,max_del_t,initial_del_t;
        type0 eq_ratio;
        
        // the atomic vectors indices
        int c_n,c_d_n;
        // dofs
        int dof_tot,dof_lcl;
        
        //predictor for c
        type0* y_0;

        
        // stuff for solution
        type0* a;
        type0* g0;
        type0* g;
        type0* h;
        type0* c0;
        type0 beta;
        
        // stuff for book keeping
        type0* t;
        type0** y;
        type0* dy;
        
        // stuff for coefficients
        type0 alpha_dy_0;
        type0* alpha_y;
        type0 dalpha_dy_0;
        type0* dalpha_y;
        
        
        //stuff for error calculation
        type0 err,err_prefac;
        type0 lwr_alpha_dy,lwr_alpha_y;
        type0* lwr_alpha;
        type0* e_n;
        
        
        type0 solve(type0,int);
        int interpolate(type0,int);
        void ord_dt(type0&,int&,int);

        
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