#ifdef Clock_Style
    ClockStyle(Clock_mbdf,mbdf)
#else
#ifndef __MAPP__clock_mbdf__
#define __MAPP__clock_mbdf__
#include "clock.h"
namespace MAPP_NS {
    class Clock_mbdf :public Clock
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
        type0* alpha_y;
        type0* dalpha_y;
        
        //stuff for error calculation
        type0 err,err_prefac;
        type0* alph_err;
        
        //number of constant steps
        int const_stps;
        
        type0 solve(type0,int);
        int interpolate(type0,int);
        void ord_dt(type0&,int&,int);
        type0 fac(int);
        type0 step_size(type0,int);
        type0 err_est(int);

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