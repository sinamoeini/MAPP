#ifdef Clock_Style
    ClockStyle(Clock_bdf,bdf)
#else
#ifndef __MAPP__clock_bdf__
#define __MAPP__clock_bdf__

#include <stdio.h>
#include "clock.h"
namespace MAPP_NS {
    class Clock_bdf :public Clock
    {
    private:
    protected:
        // inputs by user
        int max_iter,no_steps,max_order;
        TYPE0 min_gamma,gamma_red,slope;
        TYPE0 m_tol,a_tol,e_tol;
        TYPE0 min_del_t,max_del_t,initial_del_t;
        TYPE0 eq_ratio;
        
        // the atomic vectors indices
        int c_n,c_d_n;
        // dofs
        int dof_tot,dof_lcl;
        
        //predictor for c
        TYPE0* y_0;

        
        // stuff for solution
        TYPE0* a;
        TYPE0* g0;
        TYPE0* g;
        TYPE0* h;
        TYPE0* c0;
        TYPE0 beta;
        
        // stuff for book keeping
        TYPE0* t;
        TYPE0** y;
        TYPE0* dy;
        
        // stuff for coefficients
        TYPE0 alpha_dy_0;
        TYPE0* alpha_y;
        TYPE0 dalpha_dy_0;
        TYPE0* dalpha_y;
        
        
        //stuff for error calculation
        TYPE0 err,err_prefac;
        TYPE0 lwr_alpha_dy,lwr_alpha_y;
        TYPE0* lwr_alpha;
        TYPE0* e_n;
        
        
        TYPE0 solve(TYPE0,int);
        int interpolate(TYPE0,int);
        void ord_dt(TYPE0&,int&,int);

        
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