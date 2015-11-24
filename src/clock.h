#ifndef __MAPP__clock__
#define __MAPP__clock__
#include "init.h"
#include "ff.h"
#include "neighbor_dmd.h"
#include "thermo_dynamics.h"
#include "min.h"
namespace MAPP_NS
{
    enum
    {
        LS_BT,
        LS_GS,
        LS_BRENT
    };
    
    class Clock : protected InitPtrs
    {
    private:
    protected:
        VecLst* vecs_comm;
        ForceFieldDMD* forcefield_dmd;
        Neighbor_dmd* neighbor_dmd;

        
        int fe_idx;
        int stress_idx;
        int time_idx;
        type0*& nrgy_strss;
        type0 tot_t,a_tol;
        int c_dim;
        
        // dofs
        int dof_tot,dof_lcl;
        type0 old_skin;
        void rectify(type0*);
        
        virtual void print_stats();
        
        int max_step;
        type0 min_del_t,initial_del_t,max_t;
        ThermoDynamics* thermo;
        
        Min* min;
        int nmin;
        type0 f_tol;
        type0 init_f_norm;
    public:
        Clock(MAPP *);
        virtual ~Clock();
        virtual void run()=0;
        virtual void init();
        virtual void fin();
        void coef(int,char**);
    };

    class ClockImplicit: public Clock
    {
    private:
        type0 epsilon;
        type0 golden;
    protected:
        int max_order;
        // variables & vectors for implicit integrator
        type0 beta;
        Vec<type0>** vecs_0;
        type0* a;
        type0* y_0;
        type0* y_1;
        type0* g;
        type0* h;
        type0* g0;
        type0* c0;
        
        // parameters for implicit integrator
        type0 err_prefac;
        
        //parameters for the solver
        type0 prev_val;
        int pre_cond;
        int max_iter;
        type0 m_tol,min_gamma,slope,gamma_red;
        
        int ls_mode;

        inline type0 cost_func(type0);
        int line_search_gs(type0&,type0&,type0);
        int line_search_bt(type0&,type0&,type0);
        int line_search_brent(type0&,type0&,type0);
        void solve_n_err(type0&,type0&);

        
        virtual void allocate();
        virtual void deallocate();
        
        virtual void init();
        virtual void fin();
        void print_stats();
        int max_succ_q;
        type0 max_succ_dt;
        int solve_rej;
        int solve_acc;
        int intg_rej;
        int intp_rej;

        void reset();
        int test(type0,type0,type0);

    public:
        ClockImplicit(MAPP*);
        virtual ~ClockImplicit();
    };
    
    
    class ClockExplicit: public Clock
    {
    private:
    protected:
        void print_stats();
        type0 max_succ_dt;
        int intg_rej;
        int intp_rej;
    public:
        ClockExplicit(MAPP*);
        ~ClockExplicit();
        virtual void run()=0;
        virtual void init();
        virtual void fin();
    };
}

#endif
