#ifndef __MAPP__clock__
#define __MAPP__clock__
#include "init.h"
#include "atoms.h"
namespace MAPP_NS {
    enum
    {
        LS_BT,
        LS_GS
    };
    
    class Clock : protected InitPtrs
    {
    private:
        int ns_alloc;
        type0 epsilon;
        type0 golden;
    protected:
        VecLst* vecs_comm;

        ThermoDynamics* thermo;
        int fe_idx;
        int stress_idx;
        int time_idx;
        int cdof_n;
        int ls_mode;
        char* cdof;
        type0* nrgy_strss;
        type0 tot_t;

        
        // variables & vectors for implicit integrator
        type0 beta;
        type0* a;
        type0* y_0;
        type0* g;
        type0* h;
        type0* g0;
        type0* c0;
        
        // parameters for implicit integrator
        type0 a_tol,err_prefac;
        
        //parameters for the solver
        int max_iter;
        type0 m_tol,min_gamma,slope,gamma_red;
        
        // the atomic vectors indices
        int c_n,c_d_n;
        // dofs
        int dof_tot,dof_lcl;
        
        inline type0 cost_func(type0);
        int line_search_gs(type0&,type0&,type0);
        int line_search_bt(type0&,type0&,type0);
        //int test(type0,type0,type0);
        
        void rectify(type0*);
        void solve_n_err(type0&,type0&);
        type0 prev_val;
        int pre_cond;
    public:
        Clock(MAPP *);
        virtual ~Clock();
        void errors();
        virtual void run()=0;
        virtual void init()=0;
        virtual void fin()=0;
        



    };
    
    
}

#endif
