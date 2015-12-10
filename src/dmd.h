#ifndef __MAPP__dmd__
#define __MAPP__dmd__
#include "init.h"
#include "ff.h"
#include "neighbor_dmd.h"
#include "thermo_dynamics.h"
#include "min.h"
#include "ls_styles.h"
namespace MAPP_NS
{
    enum
    {
        NO_FLAG,
        F_FLAG,
        CD_FLAG,
        STEP_FLAG
    };
    
    class DMD : protected InitPtrs
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
        int ncs;
        type0 nc_dofs;
        type0 old_skin;
        void rectify(type0*);
        
        virtual void print_stats();
        
        int max_step;
        type0 min_del_t,max_t;
        ThermoDynamics* thermo;
        
        Min* min;
        int nmin;
        type0 f_tol;
        type0 f_norm0;
        type0 c_d_norm;
        int istep;
        type0 cd_tol;
        int step_tol;
        int dstep;
        int min_flag;
        void reset();
        type0 calc_nc_dofs();
        bool decide_min(int&,type0&);
        void do_min();

        type0 max_succ_dt;
        virtual void run()=0;
    public:
        DMD(MAPP *);
        virtual ~DMD();
        void run(type0);
        virtual void init();
        virtual void fin();
        void coef(int,char**);
        void dmd_min(int,char**);
    };

    class DMDImplicit: public DMD
    {
    private:
        type0 g_h;
        type0 h_h;
        
        LineSearch<DMDImplicit>* ls_dmd;
    protected:
        int max_order;
        // variables & vectors for implicit integrator
        type0 beta;
        Vec<type0>** vecs_0;
        type0* a;
        type0* y_0;
        type0* g;
        type0* h;
        type0* g0;
        type0* c0;
        type0* c1;
        
        // parameters for implicit integrator
        type0 err_prefac;

        int max_iter;
        int iter_dcr_cntr,iter_dcr_thrsh;
        type0 m_tol,max_a;
        

        int solve_n_err(type0&,type0&);
        int const_stps;


        void print_stats();
        int max_succ_q;
        int solve_rej;
        int solve_acc;
        int intg_rej;
        int intp_rej;
        

        int test();


        void reset();
        void ord_dt(type0,type0&,int&);
        void fail_stp_adj(type0,type0,type0&,int&);
        
        // allocation and deallocation functions
        // specific to each style
        // pure virtual
        virtual void allocate()=0;
        virtual void deallocate()=0;
        
        // functions to be accessed by run
        // specific to each style
        // pure virtual
        virtual void restart(type0&,int&)=0;
        virtual void store_vecs(type0)=0;
        virtual void interpolate(type0&,int&)=0;
        virtual void ord_dt(type0,type0,int,type0&,int&)=0;
        
        void run();
    public:
        DMDImplicit(MAPP*);
        virtual ~DMDImplicit();
        
        void init();
        void fin();
        
        type0 F(type0);
        type0 dF(type0,type0&);
        void ls_prep(type0&,type0&,type0&);
        void F_reset();
    };
    
    
    class DMDExplicit: public DMD
    {
    private:
    protected:
        void reset();
        void print_stats();
        int intg_rej;
        int intp_rej;
        virtual type0 est_dt()=0;
               
    public:
        DMDExplicit(MAPP*);
        ~DMDExplicit();
        virtual void run()=0;
        virtual void init();
        virtual void fin();
    };
}

#endif
