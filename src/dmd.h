#ifndef __MAPP__dmd__
#define __MAPP__dmd__
#include "init.h"
#include "ff.h"
#include "neighbor_dmd.h"
#include "thermo_dynamics.h"
#include "min.h"
#include "gmres.h"
namespace MAPP_NS
{
    enum
    {
        NO_FLAG,
        F_FLAG,
        CD_FLAG,
        C_FLAG,
        STEP_FLAG
    };
    
    class DMD : protected InitPtrs
    {
    private:
    protected:
        VecLst* vecs_comm;
        ForceFieldDMD* forcefield_dmd;
        Neighbor_dmd* neighbor_dmd;

        //int msd_idx;
        int fe_idx;
        int stress_idx;
        int time_idx;
        type0*& nrgy_strss;
        type0 eps,eps_sqr;
        int c_dim;
        
        // dofs
        int ncs;
        type0 nc_dofs;
        type0 old_skin;
        void rectify(type0*);
        
        virtual void print_stats();
        
        int max_step;
        type0 t_cur,a_tol,dt_min,t_fin;
        type0 inf;
        ThermoDynamics* thermo;
        
        Min* min;
        int nmin;
        type0 f_tol,delta_c;
        type0 f_norm0;
        type0 c_d_norm,max_c_d;
        int istep;
        type0 cd_tol,c_tol;
        int step_tol;
        int dstep;
        int min_flag;
        void reset();
        type0 calc_nc_dofs();
        bool decide_min(int&,type0&);
        void do_min();

        type0 max_succ_dt;
        virtual void run()=0;
        type0 vac_msd();
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
        
    protected:
        type0 dt,dt_p,dt_new;
        int dq,q,q_p,q_max;
        int intp_failure;
        
        int const_dt,const_q;
        
        // variables & vectors for implicit integrator
        type0 beta;
        type0 beta_inv;
        
        Vec<type0>** vecs_0;
        type0* a;
        type0* y_0;
        type0* F;
        type0* del_c;
        type0* c0;

        
        // parameters for implicit integrator
        type0 err,err_fac;

        int max_iter;

        bool solve_non_lin();
        
        


        void print_stats();
        int max_succ_q;
        int solve_rej;
        int solve_acc;
        int intg_rej;
        int intp_rej;
        



        void reset();
        void ord_dt();
        void intg_fail();
        void nonl_fail();
        
        // allocation and deallocation functions
        // specific to each style
        // pure virtual
        virtual void allocate()=0;
        virtual void deallocate()=0;
        
        // functions to be accessed by run
        // specific to each style
        // pure virtual
        virtual void restart()=0;
        virtual void update_for_next()=0;
        virtual bool interpolate()=0;
        virtual void interpolate_fail()=0;
        virtual void err_fac_calc()=0;
        virtual void ord_dt(type0&)=0;
        virtual void err_calc()=0;
        virtual void start()=0;
        
        void run();
        
        GMRES<type0,ForceFieldDMD>* gmres;
    public:
        DMDImplicit(MAPP*);
        virtual ~DMDImplicit();
        
        void init();
        void fin();
        

    };
    
    
    class DMDExplicit: public DMD
    {
    private:
    protected:
        void reset();
        void print_stats();
        int intg_rej;
        int intp_rej;
        
        // allocation and deallocation functions
        // specific to each style
        // pure virtual
        virtual void allocate()=0;
        virtual void deallocate()=0;
        
        virtual void restart(type0&,int&)=0;
        virtual void store_vecs(type0)=0;
        virtual void ord_dt(type0&,type0)=0;
        virtual void interpolate_n_err(type0&,type0&)=0;
        
        void fail_stp_adj(type0,type0&);
        void run();
    public:
        DMDExplicit(MAPP*);
        ~DMDExplicit();
        void init();
        void fin();
    };
}

#endif
