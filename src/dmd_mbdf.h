#ifdef DMD_Style
    DMDStyle(DMD_mbdf,mbdf)
#else
#ifndef __MAPP__dmd_mbdf__
#define __MAPP__dmd_mbdf__
#include "dmd.h"
namespace MAPP_NS {
    class DMD_mbdf :public DMDImplicit
    {
    private:
    protected:
        // stuff for book keeping
        type0* t;
        type0** y;
        type0* dy;
        Vec<type0>** vecs_1;
        void reset();

        // stuff for coefficients
        type0* alpha_y;
        type0* dalpha_y;
        
        //stuff for error calculation
        type0* alph_err;
        
        void interpolate(int&,int&,type0&,int&);
        void ord_dt(int&,int&,type0&,int&);
        type0 err_est(int,type0);
        
        void allocate();
        void deallocate();
        inline void fail_stp_adj(int&,int&,type0,type0&,int&);
        inline void init_stp_adj(type0&);
        inline type0 precond_rat_adj(type0);
        type0 est_dt();
    public:
        DMD_mbdf(MAPP *,int,char**);
        ~DMD_mbdf();
        void run();
        void init();
        void fin();
    };
}

#endif 
#endif