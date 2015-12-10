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
        bool init_phase;
        type0 err_est(int,type0);
    protected:
        // stuff for book keeping
        type0* t;
        type0** y;
        type0* dy;
        Vec<type0>** vecs_1;

        // stuff for coefficients
        type0* alpha_y;
        type0* dalpha_y;
        
        //stuff for error calculation
        type0* alph_err;
        
        void allocate();
        void deallocate();
        
        //functions to be accessed by run()
        void restart(type0&,int&);
        void store_vecs(type0);        
        void interpolate(type0&,int&);
        void ord_dt(type0,type0,int,type0&,int&);
    public:
        DMD_mbdf(MAPP*,int,char**);
        ~DMD_mbdf();

    };
}

#endif 
#endif