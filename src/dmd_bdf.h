#ifdef DMD_Style
    DMDStyle(DMD_bdf,bdf)
#else
#ifndef __MAPP__dmd_bdf__
#define __MAPP__dmd_bdf__
#include "dmd.h"
namespace MAPP_NS
{
    class DMD_bdf :public DMDImplicit
    {
    private:
        
        void ratio_calc(int,type0,type0&,type0&);
    protected:
        // stuff for book keeping
        type0* t;
        type0** y;
        type0* dy;
        type0* e_n;
        Vec<type0>** vecs_1;
        
        // stuff for coefficients
        type0 alpha_dy;
        type0* alpha_y;
        type0 dalpha_dy;
        type0* dalpha_y;
        
        void allocate();
        void deallocate();
        
        //functions to be accessed by run()
        void reset();
        type0 est_dt();
        void interpolate(type0&,int&);
        void ord_dt(type0,type0&,int&);
        void fail_stp_adj(type0,type0,type0&,int&);        
    public:
        DMD_bdf(MAPP*,int,char**);
        ~DMD_bdf();
        void run();
        void init();
        void fin();
    };
}

#endif 
#endif 

