#ifdef DMD_Style
    DMDStyle(DMD_fe,fe)
#else
#ifndef __MAPP__dmd_fe__
#define __MAPP__dmd_fe__
#include "dmd.h"
namespace MAPP_NS
{
    class DMD_fe :public DMDExplicit
    {
    private:
    protected:
        type0* y;
        type0* dy;
        Vec<type0>** vecs_1;
        void reset();
        
        void ord_dt(type0&,type0);
        void interpolate_n_err(type0&,type0&);
        
        void allocate();
        void deallocate();
        inline void init_stp_adj(type0&);
        void restart(type0&,int&);
        void store_vecs(type0);
        type0 est_dt();
        void fail_stp_adj(type0,type0&); 
    public:
        DMD_fe(MAPP *,int,char**);
        ~DMD_fe();
        void run();
        void init();
        void fin();
        
    };
}
#endif
#endif 
