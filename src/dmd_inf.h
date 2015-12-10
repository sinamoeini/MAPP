#ifdef DMD_Style
    DMDStyle(DMD_inf,inf)
#else
#ifndef __MAPP__dmd_inf__
#define __MAPP__dmd_inf__
#include "dmd.h"

namespace MAPP_NS
{
    enum
    {
        LS_BT,
        LS_GS,
        LS_BRENT
    };
    enum
    {
        MIN_MOD_EN_DC_PROJ,
        MIN_MOD_FORM_EN,
        MIN_MOD_DC_NORM
    };
    
    
    class DMD_inf :public DMDImplicit
    {
    private:
        type0 epsilon,epsilon_sqrt,golden,max_a,h_norm,prev_val,min_gamma,slope,gamma_red;
        int mod,max_step;
        type0* g_orig;
        type0* c1;
    protected:
        
        int old_comm_mode;
        type0 old_skin;
        void allocate();
        void deallocate();

        inline void ls_prep();
        inline type0 cost_func(type0);
        inline type0 dcost_func(type0,type0&);
        int line_search_gs(type0&,type0&,type0);
        int line_search_bt(type0&,type0&,type0);
        int line_search_brent(type0&,type0&,type0);
        int test0();
        int test1(type0,type0,type0);
        int ls_mode;
                
        void restart(type0&,int&){};
        void store_vecs(type0){};
        void interpolate(type0&,int&){};
        void ord_dt(type0,type0,int,type0&,int&){};
    public:
        DMD_inf(MAPP *,int,char**);
        ~DMD_inf();

    };
}

#endif
#endif

/*
 
*/
