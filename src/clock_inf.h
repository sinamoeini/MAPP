#ifdef Clock_Style
    ClockStyle(Clock_inf,inf)
#else
#ifndef __MAPP__clock_inf__
#define __MAPP__clock_inf__
#include "clock.h"

namespace MAPP_NS
{
    
    enum
    {
        MIN_MOD_EN_DC_PROJ,
        MIN_MOD_FORM_EN,
        MIN_MOD_DC_NORM
    };
    
    
    class Clock_inf :public ClockImplicit
    {
    private:
        type0 epsilon,epsilon_sqrt,golden,max_a,h_norm;
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
        
        
    public:
        Clock_inf(MAPP *,int,char**);
        ~Clock_inf();
        void run();
        void init();
        void fin();
    };
}

#endif
#endif

/*
 
*/
