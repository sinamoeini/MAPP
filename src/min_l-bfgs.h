#ifdef Min_Style
    MinStyle(Min_lbfgs,l-bfgs)
#else
#ifndef __MAPP__min_lbfgs__
#define __MAPP__min_lbfgs__
#include "min.h"
namespace MAPP_NS {
    class Min_lbfgs :public Min
    {
    private:
    protected:
        int f_prev_n;
        
        int m_it;
        int* s_list;
        int* y_list;
        
        type0* rho;
        type0* alpha;
        
        type0*** H_y;
        type0*** H_s;
        type0** f_H_prev;
    public:
        Min_lbfgs(MAPP *,int,char**);
        ~Min_lbfgs();
        void run();
        void init();
        void fin();

    };
}

#endif
#endif
