#ifdef Min_Style
    MinStyle(Min_cg,cg)
#else
#ifndef __MAPP__min_cg__
#define __MAPP__min_cg__

#include <iostream>
#include "min.h"
namespace MAPP_NS {
    class Min_cg : public Min {
    private:
    protected:
        int f_prev_n;
        int x_prev_n;
        int x_dim;
        int f_n;
        int h_n;
        
        
        type0** H_prev;
        type0** B_prev;
        type0** f_H_prev;
        type0** f_H;
        type0** h_H;
    public:
        Min_cg(MAPP *,int,char**);
        ~Min_cg();
        void run();
        void init();
        void fin();
        
    };

}
#endif
#endif
