#ifdef Min_Style
    MinStyle(Min_lbfgs,l-bfgs)
#else
#ifndef __MAPP__min_lbfgs__
#define __MAPP__min_lbfgs__
#include "min.h"
#include "vec.h"
namespace MAPP_NS {
    class Min_lbfgs :public Min
    {
    private:
    protected:
        
        int m_it;

        
        type0* rho;
        type0* alpha;
        
        VecTens<type0>* s;
        VecTens<type0>* y;

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
