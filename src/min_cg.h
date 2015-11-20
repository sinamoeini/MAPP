#ifdef Min_Style
    MinStyle(Min_cg,cg)
#else
#ifndef __MAPP__min_cg__
#define __MAPP__min_cg__
#include "min.h"
#include "vec.h"
namespace MAPP_NS {
    class Min_cg : public Min {
    private:
    protected:
        vvec<type0> h;
        vvec<type0> x;
        vvec<type0> x0;
        vvec<type0> f;
        vvec<type0> f0;
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
