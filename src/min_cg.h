#ifdef Min_Style
    MinStyle(Min_cg,cg)
#else
#ifndef __MAPP__min_cg__
#define __MAPP__min_cg__
#include "min.h"
namespace MAPP_NS {
    class Min_cg : public Min {
    private:
    protected:
        int f_prev_n;        
        
        type0** f_H_prev;
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
