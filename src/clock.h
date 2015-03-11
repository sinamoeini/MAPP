#ifndef __MAPP__clock__
#define __MAPP__clock__
#include "init.h"
#include "atoms.h"
namespace MAPP_NS {
    class Clock : protected InitPtrs
    {
    private:
        int ns_alloc;
    protected:
        VecLst* vecs_comm;

        ThermoDynamics* thermo;
        int fe_idx;
        int stress_idx;
        int time_idx;
        int cdof_n;
        char* cdof;
        type0* nrgy_strss;
    public:
        Clock(MAPP *);
        virtual ~Clock();
        void errors();
        virtual void run()=0;
        virtual void init()=0;
        virtual void fin()=0;
        void rectify(type0*);
        
    };
    
    
}

#endif
