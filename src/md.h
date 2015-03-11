#ifndef __MAPP__md__
#define __MAPP__md__
#include "init.h"
#include "atoms.h"
#include "thermo_dynamics.h"
namespace MAPP_NS {
    
    class MD : protected InitPtrs
    {
    private:
        int ns_alloc;
    protected:
        
        class ThermoDynamics* thermo;
        int pe_idx;
        int ke_idx;
        int temp_idx;
        int stress_idx;
        type0* nrgy_strss;
        VecLst* vecs_comm;
    public:
        
        MD(MAPP*);
        virtual ~MD()=0;
        virtual void init()=0;
        virtual void fin()=0;
        virtual void run(int)=0;

        type0 dt,boltz;
    };
}

#endif
