#ifndef __MAPP__md__
#define __MAPP__md__
#include "type_def.h"
namespace MAPP_NS
{
    class MD 
    {
    private:
    protected:
        int pe_idx;
        int ke_idx;
        int temp_idx;
        int stress_idx;
        type0*& nrgy_strss;
        class VecLst* vecs_comm;
        class ThermoDynamics* thermo;
    public:
        
        MD();
        virtual ~MD()=0;
        virtual void init()=0;
        virtual void fin()=0;
        virtual void run(int)=0;

        type0 dt,boltz,hplanck;
    };
}
#endif
