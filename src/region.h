#ifndef __MAPP__region__
#define __MAPP__region__
#include "init.h"

namespace MAPP_NS
{
    class Region:protected InitPtrs
    {
    private:
    protected:
    public:
        Region(MAPP*);
        virtual ~Region();
        virtual int belong(type0*)=0;

    };
}
#endif
