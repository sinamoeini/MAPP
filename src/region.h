#ifndef __MAPP__region__
#define __MAPP__region__
#include <iostream>
#include "init.h"
using namespace std;
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
