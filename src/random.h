#ifndef __MAPP__rand_engine__
#define __MAPP__rand_engine__
#include "type_def.h"
namespace MAPP_NS
{
    class Random
    {
    private:
        type0 reserved;
        int gauss_chk;
        int seed;
    protected:
    public:
        Random(int);
        ~Random();
        type0 uniform();
        type0 gaussian();
        type0 gaussian(type0,type0);
    };
}
#endif 
