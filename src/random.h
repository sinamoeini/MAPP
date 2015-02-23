/*--------------------------------------------
 Created by Sina on 07/30/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__rand_engine__
#define __MAPP__rand_engine__

#include <iostream>
#include "init.h"
namespace MAPP_NS {
    class Random : protected InitPtrs{
    private:
        type0 reserved;
        int gauss_chk;
        int seed;
    protected:
    public:
        Random(MAPP *,int);
        ~Random();
        type0 uniform();
        type0 gaussian();
        type0 gaussian(type0,type0);
    };
}
#endif 

