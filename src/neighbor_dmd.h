/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__neighbor_dmd__
#define __MAPP__neighbor_dmd__
#include "neighbor.h"
namespace MAPP_NS
{
    class Neighbor_dmd:public Neighbor
    {
    private:
    protected:
    public:
        Neighbor_dmd(MAPP*);
        ~Neighbor_dmd();
        
        void create_list(bool);
        void create_2nd_list();
        void init();
        void fin();
        
    };
}
#endif
