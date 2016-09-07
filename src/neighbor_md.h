#ifndef __MAPP__neighbor_md__
#define __MAPP__neighbor_md__
#include "neighbor.h"
namespace MAPP_NS
{
    class Neighbor_md:public Neighbor
    {
    private:
    protected:
    public:
        Neighbor_md();
        ~Neighbor_md();
        
        void create_list(bool);
        void init();
        void fin();
    };
}
#endif
