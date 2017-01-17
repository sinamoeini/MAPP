#ifndef __MAPP__neighbor_cna__
#define __MAPP__neighbor_cna__
#include "neighbor.h"
namespace MAPP_NS
{
    class NeighborCNA:public Neighbor
    {
    private:
        Vec<atom_type>*& type_vec;
        type0**& cut_sk_sq;
    protected:
    public:
        NeighborCNA(class Atoms*&,Vec<atom_type>*&,type0**&);
        ~NeighborCNA();
        
        void create_list(bool);
        void init();
        void fin();
    };
}
#endif
