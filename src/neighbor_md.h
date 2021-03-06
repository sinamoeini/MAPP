#ifndef __MAPP__neighbor_md__
#define __MAPP__neighbor_md__
#include "neighbor.h"
namespace MAPP_NS
{
    class NeighborMD:public Neighbor
    {
    private:
        Vec<atom_type>*& type_vec;
        type0**& cut_sk_sq;
    protected:
    public:
        NeighborMD(class Atoms*&,Vec<atom_type>*&,type0**&);
        ~NeighborMD();
        
        void create_list(bool);
        void init();
        void fin();
    };
}
#endif
