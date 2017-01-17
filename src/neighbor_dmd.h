#ifndef __MAPP__neighbor_dmd__
#define __MAPP__neighbor_dmd__
#include "neighbor.h"
namespace MAPP_NS
{
    class NeighborDMD:public Neighbor
    {
    private:
        Vec<atom_type>*& ctype_vec;
        Vec<type0>*& c_vec;
        type0**& cut_sk_sq;
        type0*& rsq_crd;
    protected:
    public:
        NeighborDMD(class Atoms*&,Vec<atom_type>*&,Vec<type0>*& c_vec,type0**&,type0*&);
        ~NeighborDMD();
        
        void create_list(bool);
        void create_2nd_list();
        void init();
        void fin();
        
    };
}
#endif
