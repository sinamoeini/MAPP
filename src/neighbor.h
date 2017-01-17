#ifndef __MAPPneighbor
#define __MAPPneighbor
#include "global.h"
namespace MAPP_NS
{
    template<typename> class Vec;
    class Neighbor
    {
    private:
    protected:
        
        class Cell;
        int no_neigh_lists;
        Cell* cell;
        class Atoms*& atoms;
    public:
        bool pair_wise;
        int** neighbor_list;
        int* neighbor_list_size;
        int neighbor_list_size_size;
        int no_pairs;
        
        int** neighbor_list_2nd;
        int* neighbor_list_size_2nd;
        int neighbor_list_size_size_2nd;
        
        Neighbor(class Atoms*&);
        virtual ~Neighbor();
        
        virtual void create_list(bool)=0;
        virtual void init();
        virtual void fin();
        
        void print_stats();
        
        void mark_redndnt_ph(byte*);
        void rename_atoms(int*);
    };
}

namespace MAPP_NS
{
    class Neighbor::Cell
    {
    private:
        /*----------------------*/
        const int m;
        int& natms;
        int& natms_ph;
        
        class Atoms*& atoms;
        
        Vec<type0>*& x;
        type0 (&s_lo)[__dim__];
        type0 (&s_hi)[__dim__];
        type0 (&cut_s)[__dim__];
        /*----------------------*/
        
        //size: dim
        type0 cell_size[__dim__];
        int cell_denom[__dim__];
        int ncells_per_dim[__dim__];
        
        //size: ncells
        int ncells;
        int* head_atm;
        
        int nneighs;
        int* rel_neigh_lst;
        
        int* cell_vec;
        int* next_vec;
        
        /*----------------------*/
        
        int ineigh;
        int icell,jcell;
        
        
        void find_cell_no(type0*&,int&);
        void box_setup();
        
    protected:
    public:
        Cell(int,class Atoms*&);
        ~Cell();
        
        void create(bool);
        void destroy();
        
        void nxt_i();
        void nxt_j();
        
        int iatm,jatm;
        type0* ix;
        type0* jx;
        type0 rsq;
    };
}
#endif

