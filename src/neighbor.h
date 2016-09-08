#ifndef __MAPP__neighbor__
#define __MAPP__neighbor__
#include "global.h"
namespace MAPP_NS
{
    template<typename> class Vec;
    class Neighbor
    {
    private:
    protected:
        
        class Cell
        {
        private:
            /*----------------------*/
            const int m;
            int& natms;
            int& natms_ph;
            
            Vec<type0>*& x;
            type0*& s_lo;
            type0*& s_hi;
            type0*& cut_s;
            /*----------------------*/
            
            //size: dim
            type0 cell_size[dimension];
            int cell_denom[dimension];
            int ncells_per_dim[dimension];
            
            //size: ncells
            int ncells;
            int* head_atm;
            
            int nneighs;
            int* rel_neigh_lst;
            
            int* cell_vec;
            int* next_vec;
            
            /*----------------------*/
            
            int curr_cell;
            
            int ineigh;
            int icell,jcell;
            
            
            void find_cell_no(type0*&,int&);
            void box_setup();

        protected:
        public:
            Cell(int,type0*&);
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
        
        int no_neigh_lists;
        Cell* cell;
        
    public:
        bool pair_wise;
        int** neighbor_list;
        int* neighbor_list_size;
        int neighbor_list_size_size;
        int no_pairs;
        
        int** neighbor_list_2nd;
        int* neighbor_list_size_2nd;
        int neighbor_list_size_size_2nd;
        
        Neighbor();
        virtual ~Neighbor();
        
        virtual void create_list(bool)=0;
        virtual void init();
        virtual void fin();
        
        void print_stats();
        
        void mark_redndnt_ph(byte*);
        void rename_atoms(int*);
    };
}
#endif

