/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__neighbor__
#define __MAPP__neighbor__
#include "init.h"
#define DEBUG_NEIGH
namespace MAPP_NS
{
    class Neighbor:protected InitPtrs
    {
    private:
    protected:
        
        class Cell
        {
        private:
#ifdef DEBUG_NEIGH
            Error*& error;
#endif
            /*----------------------*/
            const int dim;
            const int m;
            int& natms;
            int& natms_ph;
            
            Vec<type0>*& x;
            type0*& s_lo;
            type0*& s_hi;
            type0*& cut_s;
            /*----------------------*/
            
            //size: dim
            type0* cell_size;
            int* cell_denom;
            int* ncells_per_dim;
            
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
            Cell(MAPP*,int,type0*&);
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
        
        Neighbor(MAPP*);
        ~Neighbor();
        
        virtual void create_list(bool)=0;
        virtual void init();
        virtual void fin();
        
        void print_stats();
    };
}
#endif

