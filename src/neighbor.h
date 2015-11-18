/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__neighbor__
#define __MAPP__neighbor__
#include "init.h"
namespace MAPP_NS
{
    class Neighbor:protected InitPtrs
    {
    private:
    protected:
        int* bin_denom_list;
        type0* bin_size;
        
        int* tot_bin_grid;
        int tot_bin;
        
        int** bin_neigh_list;
        int* bin_neigh_list_size;
        
        int* first_atom_bin;
        int first_atom_bin_size;
        
        int* atm_bin;
        int atm_bin_size;
        
        int* next_atm;
        int next_atm_size;

        int no_neigh_lists;
        int dim;
        
        void create_bin_list();
        void bin_atoms();
        void find_bin_no(int,int,int*&,int,int*,int*);
    public:
        int pair_wise;
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
        virtual void init()=0;
        virtual void fin();
        
        void print_stats();
    };
}
#endif

