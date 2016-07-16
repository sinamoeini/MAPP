/*--------------------------------------------
 Created by Sina on 06/29/16.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__gcmc__
#define __MAPP__gcmc__
#include "init.h"
namespace MAPP_NS
{
    enum{INS_MODE,DEL_MODE};
    /*--------------------------------------------
     allocation for this constructor has 3 levels:
     
     0. buffers that have constant size during the
     life of this object. They are allocated in
     constructor and destroyed destructor:
     
     cut_s;
     s_lo_ph;
     s_hi_ph;
     cell_size;
     ncells_per_dim;
     cell_denom;
     icell_coord;
     jcell_coord;
     nimages_per_dim; 
     *nimages_per_dim;
     ins_s_trials; (!! NOT *ins_s_trials)
     rel_neigh_lst_coord;
     rel_neigh_lst;
     
     1. buffers whose size are dependent on box 
     dimensions and domain dimensions:
     
     head_atm;
     *ins_s_trials;
     ins_cell;
     ins_cell_coord;
     ins_buff;
     ins_cell;
     del_lst;
     
     2. buffers whose size are decided on fly
     
     del_ids;
     
     --------------------------------------------*/
    
    
    class GCMC:protected InitPtrs
    {
    private:
        enum{LCL_MODE,PH_MODE};
        int atm_mode;
        const int x_dim;
        const int dim;
        const int m;
        
        
        md_type gas_type;
        int igas,gas_id,ngas,ngas_before,tot_ngas;
        type0 vol;
        //constants
        type0 gas_mass,beta,kbT,T,mu,lambda,sigma,z_fac;
        
        int& natms;
        int& natms_ph;
        type0 cut;
        type0**& cut_sq;
        type0*& s_lo;
        type0*& s_hi;
        
        // size dim
        type0* vel_buff;
        type0* cut_s;
        type0* s_lo_ph;
        type0* s_hi_ph;
        type0* cell_size;
        int* ncells_per_dim;
        int* cell_denom;
        int* icell_coord;
        int* jcell_coord;
        int** nimages_per_dim;
        type0** ins_s_trials;
        int* rel_neigh_lst_coord;
        int* rel_neigh_lst;
        
        
        int* head_atm;
        //type0* *ins_s_trials;
        int* ins_cell;
        int* ins_cell_coord;
        type0* ins_buff;
        int* del_lst;
        
        int* del_ids;
        int del_ids_sz,del_ids_cpcty;
        int max_id;
        void add_del_id();
        int get_new_id();
        
        // self interaction
        int iself;
        
        
        int ineigh;
        int nneighs;
        
        
        int icell,jcell;
        int ncells;
        
    
        Vec<int>* cell_vec_p;
        Vec<int>* next_vec_p;
        
        void find_cell_no(type0*&,int&);
        bool lcl(type0*&);
        
        void next_iatm_ins();
        void next_iatm_del();
        void next_jcell_all();
        void next_jcell_lcl();
        void next_jatm_all();
        void next_jatm_self();
        void next_jatm_lcl();
        void (GCMC::*next_jatm_p)();
        class Random* random;
        
        void refresh();
        
        void ins_attmpt();
        void ins_succ();

        void del_attmpt();
        void del_succ();
        
        
        int itrial_atm;
        int ntrial_atms;
        void box_setup();
        
    protected:
    public:
        GCMC(MAPP*,int,dmd_type,type0,type0,int);
        ~GCMC();
        
        void init();
        void fin();
        
        
        void xchng(bool,int);
        

        void next_iatm();
        void next_jatm();
        
        bool decide(type0&);

        bool first_atm_lcl;
        int iatm;
        md_type& itype;
        type0* ix;
        type0* jx;
        int jatm;
        md_type jtype;
        type0 rsq;
        int xchng_mode;
        int dof_diff;
    };
    
    
    
    
    
    
    
    
    
    
    
    

}
#endif

