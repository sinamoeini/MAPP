/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "neighbor.h"
#include "ff.h"
#include "memory.h"
#include "timer.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Neighbor::Neighbor(MAPP* mapp):InitPtrs(mapp)
{
    pair_wise=true;
    dim=atoms->dimension;
    
    CREATE1D(tot_bin_grid,dim);
    CREATE1D(bin_size,dim);
    CREATE1D(bin_denom_list,dim);

    neighbor_list_size_size=0;
    tot_bin=0;
    first_atom_bin_size=0;
    next_atm_size=0;
    atm_bin_size=0;
    no_neigh_lists=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Neighbor::~Neighbor()
{
    if(dim)
    {
        delete [] tot_bin_grid;
        delete [] bin_size;
        delete [] bin_denom_list;
    }
}
/*--------------------------------------------
 create the bin neighbor list
 --------------------------------------------*/
void Neighbor::create_bin_list()
{
    type0* max_cut_s=atoms->max_cut_s;
    type0* s_lo=atoms->s_lo;
    type0* s_hi=atoms->s_hi;
    int dim=atoms->dimension;
    if(tot_bin)
    {
        for(int i=0;i<tot_bin;i++)
            if(bin_neigh_list_size[i])
                delete [] bin_neigh_list[i];
        
        delete [] bin_neigh_list;
        delete [] bin_neigh_list_size;
    }
    tot_bin=1;
    for(int i=0;i<dim;i++)
    {
        bin_size[i]=max_cut_s[i];
        tot_bin_grid[i]=static_cast<int>
        ((2.0*max_cut_s[i]+s_hi[i]-s_lo[i])/bin_size[i])+1;
        tot_bin*=tot_bin_grid[i];
    }
    
    CREATE1D(bin_neigh_list,tot_bin);
    CREATE1D(bin_neigh_list_size,tot_bin);
    for(int i=0;i<tot_bin;i++)
        bin_neigh_list_size[i]=0;
    for(int i=0;i<dim;i++)
    {
        int no=1;
        for(int j=0;j<i;j++)
            no*=tot_bin_grid[j];
        bin_denom_list[i]=no;
    }
    int* ibin_loc;
    CREATE1D(ibin_loc,dim);
    int* list_ii;
    CREATE1D(list_ii,dim);
    int tmp;
    for(int ibin=0;ibin<tot_bin;ibin++)
    {
        tmp=ibin;
        for(int i=dim-1;i>-1;i--)
        {
            ibin_loc[i]=tmp/bin_denom_list[i];
            tmp-=ibin_loc[i]*bin_denom_list[i];
        }

        for(int i=0;i<dim;i++)
            list_ii[i]=0;
        
        find_bin_no(dim,0,list_ii
        ,ibin,ibin_loc,bin_denom_list);
        
    }
    delete [] ibin_loc;
    delete [] list_ii;

 }
/*--------------------------------------------
 find the bin no based
 --------------------------------------------*/
void Neighbor::find_bin_no(int dim,int pos,int*& list_ii
,int ibin,int* ibin_loc,int* list)
{
    if(dim>0)
    {
        for(int i=0;i<atoms->dimension;i++)
        {
            list_ii[pos]=-1+i;
            find_bin_no(dim-1,pos+1,list_ii,ibin,ibin_loc,list);
        }
    }
    else
    {
        int no=0;
        int tmp;
        int chk=0;
        for(int i=0;i<pos;i++)
        {
            tmp=ibin_loc[i]+list_ii[i];
            if(tmp<0 || tmp>(tot_bin_grid[i]-1))
                chk=1;
            no+=tmp*list[i];
        }
        if(chk==0)
        {
            GROW(bin_neigh_list[ibin]
            ,bin_neigh_list_size[ibin]
            ,bin_neigh_list_size[ibin]+1);
            bin_neigh_list[ibin][bin_neigh_list_size[ibin]]=no;
            bin_neigh_list_size[ibin]++;
        }
    }    
}
/*--------------------------------------------
 bin the atoms
 --------------------------------------------*/
void Neighbor::bin_atoms()
{
    if(first_atom_bin_size)
        delete [] first_atom_bin;
    CREATE1D(first_atom_bin,tot_bin);
    first_atom_bin_size=tot_bin;
    for(int i=0;i<first_atom_bin_size;i++)
        first_atom_bin[i]=-1;
    
    int tot_natms=atoms->natms+atoms->natms_ph;
    
    if(next_atm_size)
        delete [] next_atm;
    CREATE1D(next_atm,tot_natms);
    next_atm_size=tot_natms;
    
    if(atm_bin_size)
        delete [] atm_bin;
    CREATE1D(atm_bin,atoms->natms);
    atm_bin_size=atoms->natms;
    
    int bin;
    type0* s=mapp->x->begin();
    int s_dim=mapp->x->dim;
    type0* max_cut_s=atoms->max_cut_s;
    type0* s_lo=atoms->s_lo;
    
    for(int i=tot_natms-1;i>-1;i--)
    {
        bin=0;
        for(int j=0;j<dim;j++)
        {
            bin+=static_cast<int>
            ((s[s_dim*i+j]+max_cut_s[j]-s_lo[j])/bin_size[j])
            *bin_denom_list[j];
        }
        
        if(i<atoms->natms)
            atm_bin[i]=bin;
        next_atm[i]=first_atom_bin[bin];
        first_atom_bin[bin]=i;
    }
    atoms->s2x(tot_natms);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::print_stats()
{
    if(atoms->my_p==0)
    {
        fprintf(output,"neighbor stats:\n");
        fprintf(output,"no of neigh lists "
        "generated: %d\n",no_neigh_lists);
    }
    no_neigh_lists=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::fin()
{    
    if(tot_bin)
    {
        for(int i=0;i<tot_bin;i++)
            if(bin_neigh_list_size[i])
                delete [] bin_neigh_list[i];
        
        delete [] bin_neigh_list;
        delete [] bin_neigh_list_size;
    }
    tot_bin=0;
    
    if(first_atom_bin_size)
        delete [] first_atom_bin;
    first_atom_bin_size=0;
    
    if(next_atm_size)
        delete [] next_atm;
    next_atm_size=0;
    
    if(atm_bin_size)
        delete [] atm_bin;
    atm_bin_size=0;
}

