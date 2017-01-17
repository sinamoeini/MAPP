/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "neighbor.h"
#include "ff.h"
#include "memory.h"
#include "timer.h"
#include "atoms.h"
#include "comm.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Neighbor::Neighbor(Atoms*& atoms_):
atoms(atoms_),
pair_wise(true)
{
    neighbor_list_size_size=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Neighbor::~Neighbor()
{
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
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::init()
{
    no_neigh_lists=0;
    cell=new Cell(1,atoms);
    create_list(true);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::fin()
{
    delete cell;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::mark_redndnt_ph(byte* mark)
{
    int natms=atoms->natms;
    int natms_ph=atoms->natms_ph;
    memset(mark,'0',natms_ph);
    for(int iatm=0;iatm<neighbor_list_size_size;iatm++)
        for(int j=0;j<neighbor_list_size[iatm];j++)
            if(neighbor_list[iatm][j]>=natms)
                mark[neighbor_list[iatm][j]-natms]='1';
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::rename_atoms(int* old_2_new)
{
    for(int iatm=0;iatm<neighbor_list_size_size;iatm++)
        for(int j=0;j<neighbor_list_size[iatm];j++)
            neighbor_list[iatm][j]=old_2_new[neighbor_list[iatm][j]];
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Neighbor::Cell::Cell(int m_,Atoms*& __atoms):
m(m_),
atoms(__atoms),
natms(__atoms->natms),
natms_ph(__atoms->natms_ph),
x(__atoms->x),
cut_s(__atoms->max_cut_s),
s_lo(__atoms->s_lo),
s_hi(__atoms->s_hi)
{
    ncells=0;
    head_atm=NULL;
    
    nneighs=0;
    rel_neigh_lst=NULL;
    
    next_vec=cell_vec=NULL;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Neighbor::Cell::~Cell()
{
    delete [] next_vec;
    delete [] cell_vec;
    
    delete [] rel_neigh_lst;
    
    delete [] head_atm;

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::Cell::box_setup()
{
    int ncells_=1;
    for(int i=0;i<__dim__;i++)
    {
        cell_size[i]=cut_s[i]/static_cast<type0>(m);
        ncells_per_dim[i]=static_cast<int>
        ((s_hi[i]-s_lo[i])/cell_size[i])+1+2*m;
        cell_denom[i]=ncells_;
        ncells_*=ncells_per_dim[i];
    }
    
    int countr[__dim__];
    for(int i=0;i<__dim__;i++)
        countr[i]=-m;
    int max_no_neighs=1;
    for(int i=0;i<__dim__;i++)
        max_no_neighs*=2*m+1;
    
    delete [] rel_neigh_lst;
    nneighs=0;
    rel_neigh_lst=new int[max_no_neighs];
    
    int sum,tmp;
    int rc_sq=m*m;
    for(int i=0;i<max_no_neighs;i++)
    {
        sum=0;
        for(int j=0;j<__dim__;j++)
        {
            sum+=countr[j]*countr[j];
            if(countr[j]!=0)
                sum+=1-2*std::abs(countr[j]);
        }
        
        if(sum<rc_sq)
        {
            tmp=0;
            for(int j=0;j<__dim__;j++)
                tmp+=cell_denom[j]*countr[j];
            
            rel_neigh_lst[nneighs++]=tmp;
        }
        
        countr[0]++;
        for(int j=0;j<__dim__-1;j++)
            if(countr[j]==m+1)
            {
                countr[j]=-m;
                countr[j+1]++;
            }
    }
    
    delete [] head_atm;
    head_atm=new int[ncells_];
    ncells=ncells_;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Neighbor::Cell::find_cell_no(type0*& s,int& cell_no)
{
    /*
    cell_no=0;
    auto icell_coord=
    [&s,&cell_no,this](int i)->void
    {
        if(s[i]<s_lo[i])
        {
            //0<=y y<m
            cell_no+=cell_denom[i]*(MIN(static_cast<int>((s[i]-(s_lo[i]-cut_s[i]))/cell_size[i]),m-1));
        }
        else if(s_hi[i]<=s[i])
        {
            //y<m+n+1
            cell_no+=cell_denom[i]*(MIN(MAX(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),ncells_per_dim[i]-2*m-1),ncells_per_dim[i]-1-m)+m);
        }
        else
        {
            //m<=y y<=m+n+1
            cell_no+=cell_denom[i]*(MIN(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),ncells_per_dim[i]-2*m-1)+m);
        }
    };
    XMatrixVector::UnrolledLoop[__dim__]::Do(icell_coord);
     */
    
    cell_no=0;
    for(int i=0;i<__dim__;i++)
    {
        if(s[i]<s_lo[i])
        {
            //0<=y y<m
            cell_no+=cell_denom[i]*(MIN(static_cast<int>((s[i]-(s_lo[i]-cut_s[i]))/cell_size[i]),m-1));
        }
        else if(s_hi[i]<=s[i])
        {
            //y<m+n+1
            cell_no+=cell_denom[i]*(MIN(MAX(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),ncells_per_dim[i]-2*m-1),ncells_per_dim[i]-1-m)+m);
        }
        else
        {
            //m<=y y<=m+n+1
            cell_no+=cell_denom[i]*(MIN(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),ncells_per_dim[i]-2*m-1)+m);
        }
    }
     
}
/*--------------------------------------------
 construct the cell list
 --------------------------------------------*/
void Neighbor::Cell::create(bool box_chng)
{
    if(box_chng)
        box_setup();
    
    int nall=natms+natms_ph;
    if(nall)
    {
        next_vec=new int[nall];
        cell_vec=new int[nall];
    }
    
    for(int i=0;i<ncells;i++)
        head_atm[i]=-1;
    int s_dim=x->dim;
    type0* s=x->begin()+(nall-1)*s_dim;
    

    
    for(int i=nall-1;i>-1;i--,s-=s_dim)
    {
        find_cell_no(s,cell_vec[i]);
        next_vec[i]=head_atm[cell_vec[i]];
        head_atm[cell_vec[i]]=i;
    }

    if(!natms)
    {
        iatm=-1;
        destroy();
        return;
    }
    
    iatm=-1;
    atoms->s2x_all();
    nxt_i();
}
/*--------------------------------------------
 destroy the cell list
 --------------------------------------------*/
void Neighbor::Cell::destroy()
{
    delete [] next_vec;
    delete [] cell_vec;
    next_vec=cell_vec=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::Cell::nxt_i()
{
    iatm++;
    if(iatm==natms)
    {
        iatm=-1;
        destroy();
        return;
    }
    
    ix=x->begin()+iatm*x->dim;
    
    icell=cell_vec[iatm];
    ineigh=-1;
    jatm=-1;
    while(jatm==-1 && ++ineigh<nneighs)
    {
        jcell=icell+rel_neigh_lst[ineigh];
        jatm=head_atm[jcell];
    }
    
    if(jatm==-1)
        return;
    if(jatm==iatm)
    {
        nxt_j();
        return;
    }
    jx=x->begin()+jatm*x->dim;
    rsq=XMatrixVector::rsq<__dim__>(ix,jx);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Neighbor::Cell::nxt_j()
{
    
    jatm=next_vec[jatm];
    while(jatm==-1 && ++ineigh<nneighs)
    {
        jcell=icell+rel_neigh_lst[ineigh];
        jatm=head_atm[jcell];
    }
    if(jatm==-1)
        return;
    if(jatm==iatm)
    {
        nxt_j();
        return;
    }
    jx=x->begin()+jatm*x->dim;
    rsq=XMatrixVector::rsq<__dim__>(ix,jx);
}
