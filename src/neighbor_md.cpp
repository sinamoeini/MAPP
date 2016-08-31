/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "neighbor_md.h"
#include "ff.h"
#include "memory.h"
#include "timer.h"
#ifdef DEBUG_NEIGH
#include "atom_types.h"
#endif
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Neighbor_md::Neighbor_md(MAPP* mapp):Neighbor(mapp)
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Neighbor_md::~Neighbor_md()
{
}
/*--------------------------------------------
 initiation before MD
 --------------------------------------------*/
void Neighbor_md::init()
{
    Neighbor::init();
}
/*--------------------------------------------
 finalize MD
 --------------------------------------------*/
void Neighbor_md::fin()
{
    if(neighbor_list_size_size)
    {
        for(int i=0;i<neighbor_list_size_size;i++)
            if(neighbor_list_size[i])
                delete [] neighbor_list[i];
        delete [] neighbor_list;
        delete [] neighbor_list_size;
    }
    neighbor_list_size_size=0;

    Neighbor::fin();
}
/*--------------------------------------------
 create the neighbopr list
 s_or_x: 0 for x, 1 for s
 --------------------------------------------*/
void Neighbor_md::create_list(bool box_change)
{
    timer->start(NEIGH_TIME_mode);
    
    cell->create(box_change);
    atoms->s2x_all();
    
    if(neighbor_list_size_size)
    {
        for(int i=0;i<neighbor_list_size_size;i++)
            if(neighbor_list_size[i])
                delete [] neighbor_list[i];
        delete [] neighbor_list;
        delete [] neighbor_list_size;
    }

    neighbor_list_size_size=atoms->natms;
    
    CREATE1D(neighbor_list,neighbor_list_size_size);
    CREATE1D(neighbor_list_size,neighbor_list_size_size);
    for(int i=0;i<neighbor_list_size_size;i++)
        neighbor_list_size[i]=0;

    
    md_type* type=mapp->type->begin();
    
    type0& rsq=cell->rsq;
    int& iatm=cell->iatm;
    int& jatm=cell->jatm;
    
    type0** cut_sk_sq=forcefield->cut_sk_sq;
    
    int* tmp_neigh_list=NULL;
    CREATE1D(tmp_neigh_list,neighbor_list_size_size+atoms->natms_ph);

    no_pairs=0;
    if(pair_wise)
    {
        for(;iatm!=-1;cell->nxt_i())
        {
            for(;jatm!=-1;cell->nxt_j())
            {
                if(jatm<=iatm) continue;
                if(rsq>=cut_sk_sq[type[iatm]][type[jatm]]) continue;
                tmp_neigh_list[neighbor_list_size[iatm]]=jatm;
                neighbor_list_size[iatm]++;

            }
            if(neighbor_list_size[iatm])
            {
                CREATE1D(neighbor_list[iatm],neighbor_list_size[iatm]);
                memcpy(neighbor_list[iatm],tmp_neigh_list,neighbor_list_size[iatm]*sizeof(int));
                no_pairs+=neighbor_list_size[iatm];
            }
            
        }
    }
    else
    {
        for(;iatm!=-1;cell->nxt_i())
        {
            for(;jatm!=-1;cell->nxt_j())
            {
                if(jatm==iatm) continue;
                if(rsq>=cut_sk_sq[type[iatm]][type[jatm]]) continue;
                tmp_neigh_list[neighbor_list_size[iatm]]=jatm;
                neighbor_list_size[iatm]++;
            }
            if(neighbor_list_size[iatm])
            {
                CREATE1D(neighbor_list[iatm],neighbor_list_size[iatm]);
                memcpy(neighbor_list[iatm],tmp_neigh_list,neighbor_list_size[iatm]*sizeof(int));
            }
        }
    }
    
    delete [] tmp_neigh_list;

    no_neigh_lists++;
    timer->stop(NEIGH_TIME_mode);
}
