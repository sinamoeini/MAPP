/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "neighbor_md.h"
#include "ff.h"
#include "memory.h"
#include "timer.h"
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
    atoms->s2x(atoms->natms+atoms->natms_ph);
    
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
    
    type0* x=mapp->x->begin();
    md_type* type=mapp->type->begin();
    
    int x_dim=mapp->x->dim;
    int dim=atoms->dimension;
    type0 rsq;
    int icomp,jcomp;
    int iatm;
    int& jatm=cell->jatm;
    
    type0** cut_sk_sq=forcefield->cut_sk_sq;
    
    int* tmp_neigh_list=NULL;
    CREATE1D(tmp_neigh_list,neighbor_list_size_size+atoms->natms_ph);
    
    no_pairs=0;
    if(pair_wise)
    {
        for(iatm=0;iatm<atoms->natms;iatm++,cell->nxt_i())
        {
            icomp=x_dim*iatm;
            for(;jatm!=-1;cell->nxt_j())
            {
                if(jatm>iatm)
                {
                    jcomp=x_dim*jatm;
                    
                    rsq=0.0;
                    for(int idim=0;idim<dim;idim++)
                        rsq+=(x[icomp+idim]-x[jcomp+idim])
                        *(x[icomp+idim]-x[jcomp+idim]);
                    
                    if(rsq<cut_sk_sq[type[iatm]][type[jatm]])
                    {
                        tmp_neigh_list[neighbor_list_size[iatm]]=jatm;
                        neighbor_list_size[iatm]++;
                    }
                }
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
        for(iatm=0;iatm<atoms->natms;iatm++,cell->nxt_i())
        {
            icomp=x_dim*iatm;
            for(;jatm!=-1;cell->nxt_j())
            {
                if(jatm!=iatm)
                {
                    jcomp=x_dim*jatm;
                    
                    rsq=0.0;
                    for(int idim=0;idim<dim;idim++)
                        rsq+=(x[icomp+idim]-x[jcomp+idim])
                        *(x[icomp+idim]-x[jcomp+idim]);
                    
                    if(rsq<cut_sk_sq[type[iatm]][type[jatm]])
                    {
                        tmp_neigh_list[neighbor_list_size[iatm]]=jatm;
                        neighbor_list_size[iatm]++;
                    }
                }
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
