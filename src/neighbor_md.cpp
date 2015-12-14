/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "neighbor_md.h"
#include "ff.h"
#include "memory.h"
#include "timer.h"
//#define ROW
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
}
/*--------------------------------------------
 finalize MD
 --------------------------------------------*/
void Neighbor_md::fin()
{
    if(neighbor_list_size_size)
    {
#ifdef ROW
        delete [] *neighbor_list;
        delete [] neighbor_list;
#else
        for(int i=0;i<neighbor_list_size_size;i++)
            if(neighbor_list_size[i])
                delete [] neighbor_list[i];
        delete [] neighbor_list;
#endif
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
#ifdef ROW

    timer->start(NEIGH_TIME_mode);
    
    if(box_change) create_bin_list();
    bin_atoms();
    
    
    if(neighbor_list_size_size)
    {
        delete [] *neighbor_list;
        delete [] neighbor_list;
        delete [] neighbor_list_size;
    }
    
    neighbor_list_size_size=atoms->natms;
    CREATE1D(neighbor_list_size,neighbor_list_size_size);
    CREATE1D(neighbor_list,neighbor_list_size_size);
    for(int i=0;i<neighbor_list_size_size;i++)
        neighbor_list_size[i]=0;
    
    int* tmp_neigh_list;
    int tmp_neigh_list_size=neighbor_list_size_size*(atoms->natms+atoms->natms_ph);
    CREATE1D(tmp_neigh_list,tmp_neigh_list_size);
    int* tmp_neigh_list_=tmp_neigh_list;
    
    
    type0* x=mapp->x->begin();
    md_type* type=mapp->type->begin();
    
    int x_dim=mapp->x->dim;
    int dim=atoms->dimension;
    type0 rsq;
    int icomp,jcomp;
    int ibin,jbin;
    int iatm,jatm;
    type0* cut_sk_sq=forcefield->cut_sk_sq;
    
    
    
    no_pairs=0;
    int no_pairs_=0;
    if(pair_wise)
    {
        for(iatm=0;iatm<atoms->natms;iatm++)
        {
            icomp=iatm*x_dim;
            ibin=atm_bin[iatm];
            for(int j=0;j<bin_neigh_list_size[ibin];j++)
            {
                jbin=bin_neigh_list[ibin][j];
                jatm=first_atom_bin[jbin];
                while(jatm!=-1)
                {
                    if(jatm>iatm)
                    {
                        jcomp=x_dim*jatm;
                        
                        rsq=0.0;
                        for(int idim=0;idim<dim;idim++)
                            rsq+=(x[icomp+idim]-x[jcomp+idim])*(x[icomp+idim]-x[jcomp+idim]);
                        
                        if(rsq<cut_sk_sq[COMP(type[iatm],type[jatm])])
                        {

                            *tmp_neigh_list_=jatm;
                            neighbor_list_size[iatm]++;
                            no_pairs++;
                            tmp_neigh_list_++;
                        }
                    }
                    jatm=next_atm[jatm];
                }
            }
        }
    }
    else
    {
        for(iatm=0;iatm<atoms->natms;iatm++)
        {
            icomp=iatm*x_dim;
            ibin=atm_bin[iatm];
            for(int j=0;j<bin_neigh_list_size[ibin];j++)
            {
                jbin=bin_neigh_list[ibin][j];
                jatm=first_atom_bin[jbin];
                while(jatm!=-1)
                {
                    if(jatm>iatm)
                    {
                        jcomp=x_dim*jatm;
                        
                        rsq=0.0;
                        for(int idim=0;idim<dim;idim++)
                            rsq+=(x[icomp+idim]-x[jcomp+idim])*(x[icomp+idim]-x[jcomp+idim]);
                        
                        if(rsq<cut_sk_sq[COMP(type[iatm],type[jatm])])
                        {
                            *tmp_neigh_list_=jatm;
                            neighbor_list_size[iatm]++;
                            no_pairs++;
                            tmp_neigh_list_++;
                            if(jatm>iatm)
                                no_pairs_++;
                        }
                    }
                    jatm=next_atm[jatm];
                }
            }
        }
    }
    
    CREATE1D(neighbor_list,neighbor_list_size_size);
    CREATE1D(*neighbor_list,no_pairs);
    memcpy(*neighbor_list,tmp_neigh_list,no_pairs*sizeof(int));
    if(tmp_neigh_list_size)
        delete [] tmp_neigh_list;
    for(int iatm=1;iatm<neighbor_list_size_size;iatm++)
        neighbor_list[iatm]=neighbor_list[iatm-1]+neighbor_list_size[iatm-1];
    
    if(!pair_wise) no_pairs=no_pairs_;
    
    if(atm_bin_size) delete [] atm_bin;
    
    atm_bin_size=0;
    no_neigh_lists++;
    
    timer->stop(NEIGH_TIME_mode);
#else
    timer->start(NEIGH_TIME_mode);
    
    if(box_change)
        create_bin_list();
    
    bin_atoms();
    
    
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
    int ibin,jbin;
    int iatm,jatm;
    type0* cut_sk_sq=forcefield->cut_sk_sq;
    
    int* tmp_neigh_list=NULL;
    CREATE1D(tmp_neigh_list,neighbor_list_size_size+atoms->natms_ph);
    
    no_pairs=0;
    if(pair_wise)
    {
        for(iatm=0;iatm<atoms->natms;iatm++)
        {
            icomp=x_dim*iatm;
            ibin=atm_bin[iatm];
            
            for(int j=0;j<bin_neigh_list_size[ibin];j++)
            {
                jbin=bin_neigh_list[ibin][j];
                jatm=first_atom_bin[jbin];
                while(jatm!=-1)
                {
                    if(jatm>iatm)
                    {
                        
                        jcomp=x_dim*jatm;
                        
                        rsq=0.0;
                        for(int idim=0;idim<dim;idim++)
                            rsq+=(x[icomp+idim]-x[jcomp+idim])
                            *(x[icomp+idim]-x[jcomp+idim]);
                       
                        if(rsq<cut_sk_sq[COMP(type[iatm],type[jatm])])
                        {
                            tmp_neigh_list[neighbor_list_size[iatm]]=jatm;
                            neighbor_list_size[iatm]++;
                        }
                    }
                    jatm=next_atm[jatm];
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
        for(iatm=0;iatm<atoms->natms;iatm++)
        {
            icomp=x_dim*iatm;
            ibin=atm_bin[iatm];
            
            for(int j=0;j<bin_neigh_list_size[ibin];j++)
            {
                jbin=bin_neigh_list[ibin][j];
                jatm=first_atom_bin[jbin];
                while(jatm!=-1)
                {
                    if(jatm!=iatm)
                    {
                        jcomp=x_dim*jatm;
                        
                        rsq=0.0;
                        for(int idim=0;idim<dim;idim++)
                            rsq+=(x[icomp+idim]-x[jcomp+idim])
                            *(x[icomp+idim]-x[jcomp+idim]);
                        
                        if(rsq<cut_sk_sq[COMP(type[iatm],type[jatm])])
                        {
                            tmp_neigh_list[neighbor_list_size[iatm]]=jatm;
                            neighbor_list_size[iatm]++;
                            
                            if(jatm>iatm) no_pairs++;
                        }
                    }
                    
                    jatm=next_atm[jatm];
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
    if(atm_bin_size)
        delete [] atm_bin;
    atm_bin_size=0;
    
    no_neigh_lists++;
    timer->stop(NEIGH_TIME_mode);
#endif
}
