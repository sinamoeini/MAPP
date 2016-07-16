/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "neighbor_dmd.h"
#include "ff.h"
#include "memory.h"
#include "timer.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Neighbor_dmd::Neighbor_dmd(MAPP* mapp):Neighbor(mapp)
{
    neighbor_list_size_size_2nd=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Neighbor_dmd::~Neighbor_dmd()
{
}
/*--------------------------------------------
 initiation before DMD
 --------------------------------------------*/
void Neighbor_dmd::init()
{
    Neighbor::init();
}
/*--------------------------------------------
 finalize DMD
 --------------------------------------------*/
void Neighbor_dmd::fin()
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

    if(neighbor_list_size_size_2nd)
    {
        for(int i=0;i<neighbor_list_size_size_2nd;i++)
            if(neighbor_list_size_2nd[i])
                delete [] neighbor_list_2nd[i];
        
        delete [] neighbor_list_2nd;
        delete [] neighbor_list_size_2nd;
    }
    neighbor_list_size_size_2nd=0;
    
    Neighbor::fin();
}
/*--------------------------------------------
 create the neighbopr list
 s_or_x: 0 for x, 1 for s
 --------------------------------------------*/
void Neighbor_dmd::create_list(bool box_change)
{
    timer->start(NEIGH_TIME_mode);
        
    cell->create(box_change);
    atoms->s2x(atoms->natms+atoms->natms_ph);
    
    if(neighbor_list_size_size)
    {
        for(int i=0;i<neighbor_list_size_size;i++)
            if(neighbor_list_size[i])
                delete [] neighbor_list[i];

        delete [] neighbor_list_size;
        delete [] neighbor_list;
    }
    
    neighbor_list_size_size=atoms->natms;
   
    neighbor_list=CREATE1D(neighbor_list,neighbor_list_size_size);
    CREATE1D(neighbor_list_size,neighbor_list_size_size);
    for(int i=0;i<neighbor_list_size_size;i++)
        neighbor_list_size[i]=0;
    
    type0* x=mapp->x->begin();
    type0* c=mapp->c->begin();
    dmd_type* ctype=mapp->ctype->begin();
    
    int x_dim=mapp->x->dim;
    int c_dim=mapp->c->dim;
    
    int dim=atoms->dimension;
    type0 rsq;
    int icomp,jcomp;
    int iatm;
    int& jatm=cell->jatm;
    int ic,jc;
    int rsq_chk;
    type0** cut_sk_sq=forcefield->cut_sk_sq;
    
    int* tmp_neigh_list;
    int tmp_neigh_list_size=1024;
    int tmp_neigh_list_grow=50;
    CREATE1D(tmp_neigh_list,tmp_neigh_list_size);
    
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
                    
                    rsq_chk=0;
                    ic=0;
                    while(rsq_chk==0 && ic<c_dim)
                    {
                        jc=0;
                        while(rsq_chk==0 && jc<c_dim)
                        {
                            if(c[iatm*c_dim+ic]>=0.0 &&
                               c[jatm*c_dim+jc]>=0.0 &&
                               rsq<cut_sk_sq[ctype[iatm*c_dim+ic]][ctype[jatm*c_dim+jc]])
                                rsq_chk=1;
                            jc++;
                        }
                        ic++;
                    }
                    
                    if(rsq_chk)
                    {
                        if(neighbor_list_size[iatm]+1>tmp_neigh_list_size)
                        {
                            GROW(tmp_neigh_list, tmp_neigh_list_size,tmp_neigh_list_size+tmp_neigh_list_grow);
                            tmp_neigh_list_size+=tmp_neigh_list_grow;
                        }
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
        for(iatm=0;iatm<atoms->natms;iatm++)
        {
            icomp=x_dim*iatm;
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
                    
                    rsq_chk=0;
                    ic=0;
                    while(rsq_chk==0 && ic<c_dim)
                    {
                        jc=0;
                        while(rsq_chk==0 && jc<c_dim)
                        {
                            if(c[iatm*c_dim+ic]>=0.0 &&
                               c[jatm*c_dim+jc]>=0.0 &&
                               rsq<cut_sk_sq[ctype[iatm*c_dim+ic]][ctype[jatm*c_dim+jc]])
                                rsq_chk=1;
                            jc++;
                        }
                        ic++;
                    }
                    
                    if(rsq_chk)
                    {
                        if(neighbor_list_size[iatm]+1>tmp_neigh_list_size)
                        {
                            GROW(tmp_neigh_list, tmp_neigh_list_size,tmp_neigh_list_size+tmp_neigh_list_grow);
                            tmp_neigh_list_size+=tmp_neigh_list_grow;
                        }
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
/*--------------------------------------------
 create the neighbopr list
 s_or_x: 0 for x, 1 for s
 --------------------------------------------*/
void Neighbor_dmd::create_2nd_list()
{
    timer->start(NEIGH_TIME_mode);

    type0* rsq_crd=forcefield->rsq_crd;
    
    type0 dx0,dx1,dx2,rsq;
    int iatm,jatm,icomp,jcomp,iicomp,jjcomp;
    int natms=atoms->natms;

    
    type0* x=mapp->x->begin();
    type0* c=mapp->c->begin();
    int c_dim=mapp->c->dim;
    dmd_type* type=mapp->ctype->begin();

    
    
    int** tmp_neigh_list;
    int* tmp_neigh_list_size;
    CREATE1D(tmp_neigh_list_size,c_dim);
    CREATE1D(tmp_neigh_list,c_dim);
    for(int itype=0;itype<c_dim;itype++)
    {
        CREATE1D(tmp_neigh_list[itype],1024);
        tmp_neigh_list_size[itype]=1024;
    }

    
    neighbor_list_size_size_2nd=natms*c_dim;
    CREATE1D(neighbor_list_2nd,neighbor_list_size_size_2nd);
    CREATE1D(neighbor_list_size_2nd,neighbor_list_size_size_2nd);
    
    for(int i=0;i<neighbor_list_size_size_2nd;i++)
        neighbor_list_size_2nd[i]=0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+c_dim)*iatm;
        iicomp=c_dim*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+c_dim)*jatm;
            jjcomp=c_dim*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            for(int ic_dim=0;ic_dim<c_dim;ic_dim++)
            {
                for(int jc_dim=0;jc_dim<c_dim;jc_dim++)
                {
                    if(type[iicomp+ic_dim]==type[jjcomp+jc_dim]
                       && c[iicomp+ic_dim]>=0.0
                       && c[jjcomp+jc_dim]>=0.0
                       && rsq<rsq_crd[type[iicomp+ic_dim]])
                    {
                        if(neighbor_list_size_2nd[iicomp+ic_dim]+1>tmp_neigh_list_size[ic_dim])
                        {
                            GROW(tmp_neigh_list[ic_dim],tmp_neigh_list_size[ic_dim],tmp_neigh_list_size[ic_dim]+1);
                            tmp_neigh_list_size[ic_dim]++;
                        }
                        tmp_neigh_list[ic_dim][neighbor_list_size_2nd[iicomp+ic_dim]]=jjcomp+jc_dim;
                        neighbor_list_size_2nd[iicomp+ic_dim]++;
                    }
                }
            }
        }
        
        for(int ic_dim=0;ic_dim<c_dim;ic_dim++)
        {
            CREATE1D(neighbor_list_2nd[iicomp+ic_dim],neighbor_list_size_2nd[iicomp+ic_dim]);
            memcpy(neighbor_list_2nd[iicomp+ic_dim],tmp_neigh_list[ic_dim],neighbor_list_size_2nd[iicomp+ic_dim]*sizeof(int));
        }
    }
    
    for(int ic_dim=0;ic_dim<c_dim;ic_dim++)
        if(tmp_neigh_list_size[ic_dim])
            delete [] tmp_neigh_list[ic_dim];
    if(c_dim)
    {
        delete [] tmp_neigh_list_size;
        delete [] tmp_neigh_list;
    }
    
    no_neigh_lists++;
    timer->stop(NEIGH_TIME_mode);
}
