#include "neighbor_cna.h"
#include "atoms.h"
#include "timer.h"
#include "MAPP.h"
#include "memory.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
NeighborCNA::NeighborCNA(Atoms*& atoms_,
Vec<atom_type>*& type_vec_,type0**& cut_sk_sq_):
type_vec(type_vec_),
cut_sk_sq(cut_sk_sq_),
Neighbor(atoms_)
{
    pair_wise=false;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
NeighborCNA::~NeighborCNA()
{
}
/*--------------------------------------------
 initiation before MD
 --------------------------------------------*/
void NeighborCNA::init()
{
    Neighbor::init();
}
/*--------------------------------------------
 finalize MD
 --------------------------------------------*/
void NeighborCNA::fin()
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
void NeighborCNA::create_list(bool box_change)
{
    mapp->timer->start(NEIGH_TIME_mode);
    cell->create(box_change);
    
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

    
    atom_type* type=type_vec->begin();
    
    type0& rsq=cell->rsq;
    int& iatm=cell->iatm;
    int& jatm=cell->jatm;
    
    
    int* tmp_neigh_list=NULL;
    CREATE1D(tmp_neigh_list,neighbor_list_size_size+atoms->natms_ph);
    type0* rsq_s=NULL;
    CREATE1D(rsq_s,neighbor_list_size_size+atoms->natms_ph);
    int* dummy=NULL;
    CREATE1D(dummy,neighbor_list_size_size+atoms->natms_ph);
    

    no_pairs=0;
    for(;iatm!=-1;cell->nxt_i())
    {
        for(;jatm!=-1;cell->nxt_j())
        {
            if(jatm==iatm) continue;
            if(rsq>=cut_sk_sq[type[iatm]][type[jatm]]) continue;
            tmp_neigh_list[neighbor_list_size[iatm]]=jatm;
            rsq_s[neighbor_list_size[iatm]]=rsq;
            neighbor_list_size[iatm]++;
        }
        
        for(int i=0;i<neighbor_list_size[iatm];i++) dummy[i]=i;
        
        XMath::quicksort(dummy,dummy+neighbor_list_size[iatm]
        ,[&rsq_s](int* i,int* j){return (rsq_s[*i]<rsq_s[*j]);}
        ,[](int* i,int* j){std::swap(*i,*j);});

        
        CREATE1D(neighbor_list[iatm],neighbor_list_size[iatm]);
        for(int i=0;i<neighbor_list_size[iatm];i++)
            neighbor_list[iatm][i]=tmp_neigh_list[dummy[i]];
        
        
        
    }
    
    delete [] dummy;
    delete [] tmp_neigh_list;
    delete [] rsq_s;
    no_neigh_lists++;

    mapp->timer->stop(NEIGH_TIME_mode);
}
























































