/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <mpi.h>
#include "gcmc.h"
#include "memory.h"
#include "random.h"
#include "atom_types.h"
#include "neighbor.h"
#include "ff.h"
#include "md.h"
#include "MAPP.h"
#include "atoms.h"
#include "comm.h"
#include "dynamic.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
GCMC::GCMC(Dynamic*& dynamic_,atom_type gas_type_,type0 mu_,type0 T_,int seed):
world(comm->world),
dynamic(dynamic_),
gas_type(gas_type_),
mu(mu_),
T(T_),
cut_sq(forcefield->cut_sq),
s_hi(comm->s_hi),
s_lo(comm->s_lo),
natms(atoms->natms),
natms_ph(atoms->natms_ph),
itype(gas_type)
{
    ff=dynamic_cast<ForceFieldMD*>(forcefield);
    if(!ff)
        Error::abort("gcmc requires an md type forcefield");
    if(atoms->x->dim!=__dim__)
        Error::abort("x dimension must be %d for gcmc",__dim__);
    
    random=new Random(seed);
    s_trials=new type0*[__dim__];
    *s_trials=NULL;
    del_ids=NULL;
    del_ids_sz=del_ids_cpcty=0;
    vars=lcl_vars=NULL;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
GCMC::~GCMC()
{
    delete [] del_ids;
    delete [] s_trials;
    delete random;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::add_del_id(int* new_ids,int no)
{
    if(del_ids_sz+no>del_ids_cpcty)
    {
        int* del_ids_=new int[del_ids_sz+no];
        memcpy(del_ids_,del_ids,del_ids_sz*sizeof(int));
        del_ids_cpcty=del_ids_sz+no;
        delete [] del_ids;
        del_ids=del_ids_;
    }
    memcpy(del_ids+del_ids_sz,new_ids,sizeof(int)*no);
    del_ids_sz+=no;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int GCMC::get_new_id()
{
    if(del_ids_sz)
    {
        del_ids_sz--;
        return del_ids[del_ids_sz];
    }
    else
    {
        max_id++;
        return max_id;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::init()
{
    cut=forcefield->cut[itype][0];
    for(int i=1;i<mapp->atom_types->no_types;i++)
        cut=MAX(cut,forcefield->cut[itype][i]);
    
    gas_mass=mapp->atom_types->mass[gas_type];
    kbT=md->boltz*T;
    beta=1.0/kbT;
    lambda=md->hplanck/sqrt(2.0*M_PI*kbT*gas_mass);
    sigma=sqrt(kbT/gas_mass);
    z_fac=1.0;
    for(int i=0;i<__dim__;i++) z_fac/=lambda;
    z_fac*=exp(beta*mu);
    vol=1.0;
    for(int i=0;i<__dim__;i++)vol*=atoms->H[i][i];
    
    int max_id_=0;
    int* id=atoms->id->begin();
    for(int i=0;i<natms;i++)
        max_id_=MAX(id[i],max_id_);
    MPI_Allreduce(&max_id_,&max_id,1,MPI_INT,MPI_MAX,world);
    for(int i=0;i<del_ids_sz;i++)
        max_id=MAX(max_id,del_ids[i]);
        
    ngas=0;
    atom_type* type=mapp->type->begin();
    for(int i=0;i<natms;i++)
        if(type[i]==gas_type) ngas++;
    MPI_Allreduce(&ngas,&tot_ngas,1,MPI_INT,MPI_SUM,world);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::fin()
{
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::box_setup()
{
    int sz=0;
    max_ntrial_atms=1;
    for(int i=0;i<__dim__;i++)
    {
        type0 tmp=0.0;
        for(int j=i;j<__dim__;j++)
            tmp+=atoms->B[j][i]*atoms->B[j][i];
        cut_s[i]=sqrt(tmp)*cut;
        

        
        s_lo_ph[i]=s_lo[i]-cut_s[i];
        s_hi_ph[i]=s_hi[i]+cut_s[i];
        nimages_per_dim[i][0]=static_cast<int>(floor(s_hi_ph[i]));
        nimages_per_dim[i][1]=-static_cast<int>(floor(s_lo_ph[i]));
        max_ntrial_atms*=1+nimages_per_dim[i][0]+nimages_per_dim[i][1];
        sz+=1+nimages_per_dim[i][0]+nimages_per_dim[i][1];

    }
    
    *s_trials=new type0[sz];
    for(int i=1;i<__dim__;i++)
        s_trials[i]=s_trials[i-1]+1+nimages_per_dim[i-1][0]+nimages_per_dim[i-1][1];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::box_dismantle()
{
    delete [] *s_trials;
    *s_trials=NULL;
}



