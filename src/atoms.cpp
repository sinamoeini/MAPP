#include "atoms.h"
#include "comm.h"
#include "xmath.h"
#include "error.h"
using namespace MAPP_NS;
/*---------------------------------------------------------------------------
      ___   _____   _____       ___  ___   _____
     /   | |_   _| /  _  \     /   |/   | /  ___/
    / /| |   | |   | | | |    / /|   /| | | |___
   / / | |   | |   | | | |   / / |__/ | | \___  \
  / /  | |   | |   | |_| |  / /       | |  ___| |
 /_/   |_|   |_|   \_____/ /_/        |_| /_____/
 
 ---------------------------------------------------------------------------*/
Atoms::Atoms(Communication* comm):
nvecs(0),
world(comm->world),
my_p(comm->my_p),
tot_p(comm->tot_p),
s_hi(comm->s_hi),
s_lo(comm->s_lo)
{
    vecs=NULL;
    

    tot_natms=0;
    natms=0;
    natms_ph=0;
    
    x=new Vec<type0>(this,__dim__);
    id= new Vec<int>(this,1);

    
    
    for(int idim=0;idim<__dim__;idim++)
    {
        max_cut_s[idim]=0.0;
        for(int jdim=0;jdim<__dim__;jdim++)
            H[idim][jdim]=B[idim][jdim]=0.0;
    }
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Atoms::~Atoms()
{


    while(nvecs)
        delete vecs[0];
}
/*--------------------------------------------
 add a new vec
 --------------------------------------------*/
void Atoms::add_vec(vec* v)
{
    if(v->name)
    {
        for(int ivec=0;ivec<nvecs;ivec++)
            if(vecs[ivec]->name && strcmp(v->name,vecs[ivec]->name)==0)
                Error::abort("vector %s already exist",v->name);
    }
    
    vec** vecs_=new vec*[nvecs+1];
    memcpy(vecs_,vecs,nvecs*sizeof(vec*));
    delete [] vecs;
    vecs=vecs_;
    vecs[nvecs]=v;
    nvecs++;
}
/*--------------------------------------------
 add a new vec with name
 --------------------------------------------*/
vec* Atoms::find_vec(const char* name)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        if(vecs[ivec]->name && strcmp(name,vecs[ivec]->name)==0)
            return vecs[ivec];
    Error::abort("vector %s does not exist",name);
    return NULL;
}
/*--------------------------------------------
 remove a vector
 --------------------------------------------*/
void Atoms::del_vec(vec* v)
{
    if(nvecs-1)
    {
        vec** vecs_=new vec*[nvecs-1];
        int ivec=0;
        for(int jvec=0;jvec<nvecs;jvec++)
            if(vecs[jvec]!=v)
                vecs_[ivec++]=vecs[jvec];
        
        delete [] vecs;
        vecs=vecs_;
        nvecs--;
    }
    else
    {
        delete [] vecs;
        vecs=NULL;
        nvecs=0;
    }
}

/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s(int no)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
#ifdef UNROLL
    for(int i=0;i<no;i++)
        XMatrixVector::x2s<__dim__>(x_vec+i*x_dim,B);
#else
    for(int i=0;i<no;i++,x_vec+=x_dim)
        for(int idim=0;idim<__dim__;idim++)
        {
            x_vec[idim]=x_vec[idim]*B[idim][idim];
            for(int jdim=idim+1;jdim<__dim__;jdim++)
                x_vec[idim]+=x_vec[jdim]*B[jdim][idim];
            while(x_vec[idim]<0.0)
                x_vec[idim]++;
            while(x_vec[idim]>=1.0)
                x_vec[idim]--;
        }
#endif
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x(int no)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
#ifdef UNROLL
    for(int i=0;i<no;i++)
        XMatrixVector::s2x<__dim__>(x_vec+i*x_dim,H);
#else
    for(int i=0;i<no;i++,x_vec+=x_dim)
        for(int idim=0;idim<__dim__;idim++)
        {
            x_vec[idim]=x_vec[idim]*H[idim][idim];
            for(int jdim=idim+1;jdim<__dim__;jdim++)
                x_vec[idim]+=x_vec[jdim]*H[jdim][idim];
        }
#endif
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s_lcl()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
#ifdef UNROLL
    for(int i=0;i<natms;i++)
        XMatrixVector::x2s<__dim__>(x_vec+i*x_dim,B);
#else
    for(int i=0;i<natms;i++,x_vec+=x_dim)
        for(int idim=0;idim<__dim__;idim++)
        {
            x_vec[idim]=x_vec[idim]*B[idim][idim];
            for(int jdim=idim+1;jdim<__dim__;jdim++)
                x_vec[idim]+=x_vec[jdim]*B[jdim][idim];
            
            while(x_vec[idim]<0.0)
                x_vec[idim]++;
            while(x_vec[idim]>=1.0)
                x_vec[idim]--;
        }
#endif
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x_lcl()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
#ifdef UNROLL
    for(int i=0;i<natms;i++)
        XMatrixVector::s2x<__dim__>(x_vec+i*x_dim,H);
#else
    for(int i=0;i<natms;i++,x_vec+=x_dim)
        for(int idim=0;idim<__dim__;idim++)
        {
            x_vec[idim]=x_vec[idim]*H[idim][idim];
            for(int jdim=idim+1;jdim<__dim__;jdim++)
                x_vec[idim]+=x_vec[jdim]*H[jdim][idim];
        }
#endif
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s_all()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    int nall=natms+natms_ph;
#ifdef UNROLL
    for(int i=0;i<nall;i++)
        XMatrixVector::x2s<__dim__>(x_vec+i*x_dim,B);
#else
    for(int i=0;i<nall;i++,x_vec+=x_dim)
        for(int idim=0;idim<__dim__;idim++)
        {
            x_vec[idim]=x_vec[idim]*B[idim][idim];
            for(int jdim=idim+1;jdim<__dim__;jdim++)
                x_vec[idim]+=x_vec[jdim]*B[jdim][idim];
            
            while(x_vec[idim]<0.0)
                x_vec[idim]++;
            while(x_vec[idim]>=1.0)
                x_vec[idim]--;
        }
#endif
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x_all()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    int nall=natms+natms_ph;
#ifdef UNROLL
    for(int i=0;i<nall;i++)
        XMatrixVector::s2x<__dim__>(x_vec+i*x_dim,H);
#else
    for(int i=0;i<nall;i++,x_vec+=x_dim)
        for(int idim=0;idim<__dim__;idim++)
        {
            x_vec[idim]=x_vec[idim]*H[idim][idim];
            for(int jdim=idim+1;jdim<__dim__;jdim++)
                x_vec[idim]+=x_vec[jdim]*H[jdim][idim];
        }
#endif
}
/*--------------------------------------------
 insert a number of atoms
 --------------------------------------------*/
void Atoms::insert(byte* buff,vec** vecs_,int nvecs_,int natms_)
{
    int stride=0;
    for(int ivec=0;ivec<nvecs_;ivec++) stride+=vecs_[ivec]->byte_sz;
    

    for(int ivec=0;ivec<nvecs_;ivec++)
        vecs_[ivec]->pst(buff,stride,natms_);
    
    natms+=natms_;
}
/*--------------------------------------------
 make room for some local atoms and phantom 
 atoms; this is used for grand canocical monte 
 carlo, when a successfull insertion trial has
 occured. using this function we make room for 
 the new entries. The new entries are insrerted 
 manually by GCMC
 
 *** we might need a better name for this 
 function
 --------------------------------------------*/
void Atoms::add()
{
    for(int i=0;i<nvecs;i++)
        vecs[i]->add();
    natms++;
}
/*--------------------------------------------
 delete some local atoms and phantom atoms; 
 this is used for grand canocical monte carlo, 
 when a successfull deletion trial has occured.
 it takes a the list of local atoms and phantoms
 
 !! it is assumed that both lists are ascending
 
 *** we might need a better name for this
 function
 --------------------------------------------*/
void Atoms::del(int& del_idx)
{
    for(int i=0;i<nvecs;i++)
        vecs[i]->del(del_idx);
    natms--;
}
/*--------------------------------------------
 restart
 --------------------------------------------*/
void Atoms::restart()
{
    tot_natms=0;
    natms=0;
    natms_ph=0;
    while(nvecs)
        delete vecs[0];
}


