#ifndef __MAPP__atoms__
#define __MAPP__atoms__
#define UNROLL
#define BLOCK_COMM
#include <typeinfo>
#include <cmath>
#include <mpi.h>
#include "macros.h"
#include "global.h"

enum {snd_to_bhnd,snd_to_frnt};
enum {rcv_fm_frnt,rcv_fm_bhnd};
/*-----------------------
 _     _   _____   _____  
| |   / / | ____| /  ___| 
| |  / /  | |__   | |     
| | / /   |  __|  | |     
| |/ /    | |___  | |___  
|___/     |_____| \_____| 
 -----------------------*/
namespace MAPP_NS
{
    class vec
    {
    private:
    protected:
    public:
        int dim;
        int orig_dim;
        int ncmprsd_dims;
        int byte_sz;
        unsigned int vec_sz;
        unsigned int dump_vec_sz;
        unsigned int vec_cpcty;
        
        
        char* name;
        char* print_format;
        vec(){}
        virtual ~vec(){};
        
        virtual void gather_dump(vec*)=0;
        virtual void gather_dump()=0;
        virtual void del_dump()=0;
        virtual void print_dump(FILE*,int)=0;
        virtual void* end()=0;
        
        virtual void reserve(int)=0;
        virtual void resize(int)=0;
        virtual void shrink_to_fit()=0;
        virtual void rearrange(int*,int*,int,int)=0;
        virtual void pst_to(byte*&,int)=0;
        
        virtual void del(int&)=0;
        virtual void add()=0;
        
        virtual void pop_out(byte*&,int)=0;
        virtual void pop_in(byte*&)=0;
        
        virtual void cpy(byte*&,int)=0;
        virtual void pst(byte*&,int,int)=0;
        virtual void cpy_pst(int)=0;
        
        virtual void cpy(byte*&,int*,int)=0;
        virtual void pst(byte*&,int)=0;
        virtual void cpy_pst(int*,int)=0;
        virtual void assign_print_format(const char*)=0;
    };
}
/*-----------------------
 _     _   _____   _____  
| |   / / | ____| /  ___| 
| |  / /  | |__   | |     
| | / /   |  __|  | |     
| |/ /    | |___  | |___  
|___/     |_____| \_____| 
 -----------------------*/
namespace MAPP_NS
{
    template<typename T>
    class Vec: public vec
    {
    private:
        T* vec;
        T* dump_vec;
        T def_val;
        class Atoms* atoms;
    protected:
    public:
        Vec(Atoms*,int);
        Vec(Atoms*,int,const char*);
        Vec(Vec&);
        Vec(const Vec&&);
        ~Vec();
        
        void gather_dump(class vec*);
        void gather_dump(Vec<atom_type>*);
        void gather_dump();
        void del_dump();
        void print_dump(FILE*,int);
        void assign_print_format();
        void assign_print_format(const char*);
        void change_dimension(int);
        void change_dimension(T,int,int);
        
        T*& begin();
        T*& begin_dump();
        void* end();
        void reserve(int);
        void resize(int);
        void shrink_to_fit();
        // rearrange
        void rearrange(int*,int*,int,int);
        void pst_to(byte*&,int);
        
        // these two come in pair for gcmc
        void del(int&);
        void add();
        
        // these two come in pair for exchange atoms one by one
        void pop_out(byte*&,int);
        void pop_in(byte*&);
        
        // these two come in pair for sending & recieving phantom atoms one by one
        void cpy(byte*&,int);
        void pst(byte*&,int,int);
        void cpy_pst(int);
        
        // these two come in pair for sending & recieving phantom atoms collectively
        void cpy(byte*&,int*,int);
        void pst(byte*&,int);
        void cpy_pst(int*,int);
        
        void operator() (Vec<T>*&);
        void operator() (Vec<T>&);
        T& operator () (const int,const int);
        T& operator [] (const int);
        const T& operator () (const int iatm,const int idim) const;
        const T& operator [] (const int i) const;
        
        void compress(Vec<atom_type>*&,int,T);
        void decompress_dump(Vec<atom_type>*&);
        
    };
}
/*-------------------------------------------------
      ___   _____   _____       ___  ___   _____
     /   | |_   _| /  _  \     /   |/   | /  ___/
    / /| |   | |   | | | |    / /|   /| | | |___
   / / | |   | |   | | | |   / / |__/ | | \___  \
  / /  | |   | |   | |_| |  / /       | |  ___| |
 /_/   |_|   |_|   \_____/ /_/        |_| /_____/
 -------------------------------------------------*/
namespace MAPP_NS
{
    class Atoms
    {
    private:
    protected:
    public:
        //keeep these
        int natms;
        int natms_ph;
        int tot_natms;
        type0 H[__dim__][__dim__];
        type0 B[__dim__][__dim__];
        Vec<type0>* x;
        Vec<int>* id;
        vec** vecs;
        int nvecs;
        void add_vec(vec*);
        vec* find_vec(const char*);
        void del_vec(vec*);
        void x2s(int);
        void s2x(int);
        void x2s_lcl();
        void s2x_lcl();
        void x2s_all();
        void s2x_all();
        
        void insert(byte*,vec**,int,int);
        void add();
        void del(int&);
        Atoms(class Communication*);
        ~Atoms();
        
        //move these to comm
        const int tot_p;
        const int my_p;
        MPI_Comm& world;
        type0 (&s_lo)[__dim__];
        type0 (&s_hi)[__dim__];

        
        //dont know what to do with this
        type0 max_cut;
        type0 max_cut_s[__dim__];
        void restart();        
    };
}
using namespace MAPP_NS;
/*-----------------------
 _     _   _____   _____  
| |   / / | ____| /  ___| 
| |  / /  | |__   | |     
| | / /   |  __|  | |     
| |/ /    | |___  | |___  
|___/     |_____| \_____| 
 
 -----------------------*/
/*--------------------------------------------
 constructor
 --------------------------------------------*/
template<typename T>
Vec<T>::Vec(Atoms* atoms_,int dim_):
atoms(atoms_),
vec()
{
    def_val=0;
    ncmprsd_dims=0;
    dim=dim_;
    orig_dim=dim;
    byte_sz=sizeof(T)*dim;
    
    vec_sz=vec_cpcty=atoms->natms;
    dump_vec_sz=0;
    
    vec=NULL;
    if(vec_sz) vec=new T[vec_sz*dim];
    
    dump_vec=NULL;
    
    name=NULL;
    assign_print_format();
    
    atoms->add_vec(this);
}
/*--------------------------------------------
 constructor with name;
 --------------------------------------------*/
template<typename T>
Vec<T>::Vec(Atoms* atoms_,int dim_,const char* name_):
atoms(atoms_),
vec()
{
    def_val=0;
    ncmprsd_dims=0;
    dim=dim_;
    orig_dim=dim;
    byte_sz=sizeof(T)*dim;
    
    vec_sz=vec_cpcty=atoms->natms;
    dump_vec_sz=0;
    
    vec=NULL;
    if(vec_sz) vec=new T[vec_sz*dim];
        
    dump_vec=NULL;
    
    name=NULL;
    assign_print_format();
    
    int length=static_cast<int>(strlen(name_))+1;
    name=new char[length];
    memcpy(name,name_,length*sizeof(char));
    
    atoms->add_vec(this);
}
/*--------------------------------------------
 copy costructor
 --------------------------------------------*/
template<typename T>
Vec<T>::Vec(Vec<T>& other)
{
    atoms=other.atoms;
    
    def_val=other.def_val;
    ncmprsd_dims=other.ncmprsd_dims;
    dim=other.dim;
    orig_dim=other.orig_dim;
    byte_sz=other.byte_sz;
    
    
    vec_sz=other.vec_sz;
    vec_cpcty=other.vec_cpcty;
    dump_vec_sz=other.dump_vec_sz;
    
    vec=NULL;
    if(vec_cpcty)
    {
        vec=new T[vec_cpcty*dim];
        memcpy(vec,other.vec,vec_sz*byte_sz);
    }
    
    dump_vec=NULL;
    if(dump_vec_sz)
    {
        dump_vec=new T[dump_vec_sz*dim];
        memcpy(dump_vec,other.dump_vec,dump_vec_sz*byte_sz);
    }
    
    name=NULL;
    if(other.name)
    {
        size_t name_sz=strlen(other.name);
        size_t cpy_sz=strlen("_copy");
        name=new char[name_sz+cpy_sz+1];
        memcpy(name,other.name,name_sz*sizeof(char));
        memcpy(name+name_sz,"_copy",(cpy_sz+1)*sizeof(char));
    }
        
    size_t print_sz=strlen(other.print_format)+1;
    print_format=new char[print_sz];
    memcpy(print_format,other.print_format,print_sz*sizeof(char));
    
    atoms->add_vec(this);
}
/*--------------------------------------------
 move costructor
 --------------------------------------------*/
template<typename T>
Vec<T>::Vec(const Vec<T>&& other)
{
    atoms=other.atoms;
    
    def_val=other.def_val;
    ncmprsd_dims=other.ncmprsd_dims;
    dim=other.dim;
    orig_dim=other.orig_dim;
    byte_sz=other.byte_sz;
    
    
    vec_sz=other.vec_sz;
    other.vec_sz=0;
    vec_cpcty=other.vec_cpcty;
    other.vec_cpcty=0;
    dump_vec_sz=other.dump_vec_sz;
    other.dump_vec_sz=0;
    vec=other.vec;
    other.vec=NULL;
    dump_vec=other.dump_vec;
    other.dump_vec=NULL;
    name=other.name;
    other.name=NULL;
    print_format=other.print_format;
    other.print_format=NULL;
    
    atoms->add_vec(this);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<typename T>
inline Vec<T>::~Vec()
{
    delete [] vec;
    delete [] dump_vec;
    delete [] name;
    delete [] print_format;
    atoms->del_vec(this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void Vec<T>::change_dimension(int dim_)
{
    if(dim_==dim)
        return;
    
    int min_dim=MIN(dim_,dim);

    
    T* new_vec=NULL;
    if(vec_cpcty) new_vec=new T[vec_cpcty*dim_];
    
    for(int i=0;i<vec_sz;i++)
        for(int j=0;j<min_dim;j++)
            new_vec[i*dim_+j]=vec[i*dim+j];
    
    delete [] vec;
    
    vec=new_vec;
    int del_dim=dim_-dim;
    dim+=del_dim;
    orig_dim+=del_dim;
    byte_sz=sizeof(T)*dim;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void Vec<T>::change_dimension(T def_val_,int ncmprsd_dims_,int del_dim)
{
    int dim_=dim-del_dim;
    ncmprsd_dims=ncmprsd_dims_;
    def_val=def_val_;
    
    if(dim_==dim)
        return;
    
    
    int min_dim=MIN(dim_,dim);

    
    T* new_vec=NULL;
    if(vec_cpcty) new_vec=new T[vec_cpcty*dim_];
    
    for(int i=0;i<vec_sz;i++)
        for(int j=0;j<min_dim;j++)
            new_vec[i*dim_+j]=vec[i*dim+j];
    
    
    delete [] vec;
    
    vec=new_vec;
    dim=dim_;
    byte_sz=sizeof(T)*dim;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::resize(int natms)
{
    if(vec_cpcty>=natms)
    {
        vec_sz=natms;
        return;
    }
    
    
    T* new_vec=NULL;
    if(natms) new_vec=new T[natms*dim];
    
    memcpy(new_vec,vec,byte_sz*vec_sz);
    
    delete [] vec;
    vec_sz=natms;
    vec_cpcty=natms;
    vec=new_vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::shrink_to_fit()
{
    if(vec_cpcty==vec_sz)
        return;
    
    T* new_vec=NULL;
    if(vec_sz) new_vec=new T[vec_sz*dim];
    
    memcpy(new_vec,vec,byte_sz*vec_sz);
    
    delete [] vec;
    vec_cpcty=vec_sz;
    vec=new_vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::reserve(int xtra_atms)
{
    if(vec_cpcty>=vec_sz+xtra_atms)
        return;
    
    T* new_vec=new T[(vec_sz+xtra_atms)*dim];
    
    memcpy(new_vec,vec,byte_sz*vec_sz);
    
    delete [] vec;
    vec_cpcty=vec_sz+xtra_atms;
    vec=new_vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::rearrange(int* old_pos
,int* new_pos,int sz,int new_vec_cpcty)
{
    T* new_vec=NULL;
    if(new_vec_cpcty) new_vec=new T[new_vec_cpcty*dim];
    
    for(int i=0;i<sz;i++)
        memcpy(new_vec+new_pos[i]*dim,vec+old_pos[i]*dim,byte_sz);
    
    delete [] vec;
    vec=new_vec;
    vec_cpcty=new_vec_cpcty;
    vec_sz=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void Vec<T>::decompress_dump(Vec<atom_type>*& map)
{
    T* dump_vec_=NULL;
    if(dump_vec_sz) dump_vec_=new T[orig_dim*dump_vec_sz];

    T* v0=dump_vec;
    T* v1=dump_vec_;
    for(int i=0;i<dump_vec_sz;i++,v0+=dim,v1+=orig_dim)
        memcpy(v1,v0,byte_sz);
    delete [] dump_vec;
    dump_vec=dump_vec_;
    int map_dim=map->dim;
    
    int xtndd_dim=orig_dim-dim+map_dim;
    v0=dump_vec+(dim-map_dim);
    atom_type* map_vec=map->begin();
    for(int i=0;i<dump_vec_sz;i++,v0+=dim,v0+=orig_dim,map_vec+=map_dim)
    {
        for(int j=0;j<xtndd_dim;j++)
            v0[j]=def_val;
        for(int j=0;j<map_dim;j++)
            v0[map_vec[j]]=v0[j];
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void Vec<T>::compress(Vec<atom_type>*& map,int ncmprsd_dims_,T def_val_)
{
    def_val=def_val_;
    
    if(map)
    {
        int map_dim=map->dim;
        if(map_dim==ncmprsd_dims_)
            return;
        
        T* vec_=vec+dim-ncmprsd_dims_;
        atom_type* map_vec=map->begin();
        
        for(int i=0;i<vec_sz;i++,map_vec+=map_dim,vec_+=dim)
            for(int j=0;j<map_dim;j++)
                vec_[j]=vec_[map[j]];
        int dim_reduc=ncmprsd_dims_-map_dim;
        change_dimension(dim-dim_reduc);
        orig_dim+=dim_reduc;
        return;
    }
    
    int map_dim=ncmprsd_dims_;
    map=new Vec<atom_type>(atoms,map_dim);
    
    int max_dim_lcl=0,max_dim_tmp,dim_reduc;
    T* vec_=vec+dim-map_dim;
    for(int i=0;i<vec_sz;i++,vec_+=dim)
    {
        max_dim_tmp=0;
        for(int j=0;j<map_dim;j++)
            if(vec_[j]!=def_val)
                max_dim_tmp++;
        max_dim_lcl=MAX(max_dim_tmp,max_dim_lcl);
    }
    MPI_Allreduce(&max_dim_lcl,&dim_reduc,1,MPI_INT,MPI_MAX,atoms->world);
    dim_reduc=map_dim-dim_reduc;
    
    if(dim_reduc==0)
    {
        atom_type* map_vec=map->begin();
        for(int i=0;i<vec_sz;i++,map_vec+=map_dim)
            for(atom_type j=0;j<map_dim;j++)
                map_vec[j]=j;
        return;
    }
    
    vec_=vec+dim-map_dim;
    atom_type* map_vec=map->begin();
    for(int i=0,icurs;i<vec_sz;i++,map_vec+=map_dim,vec_+=dim)
    {
        icurs=0;
        for(int j=0;j<map_dim;j++)
            map_vec[j]=j;
        for(int j=0;j<map_dim;j++)
            if(vec_[j]!=def_val)
            {
                std::swap(vec_[j],vec_[icurs]);
                std::swap(map_vec[j],map_vec[icurs]);
                icurs++;
            }
    }
    
    change_dimension(dim-dim_reduc);
    orig_dim+=dim_reduc;
    
    map->change_dimension(map_dim-dim_reduc);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::pst_to(byte*& buff,int iatm)
{
    memcpy(&vec[iatm*dim],buff,byte_sz);
    buff+=byte_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline T& Vec<T>::operator[](const int i)
{
    return vec[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline T& Vec<T>::operator()(const int iatm,const int idim)
{
    return vec[iatm*dim+idim];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
const T& Vec<T>::operator () (const int iatm,const int idim) const
{
    return vec[iatm*dim+idim];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
const T& Vec<T>::operator [] (const int i) const
{
    return vec[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void* Vec<T>::end()
{
    return &vec[vec_sz*dim];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline T*& Vec<T>::begin()
{
    return vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline T*& Vec<T>::begin_dump()
{
    return dump_vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::operator()(Vec<T>*& v)
{
    memcpy(vec,v->vec,byte_sz*vec_sz);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::operator()(Vec<T>& v)
{
    memcpy(vec,v.vec,byte_sz*vec_sz);
}
/*------------------------------------
 _____   _   _       ___  ___   _____  
|  _  \ | | | |     /   |/   | |  _  \ 
| | | | | | | |    / /|   /| | | |_| | 
| | | | | | | |   / / |__/ | | |  ___/ 
| |_| | | |_| |  / /       | | | |     
|_____/ \_____/ /_/        |_| |_|
 
 ------------------------------------*/
template<typename T>
void Vec<T>::assign_print_format(const char* form)
{
    size_t sz=strlen(form)+1;
    delete [] print_format;
    print_format=new char[sz];
    memcpy(print_format,form,sz*sizeof(char));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void Vec<T>::assign_print_format()
{
    print_format=NULL;
    if(std::is_same<T,char>::value)
        assign_print_format("%c ");
    else if(std::is_same<T,unsigned char>::value)
        assign_print_format("%c ");
    else if(std::is_same<T,short int>::value)
        assign_print_format("%d ");
    else if(std::is_same<T,unsigned short int>::value)
        assign_print_format("%d ");
    else if(std::is_same<T,int>::value)
        assign_print_format("%d ");
    else if(std::is_same<T,unsigned int>::value)
        assign_print_format("%d ");
    else if(std::is_same<T,long int>::value)
        assign_print_format("%ld ");
    else if(std::is_same<T,unsigned long int>::value)
        assign_print_format("%ld ");
    else if(std::is_same<T,long long int>::value)
        assign_print_format("%lld ");
    else if(std::is_same<T,unsigned long long int>::value)
        assign_print_format("%lld ");
    else if(std::is_same<T,float>::value)
        assign_print_format("%f ");
    else if(std::is_same<T,double>::value)
        assign_print_format("%lf ");
    else if(std::is_same<T,long double>::value)
        assign_print_format("%Lf ");
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void Vec<T>::gather_dump(class vec* map)
{
    Vec<atom_type>* map_atom_type
    =dynamic_cast<Vec<atom_type>*>(map);
    if(map_atom_type)
    {
        gather_dump(map_atom_type);
        return;
    }
}
/*--------------------------------------------
 
                            ncmprsd_dims
                          ______/\______
                         |              |
                 orig_dim
  __________________/\__________________
 |                                      |
 x-----------------------x--------------x
 
 
 after transformation:
 
 
                    dim-orig_dim+ncmprsd_dims
                          ____/\____
                         |          |
                   dim
  __________________/\______________
 |                                  |
 x-----------------------x----------x
 |___________  __________|
             \/
 
     orig_dim-ncmprsd_dims
 
 --------------------------------------------*/
template<typename T>
void Vec<T>::gather_dump(Vec<atom_type>* map)
{
    int natms=atoms->natms;
    T* vec_=vec;
    int byte_sz_=byte_sz;
    
    if(ncmprsd_dims)
    {
        atom_type* map_vec=map->begin();
        int map_dim=map->dim;
        int ncmp=(dim-orig_dim+ncmprsd_dims)/map_dim;
        int map_dim_=ncmprsd_dims/ncmp;
        int strt_dim=orig_dim-ncmprsd_dims;
        
        vec_=new T[natms*orig_dim];
        byte_sz_=orig_dim*sizeof(T);
        
        int ipos=0;
        int ipos_=0;
        int map_ipos=0;
        for(int iatm=0;iatm<natms;iatm++)
        {
            for(int idim=0;idim<strt_dim;idim++)
                vec_[ipos_++]=vec[ipos++];
            
            for(int idim=0;idim<ncmprsd_dims;idim++)
                vec_[ipos_+idim]=def_val;
            
            for(int icmp=0;icmp<ncmp;icmp++)
            {
                for(int idim=0;idim<map_dim;idim++)
                    vec_[ipos_+map_vec[map_ipos+idim]]=vec[ipos+idim];
                
                ipos+=map_dim;
                ipos_+=map_dim_;
            }
            
            map_ipos+=map_dim;
        }
    }
    
    int my_p=atoms->my_p;
    int tot_p=atoms->tot_p;
    int tot_natms=atoms->tot_natms;
    int* rcv_size=NULL;
    
    if(my_p!=0)
        MPI_Send(&natms,1,MPI_INT,0,my_p,atoms->world);
    
    if(my_p==0)
    {
        rcv_size=new int[tot_p];
        for(int ip=1;ip<tot_p;ip++)
            MPI_Recv(&rcv_size[ip],1
            ,MPI_INT,ip,ip,atoms->world,MPI_STATUS_IGNORE);

        dump_vec=new T[orig_dim*tot_natms];
        dump_vec_sz=tot_natms;
        memcpy(dump_vec,vec_,natms*byte_sz_);
    }
    
    if(my_p!=0)
        MPI_Send(vec_,natms*byte_sz_,MPI_BYTE,0,my_p,atoms->world);
    
    if(my_p==0)
    {
        int tot_atoms=natms;
        for(int ip=1;ip<tot_p;ip++)
        {
            MPI_Recv(&dump_vec[tot_atoms*orig_dim]
            ,(byte_sz_*rcv_size[ip]),MPI_BYTE
            ,ip,ip,atoms->world,MPI_STATUS_IGNORE);
            tot_atoms+=rcv_size[ip];
        }
        
        delete [] rcv_size;
    }
    
    if(ncmprsd_dims && natms)
        delete [] vec_;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void Vec<T>::gather_dump()
{
    int natms=atoms->natms;
    int my_p=atoms->my_p;
    int tot_p=atoms->tot_p;
    int tot_natms=atoms->tot_natms;
    int* rcv_size=NULL;
    
    if(my_p!=0)
        MPI_Send(&natms,1,MPI_INT,0,my_p,atoms->world);
    
    if(my_p==0)
    {
        rcv_size=new int[tot_p];
        for(int ip=1;ip<tot_p;ip++)
            MPI_Recv(&rcv_size[ip],1
            ,MPI_INT,ip,ip,atoms->world,MPI_STATUS_IGNORE);
        
        dump_vec=new T[orig_dim*tot_natms];
        dump_vec_sz=tot_natms;
        memcpy(dump_vec,vec,natms*byte_sz);
    }
    
    if(my_p!=0)
        MPI_Send(vec,natms*byte_sz,MPI_BYTE,0,my_p,atoms->world);
    
    if(my_p==0)
    {
        int tot_atoms=natms;
        for(int ip=1;ip<tot_p;ip++)
        {
            MPI_Recv(&dump_vec[tot_atoms*dim]
            ,(byte_sz*rcv_size[ip]),MPI_BYTE
            ,ip,ip,atoms->world,MPI_STATUS_IGNORE);
            tot_atoms+=rcv_size[ip];
        }
        
        delete [] rcv_size;
    }
}
/*--------------------------------------------
 delete dump
 --------------------------------------------*/
template<typename T>
void Vec<T>::del_dump()
{
    delete [] dump_vec;
    dump_vec_sz=0;
    dump_vec=NULL;
}
/*--------------------------------------------
 delete dump
 --------------------------------------------*/
template<typename T>
void Vec<T>::print_dump(FILE* fp,int iatm)
{
    for(int i=0;i<orig_dim;i++)
        fprintf(fp,print_format,dump_vec[iatm*orig_dim+i]);
}
/*-----------------------------------------------------------------
  _____   _____       ___  ___   _____
 /  ___| /  ___|     /   |/   | /  ___|
 | |     | |        / /|   /| | | |
 | |  _  | |       / / |__/ | | | |
 | |_| | | |___   / /       | | | |___
 \_____/ \_____| /_/        |_| \_____|
 -----------------------------------------------------------------*/
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
template<typename T>
inline void Vec<T>::add()
{
    resize(vec_sz+1);
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
template<typename T>
inline void Vec<T>::del(int &del_idx)
{
    vec_sz--;
    if(del_idx==vec_sz) return;
    memcpy(vec+del_idx*dim,vec+vec_sz*dim,byte_sz);
}
/*-----------------------------------------------------------------
 _____  __    __  _____   _   _       ___   __   _   _____   _____  
| ____| \ \  / / /  ___| | | | |     /   | |  \ | | /  ___| | ____| 
| |__    \ \/ /  | |     | |_| |    / /| | |   \| | | |     | |__   
|  __|    }  {   | |     |  _  |   / / | | | |\   | | |  _  |  __|  
| |___   / /\ \  | |___  | | | |  / /  | | | | \  | | |_| | | |___  
|_____| /_/  \_\ \_____| |_| |_| /_/   |_| |_|  \_| \_____/ |_____|
 
 -----------------------------------------------------------------*/
/*--------------------------------------------
 1. removes atom iatm and copy it, using buff
 as destenation
 2. replaces iatm with last atom
 3. moves buff to the end of copied data
 4. decreases the size
 ** used during the EXCHANGE of atoms
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::pop_out(byte*& buff,int iatm)
{
    memcpy(buff,&vec[iatm*dim],byte_sz);
    memcpy(&vec[iatm*dim],&vec[(vec_sz-1)*dim],byte_sz);
    buff+=byte_sz;
    vec_sz--;
}
/*--------------------------------------------
 1. adds one atom to the bottom of the vector
 using buff as source
 2. moves buff to the end of copied data
 3. increases the size
 !! does NOT check for capicity
  ** used during the EXCHANGE of atoms
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::pop_in(byte*& buff)
{
    memcpy(&vec[vec_sz*dim],buff,byte_sz);
    buff+=byte_sz;
    vec_sz++;
}
/*------------------------------------------------------------------------------------------------
 _____   _   _       ___   __   _   _____   _____       ___  ___        _       _   _____   _____
|  _  \ | | | |     /   | |  \ | | |_   _| /  _  \     /   |/   |      | |     | | /  ___/ |_   _| 
| |_| | | |_| |    / /| | |   \| |   | |   | | | |    / /|   /| |      | |     | | | |___    | |   
|  ___/ |  _  |   / / | | | |\   |   | |   | | | |   / / |__/ | |      | |     | | \___  \   | |   
| |     | | | |  / /  | | | | \  |   | |   | |_| |  / /       | |      | |___  | |  ___| |   | |   
|_|     |_| |_| /_/   |_| |_|  \_|   |_|   \_____/ /_/        |_|      |_____| |_| /_____/   |_|   
 
 ------------------------------------------------------------------------------------------------*/
/*--------------------------------------------
 1. copies atom iatm, using buff as 
 destenation
 2. moves buff to the end of copied data
 3. does NOT increase/decrease the size
 ** used during the constructing the PHANTOM
 atoms
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::cpy(byte*& buff,int iatm)
{
    memcpy(buff,&vec[iatm*dim],byte_sz);
    buff+=byte_sz;
}
/*--------------------------------------------
 1. adds xtra_natms number of atom to the 
 bottom of the vector using buff as source
 2. moves buff to the end of 1st atom copied
 3. increases the size by xtra_natms
 !! DOES check for capicity
 ** used during the constructing the PHANTOM
 atoms
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::pst(byte*& buff,int stride,int xtra_natms)
{
    reserve(xtra_natms);
    T* vec_buff=&vec[vec_sz*dim];
    for(int i=0;i<xtra_natms;i++)
    {
        memcpy(vec_buff,buff,byte_sz);
        buff+=stride;
        vec_buff+=dim;
    }
    buff+=byte_sz-stride*xtra_natms;
    vec_sz+=xtra_natms;
}
/*--------------------------------------------
 1. copies from atom iatm to the bottom of vector
 2. increases the size by one
 !! DOES check for capicity
 ** used during the constructing the PHANTOM
 atoms
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::cpy_pst(int iatm)
{
    reserve(1);
    memcpy(&vec[vec_sz*dim],&vec[iatm*dim],byte_sz);
    vec_sz++;
}
/*-----------------------------------------------
 _   _   _____   _____       ___   _____   _____
| | | | |  _  \ |  _  \     /   | |_   _| | ____| 
| | | | | |_| | | | | |    / /| |   | |   | |__   
| | | | |  ___/ | | | |   / / | |   | |   |  __|  
| |_| | | |     | |_| |  / /  | |   | |   | |___  
\_____/ |_|     |_____/ /_/   |_|   |_|   |_____| 
 
 -----------------------------------------------*/
/*--------------------------------------------
 1. copies no atoms, given by lst,
 using buff as destenation
 2. moves buff to the end of copied data
 3. does NOT increase/decrease the size
 ** used during the UPDATING the phantom atoms
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::cpy(byte*& buff,int* lst,int no)
{
    for(int i=0;i<no;i++)
    {
        memcpy(buff,&vec[lst[i]*dim],byte_sz);
        buff+=byte_sz;
    }
}
/*--------------------------------------------
 1. adds xtra_natms number of atom to the
 bottom of the vector using buff as source
 2. moves buff by specific stride
 3. increases the size by xtra_natms
 !! does NOT check for capicity
 ** used during the UPDATING the phantom atoms
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::pst(byte*& buff,int xtra_natms)
{
    memcpy(&vec[vec_sz*dim],buff,byte_sz*xtra_natms);
    buff+=byte_sz*xtra_natms;
    vec_sz+=xtra_natms;
}
/*--------------------------------------------
 1. copies no atoms, given by lst, to the
 bottom of vector
 2. increases the size by no
 !! does NOT check for capicity
 ** used during the UPDATING the phantom atoms
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::cpy_pst(int* lst,int no)
{
    T* vec_buff=&vec[vec_sz*dim];
    for(int i=0;i<no;i++)
    {
        memcpy(vec_buff,&vec[dim*lst[i]],byte_sz);
        vec_buff+=dim;
    }
    vec_sz+=no;
}

#endif
