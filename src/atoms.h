#ifndef __MAPP__atoms__
#define __MAPP__atoms__
#include <typeinfo>
#include <cmath>
#include <mpi.h>
#include "type_def.h"
#include "macros.h"
enum {snd_to_bhnd,snd_to_frnt};
enum {rcv_fm_frnt,rcv_fm_bhnd};
enum {bhnd,frnt};

/*-----------------------
 _     _   _____   _____  
| |   / / | ____| /  ___| 
| |  / /  | |__   | |     
| | / /   |  __|  | |     
| |/ /    | |___  | |___  
|___/     |_____| \_____| 
 
 -----------------------*/
namespace MAPP_NS {
    class Atoms;

    class vec
    {
    private:
    protected:
    public:
        int dim;
        int byte_sz;
        unsigned int vec_sz;
        unsigned int dump_vec_sz;
        unsigned int vec_cpcty;
        int orig_dim;
        int cmpctd_dims;
        
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
        T dump_def_val;
        Atoms* atoms;
    protected:
    public:
        Vec(Atoms*,int);
        Vec(Atoms*,int,const char*);
        Vec(Vec&);
        Vec(const Vec&&);
        ~Vec();
        
        void gather_dump(class vec*);
        void gather_dump(Vec<dmd_type>*);
        void gather_dump();
        void del_dump();
        void print_dump(FILE*,int);
        void assign_print_format();
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
    class ForceField;
    class Neighbor;
    class MAPP;
    class VecLst;
    class Error;
    class Timer;
    class Atoms
    {
        class Communincation;
        class Xchng;
        class Swap;
    private:
        unsigned long xchng_id;
        Error*& error;
        ForceField*& forcefield;
        Neighbor*& neighbor;
        Timer*& timer;
        FILE*& output;
        
        
        
        bool box_chng;
        bool grid_established;
        
        
        type0 skin;
        type0 max_cut;
        VecLst* vec_list;
        
        Vec<type0>* x0;
        Vec<int>* id_arch;
        
        void re_arrange(vec**,int);
    protected:
    public:
        
        Communincation* comm;
        Xchng* xchng;
        Swap* swap;
        
        MPI_Comm& world;
        vec** vecs;
        int nvecs;
        void add_vec(vec*);
        void add_vec(vec*,const char*);
        vec* find_vec(const char*);
        void del_vec(vec*);

        Vec<type0>* x;
        Vec<int>* id;
        
        const int dimension;
        type0** H;
        type0** B;
        
        int natms;
        int natms_ph;
        int tot_natms;

        int tot_p;
        int my_p;
        
        
        Atoms(MAPP*,int);
        ~Atoms();


        void set_skin(type0);
        type0 get_skin();
        
        type0* max_cut_s;
        type0* s_lo;
        type0* s_hi;
        
        void x2s(int);
        void s2x(int);
        void restart();
        void auto_grid();
        void man_grid(int*);
        void insert(byte*,vec**,int,int);
        void add();
        void del(int&);
        void init_xchng();
        void fin_xchng();
        
        void init(VecLst*,bool);
        void fin();
        void update(vec**,int);
        void update(vec*);
        void reset();
        bool xchng_chk(unsigned long&);
        
    };
}
/*------------------------------------------------------------------------------------------------------------------------
 _____   _____       ___  ___       ___  ___   _   _   __   _   _   __   _   _____       ___   _____   _   _____   __   _  
/  ___| /  _  \     /   |/   |     /   |/   | | | | | |  \ | | | | |  \ | | /  ___|     /   | |_   _| | | /  _  \ |  \ | | 
| |     | | | |    / /|   /| |    / /|   /| | | | | | |   \| | | | |   \| | | |        / /| |   | |   | | | | | | |   \| | 
| |     | | | |   / / |__/ | |   / / |__/ | | | | | | | |\   | | | | |\   | | |       / / | |   | |   | | | | | | | |\   | 
| |___  | |_| |  / /       | |  / /       | | | |_| | | | \  | | | | | \  | | |___   / /  | |   | |   | | | |_| | | | \  | 
\_____| \_____/ /_/        |_| /_/        |_| \_____/ |_|  \_| |_| |_|  \_| \_____| /_/   |_|   |_|   |_| \_____/ |_|  \_|
 ------------------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class Atoms::Communincation
    {
    private:
        MPI_Comm& world;
        
        
        const int dimension;
        /*
         * dimension of the domain
         * set in: constructor
         */
        
        int tot_n;
        /*
         * total number of nodes
         * set in: constructor
         */
        
        int my_n;
        /*
         * my node id
         * set in: constructor
         */
        
        int* p_per_n;
        /*
         * number of processors per each node
         * 1st input: node number
         * allocated in: constructor
         * deallocated in: destructor
         * set in: constructor
         */
        
        int** n_p_grid;
        /*
         * the processor ids in each node
         * 1st input: node number
         * 2nd input: ith processor in that node
         * 0 <= i < p_per_n[node]
         * n_p_grid[i][j] id of jth proc in ith node
         * allocated in: constructor
         * deallocated in: destructor
         * set in: constructor
         */
        
        
    protected:
    public:
        
        int tot_p;
        /*
         * total number of procs. in this simulation
         * set in: constructor
         */
        int my_p;
        /*
         * the id of my proc.
         * set in: constructor
         */
        
        int* tot_p_grid;
        /*
         * tot_p_grid[idm]: number of processors in idim dimension
         * allocated in: constructor
         * deallocated in: destructor
         * set in: auto_grid() / man_grid ()
         */
        
        int* my_loc;
        /*
         * the location of my proc. in proc. grid (domain)
         * my_loc[idim] my location in dimension idim
         * allocated in: constructor
         * deallocated in: destructor
         * set in: auto_grid() / man_grid ()
         */
        
        int** neigh_p;
        /*
         * the number of my neighbor processors
         * 1st input: dimension
         * 2nd input: direction; 0 for behind me, 1 for infront of me
         * neigh_p[idim][idir]: id of my neighbor proc. in my idim and idir direction
         * allocated in: constructor
         * deallocated in: destructor
         * set in: auto_grid() / man_grid ()
         */
        
        type0* s_lo;
        type0* s_hi;
        
        Communincation(Atoms*);
        ~Communincation();
        
        void auto_grid(type0**);
        /*
         * does an automatic griding of procs.
         * gets H as an input
         */
        
        void man_grid(int*);
        /*
         * does a manual griding of procs.
         * gets the grid as input
         */
        
    };
}
/*----------------------------------------------
 _     _   _____   _____   _       _____   _____  
| |   / / | ____| /  ___| | |     /  ___/ |_   _| 
| |  / /  | |__   | |     | |     | |___    | |   
| | / /   |  __|  | |     | |     \___  \   | |   
| |/ /    | |___  | |___  | |___   ___| |   | |   
|___/     |_____| \_____| |_____| /_____/   |_|
 
 ----------------------------------------------*/

namespace MAPP_NS
{
    class VecLst
    {
    private:
    protected:
    public:
        vec** xchng_vecs;
        int nxchng_vecs;
        int tot_xchng_byte_sz;
        vec** updt_vecs;
        int nupdt_vecs;
        int tot_updt_byte_sz;
        
        vec** arch_vecs;
        int narch_vecs;
        int tot_arch_byte_sz;
        
        VecLst(Atoms*);
       
        ~VecLst();
        
        void add_updt(vec*);
        void del_updt(vec*);
        void add_xchng(vec*);
        void del_xchng(vec*);
        void add_arch(vec*);
        void del_arch(vec*);
        
        void add_updt(vec**,int);
        void del_updt(vec**,int);
        void add_xchng(vec**,int);
        void del_xchng(vec**,int);
        void add_arch(vec**,int);
        void del_arch(vec**,int);
        
    };
}
/*---------------------------------------
__    __  _____   _   _   __   _   _____  
\ \  / / /  ___| | | | | |  \ | | /  ___| 
 \ \/ /  | |     | |_| | |   \| | | |     
  }  {   | |     |  _  | | |\   | | |  _
 / /\ \  | |___  | | | | | | \  | | |_| | 
/_/  \_\ \_____| |_| |_| |_|  \_| \_____/

 ---------------------------------------*/
namespace MAPP_NS
{
    class Atoms::Xchng
    {
    private:
        MPI_Comm& world;
        unsigned long& xchng_id;
        
        int** neigh_p;
        int my_p;
        
        int buff_grw;
        byte** snd_buff;
        int* snd_buff_sz;
        int* snd_buff_cpcty;
        
        byte* rcv_buff;
        int rcv_buff_sz;
        int rcv_buff_cpcty;
        
        type0* s_lo;
        type0* s_hi;
        void load(int&,int);
        void load(byte*&,int&);
        int xchng_buff(int,int);
        
        vec** xchng_vecs;
        int nxchng_vecs;
        int tot_xchng_sz;
        
        const int dimension;
        
        Vec<type0>*& x;
    protected:
    public:        
        Xchng(Atoms*,vec**,int);
        ~Xchng();
        void full_xchng();
    };
}
/*--------------------------------------
 _____   _          __      ___   _____  
/  ___/ | |        / /     /   | |  _  \ 
| |___  | |  __   / /     / /| | | |_| | 
\___  \ | | /  | / /     / / | | |  ___/ 
 ___| | | |/   |/ /     / /  | | | |     
/_____/ |___/|___/     /_/   |_| |_|
 
 --------------------------------------*/
namespace MAPP_NS
{
    class Atoms::Swap
    {
    private:
        Atoms* atoms;
        MPI_Comm& world;
        int**& neigh_p;
        type0**& H;
        type0**& B;
        type0*& s_lo;
        type0*& s_hi;
        type0*& max_cut_s;
        int my_p;
        
        
        int tot_ncomms;
        int* strt;
        int* ncomms;
        int* pbc_correction;
        type0* s_bnd;
        
        int** snd_atms_lst;
        int* snd_atms_lst_sz;
        int* snd_atms_lst_cpcty;
        int max_snd_atms_lst_sz;
        int snd_atms_lst_grw;

        int* rcv_atms_lst_sz;
        int max_rcv_atms_lst_sz;

        byte* snd_buff;        
        int snd_buff_sz;
        int snd_buff_cpcty;
        int snd_buff_grw;

        byte* rcv_buff;
        int rcv_buff_sz;
        int rcv_buff_cpcty;
        int rcv_buff_grw;
        

        
        vec** swap_vecs;
        int nswap_vecs;
        int tot_swap_vecs_sz;
        
        Vec<type0>*& x;
        
        const int dimension;
        void add_to_snd_lst(int&,int&);
        void reserve_snd_buff(int);
        void reserve_rcv_buff(int);
        
        class LoadUnLoadUpdate;
        class LoadUnLoadUpdateComm;
        class LoadUnLoadUpdateSelfComm;
        LoadUnLoadUpdate** comm_manager;
        
    protected:
    public:
        Swap(Atoms*,vec**,int);
        ~Swap();
        void reset();
        void update(vec*,bool);
        void update_(vec*,bool);
        void update(vec**,int,bool);
        void update_(vec**,int,bool);
        void list();
        void eliminate_redundancy();
    };

    
    class Atoms::Swap::LoadUnLoadUpdate
    {
    private:
    protected:
    public:
        LoadUnLoadUpdate(){};
        virtual ~LoadUnLoadUpdate(){};
        virtual void load(int&)=0;
        virtual void unload(int&,int&,int&)=0;
        virtual void load_unload(int&,int&,int&)=0;
        virtual void update_mult(int&,int&,int&,vec**&,int&,int&)=0;
        virtual void update_sing(int&,int&,int&,vec*&)=0;
        virtual void xchng_buff(int&,int&,byte*&,int&,int&,byte*&)=0;
    };
    
    class Atoms::Swap::LoadUnLoadUpdateComm : public Atoms::Swap::LoadUnLoadUpdate
    {
    private:
        MPI_Comm& world;
        
        int**& snd_atms_lst;
        int*& snd_atms_lst_sz;

        int*& rcv_atms_lst_sz;

        byte*& snd_buff;
        int& snd_buff_sz;
        int& snd_buff_cpcty;
        int& snd_buff_grw;
        
        byte*& rcv_buff;
        int& rcv_buff_sz;
        int& rcv_buff_cpcty;
        int& rcv_buff_grw;
        
        vec**& swap_vecs;
        int& nswap_vecs;
        int& tot_swap_vecs_sz;
    protected:
    public:
        LoadUnLoadUpdateComm(Swap*);
        void load(int&);
        void unload(int&,int&,int&);
        void load_unload(int&,int&,int&);
        void update_mult(int&,int&,int&,vec**&,int&,int&);
        void update_sing(int&,int&,int&,vec*&);
        void xchng_buff(int&,int&,byte*&,int&,int&,byte*&);
    };
    
    class Atoms::Swap::LoadUnLoadUpdateSelfComm : public Atoms::Swap::LoadUnLoadUpdate
    {
    private:
        int**& snd_atms_lst;
        int*& snd_atms_lst_sz;
        
        int*& rcv_atms_lst_sz;
        
        vec**& swap_vecs;
        int& nswap_vecs;
    protected:
    public:
        LoadUnLoadUpdateSelfComm(Swap*);
        void load(int&);
        void unload(int&,int&,int&);
        void load_unload(int&,int&,int&);
        void update_mult(int&,int&,int&,vec**&,int&,int&);
        void update_sing(int&,int&,int&,vec*&);
        void xchng_buff(int&,int&,byte*&,int&,int&,byte*&);
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
    dim=dim_;
    byte_sz=sizeof(T)*dim;
    vec_sz=vec_cpcty=atoms->natms;
    vec=new T[vec_sz*dim];
    atoms->add_vec(this);
    orig_dim=dim;
    cmpctd_dims=0;
    dump_vec_sz=0;
    dump_def_val=0;
    name=NULL;
    assign_print_format();
}
/*--------------------------------------------
 constructor with name;
 --------------------------------------------*/
template<typename T>
Vec<T>::Vec(Atoms* atoms_,int dim_,const char* name_):
atoms(atoms_),
vec()
{
    dim=dim_;
    byte_sz=sizeof(T)*dim;
    vec_sz=vec_cpcty=atoms->natms;
    vec=new T[vec_sz*dim];
    atoms->add_vec(this);
    orig_dim=dim;
    cmpctd_dims=0;
    dump_vec_sz=0;
    dump_def_val=0;
    name=NULL;
    assign_print_format();
    
    int length=static_cast<int>(strlen(name_))+1;
    name=new char[length];
    memcpy(name,name_,length*sizeof(char));
}
/*--------------------------------------------
 copy costructor
 --------------------------------------------*/
template<typename T>
Vec<T>::Vec(Vec<T>& old)
{
    atoms=old.atoms;
    
    dim=old.dim;
    byte_sz=old.byte_sz;
    vec_sz=old.vec_sz;
    vec_cpcty=old.vec_cpcty;
    vec=old.vec;
    orig_dim=old.orig_dim;
    cmpctd_dims=old.cmpctd_dims;
    dump_vec_sz=old.dump_vec_sz;
    name=old.name;
    print_format=old.print_format;
}
/*--------------------------------------------
 move costructor
 --------------------------------------------*/
template<typename T>
Vec<T>::Vec(const Vec<T>&& old)
{
    atoms=old.atoms;
    
    dim=old.dim;
    byte_sz=old.byte_sz;
    vec_sz=old.vec_sz;
    vec_cpcty=old.vec_cpcty;
    vec=old.vec;
    orig_dim=old.orig_dim;
    cmpctd_dims=old.cmpctd_dims;
    dump_vec_sz=old.dump_vec_sz;
    name=old.name;
    print_format=old.print_format;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<typename T>
inline Vec<T>::~Vec()
{
    if(vec_cpcty)
        delete [] vec;
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

    
    T* new_vec=new T[vec_cpcty*dim_];
    
    for(int i=0;i<vec_sz;i++)
        for(int j=0;j<min_dim;j++)
            new_vec[i*dim_+j]=vec[i*dim+j];
    
    if(vec_cpcty)
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
void Vec<T>::change_dimension(T dump_def_val_,int cmpctd_dims_,int del_dim)
{
    int dim_=dim-del_dim;
    cmpctd_dims=cmpctd_dims_;
    dump_def_val=dump_def_val_;
    
    if(dim_==dim)
        return;
    
    
    int min_dim=MIN(dim_,dim);

    
    T* new_vec=new T[vec_cpcty*dim_];
    
    for(int i=0;i<vec_sz;i++)
        for(int j=0;j<min_dim;j++)
            new_vec[i*dim_+j]=vec[i*dim+j];
    
    if(vec_cpcty)
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
    
    T* vec_tmp=new T[natms*dim];
    
    memcpy(vec_tmp,vec,byte_sz*vec_sz);
    
    if(vec_cpcty) delete [] vec;
    
    vec_sz=natms;
    vec_cpcty=natms;
    vec=vec_tmp;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::shrink_to_fit()
{
    if(vec_cpcty==vec_sz)
        return;
    
    T* vec_tmp=new T[vec_sz*dim];
    
    memcpy(vec_tmp,vec,byte_sz*vec_sz);
    
    if(vec_cpcty)
        delete [] vec;
    vec_cpcty=vec_sz;
    vec=vec_tmp;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::reserve(int xtra_atms)
{
    if(vec_cpcty>=vec_sz+xtra_atms)
        return;
    
    T* vec_tmp=new T[(vec_sz+xtra_atms)*dim];
    
    memcpy(vec_tmp,vec,byte_sz*vec_sz);
    
    if(vec_cpcty)
        delete [] vec;
    
    vec_cpcty=vec_sz+xtra_atms;
    vec=vec_tmp;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
inline void Vec<T>::rearrange(int* old_pos
,int* new_pos,int sz,int new_vec_cpcty)
{
    T* vec_tmp=new T[new_vec_cpcty*dim];
    
    for(int i=0;i<sz;i++)
        memcpy(vec_tmp+new_pos[i]*dim,vec+old_pos[i]*dim,byte_sz);
    
    if(vec_cpcty)
        delete [] vec;
    vec=vec_tmp;
    vec_cpcty=new_vec_cpcty;
    vec_sz=0;
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
void Vec<T>::assign_print_format()
{
    print_format=NULL;
    int lngth;
    if(typeid(T)==typeid(char))
    {
        lngth=static_cast<int>(strlen((char*)"%c "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%c ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(unsigned char))
    {
        lngth=static_cast<int>(strlen((char*)"%c "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%c ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(short int))
    {
        lngth=static_cast<int>(strlen((char*)"%d "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%d ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(unsigned short int))
    {
        lngth=static_cast<int>(strlen((char*)"%d "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%d ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(int))
    {
        lngth=static_cast<int>(strlen((char*)"%d "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%d ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(unsigned int))
    {
        lngth=static_cast<int>(strlen((char*)"%d "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%d ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(long int))
    {
        lngth=static_cast<int>(strlen((char*)"%ld "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%ld ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(unsigned long int))
    {
        lngth=static_cast<int>(strlen((char*)"%ld "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%ld ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(long long int))
    {
        lngth=static_cast<int>(strlen((char*)"%lld "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%lld ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(unsigned long long int))
    {
        lngth=static_cast<int>(strlen((char*)"%lld "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%lld ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(float))
    {
        lngth=static_cast<int>(strlen((char*)"%f "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%f ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(double))
    {
        lngth=static_cast<int>(strlen((char*)"%lf "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%lf ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(long double))
    {
        lngth=static_cast<int>(strlen((char*)"%Lf "))+1;
        print_format=new char[lngth];
        memcpy(print_format,(char*)"%Lf ",lngth*sizeof(char));
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void Vec<T>::gather_dump(class vec* map)
{
    Vec<dmd_type>* map_dmd_type
    =dynamic_cast<Vec<dmd_type>*>(map);
    if(map_dmd_type!=NULL)
    {
        gather_dump(map_dmd_type);
        return;
    }
}
/*--------------------------------------------
 
                            cmpctd_dims
                          ______/\______
                         |              |
                 orig_dim
  __________________/\__________________
 |                                      |
 x-----------------------x--------------x
 
 
 after transformation:
 
 
                    dim-orig_dim+cmpctd_dims
                          ____/\____
                         |          |
                   dim
  __________________/\______________
 |                                  |
 x-----------------------x----------x
 |___________  __________|
             \/
 
     orig_dim-cmpctd_dims
 
 --------------------------------------------*/
template<typename T>
void Vec<T>::gather_dump(Vec<dmd_type>* map)
{
    int natms=atoms->natms;
    T* vec_=vec;
    int byte_sz_=byte_sz;
    
    if(cmpctd_dims)
    {
        dmd_type* map_vec=map->begin();
        int map_dim=map->dim;
        int ncmp=(dim-orig_dim+cmpctd_dims)/map_dim;
        int map_dim_=cmpctd_dims/ncmp;
        int strt_dim=orig_dim-cmpctd_dims;
        
        vec_=new T[natms*orig_dim];
        byte_sz_=orig_dim*sizeof(T);
        
        int ipos=0;
        int ipos_=0;
        int map_ipos=0;
        for(int iatm=0;iatm<natms;iatm++)
        {
            for(int idim=0;idim<strt_dim;idim++)
                vec_[ipos_++]=vec[ipos++];
            
            for(int idim=0;idim<cmpctd_dims;idim++)
                vec_[ipos_+idim]=dump_def_val;
            
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
    MPI_Comm& world=atoms->world;
    
    if(my_p!=0)
        MPI_Send(&natms,1,MPI_INT,0,my_p,world);
    
    if(my_p==0)
    {
        rcv_size=new int[tot_p];
        MPI_Status status;
        for(int ip=1;ip<tot_p;ip++)
            MPI_Recv(&rcv_size[ip],1
            ,MPI_INT,ip,ip,world,&status);

        dump_vec=new T[orig_dim*tot_natms];
        dump_vec_sz=tot_natms;
        memcpy(dump_vec,vec_,natms*byte_sz_);
    }
    
    if(my_p!=0)
        MPI_Send(vec_,natms*byte_sz_,MPI_BYTE,0,my_p,world);
    
    if(my_p==0)
    {
        MPI_Status status;
        int tot_atoms=natms;
        for(int ip=1;ip<tot_p;ip++)
        {
            MPI_Recv(&dump_vec[tot_atoms*orig_dim]
            ,(byte_sz_*rcv_size[ip]),MPI_BYTE
            ,ip,ip,world,&status);
            tot_atoms+=rcv_size[ip];
        }
        
        delete [] rcv_size;
    }
    
    if(cmpctd_dims && natms)
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
    MPI_Comm& world=atoms->world;
    
    if(my_p!=0)
        MPI_Send(&natms,1,MPI_INT,0,my_p,world);
    
    if(my_p==0)
    {
        rcv_size=new int[tot_p];
        MPI_Status status;
        for(int ip=1;ip<tot_p;ip++)
            MPI_Recv(&rcv_size[ip],1
                     ,MPI_INT,ip,ip,world,&status);
        
        dump_vec=new T[orig_dim*tot_natms];
        dump_vec_sz=tot_natms;
        memcpy(dump_vec,vec,natms*byte_sz);
    }
    
    if(my_p!=0)
        MPI_Send(vec,natms*byte_sz,MPI_BYTE,0,my_p,world);
    
    if(my_p==0)
    {
        MPI_Status status;
        int tot_atoms=natms;
        for(int ip=1;ip<tot_p;ip++)
        {
            MPI_Recv(&dump_vec[tot_atoms*dim]
            ,(byte_sz*rcv_size[ip]),MPI_BYTE
            ,ip,ip,world,&status);
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
    if(atoms->my_p==0 && dump_vec_sz)
        delete [] dump_vec;
    dump_vec_sz=0;
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
