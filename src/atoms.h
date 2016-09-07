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
        Atoms<MAPP_NS::dimension>* atoms;
    protected:
    public:
        Vec(Atoms<MAPP_NS::dimension>*,int);
        Vec(Atoms<MAPP_NS::dimension>*,int,const char*);
        Vec(Vec&);
        Vec(const Vec&&);
        ~Vec();
        
        void gather_dump(class vec*);
        void gather_dump(Vec<dmd_type>*);
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
        
        void compress(Vec<dmd_type>*&,int,T);
        void decompress_dump(Vec<dmd_type>*&);

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
    class VecLst;
    
    template<const int>
    class Atoms
    {
        class Communincation;
        class Xchng;
        class Swap;
    private:
        unsigned long xchng_id;
        
        
        
        bool box_chng;
        bool grid_established;
        
        
        type0 skin;
        type0 max_cut;
        VecLst* vec_list;
        
        Vec<type0>* x0;
        Vec<int>* id_arch;
        
        void store_x0();
        bool decide();
        
        void store_arch_vecs();
        void restore_arch_vecs_();
        void restore_arch_vecs();
        void reorder_vecs(VecLst*);
    protected:
    public:
        
        Communincation* comm;
        Xchng* xchng;
        Swap* swap;
        
        vec** vecs;
        int nvecs;
        void add_vec(vec*);
        void add_vec(vec*,const char*);
        vec* find_vec(const char*);
        void del_vec(vec*);

        Vec<type0>* x;
        Vec<int>* id;
        
        type0** H;
        type0** B;
        
        int natms;
        int natms_ph;
        int tot_natms;

        int tot_p;
        int my_p;
        
        
        Atoms();
        ~Atoms();


        void set_skin(type0);
        type0 get_skin();
        
        type0* max_cut_s;
        type0* s_lo;
        type0* s_hi;
        
        void x2s(int);
        void s2x(int);
        void x2s_lcl();
        void s2x_lcl();
        void x2s_all();
        void s2x_all();
        
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
    template<const int _dimension_>
    class Atoms<_dimension_>::Communincation
    {
    private:
        
        
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
        
        int& tot_p;
        /*
         * total number of procs. in this simulation
         * set in: constructor
         */
        int& my_p;
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
        
        type0*& s_lo;
        type0*& s_hi;
        
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
        
        VecLst(Atoms<MAPP_NS::dimension>*);
       
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
    template<const int _dimension_>
    class Atoms<_dimension_>::Swap
    {
    private:
        /*things that reference cannot be removed*/
        int& natms;
        int& natms_ph;
        
        /*things that reference can be removed*/
        type0**& H;
        type0**& B;
        type0*& max_cut_s;
        int& my_p;
        int**& neigh_p;
        type0*& s_lo;
        type0*& s_hi;
        type0& max_cut;
        Vec<type0>*& x;
        int& nvecs;
        vec**& vecs;
        int& nxchng_vecs;
        
        
        int tot_ncomms;
        int ncomms[2*_dimension_];
        bool pbc_correction[2*_dimension_];
        type0 s_bnd[2*_dimension_];
        
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
        
        void add_to_snd_lst(int&,int&);
        void reserve_snd_buff(int);
        void reserve_rcv_buff(int);
        
        class LoadUnLoadUpdate;
        class LoadUnLoadUpdateComm;
        class LoadUnLoadUpdateSelfComm;
        LoadUnLoadUpdate* comm_manager[_dimension_];
        
    protected:
    public:
        Swap(Atoms*,vec**,int);
        ~Swap();
        void reset();
        void update(vec*,bool);
        void update(vec**,int,bool);
        void list();
        void rm_rdndncy();
    };

    template<const int _dimension_>
    class Atoms<_dimension_>::Swap::LoadUnLoadUpdate
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
    
    
    template<const int _dimension_>
    class Atoms<_dimension_>::Swap::LoadUnLoadUpdateComm : public Atoms::Swap::LoadUnLoadUpdate
    {
    private:
        
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
#ifndef BLOCK_COMM
        MPI_Request request[2];
#endif
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
    
    template<const int _dimension_>
    class Atoms<_dimension_>::Swap::LoadUnLoadUpdateSelfComm : public Atoms::Swap::LoadUnLoadUpdate
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
    template<const int _dimension_>
    class Atoms<_dimension_>::Xchng
    {
    private:
        
        /*things that reference cannot be removed*/
        int& natms;
        unsigned long& xchng_id;
        
        /*things that reference cannot be removed*/
        int& my_p;
        int**& neigh_p;
        type0*& s_lo;
        type0*& s_hi;
        Vec<type0>*& x;
        
        int buff_grw;
        byte* snd_buff[2];
        int snd_buff_sz[2];
        int snd_buff_cpcty[2];
        
        byte* rcv_buff;
        int rcv_buff_sz;
        int rcv_buff_cpcty;
        
        void load(int&,int);
        void load(byte*&,int&);
        int xchng_buff(int,int);
        
        vec** xchng_vecs;
        int nxchng_vecs;
        int tot_xchng_sz;
        
#ifndef BLOCK_COMM
        MPI_Request request[2];
#endif
    protected:
    public:        
        Xchng(Atoms*,vec**,int);
        ~Xchng();
        void full_xchng();
    };
}
using namespace MAPP_NS;
#include "xmath.h"
#include "error.h"
#include "timer.h"
#include "ff.h"
#include "neighbor.h"
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
Vec<T>::Vec(Atoms<MAPP_NS::dimension>* atoms_,int dim_):
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
Vec<T>::Vec(Atoms<MAPP_NS::dimension>* atoms_,int dim_,const char* name_):
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
    if(other.name!=NULL)
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
void Vec<T>::decompress_dump(Vec<dmd_type>*& map)
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
    dmd_type* map_vec=map->begin();
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
void Vec<T>::compress(Vec<dmd_type>*& map,int ncmprsd_dims_,T def_val_)
{
    def_val=def_val_;
    
    if(map)
    {
        int map_dim=map->dim;
        if(map_dim==ncmprsd_dims_)
            return;
        
        T* vec_=vec+dim-ncmprsd_dims_;
        dmd_type* map_vec=map->begin();
        
        for(int i=0;i<vec_sz;i++,map_vec+=map_dim,vec_+=dim)
            for(int j=0;j<map_dim;j++)
                vec_[j]=vec_[map[j]];
        int dim_reduc=ncmprsd_dims_-map_dim;
        change_dimension(dim-dim_reduc);
        orig_dim+=dim_reduc;
        return;
    }
    
    int map_dim=ncmprsd_dims_;
    map=new Vec<dmd_type>(atoms,map_dim);
    
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
    MPI_Allreduce(&max_dim_lcl,&dim_reduc,1,MPI_INT,MPI_MAX,world);
    dim_reduc=map_dim-dim_reduc;
    
    if(dim_reduc==0)
    {
        dmd_type* map_vec=map->begin();
        for(int i=0;i<vec_sz;i++,map_vec+=map_dim)
            for(dmd_type j=0;j<map_dim;j++)
                map_vec[j]=j;
        return;
    }
    
    vec_=vec+dim-map_dim;
    dmd_type* map_vec=map->begin();
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
    if(typeid(T)==typeid(char))
        assign_print_format("%c ");
    else if(typeid(T)==typeid(unsigned char))
        assign_print_format("%c ");
    else if(typeid(T)==typeid(short int))
        assign_print_format("%d ");
    else if(typeid(T)==typeid(unsigned short int))
        assign_print_format("%d ");
    else if(typeid(T)==typeid(int))
        assign_print_format("%d ");
    else if(typeid(T)==typeid(unsigned int))
        assign_print_format("%d ");
    else if(typeid(T)==typeid(long int))
        assign_print_format("%ld ");
    else if(typeid(T)==typeid(unsigned long int))
        assign_print_format("%ld ");
    else if(typeid(T)==typeid(long long int))
        assign_print_format("%lld ");
    else if(typeid(T)==typeid(unsigned long long int))
        assign_print_format("%lld ");
    else if(typeid(T)==typeid(float))
        assign_print_format("%f ");
    else if(typeid(T)==typeid(double))
        assign_print_format("%lf ");
    else if(typeid(T)==typeid(long double))
        assign_print_format("%Lf ");
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
void Vec<T>::gather_dump(Vec<dmd_type>* map)
{
    int natms=atoms->natms;
    T* vec_=vec;
    int byte_sz_=byte_sz;
    
    if(ncmprsd_dims)
    {
        dmd_type* map_vec=map->begin();
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
        MPI_Send(&natms,1,MPI_INT,0,my_p,world);
    
    if(my_p==0)
    {
        rcv_size=new int[tot_p];
        for(int ip=1;ip<tot_p;ip++)
            MPI_Recv(&rcv_size[ip],1
            ,MPI_INT,ip,ip,world,MPI_STATUS_IGNORE);

        dump_vec=new T[orig_dim*tot_natms];
        dump_vec_sz=tot_natms;
        memcpy(dump_vec,vec_,natms*byte_sz_);
    }
    
    if(my_p!=0)
        MPI_Send(vec_,natms*byte_sz_,MPI_BYTE,0,my_p,world);
    
    if(my_p==0)
    {
        int tot_atoms=natms;
        for(int ip=1;ip<tot_p;ip++)
        {
            MPI_Recv(&dump_vec[tot_atoms*orig_dim]
            ,(byte_sz_*rcv_size[ip]),MPI_BYTE
            ,ip,ip,world,MPI_STATUS_IGNORE);
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
        MPI_Send(&natms,1,MPI_INT,0,my_p,world);
    
    if(my_p==0)
    {
        rcv_size=new int[tot_p];
        for(int ip=1;ip<tot_p;ip++)
            MPI_Recv(&rcv_size[ip],1
            ,MPI_INT,ip,ip,world,MPI_STATUS_IGNORE);
        
        dump_vec=new T[orig_dim*tot_natms];
        dump_vec_sz=tot_natms;
        memcpy(dump_vec,vec,natms*byte_sz);
    }
    
    if(my_p!=0)
        MPI_Send(vec,natms*byte_sz,MPI_BYTE,0,my_p,world);
    
    if(my_p==0)
    {
        int tot_atoms=natms;
        for(int ip=1;ip<tot_p;ip++)
        {
            MPI_Recv(&dump_vec[tot_atoms*dim]
            ,(byte_sz*rcv_size[ip]),MPI_BYTE
            ,ip,ip,world,MPI_STATUS_IGNORE);
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
/*--------------------------------------
 _____   _          __      ___   _____
/  ___/ | |        / /     /   | |  _  \
| |___  | |  __   / /     / /| | | |_| |
\___  \ | | /  | / /     / / | | |  ___/
 ___| | | |/   |/ /     / /  | | | |
/_____/ |___/|___/     /_/   |_| |_|
 
 --------------------------------------*/
template<const int _dimension_>
Atoms<_dimension_>::Swap::Swap(Atoms<_dimension_>* atoms,vec** vecs_,int nvecs_):
natms(atoms->natms),
natms_ph(atoms->natms_ph),
H(atoms->H),
B(atoms->B),
max_cut_s(atoms->max_cut_s),

my_p(atoms->my_p),
neigh_p(atoms->comm->neigh_p),
s_lo(atoms->s_lo),
s_hi(atoms->s_hi),
max_cut(atoms->max_cut),
x(atoms->x),
nvecs(atoms->nvecs),
vecs(atoms->vecs),
nxchng_vecs(atoms->vec_list->nxchng_vecs),

swap_vecs(vecs_),
nswap_vecs(nvecs_)
{

    snd_buff=NULL;
    snd_buff_cpcty=0;
    snd_buff_grw=1024;
    
    rcv_buff=NULL;
    rcv_buff_cpcty=0;
    rcv_buff_grw=1024;
    snd_atms_lst_grw=8;
    
    tot_swap_vecs_sz=0;
    for(int ivec=0;ivec<nswap_vecs;ivec++)
        tot_swap_vecs_sz+=swap_vecs[ivec]->byte_sz;
    
    int icurs=0;
    for(int idim=0;idim<_dimension_;idim++)
    {
        if(my_p==neigh_p[idim][0] && my_p==neigh_p[idim][1])
            comm_manager[idim]=new LoadUnLoadUpdateSelfComm(this);
        else
            comm_manager[idim]=new LoadUnLoadUpdateComm(this);
        
        // snd_to_bhnd && rcv_fm_frnt
        if(atoms->comm->my_loc[idim]==atoms->comm->tot_p_grid[idim]-1)
            pbc_correction[icurs]=true;
        else
            pbc_correction[icurs]=false;
        icurs++;
        
        // snd_to_frnt && rcv_fm_bhnd
        if(atoms->comm->my_loc[idim]==0)
            pbc_correction[icurs]=true;
        else
            pbc_correction[icurs]=false;
        icurs++;
    }

    tot_ncomms=0;
    snd_atms_lst=NULL;
    snd_atms_lst_cpcty=snd_atms_lst_sz=rcv_atms_lst_sz=NULL;
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
Atoms<_dimension_>::Swap::~Swap()
{
    delete [] rcv_buff;
    delete [] snd_buff;
    
    for(int i=0;i<tot_ncomms;i++)
        delete [] snd_atms_lst[i];
    
    delete [] snd_atms_lst;
    delete [] snd_atms_lst_cpcty;
    delete [] snd_atms_lst_sz;
    delete [] rcv_atms_lst_sz;
    
    
    for(int idim=0;idim<_dimension_;idim++)
        delete comm_manager[idim];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::Swap::reset()
{
    for(int idim=0;idim<_dimension_;idim++)
    {
        type0 tmp=0.0;
        for(int jdim=idim;jdim<_dimension_;jdim++)
            tmp+=B[jdim][idim]*B[jdim][idim];
        
        max_cut_s[idim]=sqrt(tmp);
    }
    for(int idim=0;idim<_dimension_;idim++)
        max_cut_s[idim]*=max_cut;
    
    int icurs=0;
    int tot_ncomms_=0;
    for(int idim=0;idim<_dimension_;idim++)
    {
        // snd_to_bhnd && rcv_fm_frnt
        s_bnd[icurs]=s_lo[idim]+max_cut_s[idim];
        tot_ncomms_+=static_cast<int>(max_cut_s[idim]/(s_hi[idim]-s_lo[idim]))+1;
        ncomms[icurs++]=tot_ncomms_;
        
        // snd_to_frnt && rcv_fm_bhnd
        s_bnd[icurs]=s_hi[idim]-max_cut_s[idim];
        tot_ncomms_+=static_cast<int>(max_cut_s[idim]/(s_hi[idim]-s_lo[idim]))+1;
        ncomms[icurs++]=tot_ncomms_;
    }
    
    if(tot_ncomms_==tot_ncomms)
        return;
    
    for(int i=0;i<tot_ncomms;i++)
        delete [] snd_atms_lst[i];
    delete [] snd_atms_lst;
    delete [] snd_atms_lst_cpcty;
    delete [] snd_atms_lst_sz;
    delete [] rcv_atms_lst_sz;
    
    tot_ncomms=tot_ncomms_;
    snd_atms_lst=new int*[tot_ncomms];
    snd_atms_lst_cpcty=new int[tot_ncomms];
    snd_atms_lst_sz=new int[tot_ncomms];
    rcv_atms_lst_sz=new int[tot_ncomms];
    for(int i=0;i<tot_ncomms;i++)
    {
        snd_atms_lst_cpcty[i]=0;
        snd_atms_lst[i]=NULL;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::add_to_snd_lst(int& icomm,int& iatm)
{

    if(snd_atms_lst_sz[icomm]+1>snd_atms_lst_cpcty[icomm])
    {
        int* tmp_lst=new int[snd_atms_lst_sz[icomm]+1+snd_atms_lst_grw];
        memcpy(tmp_lst,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]*sizeof(int));
        delete [] snd_atms_lst[icomm];
        snd_atms_lst[icomm]=tmp_lst;
        snd_atms_lst_cpcty[icomm]=snd_atms_lst_sz[icomm]+1+snd_atms_lst_grw;
    }
    snd_atms_lst[icomm][snd_atms_lst_sz[icomm]]=iatm;
    snd_atms_lst_sz[icomm]++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::reserve_rcv_buff(int xtra)
{
    if(rcv_buff_cpcty<xtra+rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[xtra+rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=xtra+rcv_buff_sz+rcv_buff_grw;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::reserve_snd_buff(int xtra)
{
    if(snd_buff_cpcty<xtra+snd_buff_sz)
    {
        delete [] snd_buff;
        snd_buff=new byte[xtra+snd_buff_sz+snd_buff_grw];
        snd_buff_cpcty=xtra+snd_buff_sz+snd_buff_grw;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::update(vec** updt_vecs,int nupdt_vecs,bool x_xst)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    int tot_byte_sz=0;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        updt_vecs[ivec]->vec_sz=natms;
        tot_byte_sz+=updt_vecs[ivec]->byte_sz;
    }
    snd_buff_sz=rcv_buff_sz=0;
    reserve_snd_buff(tot_byte_sz*max_snd_atms_lst_sz);
    reserve_rcv_buff(tot_byte_sz*max_rcv_atms_lst_sz);
    

    int icurs=0;
    int icomm=0;
    for(int idim=0;idim<_dimension_;idim++)
    {
        for(int idir=0;idir<2;idir++)
        {
            while(icomm<ncomms[icurs])
            {
                comm_manager[idim]->update_mult(icomm
                ,neigh_p[idim][idir]
                ,neigh_p[idim][1-idir]
                ,updt_vecs,nupdt_vecs,tot_byte_sz);
                
                if(x_xst && pbc_correction[icurs])
                {
                    if(idir)
                    {
                        for(int iatm=x->vec_sz-1;iatm>x->vec_sz-rcv_atms_lst_sz[icomm]-1;iatm--)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                x_vec[iatm*x_dim+jdim]-=H[idim][jdim];
                    }
                    else
                    {
                        for(int iatm=x->vec_sz-1;iatm>x->vec_sz-rcv_atms_lst_sz[icomm]-1;iatm--)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                x_vec[iatm*x_dim+jdim]+=H[idim][jdim];
                    }
                }
                icomm++;
            }
            icurs++;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::update(vec* updt_vec,bool x_xst)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    snd_buff_sz=0;
    reserve_snd_buff(updt_vec->byte_sz*max_snd_atms_lst_sz);
    updt_vec->vec_sz=natms;
    
    int icurs=0;
    int icomm=0;
    for(int idim=0;idim<_dimension_;idim++)
    {
        for(int idir=0;idir<2;idir++)
        {
            while(icomm<ncomms[icurs])
            {
                comm_manager[idim]->update_sing(icomm
                ,neigh_p[idim][idir]
                ,neigh_p[idim][1-idir]
                ,updt_vec);
                
                if(x_xst && pbc_correction[icurs])
                {
                    if(idir)
                    {
                        for(int iatm=x->vec_sz-1;iatm>x->vec_sz-rcv_atms_lst_sz[icomm]-1;iatm--)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                x_vec[iatm*x_dim+jdim]-=H[idim][jdim];
                    }
                    else
                    {
                        for(int iatm=x->vec_sz-1;iatm>x->vec_sz-rcv_atms_lst_sz[icomm]-1;iatm--)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                x_vec[iatm*x_dim+jdim]+=H[idim][jdim];
                    }
                }
                icomm++;
            }
            icurs++;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::list()
{
    natms_ph=0;
    type0*& x_vec=x->begin();
    int x_dim=x->dim;
    int icurs=0;
    int icomm=0;
    int lo_atm,hi_atm;
    int last_atm;
    bool dir;
    type0 inc;
    max_snd_atms_lst_sz=max_rcv_atms_lst_sz=0;
    
    for(int idim=0;idim<_dimension_;idim++)
    {
        
        last_atm=x->vec_sz;
        
        inc=1.0;
        dir=true;
        for(int idir=0;idir<2;idir++)
        {
            lo_atm=0;
            hi_atm=last_atm;
            while(icomm<ncomms[icurs])
            {
                snd_buff_sz=0;
                rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm]=0;
                
                for(int iatm=lo_atm;iatm<hi_atm;iatm++)
                    if((x_vec[iatm*x_dim+idim]<s_bnd[icurs])==dir)
                        add_to_snd_lst(icomm,iatm);
                
                comm_manager[idim]->load_unload(icomm
                ,neigh_p[idim][idir]
                ,neigh_p[idim][1-idir]);
                
                lo_atm=x->vec_sz-rcv_atms_lst_sz[icomm];
                hi_atm=x->vec_sz;
                if(pbc_correction[icurs])
                    for(int iatm=lo_atm;iatm<hi_atm;iatm++)
                        x_vec[iatm*x_dim+idim]+=inc;
                
                max_snd_atms_lst_sz=MAX(max_snd_atms_lst_sz,snd_atms_lst_sz[icomm]);
                max_rcv_atms_lst_sz=MAX(max_rcv_atms_lst_sz,rcv_atms_lst_sz[icomm]);
                icomm++;
            }
            icurs++;
            inc-=2.0;
            dir=!dir;
        }
    }
    natms_ph=x->vec_sz-natms;
    for(int ivec=nxchng_vecs;ivec<nvecs;ivec++)
    {
        vecs[ivec]->vec_sz=0;
        vecs[ivec]->resize(x->vec_sz);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::rm_rdndncy()
{
    snd_buff_sz=rcv_buff_sz=0;
    reserve_rcv_buff(max_snd_atms_lst_sz);
    reserve_snd_buff(natms_ph);
    
    byte* mark=snd_buff;
    neighbor->mark_redndnt_ph(mark);
    
    
    int rcv_atms_lst_sz_;
    int snd_atms_lst_sz_=0;
    int snd_atms_lst_cpcty_=max_snd_atms_lst_sz;
    int* snd_atms_lst_=NULL;
    if(snd_atms_lst_cpcty_) snd_atms_lst_=new int[snd_atms_lst_cpcty_];
    
    int nlocomm;
    byte* mark_=mark+natms_ph;
    int icurs=2*_dimension_-1;
    int jcomm=tot_ncomms-1;
    max_snd_atms_lst_sz=max_rcv_atms_lst_sz=0;
    for(int idim=_dimension_-1;idim>-1;idim--)
    {
        for(int idir=1;idir>-1;idir--)
        {
            if(icurs!=0)
                nlocomm=ncomms[icurs-1];
            else
                nlocomm=0;
            
            while(jcomm>nlocomm-1)
            {
                mark_-=rcv_atms_lst_sz[jcomm];
                comm_manager[idim]->xchng_buff(
                neigh_p[idim][1-idir],rcv_atms_lst_sz[jcomm],mark_,
                neigh_p[idim][idir],snd_atms_lst_sz[jcomm],rcv_buff);
                
                snd_atms_lst_sz_=0;
                for(int i=0; i<snd_atms_lst_sz[jcomm];i++)
                {
                    if(rcv_buff[i]=='1')
                    {
                        if(snd_atms_lst[jcomm][i]>=natms)
                            mark[snd_atms_lst[jcomm][i]-natms]='1';
                        snd_atms_lst_[snd_atms_lst_sz_++]=snd_atms_lst[jcomm][i];
                    }
                }
                memcpy(snd_atms_lst[jcomm],snd_atms_lst_,snd_atms_lst_sz_*sizeof(int));
                snd_atms_lst_sz[jcomm]=snd_atms_lst_sz_;
                
                rcv_atms_lst_sz_=0;
                for(int i=0; i<rcv_atms_lst_sz[jcomm];i++)
                    if(mark_[i]=='1')
                        rcv_atms_lst_sz_++;
                rcv_atms_lst_sz[jcomm]=rcv_atms_lst_sz_;
                
                max_snd_atms_lst_sz=MAX(max_snd_atms_lst_sz,snd_atms_lst_sz[jcomm]);
                max_rcv_atms_lst_sz=MAX(max_rcv_atms_lst_sz,rcv_atms_lst_sz[jcomm]);
                
                jcomm--;
            }
            icurs--;
        }
    }
    
    delete [] snd_atms_lst_;
    
    int old_2_new_cpcty=natms+natms_ph;
    int* old_2_new=NULL;
    if(old_2_new_cpcty) old_2_new=new int[old_2_new_cpcty];
    
    int list_sz=0;
    int list_cpcty=natms_ph;
    int* list=NULL;
    if(list_cpcty) list=new int[list_cpcty];
    
    for(int iatm=0;iatm<natms;iatm++)
        old_2_new[iatm]=iatm;
    
    icurs=natms;
    for(int iatm=natms;iatm<natms+natms_ph;iatm++)
        if(mark[iatm-natms]=='1')
        {
            old_2_new[iatm]=icurs++;
            list[list_sz++]=iatm;
        }
    
    int new_natms_ph=list_sz;

    for(int icomm=0;icomm<tot_ncomms;icomm++)
        for(int i=0; i<snd_atms_lst_sz[icomm];i++)
            snd_atms_lst[icomm][i]=old_2_new[snd_atms_lst[icomm][i]];

    neighbor->rename_atoms(old_2_new);
    
    delete [] old_2_new;
    
    int* list_=list;
    
    int vec_sz=natms;
    while (*list_==natms+icurs)
    {
        list_++;
        vec_sz++;
        list_sz--;
    }
    
    for(int ivec=0;ivec<nswap_vecs;ivec++)
    {
        swap_vecs[ivec]->vec_sz=vec_sz;
        swap_vecs[ivec]->cpy_pst(list_,list_sz);
    }
    
    delete [] list;
    
    natms_ph=new_natms_ph;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
Atoms<_dimension_>::Swap::LoadUnLoadUpdateComm::
LoadUnLoadUpdateComm(Swap* swap):
LoadUnLoadUpdate(),

snd_atms_lst(swap->snd_atms_lst),
snd_atms_lst_sz(swap->snd_atms_lst_sz),

rcv_atms_lst_sz(swap->rcv_atms_lst_sz),

snd_buff(swap->snd_buff),
snd_buff_sz(swap->snd_buff_sz),
snd_buff_cpcty(swap->snd_buff_cpcty),
snd_buff_grw(swap->snd_buff_grw),

rcv_buff(swap->rcv_buff),
rcv_buff_sz(swap->rcv_buff_sz),
rcv_buff_cpcty(swap->rcv_buff_cpcty),
rcv_buff_grw(swap->rcv_buff_grw),

swap_vecs(swap->swap_vecs),
nswap_vecs(swap->nswap_vecs),
tot_swap_vecs_sz(swap->tot_swap_vecs_sz)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateComm::load
(int& iatm)
{
    if(snd_buff_cpcty<snd_buff_sz+tot_swap_vecs_sz)
    {
        byte* tmp_snd_buff=new byte[snd_buff_sz+tot_swap_vecs_sz+snd_buff_grw];
        memcpy(tmp_snd_buff,snd_buff,snd_buff_sz);
        delete [] snd_buff;
        snd_buff=tmp_snd_buff;
        snd_buff_cpcty=snd_buff_sz+tot_swap_vecs_sz+snd_buff_grw;
    }
    
    byte* tmp_snd_buff=snd_buff+snd_buff_sz;
    for(int ivec=0;ivec<nswap_vecs;ivec++)
        swap_vecs[ivec]->cpy(tmp_snd_buff,iatm);
    snd_buff_sz+=tot_swap_vecs_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateComm::unload
(int& icomm,int& snd_p,int& rcv_p)
{
    
    MPI_Sendrecv(&snd_atms_lst_sz[icomm],1,MPI_INT,snd_p,0,
                 &rcv_atms_lst_sz[icomm],1,MPI_INT,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
    
    rcv_buff_sz=rcv_atms_lst_sz[icomm]*tot_swap_vecs_sz;
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=rcv_atms_lst_sz[icomm]*tot_swap_vecs_sz+rcv_buff_grw;
    }

#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff,snd_buff_sz,MPI_BYTE,snd_p,0,
                 rcv_buff,rcv_buff_sz,MPI_BYTE,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_atms_lst_sz[icomm])
        MPI_Irecv(rcv_buff,rcv_buff_sz,
        MPI_BYTE,rcv_p,
        0,
        world,&request[0]);
    
    if(snd_atms_lst_sz[icomm])
        MPI_Isend(snd_buff,snd_buff_sz,
        MPI_BYTE,snd_p,
        0,
        world,&request[1]);
    
    if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm])
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    
    byte* tmp_rcv_buff0=rcv_buff;
    for(int ivec=0;ivec<nswap_vecs;ivec++)
        swap_vecs[ivec]->pst(tmp_rcv_buff0,tot_swap_vecs_sz,rcv_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateComm::load_unload
(int& icomm,int& snd_p,int& rcv_p)
{
    snd_buff_sz=snd_atms_lst_sz[icomm]*tot_swap_vecs_sz;
    if(snd_buff_cpcty<snd_buff_sz)
    {
        delete [] snd_buff;
        snd_buff=new byte[snd_buff_sz+snd_buff_grw];
        snd_buff_cpcty=snd_buff_sz+snd_buff_grw;
    }
    
    byte* tmp_snd_buff=snd_buff;
    for(int ivec=0;ivec<nswap_vecs;ivec++)
        swap_vecs[ivec]->cpy(tmp_snd_buff,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
    
    MPI_Sendrecv(&snd_atms_lst_sz[icomm],1,MPI_INT,snd_p,0,
                 &rcv_atms_lst_sz[icomm],1,MPI_INT,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
    
    rcv_buff_sz=rcv_atms_lst_sz[icomm]*tot_swap_vecs_sz;
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=rcv_atms_lst_sz[icomm]*tot_swap_vecs_sz+rcv_buff_grw;
    }

#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff,snd_buff_sz,MPI_BYTE,snd_p,0,
                 rcv_buff,rcv_buff_sz,MPI_BYTE,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_atms_lst_sz[icomm])
        MPI_Irecv(rcv_buff,rcv_buff_sz,
        MPI_BYTE,rcv_p,
        0,
        world,&request[0]);
    
    if(snd_atms_lst_sz[icomm])
        MPI_Isend(snd_buff,snd_buff_sz,
        MPI_BYTE,snd_p,
        0,
        world,&request[1]);
    
    if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm])
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    
    byte* tmp_rcv_buff=rcv_buff;
    for(int ivec=0;ivec<nswap_vecs;ivec++)
    {
        swap_vecs[ivec]->reserve(rcv_atms_lst_sz[icomm]);
        swap_vecs[ivec]->pst(tmp_rcv_buff,rcv_atms_lst_sz[icomm]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateComm::update_mult
(int& icomm,int& snd_p,int& rcv_p,vec**& vecs
,int& nvecs,int& vecs_byte_sz)
{
    byte* tmp_snd_buff=snd_buff;
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->cpy(tmp_snd_buff,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff,snd_atms_lst_sz[icomm]*vecs_byte_sz,MPI_BYTE,snd_p,0,
                 rcv_buff,rcv_atms_lst_sz[icomm]*vecs_byte_sz,MPI_BYTE,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_atms_lst_sz[icomm])
        MPI_Irecv(rcv_buff,rcv_atms_lst_sz[icomm]*vecs_byte_sz,
        MPI_BYTE,rcv_p,
        0,
        world,&request[0]);
    
    if(snd_atms_lst_sz[icomm])
        MPI_Isend(snd_buff,snd_atms_lst_sz[icomm]*vecs_byte_sz,
        MPI_BYTE,snd_p,
        0,
        world,&request[1]);
    
    if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm])
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    byte* tmp_rcv_buff=rcv_buff;
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->pst(tmp_rcv_buff,rcv_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateComm::update_sing
(int& icomm,int& snd_p,int& rcv_p,vec*& v)
{
    byte* tmp_snd_buff=snd_buff;
    v->cpy(tmp_snd_buff,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff,snd_atms_lst_sz[icomm]*v->byte_sz,MPI_BYTE,snd_p,0,
                 v->end(),rcv_atms_lst_sz[icomm]*v->byte_sz,MPI_BYTE,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_atms_lst_sz[icomm])
        MPI_Irecv(v->end(),rcv_atms_lst_sz[icomm]*v->byte_sz,
        MPI_BYTE,rcv_p,
        0,
        world,&request[0]);
    
    if(snd_atms_lst_sz[icomm])
        MPI_Isend(snd_buff,snd_atms_lst_sz[icomm]*v->byte_sz,
        MPI_BYTE,snd_p,
        0,
        world,&request[1]);
    
    if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm])
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    v->vec_sz+=rcv_atms_lst_sz[icomm];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateComm::xchng_buff
(int& snd_p,int& snd_buff_sz_,byte*& snd_buff_
,int& rcv_p,int& rcv_buff_sz_,byte*& rcv_buff_)
{
#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff_,snd_buff_sz_,MPI_BYTE,snd_p,0,
                 rcv_buff_,rcv_buff_sz_,MPI_BYTE,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_buff_sz_)
        MPI_Irecv(rcv_buff_,rcv_buff_sz_,
        MPI_BYTE,rcv_p,
        0,
        world,&request[0]);
    
    if(snd_buff_sz_)
        MPI_Isend(snd_buff_,snd_buff_sz_,
        MPI_BYTE,snd_p,
        0,
        world,&request[1]);
    
    if(rcv_buff_sz_ && snd_buff_sz_)
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_buff_sz_ && snd_buff_sz_==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_buff_sz_==0 && snd_buff_sz_)
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
Atoms<_dimension_>::Swap::LoadUnLoadUpdateSelfComm::
LoadUnLoadUpdateSelfComm(Swap* swap):
LoadUnLoadUpdate(),

snd_atms_lst(swap->snd_atms_lst),
snd_atms_lst_sz(swap->snd_atms_lst_sz),

rcv_atms_lst_sz(swap->rcv_atms_lst_sz),

swap_vecs(swap->swap_vecs),
nswap_vecs(swap->nswap_vecs)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateSelfComm::
load(int& iatm)
{
    for(int ivec=0;ivec<nswap_vecs;ivec++)
        swap_vecs[ivec]->cpy_pst(iatm);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateSelfComm::unload
(int& icomm,int& ,int& )
{
    rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateSelfComm::load_unload
(int& icomm,int&,int&)
{
    rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm];
    
    for(int ivec=0;ivec<nswap_vecs;ivec++)
    {
        swap_vecs[ivec]->reserve(rcv_atms_lst_sz[icomm]);
        swap_vecs[ivec]->cpy_pst(snd_atms_lst[icomm],rcv_atms_lst_sz[icomm]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateSelfComm::update_mult
(int& icomm,int&,int&,vec**& vecs
,int& nvecs,int&)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->cpy_pst(snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateSelfComm::update_sing
(int& icomm,int&,int&,vec*& v)
{
    v->cpy_pst(snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Swap::LoadUnLoadUpdateSelfComm::xchng_buff
(int&,int& snd_buff_sz_,byte*& snd_buff_
,int&,int& rcv_buff_sz_,byte*& rcv_buff_)
{
    memcpy(rcv_buff_,snd_buff_,rcv_buff_sz_);
}
/*---------------------------------------
__    __  _____   _   _   __   _   _____  
\ \  / / /  ___| | | | | |  \ | | /  ___| 
 \ \/ /  | |     | |_| | |   \| | | |     
  }  {   | |     |  _  | | |\   | | |  _
 / /\ \  | |___  | | | | | | \  | | |_| | 
/_/  \_\ \_____| |_| |_| |_|  \_| \_____/

 ---------------------------------------*/
template<const int _dimension_>
Atoms<_dimension_>::Xchng::Xchng(Atoms* atoms
,vec** x_vecs,int nx_vecs):
natms(atoms->natms),
xchng_id(atoms->xchng_id),

my_p(atoms->my_p),
neigh_p(atoms->comm->neigh_p),
s_lo(atoms->s_lo),
s_hi(atoms->s_hi),
x(atoms->x),

xchng_vecs(x_vecs),
nxchng_vecs(nx_vecs)
{
    buff_grw=1024;
    
    snd_buff[0]=snd_buff[1]=NULL;
    snd_buff_sz[0]=snd_buff_sz[1]=0;
    snd_buff_cpcty[0]=snd_buff_cpcty[1]=0;
    
    rcv_buff=NULL;
    rcv_buff_sz=0;
    rcv_buff_cpcty=0;
    
    tot_xchng_sz=0;
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        tot_xchng_sz+=xchng_vecs[ivec]->byte_sz;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<const int _dimension_>
Atoms<_dimension_>::Xchng::~Xchng()
{

    delete [] snd_buff[0];
    delete [] snd_buff[1];
    delete [] rcv_buff;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Xchng::load(int& iatm,int idir)
{
    if(snd_buff_cpcty[idir]<snd_buff_sz[idir]+tot_xchng_sz)
    {
        byte* tmp_buff=new byte[snd_buff_sz[idir]+tot_xchng_sz+buff_grw];
        memcpy(tmp_buff,snd_buff[idir],snd_buff_sz[idir]);
        delete [] snd_buff[idir];
        snd_buff[idir]=tmp_buff;
        snd_buff_cpcty[idir]=snd_buff_sz[idir]+tot_xchng_sz+buff_grw;
    }
    byte* tmp_buff=&snd_buff[idir][snd_buff_sz[idir]];
    
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        xchng_vecs[ivec]->pop_out(tmp_buff,iatm);
    
    snd_buff_sz[idir]+=tot_xchng_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Atoms<_dimension_>::Xchng::load(byte*& buff,int& idir)
{
    if(snd_buff_cpcty[idir]<snd_buff_sz[idir]+tot_xchng_sz)
    {
        byte* tmp_buff=new byte[snd_buff_sz[idir]+tot_xchng_sz+buff_grw];
        memcpy(tmp_buff,snd_buff[idir],snd_buff_sz[idir]);
        delete [] snd_buff[idir];
        snd_buff[idir]=tmp_buff;
        snd_buff_cpcty[idir]=snd_buff_sz[idir]+tot_xchng_sz+buff_grw;
    }
    memcpy(&snd_buff[idir][snd_buff_sz[idir]],buff,tot_xchng_sz);
    snd_buff_sz[idir]+=tot_xchng_sz;
    buff+=tot_xchng_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
int Atoms<_dimension_>::Xchng::xchng_buff(int idim,int idir)
{
    rcv_buff_sz=0;
    int max_snd_sz;
    MPI_Allreduce(&snd_buff_sz[idir],&max_snd_sz,1,MPI_INT,MPI_MAX,world);
    if(max_snd_sz==0)
        return 0;
    
    MPI_Sendrecv(&snd_buff_sz[idir],1,MPI_INT,neigh_p[idim][idir],0,
                 &rcv_buff_sz,1,MPI_INT,neigh_p[idim][1-idir],0,
                 world,MPI_STATUS_IGNORE);
    
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+buff_grw];
        rcv_buff_cpcty=rcv_buff_sz+buff_grw;
    }
    
#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff[idir],snd_buff_sz[idir],MPI_BYTE,neigh_p[idim][idir],0,
                 rcv_buff,rcv_buff_sz,MPI_BYTE,neigh_p[idim][1-idir],0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_buff_sz)
        MPI_Irecv(rcv_buff,rcv_buff_sz,
        MPI_BYTE,neigh_p[idim][1-idir],
        0,
        world,&request[0]);
    
    if(snd_buff_sz[idir])
        MPI_Isend(snd_buff[idir],snd_buff_sz[idir],
        MPI_BYTE,neigh_p[idim][idir],
        0,
        world,&request[1]);

    if(rcv_buff_sz && snd_buff_sz[idir])
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_buff_sz && snd_buff_sz[idir]==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_buff_sz==0 && snd_buff_sz[idir])
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    
    snd_buff_sz[idir]=0;
    
    /* 
     * here make sure, that we have enough space
     * for the potential atoms to be inserted in
     * xchng vectors
     */
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        xchng_vecs[ivec]->reserve(rcv_buff_sz/tot_xchng_sz);
    
    return 1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::Xchng::full_xchng()
{
    
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        xchng_vecs[ivec]->resize(natms);
    int disp;
    type0 s,ds_lo,ds_hi;
    int iatm;
    int xchng_lcl=0;
    
    for(int idim=0;idim<_dimension_;idim++)
    {
        if(my_p==neigh_p[idim][0] && my_p==neigh_p[idim][1])
            continue;
        
        disp=idim*sizeof(type0);
        
        snd_buff_sz[0]=snd_buff_sz[1]=0;
        iatm=0;
        while (iatm<x->vec_sz)
        {
            s=(*x)(iatm,idim);
            
            if(s_lo[idim]<=s && s<s_hi[idim])
                iatm++;
            else
            {
                if(s>=s_hi[idim])
                {
                    ds_hi=s-s_hi[idim];
                    ds_lo=1.0+s_lo[idim]-s;
                }
                else
                {
                    ds_hi=1.0+s-s_hi[idim];
                    ds_lo=s_lo[idim]-s;
                }
                
                if(ds_hi<ds_lo)
                    load(iatm,snd_to_frnt);
                else
                    load(iatm,snd_to_bhnd);                
              
            }
        }
        
        if(snd_buff_sz[0]+snd_buff_sz[1])
            xchng_lcl=1;
        
        byte* buff_tmp;
        for(int idir=0;idir<2;idir++)
        {
            while(xchng_buff(idim,idir))
            {
                buff_tmp=rcv_buff;
                while(rcv_buff_sz)
                {
                    s=*(type0*)(buff_tmp+disp);
                    if(s_lo[idim]<=s && s<s_hi[idim])
                        for(int ivec=0;ivec<nxchng_vecs;ivec++)
                            xchng_vecs[ivec]->pop_in(buff_tmp);
                    else
                        load(buff_tmp,idir);
                    
                    rcv_buff_sz-=tot_xchng_sz;
                }
            }
            
        }
    }
    int xchng;
    MPI_Allreduce(&xchng_lcl,&xchng,1,MPI_INT,MPI_MAX,world);
    if(xchng)
        xchng_id++;
    natms=x->vec_sz;
}
/*----------------------------------------------------------------------------------------------------------------
 _____   _____       ___  ___       ___  ___   _   _   __   _   _   _____       ___   _____   _   _____   __   _  
/  ___| /  _  \     /   |/   |     /   |/   | | | | | |  \ | | | | /  ___|     /   | |_   _| | | /  _  \ |  \ | | 
| |     | | | |    / /|   /| |    / /|   /| | | | | | |   \| | | | | |        / /| |   | |   | | | | | | |   \| | 
| |     | | | |   / / |__/ | |   / / |__/ | | | | | | | |\   | | | | |       / / | |   | |   | | | | | | | |\   | 
| |___  | |_| |  / /       | |  / /       | | | |_| | | | \  | | | | |___   / /  | |   | |   | | | |_| | | | \  | 
\_____| \_____/ /_/        |_| /_/        |_| \_____/ |_|  \_| |_| \_____| /_/   |_|   |_|   |_| \_____/ |_|  \_|
 
 ----------------------------------------------------------------------------------------------------------------*/
template<const int _dimension_>
Atoms<_dimension_>::Communincation::Communincation(Atoms* atoms):
s_lo(atoms->s_lo),
s_hi(atoms->s_hi),
my_p(atoms->my_p),
tot_p(atoms->tot_p)
{
    /* beginning of communication related parameters */
    
    MPI_Comm_rank(world,&my_p);
    MPI_Comm_size(world,&tot_p);
    
    char** names=new char*[tot_p];
    int* rank=new int[tot_p];
    char* name=new char[MPI_MAX_PROCESSOR_NAME];
    int name_sz;
    MPI_Get_processor_name(name,&name_sz);
    name_sz++;
    int name_sz_;
    for(int ip=0;ip<tot_p;ip++)
    {
        rank[ip]=ip;
        if(ip==my_p)
            name_sz_=name_sz;
        MPI_Bcast(&name_sz_,1,MPI_INT,ip,world);
        names[ip]=new char[name_sz_];
        if(ip==my_p)
            memcpy(names[ip],name,name_sz_*sizeof(char));
        MPI_Bcast(names[ip],name_sz_,MPI_CHAR,ip,world);
    }
    delete [] name;
    
    
    tot_n=0;
    p_per_n=NULL;
    int* p_per_n_;
    for(int i=0,r,last=0;i<tot_p;i++)
    {
        if(i==my_p)
            my_n=tot_n;
        for(int j=i+1;j<tot_p;j++)
            if(strcmp(names[i],names[j])==0)
            {
                if(j==my_p)
                    my_n=tot_n;
                name=names[i+1];
                names[i+1]=names[j];
                names[j]=name;
                r=rank[i+1];
                rank[i+1]=rank[j];
                rank[j]=r;
                i++;
            }
        
        p_per_n_=new int[tot_n+1];
        memcpy(p_per_n_,p_per_n,tot_n*sizeof(int));
        delete [] p_per_n;
        p_per_n=p_per_n_;
        p_per_n[tot_n++]=i+1-last;
        last=i+1;
    }
    for(int ip=0;ip<tot_p;ip++)
        delete [] names[ip];
    delete [] names;
    n_p_grid=new int*[tot_n];
    *n_p_grid=rank;
    for(int i=1;i<tot_n;i++)
        n_p_grid[i]=n_p_grid[i-1]+p_per_n[i-1];
    
    /* end of communication related parameters */
    
    tot_p_grid=new int[_dimension_];
    my_loc=new int[_dimension_];
    neigh_p=new int*[_dimension_];
    *neigh_p=new int[2*_dimension_];
    for(int idim=1;idim<_dimension_;idim++)
        neigh_p[idim]=neigh_p[idim-1]+2;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<const int _dimension_>
Atoms<_dimension_>::Communincation::~Communincation()
{
    delete [] *neigh_p;
    delete [] neigh_p;
    delete [] my_loc;
    delete [] tot_p_grid;
    
    delete [] *n_p_grid;
    delete [] n_p_grid;
    delete [] p_per_n;
}
/*--------------------------------------------
 automatically grid the domain
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::Communincation::auto_grid(type0** H)
{
    bool eq_p_per_n=true;
    for(int i=0;i<tot_n;i++)
        if(p_per_n[i]!=p_per_n[0])
            eq_p_per_n=false;
    
    int* fac_list;
    int fac_list_size;
    XMath::fac_list(tot_p,_dimension_,fac_list,fac_list_size);
    
    type0 ratio=-1.0;
    type0 tmp_ratio;
    
    if(eq_p_per_n && tot_n>1)
    {
        int prin_dimension=0;
        type0 min_area=1.0/H[0][0];
        
        for(int i=1;i<_dimension_;i++)
            if(1.0/H[i][i]<min_area)
            {
                min_area=1.0/H[i][i];
                prin_dimension=i;
            }

        int* ifac_list=fac_list;
        for(int ifac=0;ifac<fac_list_size;ifac++,ifac_list+=_dimension_)
        {
            if(ifac_list[prin_dimension]!=tot_n) continue;
                
            tmp_ratio=0.0;
            for(int i=0;i<_dimension_;i++)
                tmp_ratio+=static_cast<type0>(ifac_list[i])/H[i][i];
            
            if(tmp_ratio<ratio||ratio<0.0)
            {
                ratio=tmp_ratio;
                for(int i=0;i<_dimension_;i++)
                    tot_p_grid[i]=ifac_list[i];
            }
            
        }
        
        delete [] fac_list;
        
        
        if(_dimension_>1)
        {
            int* tmp_tot_p_grid=new int[_dimension_-1];
            int* dim_indx=new int[_dimension_-1];
            int* my_tmp_loc=new int[_dimension_-1];
            int* my_tmp_denom=new int[_dimension_-1];
            
            
            int pos=0;
            
            for(int i=0;i<_dimension_;i++)
            {
                if(i!=prin_dimension)
                {
                    tmp_tot_p_grid[pos]=tot_p_grid[i];
                    dim_indx[pos]=i;
                    pos++;
                }
            }
            
            /*
            int my_p_in_my_node=0;
            for(;p_per_n[my_p_in_my_node]==my_p;my_p_in_my_node++);
             */
            
            int my_p_in_my_node=-1;
            for(int i=0;i<p_per_n[my_n];i++)
                if(n_p_grid[my_n][i]==my_p)
                    my_p_in_my_node=i;
            
            for(int i=0;i<_dimension_-1;i++)
            {
                int no=1;
                for(int j=0;j<i;j++)
                    no*=tmp_tot_p_grid[j];
                my_tmp_denom[i]=no;
            }
            
            
            int ttmp=my_p_in_my_node;
            for(int i=_dimension_-2;i>-1;i--)
            {
                my_tmp_loc[i]=ttmp/my_tmp_denom[i];
                ttmp-=my_tmp_loc[i]*my_tmp_denom[i];
            }
            
            
            for(int i=0;i<_dimension_-1;i++)
            {
                my_loc[dim_indx[i]]=my_tmp_loc[i];
            }
            my_loc[prin_dimension]=my_n;
            
            
            
            int nop;
            int nom;
            int tmp0,tmp1;
            for(int i=0;i<_dimension_-1;i++)
            {
                
                tmp0=my_tmp_loc[i]-1;
                if(tmp0==-1)
                    tmp0=tmp_tot_p_grid[i]-1;
                
                tmp1=my_tmp_loc[i]+1;
                if(tmp1==tmp_tot_p_grid[i])
                    tmp1=0;
                
                nop=0;
                nom=0;
                
                for(int j=0;j<_dimension_-1;j++)
                {
                    if(i==j)
                    {
                        nom+=tmp0*my_tmp_denom[j];
                        nop+=tmp1*my_tmp_denom[j];
                    }
                    else
                    {
                        nom+=my_tmp_loc[j]*my_tmp_denom[j];
                        nop+=my_tmp_loc[j]*my_tmp_denom[j];
                    }
                }
                
                neigh_p[dim_indx[i]][1]=n_p_grid[my_n][nop];
                neigh_p[dim_indx[i]][0]=n_p_grid[my_n][nom];
            }
            
            tmp0=my_n-1;
            if(tmp0==-1)
                tmp0=tot_n-1;
            neigh_p[prin_dimension][0]=n_p_grid[tmp0][my_p_in_my_node];
            
            tmp1=my_n+1;
            if(tmp1==tot_n)
                tmp1=0;
            neigh_p[prin_dimension][1]=n_p_grid[tmp1][my_p_in_my_node];
            
            
            delete [] dim_indx;
            delete [] tmp_tot_p_grid;
            delete [] my_tmp_loc;
            delete [] my_tmp_denom;
            
            
        }
        else
        {
            tot_p_grid[0]=tot_p;
            
            int my_p_in_my_node=-1;
            for(int i=0;i<p_per_n[my_n];i++)
                if(n_p_grid[my_n][i]==my_p)
                    my_p_in_my_node=i;
            
            
            my_loc[0]=p_per_n[my_n]*my_n;
            my_loc[0]+=my_p_in_my_node;
            
            int tmp0p,tmp1p,tmp0n,tmp1n;
            tmp0p=my_p_in_my_node;
            tmp1p=my_p_in_my_node;
            tmp0n=my_n;
            tmp1n=my_n;
            
            tmp0p--;
            if(tmp0p==-1)
            {
                tmp0p=p_per_n[my_n]-1;
                tmp0n--;
                if(tmp0n==-1)
                    tmp0n=tot_n-1;
            }
            
            tmp1p++;
            if(tmp1p==p_per_n[my_n])
            {
                tmp1p=0;
                tmp1n++;
                if(tmp1n==tot_n)
                    tmp1n=0;
            }
            
            neigh_p[0][0]=n_p_grid[tmp0n][tmp0p];
            neigh_p[0][1]=n_p_grid[tmp1n][tmp1p];
            
        }
        
    }
    else
    {
        int* ifac_list=fac_list;
        for(int ifac=0;ifac<fac_list_size;ifac++,ifac_list+=_dimension_)
        {
            tmp_ratio=0.0;
            for(int i=0;i<_dimension_;i++)
                tmp_ratio+=static_cast<type0>(ifac_list[i])/H[i][i];
            
            if(tmp_ratio<ratio||ratio<0.0)
            {
                ratio=tmp_ratio;
                for(int i=0;i<_dimension_;i++)
                    tot_p_grid[i]=ifac_list[i];
            }
            
        }
        
        delete [] fac_list;
        
        int* list=new int[_dimension_];
        
        MPI_Comm cartesian;
        for(int i=0;i<_dimension_;i++)
            list[i]=1;
        MPI_Cart_create(world,_dimension_,tot_p_grid,list,1,&cartesian);
        MPI_Cart_get(cartesian,_dimension_,tot_p_grid,list,my_loc);
        for(int i=0;i<_dimension_;i++)
            MPI_Cart_shift(cartesian,i,1,&neigh_p[i][0],&neigh_p[i][1]);
        
        MPI_Comm_free(&cartesian);
        
        delete [] list;
    }
    
    
    
    for(int i=0;i<_dimension_;i++)
    {
        s_lo[i]=
        static_cast<type0>(my_loc[i])
        /static_cast<type0>(tot_p_grid[i]);
        
        s_hi[i]=
        static_cast<type0>(my_loc[i]+1)
        /static_cast<type0>(tot_p_grid[i]);
    }
}
/*--------------------------------------------
 manually grid the domain
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::Communincation::man_grid(int* n)
{
    for(int i=0;i<_dimension_;i++)
        tot_p_grid[i]=n[i];
    
    int* list=new int[_dimension_];
    for(int i=0;i<_dimension_;i++)
        list[i]=1;
    
    MPI_Comm cartesian;
    MPI_Cart_create(world,_dimension_,tot_p_grid,list,1,&cartesian);
    MPI_Cart_get(cartesian,_dimension_,tot_p_grid,list,my_loc);
    for(int i=0;i<_dimension_;i++)
        MPI_Cart_shift(cartesian,i,1,&neigh_p[i][0],&neigh_p[i][1]);
    
    MPI_Comm_free(&cartesian);
    
    delete [] list;
    
    for(int i=0;i<_dimension_;i++)
    {
        s_lo[i]=
        static_cast<type0>(my_loc[i])
        /static_cast<type0>(tot_p_grid[i]);
        
        s_hi[i]=
        (static_cast<type0>(my_loc[i])+1.0)
        /static_cast<type0>(tot_p_grid[i]);
    }
}
/*---------------------------------------------------------------------------
      ___   _____   _____       ___  ___   _____
     /   | |_   _| /  _  \     /   |/   | /  ___/
    / /| |   | |   | | | |    / /|   /| | | |___
   / / | |   | |   | | | |   / / |__/ | | \___  \
  / /  | |   | |   | |_| |  / /       | |  ___| |
 /_/   |_|   |_|   \_____/ /_/        |_| /_____/
 
 ---------------------------------------------------------------------------*/
template<const int _dimension_>
Atoms<_dimension_>::Atoms():
nvecs(0)
{
    vecs=NULL;
    id_arch=NULL;
    
    grid_established=false;
    box_chng=false;
    skin=0.5;
    tot_natms=0;
    natms=0;
    natms_ph=0;
    xchng_id=0;
    
    x=new Vec<type0>(this,_dimension_);
    id= new Vec<int>(this,1);
    

    max_cut_s=new type0[_dimension_];
    s_lo=new type0[_dimension_];
    s_hi=new type0[_dimension_];
    H=new type0*[_dimension_];
    B=new type0*[_dimension_];
    *H=new type0[_dimension_*_dimension_];
    *B=new type0[_dimension_*_dimension_];
    
    for(int idim=1;idim<_dimension_;idim++)
    {
        B[idim]=B[idim-1]+_dimension_;
        H[idim]=H[idim-1]+_dimension_;
    }
    
    for(int idim=0;idim<_dimension_;idim++)
    {
        max_cut_s[idim]=0.0;
        for(int jdim=0;jdim<_dimension_;jdim++)
            H[idim][jdim]=B[idim][jdim]=0.0;
    }
    
    comm=new Communincation(this);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<const int _dimension_>
Atoms<_dimension_>::~Atoms()
{
    delete comm;
    delete [] *B;
    delete [] B;
    delete [] *H;
    delete [] H;
    delete [] s_hi;
    delete [] s_lo;
    delete [] max_cut_s;

    while(nvecs)
        delete vecs[0];
}
/*--------------------------------------------
 add a new vec
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::add_vec(vec* v)
{
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
template<const int _dimension_>
void Atoms<_dimension_>::add_vec(vec* v,const char* name)
{
    // first see if the name is duplicate
    for(int ivec=0;ivec<nvecs;ivec++)
        if(strcmp(name,vecs[ivec]->name)==0)
            error->abort("vector %s already exist",name);
    
    
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
template<const int _dimension_>
vec* Atoms<_dimension_>::find_vec(const char* name)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        if(vecs[ivec]->name!=NULL && strcmp(name,vecs[ivec]->name)==0)
            return vecs[ivec];

    error->abort("vector %s does not exist",name);
    
    return NULL;
}
/*--------------------------------------------
 remove a vector
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::del_vec(vec* v)
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
 automatically grid the domain
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::auto_grid()
{
    type0 det=1.0;
    for(int idim=0;idim<_dimension_;idim++)
        det*=H[idim][idim]*H[idim][idim];
    if(det==0.0)
        error->abort("cannot perform automatic "
        "domain griding, the domain has not been established");
        
    comm->auto_grid(H);
    grid_established=true;
    
    if(my_p==0)
    {
        fprintf(output,"automatic grid performed: ");
        
        for(int i=0;i<_dimension_-1;i++)
            fprintf(output,"%dx",comm->tot_p_grid[i]);
        
        fprintf(output,"%d",comm->tot_p_grid[_dimension_-1]);
        
        fprintf(output,"\n");
    }
}
/*--------------------------------------------
 manually grid the domain
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::man_grid(int* n)
{
    comm->man_grid(n);
    xchng=new Xchng(this,vecs,nvecs);
    x2s_lcl();
    xchng->full_xchng();
    s2x_lcl();
    delete xchng;
    grid_established=true;
    
    if(my_p==0)
    {
        fprintf(output,"manual grid performed: ");
        
        for(int i=0;i<_dimension_-1;i++)
            fprintf(output,"%dx",comm->tot_p_grid[i]);
        
        fprintf(output,"%d",comm->tot_p_grid[_dimension_-1]);
        
        fprintf(output,"\n");
    }
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::x2s(int no)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
#ifdef UNROLL
    for(int i=0;i<no;i++)
        XMatrixVector::x2s<_dimension_>(x_vec+i*x_dim,B);
#else
    for(int i=0;i<no;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*B[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
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
template<const int _dimension_>
void Atoms<_dimension_>::s2x(int no)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
#ifdef UNROLL
    for(int i=0;i<no;i++)
        XMatrixVector::s2x<_dimension_>(x_vec+i*x_dim,H);
#else
    for(int i=0;i<no;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*H[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
                x_vec[idim]+=x_vec[jdim]*H[jdim][idim];
        }
#endif
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::x2s_lcl()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
#ifdef UNROLL
    for(int i=0;i<natms;i++)
        XMatrixVector::x2s<_dimension_>(x_vec+i*x_dim,B);
#else
    for(int i=0;i<natms;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*B[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
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
template<const int _dimension_>
void Atoms<_dimension_>::s2x_lcl()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
#ifdef UNROLL
    for(int i=0;i<natms;i++)
        XMatrixVector::s2x<_dimension_>(x_vec+i*x_dim,H);
#else
    for(int i=0;i<natms;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*H[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
                x_vec[idim]+=x_vec[jdim]*H[jdim][idim];
        }
#endif
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::x2s_all()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    int nall=natms+natms_ph;
#ifdef UNROLL
    for(int i=0;i<nall;i++)
        XMatrixVector::x2s<_dimension_>(x_vec+i*x_dim,B);
#else
    for(int i=0;i<nall;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*B[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
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
template<const int _dimension_>
void Atoms<_dimension_>::s2x_all()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    int nall=natms+natms_ph;
#ifdef UNROLL
    for(int i=0;i<nall;i++)
        XMatrixVector::s2x<_dimension_>(x_vec+i*x_dim,H);
#else
    for(int i=0;i<nall;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*H[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
                x_vec[idim]+=x_vec[jdim]*H[jdim][idim];
        }
#endif
}
/*--------------------------------------------
 insert a number of atoms
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::insert(byte* buff,vec** vecs_,int nvecs_,int natms_)
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
template<const int _dimension_>
void Atoms<_dimension_>::add()
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
template<const int _dimension_>
void Atoms<_dimension_>::del(int& del_idx)
{
    for(int i=0;i<nvecs;i++)
        vecs[i]->del(del_idx);
    natms--;
}
/*--------------------------------------------
 restart
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::restart()
{
    tot_natms=0;
    natms=0;
    natms_ph=0;
    while(nvecs)
        delete vecs[0];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
bool Atoms<_dimension_>::xchng_chk(unsigned long& xchng_id_)
{
    if(xchng_id_==xchng_id)
        return true;
    else
    {
        xchng_id_=xchng_id;
        return false;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::reorder_vecs(VecLst* vec_list_)
{
    int* vec_mapping=new int[nvecs];
    for(int ivec=0;ivec<nvecs;ivec++) vec_mapping[ivec]=-1;
    int vecs_mapped=0;
    for(int ivec=0;ivec<nvecs;ivec++)
    {
        int jvec=0;
        while (jvec<vec_list_->nxchng_vecs && vec_list_->xchng_vecs[jvec]!=vecs[ivec])
            jvec++;
        if(jvec!=vec_list_->nxchng_vecs)
        {
            vec_mapping[ivec]=jvec;
            vecs_mapped++;
        }
    }
    
    if(vecs_mapped!=vec_list_->nxchng_vecs)
        error->abort("not all xchng vectors were found");
    
    for(int ivec=0;ivec<nvecs;ivec++)
        if(vec_mapping[ivec]==-1)
        {
            vec_mapping[ivec]=vecs_mapped;
            vecs_mapped++;
        }
    
    vec** vecs_= new vec*[nvecs];
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs_[vec_mapping[ivec]]=vecs[ivec];
    if(nvecs)
    {
        delete [] vec_mapping;
        delete [] vecs;
    }
    vecs=vecs_;
}
/*--------------------------------------------
 init a simulation
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::init(VecLst* vec_list_,bool box_chng_)
{
    
    
    type0 det=1.0;
    for(int idim=0;idim<_dimension_;idim++)
        det*=H[idim][idim]*H[idim][idim];
    if(det==0.0)
        error->abort("cannot start the simulation "
        "domain griding, the domain is not defined");
    if(!grid_established)
        error->abort("cannot start the simulation "
        "domain grid has not been established");
    

    if(forcefield==NULL)
        error->abort("cannot start the simulation, force field has not been established");
    
    if(neighbor==NULL)
        error->abort("cannot start the simulation, neighbior list has not been established");
    
    
    timer->init();
    timer->start(COMM_TIME_mode);
    

    
    reorder_vecs(vec_list_);
    box_chng=box_chng_;
    vec_list=vec_list_;
    
    store_arch_vecs();
    
    x0=new Vec<type0>(this,_dimension_);
    


    xchng=new Xchng(this,vec_list->xchng_vecs,vec_list->nxchng_vecs);

    x2s_lcl();
    xchng->full_xchng();
    
    forcefield->init();
    max_cut=forcefield->max_cut();

    swap=new Swap(this,vec_list->updt_vecs,vec_list->nupdt_vecs);
    swap->list();
    
    neighbor->init();
    swap->rm_rdndncy();
    store_x0();
    
    timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::fin()
{
    delete xchng;
    delete swap;
    delete x0;

    restore_arch_vecs();

    for(int ivec=0;ivec<nvecs;ivec++)
    {
        vecs[ivec]->vec_sz=natms;
        vecs[ivec]->shrink_to_fit();
    }
    
    neighbor->fin();
    forcefield->fin();
    timer->fin();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::store_x0()
{
    type0* x_vec=x->begin();
    type0* x0_vec=x0->begin();
    int x_dim=x->dim;
    int last_atm=natms;
    if(box_chng) last_atm+=natms_ph;
    if(x_dim==_dimension_)
    {
        memcpy(x0_vec,x_vec,last_atm*_dimension_*sizeof(type0));
        return;
    }
    for(int iatm=0;iatm<last_atm;iatm++,x0_vec+=_dimension_,x_vec+=x_dim)
        memcpy(x0_vec,x_vec,_dimension_*sizeof(type0));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline bool Atoms<_dimension_>::decide()
{
    type0 skin_sq=0.25*skin*skin;
    type0 dx_sq;
    int succ,succ_lcl=1;
    type0* x_vec=x->begin();
    type0* x0_vec=x0->begin();
    int x_dim=x->dim;
    int last_atm=natms;
    if(box_chng) last_atm+=natms_ph;
    
    for(int iatm=0;succ_lcl && iatm<last_atm;iatm++,x0_vec+=_dimension_,x_vec+=x_dim)
    {
        dx_sq=0.0;
        for(int idim=0;idim<_dimension_;idim++)
            dx_sq+=(x0_vec[idim]-x_vec[idim])*(x0_vec[idim]-x_vec[idim]);
        if(dx_sq>skin_sq)
            succ_lcl=0;
    }
    
    MPI_Allreduce(&succ_lcl,&succ,1,MPI_INT,MPI_MIN,world);
    if(succ) return true;
    return false;
}
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::update(vec* updt_vec)
{
    update(&updt_vec,1);
}
/*--------------------------------------------
 update a number of vectors
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::update(vec** updt_vecs,int nupdt_vecs)
{
    timer->start(COMM_TIME_mode);
    bool x_xst=false;
    for(int ivec=0;x_xst==false && ivec<nupdt_vecs;ivec++)
        if(updt_vecs[ivec]==x)
            x_xst=true;
    if(x_xst==false)
    {
        if(nupdt_vecs==1)
            swap->update(updt_vecs[0],false);
        else
            swap->update(updt_vecs,nupdt_vecs,false);
        timer->stop(COMM_TIME_mode);
        return;
    }
    
    
    if(box_chng)
    {
        if(nupdt_vecs==1)
            swap->update(x,true);
        else
            swap->update(updt_vecs,nupdt_vecs,true);

        if(decide())
        {
            timer->stop(COMM_TIME_mode);
            return;
        }
        
        
        x2s_lcl();
        xchng->full_xchng();
        
        swap->reset();
        swap->list();
        
        neighbor->create_list(box_chng);

        store_x0();
    }
    else
    {
        if(decide())
        {
            if(nupdt_vecs==1)
                swap->update(x,true);
            else
                swap->update(updt_vecs,nupdt_vecs,true);
            
            timer->stop(COMM_TIME_mode);
            return;
        }
        
        x2s_lcl();
        xchng->full_xchng();
        
        swap->list();
        
        neighbor->create_list(box_chng);
        
        
        store_x0();
    }

    timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::init_xchng()
{
    timer->start(COMM_TIME_mode);

    x2s_lcl();
    xchng->full_xchng();
    swap->list();

    timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::fin_xchng()
{
    timer->start(COMM_TIME_mode);
    swap->list();
    neighbor->create_list(box_chng);
    store_x0();
    timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::set_skin(type0 skin_)
{
    skin=skin_;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
type0 Atoms<_dimension_>::get_skin()
{
    return skin;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::reset()
{
    xchng=new Xchng(this,vecs,nvecs);
    x2s_lcl();
    xchng->full_xchng();
    s2x_lcl();
    delete xchng;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::store_arch_vecs()
{
    if(vec_list->narch_vecs==0) return;
    id_arch=new Vec<int>(this,1);
    int* id_0=id->begin();
    int* id_1=id_arch->begin();
    memcpy(id_1,id_0,natms*sizeof(int));
    del_vec(id_arch);
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        del_vec(vec_list->arch_vecs[ivec]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::restore_arch_vecs()
{
    if(vec_list->narch_vecs==0) return;
    
    int byte_sz=0;
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        byte_sz+=vec_list->arch_vecs[ivec]->byte_sz;
    
    int natms_old=id_arch->vec_sz;
    int* id_old=id_arch->begin();
    int* key_old=NULL;
    if(natms_old) key_old=new int[natms_old];
    for(int i=0;i<natms_old;i++) key_old[i]=i;
    XMath::quicksort(key_old,key_old+natms_old
    ,[&id_old](int* rank_i,int* rank_j){return (id_old[*rank_i]<id_old[*rank_j]);}
    ,[](int* rank_i,int* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    int natms_new=id->vec_sz;
    int* id_new=id->begin();
    int* key_new=NULL;
    if(natms_new) key_new=new int[natms_new];
    for(int i=0;i<natms_new;i++) key_new[i]=i;
    XMath::quicksort(key_new,key_new+natms_new
    ,[&id_new](int* rank_i,int* rank_j){return (id_new[*rank_i]<id_new[*rank_j]);}
    ,[](int* rank_i,int* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    
    int nsnd=0;
    int* idx_lst_snd=NULL;
    int* id_lst_snd=NULL;
    if(natms_old)
    {
        idx_lst_snd=new int[natms_old];
        id_lst_snd=new int[natms_old];
    }
    
    int nrcv=0;
    int* idx_lst_rcv=NULL;
    int* id_lst_rcv=NULL;
    if(natms_new)
    {
        idx_lst_rcv=new int[natms_new];
        id_lst_rcv=new int[natms_new];
    }
    
    int nkeep=0;
    int* idx_keep_old=NULL;
    int* idx_keep_new=NULL;
    if(MIN(natms_old,natms_new))
    {
        idx_keep_old=new int[MIN(natms_old,natms_new)];
        idx_keep_new=new int[MIN(natms_old,natms_new)];
    }
    
    int i_old=0,i_new=0;
    for(;i_old<natms_old && i_new<natms_new;)
    {
        if(id_old[key_old[i_old]]<id_new[key_new[i_new]])
        {
            id_lst_snd[nsnd]=id_old[key_old[i_old]];
            idx_lst_snd[nsnd++]=key_old[i_old++];
        }
        else if(id_old[key_old[i_old]]>id_new[key_new[i_new]])
        {
            id_lst_rcv[nrcv]=id_new[key_new[i_new]];
            idx_lst_rcv[nrcv++]=key_new[i_new++];
        }
        else
        {
            idx_keep_old[nkeep]=key_old[i_old++];
            idx_keep_new[nkeep++]=key_new[i_new++];
        }
    }
    
    for(;i_old<natms_old;)
    {
        id_lst_snd[nsnd]=id_old[key_old[i_old]];
        idx_lst_snd[nsnd++]=key_old[i_old++];
    }
    for(;i_new<natms_new;)
    {
        id_lst_rcv[nrcv]=id_new[key_new[i_new]];
        idx_lst_rcv[nrcv++]=key_new[i_new++];
    }
    
    delete [] key_old;
    delete [] key_new;
    
    /*-------------------------------------------------*/
    
    auto shrnk_2_fit=
    [](int*& arr,const int sz,const int cpcty)->void
    {
        if(cpcty==0 || cpcty==sz)
            return;
        if(sz==0)
        {
            delete [] arr;
            arr=NULL;
        }
        
        int* _arr=new int[sz];
        memcpy(_arr,arr,sz*sizeof(int));
        delete [] arr;
        arr=_arr;
    };
    
    /*---------------------------------------------------*/
    
    shrnk_2_fit(id_lst_snd,nsnd,natms_old);
    shrnk_2_fit(idx_lst_snd,nsnd,natms_old);
    shrnk_2_fit(id_lst_rcv,nrcv,natms_new);
    shrnk_2_fit(idx_lst_rcv,nrcv,natms_new);
    shrnk_2_fit(idx_keep_old,nkeep,MIN(natms_old,natms_new));
    shrnk_2_fit(idx_keep_new,nkeep,MIN(natms_old,natms_new));

#ifdef ARCH_RESTORE_DEBUG
    int tot_nsnd;
    int tot_nrcv;
    MPI_Allreduce(&nsnd,&tot_nsnd,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&nrcv,&tot_nrcv,1,MPI_INT,MPI_SUM,world);
    if(tot_nsnd==tot_nrcv)
        error->abort("sanity check failed number of atoms to be sent(%d) "
        "is unequal to number of atoms to be recieved(%d)",tot_nsnd,tot_nrcv);
#endif
    
    auto sort_idx_by_p=
    [this] (const int* lst,int lst_sz,const int* sub_lst,int sub_lst_sz,
            int*& sub_lst_idx,int* comm_size)->void
    {
        for(int ip=0;ip<tot_p;ip++)
            comm_size[ip]=0;
        
        int max_lst_sz;
        MPI_Allreduce(&lst_sz,&max_lst_sz,1,MPI_INT,MPI_MAX,world);
        int mother_lst_sz;
        int* mother_lst=NULL;
        if(max_lst_sz) mother_lst=new int[max_lst_sz];
        
        int fnd_sz,ufnd_sz;
        
        int* ufnd=NULL;
        int* ufnd_idx=NULL;
        int* fnd_idx=NULL;
        if(sub_lst_sz)
        {
            ufnd=new int[sub_lst_sz];
            ufnd_idx=new int[sub_lst_sz];
            fnd_idx=new int[sub_lst_sz];
        }
        
        fnd_sz=0;
        ufnd_sz=sub_lst_sz;
        memcpy(ufnd,sub_lst,sub_lst_sz*sizeof(int));
        memcpy(ufnd_idx,sub_lst_idx,sub_lst_sz*sizeof(int));
        for(int imother_lst,iufnd,ufnd_tmp_sz,ip=0;ip<tot_p;ip++)
        {
            if(ip==my_p)
            {
                mother_lst_sz=lst_sz;
                memcpy(mother_lst,lst,mother_lst_sz*sizeof(int));
            }
            MPI_Bcast(&mother_lst_sz,1,MPI_INT,ip,world);
            MPI_Bcast(mother_lst,mother_lst_sz,MPI_INT,ip,world);
            
            if(ip==my_p || ufnd_sz==0 || mother_lst_sz==0) continue;
            
            ufnd_tmp_sz=0;
            for(imother_lst=0,iufnd=0;imother_lst<mother_lst_sz && iufnd<ufnd_sz;)
            {
                if(mother_lst[imother_lst]<ufnd[iufnd])
                    imother_lst++;
                else if(mother_lst[imother_lst]>ufnd[iufnd])
                {
                    ufnd_idx[ufnd_tmp_sz]=ufnd_idx[iufnd];
                    ufnd[ufnd_tmp_sz++]=ufnd[iufnd++];
                }
                else
                {
                    fnd_idx[fnd_sz++]=ufnd_idx[iufnd++];
                    imother_lst++;
                    comm_size[ip]++;
                }
            }
            
            for(;iufnd<ufnd_sz;)
            {
                ufnd_idx[ufnd_tmp_sz]=ufnd_idx[iufnd];
                ufnd[ufnd_tmp_sz++]=ufnd[iufnd++];
            }
            
            ufnd_sz=ufnd_tmp_sz;
            
        }
        
        delete [] sub_lst_idx;
        sub_lst_idx=fnd_idx;
        delete [] ufnd;
        delete [] ufnd_idx;
        delete [] mother_lst;
    };
   
    int* nsnd_comm=new int[tot_p];
    int* nrcv_comm=new int[tot_p];
    sort_idx_by_p(id_lst_rcv,nrcv,id_lst_snd,nsnd,idx_lst_snd,nsnd_comm);
    sort_idx_by_p(id_lst_snd,nsnd,id_lst_rcv,nrcv,idx_lst_rcv,nrcv_comm);
    delete [] id_lst_snd;
    delete [] id_lst_rcv;
 
#ifdef ARCH_RESTORE_DEBUG
    int comm_sz_snd=0,comm_sz_rcv=0;
    for(int i=0;i<tot_p;i++)
    {
        comm_sz_snd+=nsnd_comm[i];
        comm_sz_rcv+=nrcv_comm[i];
    }
    if(nsnd!=comm_sz_snd)
        error->abort_sing("sanity check failed just found %d out of %d atoms that should be sent");
    if(nrcv!=comm_sz_rcv)
        error->abort_sing("sanity check failed just found %d out of %d atoms that should be recieved");
#endif
    
    byte** snd_buff=new byte*[tot_p];
    *snd_buff=NULL;
    if (nsnd)
    {
        *snd_buff=new byte[byte_sz*nsnd];
        for(int ip=1;ip<tot_p;ip++)
            snd_buff[ip]=snd_buff[ip-1]+byte_sz*nsnd_comm[ip-1];
    }
    
    byte** rcv_buff=new byte*[tot_p];
    *rcv_buff=NULL;
    if(nrcv)
    {
        *rcv_buff=new byte[byte_sz*nrcv];
        for(int ip=1;ip<tot_p;ip++)
            rcv_buff[ip]=rcv_buff[ip-1]+byte_sz*nrcv_comm[ip-1];
    }
    
    byte* _snd_buff=*snd_buff;
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        for(int i=0;i<nsnd;i++)
            vec_list->arch_vecs[ivec]->cpy(_snd_buff,idx_lst_snd[i]);
    delete [] idx_lst_snd;
    
    
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        vec_list->arch_vecs[ivec]->rearrange(idx_keep_old,idx_keep_new,nkeep,natms_new);
    delete [] idx_keep_old;
    delete [] idx_keep_new;

#ifndef BLOCK_COMM
    MPI_Request request[2];
#endif
    for(int idisp=1;idisp<tot_p;idisp++)
    {
        int rcv_p=my_p-idisp;
        if(rcv_p<0) rcv_p+=tot_p;
        int snd_p=my_p+idisp;
        if(snd_p>=tot_p) snd_p-=tot_p;
#ifdef BLOCK_COMM
        MPI_Sendrecv(snd_buff[snd_p],nsnd_comm[snd_p]*byte_sz,MPI_BYTE,snd_p,0,
                     rcv_buff[rcv_p],nrcv_comm[rcv_p]*byte_sz,MPI_BYTE,rcv_p,0,
                     world,MPI_STATUS_IGNORE);
#else
        if(nrcv_comm[rcv_p])
            MPI_Irecv(rcv_buff[rcv_p],nrcv_comm[rcv_p]*byte_sz,
                      MPI_BYTE,rcv_p,
                      0,
                      world,&request[0]);
        
        if(nsnd_comm[snd_p])
            MPI_Isend(snd_buff[snd_p],nsnd_comm[snd_p]*byte_sz,
                      MPI_BYTE,snd_p,
                      0,
                      world,&request[1]);
        
        if(nrcv_comm[rcv_p] && nsnd_comm[snd_p])
            MPI_Waitall(2,request,MPI_STATUS_IGNORE);
        else if(nrcv_comm[rcv_p] && nsnd_comm[snd_p]==0)
            MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
        else if(nrcv_comm[rcv_p]==0 && nsnd_comm[snd_p])
            MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    }
    delete [] nsnd_comm;
    delete [] nrcv_comm;
    
    delete [] *snd_buff;
    delete [] snd_buff;
    
    byte* _rcv_buff=*rcv_buff;
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        for(int i=0;i<nrcv;i++)
            vec_list->arch_vecs[ivec]->cpy(_rcv_buff,idx_lst_rcv[i]);
    
    delete [] idx_lst_rcv;
    delete [] *rcv_buff;
    delete [] rcv_buff;
    
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        add_vec(vec_list->arch_vecs[ivec]);
    add_vec(id_arch);
    delete id_arch;
    id_arch=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
void Atoms<_dimension_>::restore_arch_vecs_()
{
    if(vec_list->narch_vecs==0) return;
    auto shrnk_2_fit=
    [](int*& arr,int sz,int cpcty)->void
    {
        if(cpcty==sz)
            return;
        if(sz==0)
        {
            delete [] arr;
            arr=NULL;
        }
        
        int* _arr=new int[sz];
        memcpy(_arr,arr,sz*sizeof(int));
        delete [] arr;
        arr=_arr;
    };
    
    auto sort_idx_by_p=
    [this] (const int* lst,int lst_sz,const int* sub_lst,int sub_lst_sz,
            int*& sub_lst_idx,int* comm_size)->void
    {
        for(int ip=0;ip<tot_p;ip++)
            comm_size[ip]=0;
        
        int max_lst_sz;
        MPI_Allreduce(&lst_sz,&max_lst_sz,1,MPI_INT,MPI_MAX,world);
        int* mother_lst=NULL;
        int mother_lst_sz;
        if(max_lst_sz)
            mother_lst=new int[max_lst_sz];
        
        int fnd_sz=0,ufnd_sz=sub_lst_sz,ufnd_tmp_sz;
        int* ufnd=NULL;
        int* ufnd_tmp=NULL;
        int* fnd_idx=NULL;
        int* ufnd_idx=NULL;
        int* ufnd_idx_tmp=NULL;
        if(sub_lst_sz)
        {
            ufnd=new int[sub_lst_sz];
            ufnd_tmp=new int[sub_lst_sz];
            fnd_idx=new int[sub_lst_sz];
            ufnd_idx=new int[sub_lst_sz];
            ufnd_idx_tmp=new int[sub_lst_sz];
        }
        
        memcpy(ufnd,sub_lst,sub_lst_sz*sizeof(int));
        memcpy(ufnd_idx,sub_lst_idx,sub_lst_sz*sizeof(int));
        
        
        int imother_lst,iufnd;
        for(int ip=0;ip<tot_p;ip++)
        {
            if(ip==my_p)
            {
                mother_lst_sz=lst_sz;
                memcpy(mother_lst,lst,mother_lst_sz*sizeof(int));
            }
            MPI_Bcast(&mother_lst_sz,1,MPI_INT,ip,world);
            MPI_Bcast(mother_lst,mother_lst_sz,MPI_INT,ip,world);
            
            if(ip==my_p || ufnd_sz==0 || mother_lst_sz==0) continue;
            
            ufnd_tmp_sz=0;
            for(imother_lst=0,iufnd=0;imother_lst<mother_lst_sz && iufnd<ufnd_sz;)
            {
                if(mother_lst[imother_lst]<ufnd[iufnd])
                    imother_lst++;
                else if(mother_lst[imother_lst]>ufnd[iufnd])
                {
                    ufnd_idx_tmp[ufnd_tmp_sz]=ufnd_idx[iufnd];
                    ufnd_tmp[ufnd_tmp_sz++]=ufnd[iufnd++];
                }
                else
                {
                    fnd_idx[fnd_sz++]=ufnd_idx[iufnd++];
                    imother_lst++;
                    comm_size[ip]++;
                }
            }
            
            for(;iufnd<ufnd_sz;)
            {
                ufnd_idx_tmp[ufnd_tmp_sz]=ufnd_idx[iufnd];
                ufnd_tmp[ufnd_tmp_sz++]=ufnd[iufnd++];
            }
            
            std::swap(ufnd_tmp,ufnd);
            std::swap(ufnd_idx_tmp,ufnd_idx);
            ufnd_sz=ufnd_tmp_sz;
            
        }
        
        delete [] sub_lst_idx;
        sub_lst_idx=fnd_idx;
        
        delete [] ufnd;
        delete [] ufnd_tmp;
        delete [] ufnd_idx;
        delete [] ufnd_idx_tmp;
        delete [] mother_lst;
    };
    
    class act
    {
    public:
        virtual void eq(int)=0;
        virtual void neq(int)=0;
    };
    
    class act_arch: public act
    {
        /*
         for iact to manage: it manages the arch_id
         the arch ids are catagorized into 2 groups:
         0. ids that are found here and we need the map for their new place: orig_inds,natms_to_kep
         1. ids that are not found here and therefore must be sent to
         another proc: ids_to_snd,inds_to_snd,natms_to_snd
         */
    public:
        int*& idx_lst_arch;
        int*& id_lst_snd;
        int*& idx_lst_snd;
        
        int& nkep;
        int& nsnd;
        
        int* id_lst_arch;
        
        act_arch(
        int*& _id_lst_snd,
        int*& _idx_lst_snd,
        int*& _idx_lst_arch,
        int& _nkep,
        int& _nsnd,
        int* _id_lst_arch):
        id_lst_snd(_id_lst_snd),
        idx_lst_snd(_idx_lst_snd),
        idx_lst_arch(_idx_lst_arch),
        nkep(_nkep),
        nsnd(_nsnd),
        id_lst_arch(_id_lst_arch)
        {}
        void eq(int idx)
        {
            idx_lst_arch[nkep++]=idx;
        }
        void neq(int idx)
        {
            idx_lst_snd[nsnd]=idx;
            id_lst_snd[nsnd++]=id_lst_arch[idx];
        }
    };
    
    class act_curr: public act
    {
        /*
         for jact to manage: it manages the orig_id
         the current ids are catagorized into 2 groups:
         0. ids that are found here and we need the map for their new place: arch_inds,natms_to_kep
         1. ids that are not found here and therefore must be recieved from
         another proc: ids_to_rcv,inds_to_rcv,natms_to_rcv
         */
        
    public:
        int*& idx_lst_curr;
        int*& id_lst_rcv;
        int*& idx_lst_rcv;
        
        int& nkep;
        int& nrcv;
        
        int* id_lst_curr;
        
        act_curr(
        int*& _id_lst_rcv,
        int*& _idx_lst_rcv,
        int*& _idx_lst_curr,
        int& _nkep,
        int& _nrcv,
        int* _id_lst_curr):
        id_lst_rcv(_id_lst_rcv),
        idx_lst_rcv(_idx_lst_rcv),
        idx_lst_curr(_idx_lst_curr),
        nkep(_nkep),
        nrcv(_nrcv),
        id_lst_curr(_id_lst_curr)
        {}
        void eq(int idx)
        {
            idx_lst_curr[nkep++]=idx;
        }
        void neq(int idx)
        {
            idx_lst_rcv[nrcv]=idx;
            id_lst_rcv[nrcv++]=id_lst_curr[idx];
        }
    };
    
    // input
    int* id_lst_arch=id_arch->begin();
    int* id_lst_curr=id->begin();
    int natms_arch=id_arch->vec_sz;
    int natms_curr=id->vec_sz;
    int byte_sz=0;
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        byte_sz+=vec_list->arch_vecs[ivec]->byte_sz;
    
    //output
    int* idx_lst_arch=NULL;
    int* idx_lst_snd=NULL;
    int* id_lst_snd=NULL;
    
    int* idx_lst_curr=NULL;
    int* idx_lst_rcv=NULL;
    int* id_lst_rcv=NULL;
    
    if(natms_arch)
    {
        idx_lst_arch=new int[natms_arch];
        idx_lst_snd=new int[natms_arch];
        id_lst_snd=new int[natms_arch];
    }
    
    if(natms_curr)
    {
        idx_lst_curr=new int[natms_curr];
        idx_lst_rcv=new int[natms_curr];
        id_lst_rcv=new int[natms_curr];
    }
    
    int nkep,nkep_arch=0,nkep_curr=0,nrcv=0,nsnd=0;
    
    act_arch* _act_arch_=new act_arch(id_lst_snd,idx_lst_snd,idx_lst_arch,nkep_arch,nsnd,id_lst_arch);
    act_curr* _act_curr_=new act_curr(id_lst_rcv,idx_lst_rcv,idx_lst_curr,nkep_curr,nrcv,id_lst_curr);
    act* _act_arch=_act_arch_;
    act* _act_curr=_act_curr_;
    XMath::srch_lst_lst(id_lst_arch,natms_arch,_act_arch,id_lst_curr,natms_curr,_act_curr);
    delete _act_curr_;
    delete _act_arch_;
    
    nkep=nkep_curr;
    
    shrnk_2_fit(idx_lst_arch,nkep_arch,natms_arch);
    shrnk_2_fit(idx_lst_snd,nsnd,natms_arch);
    shrnk_2_fit(id_lst_snd,nsnd,natms_arch);
    
    shrnk_2_fit(idx_lst_curr,nkep_curr,natms_curr);
    shrnk_2_fit(idx_lst_rcv,nrcv,natms_curr);
    shrnk_2_fit(id_lst_rcv,nrcv,natms_curr);
    
    
    int* nsnd_comm=new int[tot_p];
    int* nrcv_comm=new int[tot_p];
    sort_idx_by_p(id_lst_rcv,nrcv,id_lst_snd,nsnd,idx_lst_snd,nsnd_comm);
    sort_idx_by_p(id_lst_snd,nsnd,id_lst_rcv,nrcv,idx_lst_rcv,nrcv_comm);
    delete [] id_lst_snd;
    delete [] id_lst_rcv;
    
    byte** snd_buff=new byte*[tot_p];
    *snd_buff=NULL;
    if (nsnd)
    {
        *snd_buff=new byte[byte_sz*nsnd];
        for(int ip=1;ip<tot_p;ip++)
            snd_buff[ip]=snd_buff[ip-1]+byte_sz*nsnd_comm[ip-1];
    }
    
    byte** rcv_buff=new byte*[tot_p];
    *rcv_buff=NULL;
    if(nrcv)
    {
        *rcv_buff=new byte[byte_sz*nrcv];
        for(int ip=1;ip<tot_p;ip++)
            rcv_buff[ip]=rcv_buff[ip-1]+byte_sz*nrcv_comm[ip-1];
    }
    
    byte* _snd_buff=*snd_buff;
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        for(int i=0;i<nsnd;i++)
            vec_list->arch_vecs[ivec]->cpy(_snd_buff,idx_lst_snd[i]);
    
    delete [] idx_lst_snd;
    
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        vec_list->arch_vecs[ivec]->rearrange(idx_lst_arch,idx_lst_curr,nkep,natms_curr);
    
    delete [] idx_lst_arch;
    delete [] idx_lst_curr;
    
    for(int idisp=1;idisp<tot_p;idisp++)
    {
        int rcv_p=my_p-idisp;
        if(rcv_p<0) rcv_p+=tot_p;
        int snd_p=my_p+idisp;
        if(snd_p>=tot_p) snd_p-=tot_p;
        MPI_Sendrecv(snd_buff[snd_p],nsnd_comm[snd_p]*byte_sz,MPI_BYTE,snd_p,0,
                     rcv_buff[rcv_p],nrcv_comm[rcv_p]*byte_sz,MPI_BYTE,rcv_p,0,
                     world,MPI_STATUS_IGNORE);
    }
    
    delete [] *snd_buff;
    delete [] snd_buff;
    delete [] nsnd_comm;
    delete [] nrcv_comm;
    
    byte* _rcv_buff=*rcv_buff;
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        for(int i=0;i<nrcv;i++)
            vec_list->arch_vecs[ivec]->cpy(_rcv_buff,idx_lst_rcv[i]);
    
    delete [] idx_lst_rcv;
    delete [] *rcv_buff;
    delete [] rcv_buff;
    
    
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        add_vec(vec_list->arch_vecs[ivec]);
    add_vec(id_arch);
    delete id_arch;
    id_arch=NULL;
}



#endif
