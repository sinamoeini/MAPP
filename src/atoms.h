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
    template<const int> class Communincation;
    template<const int> class Exchange;
    template<const int> class Update;
    class Atoms
    {
    private:
        
        bool box_chng;
        bool grid_established;
        
        
        type0 skin;
        
        
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
        type0 max_cut;
        VecLst* vec_list;
        unsigned long xchng_id;
        
        Communincation<__dim__>* comm;
        Exchange<__dim__>* xchng;
        Update<__dim__>* updt;
        
        vec** vecs;
        int nvecs;
        void add_vec(vec*);
        vec* find_vec(const char*);
        void del_vec(vec*);

        Vec<type0>* x;
        Vec<int>* id;
        
        type0 H[__dim__][__dim__];
        type0 B[__dim__][__dim__];
        
        int natms;
        int natms_ph;
        int tot_natms;

        int tot_p;
        int my_p;
        
        
        Atoms();
        ~Atoms();


        void set_skin(type0);
        type0 get_skin();
        
        type0 max_cut_s[__dim__];
        type0 s_lo[__dim__];
        type0 s_hi[__dim__];
        
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
    class Communincation
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
        
        int tot_p_grid[_dimension_];
        /*
         * tot_p_grid[idm]: number of processors in idim dimension
         * allocated in: constructor
         * deallocated in: destructor
         * set in: auto_grid() / man_grid ()
         */
        
        int my_loc[_dimension_];
        /*
         * the location of my proc. in proc. grid (domain)
         * my_loc[idim] my location in dimension idim
         * allocated in: constructor
         * deallocated in: destructor
         * set in: auto_grid() / man_grid ()
         */
        
        int neigh_p[_dimension_][2];
        /*
         * the number of my neighbor processors
         * 1st input: dimension
         * 2nd input: direction; 0 for behind me, 1 for infront of me
         * neigh_p[idim][idir]: id of my neighbor proc. in my idim and idir direction
         * allocated in: constructor
         * deallocated in: destructor
         * set in: auto_grid() / man_grid ()
         */
        
        type0 (&s_lo)[__dim__];
        type0 (&s_hi)[__dim__];
        
        Communincation(Atoms*);
        ~Communincation();
        
        void auto_grid(type0(&)[_dimension_][_dimension_]);
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
    class Update
    {
    private:
        /*things that reference cannot be removed*/
        int& natms;
        int& natms_ph;
        
        /*things that reference can be removed*/
        type0 (&H)[_dimension_][_dimension_];
        type0 (&B)[_dimension_][_dimension_];
        type0 (&max_cut_s)[_dimension_];
        int& my_p;
        int (&neigh_p)[_dimension_][2];
        type0 (&s_lo)[_dimension_];
        type0 (&s_hi)[_dimension_];
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
        
        
        vec** updt_vecs;
        int nupdt_vecs;
        int tot_updt_vecs_sz;
        
        void add_to_snd_lst(int&,int&);
        void reserve_snd_buff(int);
        void reserve_rcv_buff(int);
        
        class LoadUnLoadUpdate;
        class LoadUnLoadUpdateComm;
        class LoadUnLoadUpdateSelfComm;
        LoadUnLoadUpdate* comm_manager[_dimension_];
        
    protected:
    public:
        Update(Atoms*,vec**,int);
        ~Update();
        void reset();
        void update(vec*,bool);
        void update(vec**,int,bool);
        void list();
        void rm_rdndncy();
    };

    template<const int _dimension_>
    class Update<_dimension_>::LoadUnLoadUpdate
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
    class Update<_dimension_>::LoadUnLoadUpdateComm : public Update<_dimension_>::LoadUnLoadUpdate
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
        
        vec**& updt_vecs;
        int& nupdt_vecs;
        int& tot_updt_vecs_sz;
#ifndef BLOCK_COMM
        MPI_Request request[2];
#endif
    protected:
    public:
        LoadUnLoadUpdateComm(Update*);
        void load(int&);
        void unload(int&,int&,int&);
        void load_unload(int&,int&,int&);
        void update_mult(int&,int&,int&,vec**&,int&,int&);
        void update_sing(int&,int&,int&,vec*&);
        void xchng_buff(int&,int&,byte*&,int&,int&,byte*&);
    };
    
    template<const int _dimension_>
    class Update<_dimension_>::LoadUnLoadUpdateSelfComm : public Update<_dimension_>::LoadUnLoadUpdate
    {
    private:
        int**& snd_atms_lst;
        int*& snd_atms_lst_sz;
        
        int*& rcv_atms_lst_sz;
        
        vec**& updt_vecs;
        int& nupdt_vecs;
    protected:
    public:
        LoadUnLoadUpdateSelfComm(Update*);
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
    class Exchange
    {
    private:
        
        /*things that reference cannot be removed*/
        int& natms;
        unsigned long& xchng_id;
        
        /*things that reference cannot be removed*/
        int& my_p;
        int (&neigh_p)[_dimension_][2];
        type0 (&s_lo)[_dimension_];
        type0 (&s_hi)[_dimension_];
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
        Exchange(Atoms*,vec**,int);
        ~Exchange();
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
    Vec<dmd_type>* map_dmd_type
    =dynamic_cast<Vec<dmd_type>*>(map);
    if(map_dmd_type)
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
Update<_dimension_>::Update(Atoms* atoms,vec** vecs_,int nvecs_):
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

updt_vecs(vecs_),
nupdt_vecs(nvecs_)
{

    snd_buff=NULL;
    snd_buff_cpcty=0;
    snd_buff_grw=1024;
    
    rcv_buff=NULL;
    rcv_buff_cpcty=0;
    rcv_buff_grw=1024;
    snd_atms_lst_grw=8;
    
    tot_updt_vecs_sz=0;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        tot_updt_vecs_sz+=updt_vecs[ivec]->byte_sz;
    
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
Update<_dimension_>::~Update()
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
void Update<_dimension_>::reset()
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
inline void Update<_dimension_>::add_to_snd_lst(int& icomm,int& iatm)
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
inline void Update<_dimension_>::reserve_rcv_buff(int xtra)
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
inline void Update<_dimension_>::reserve_snd_buff(int xtra)
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
inline void Update<_dimension_>::update(vec** updt_vecs,int nupdt_vecs,bool x_xst)
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
inline void Update<_dimension_>::update(vec* updt_vec,bool x_xst)
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
inline void Update<_dimension_>::list()
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
inline void Update<_dimension_>::rm_rdndncy()
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
    
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        updt_vecs[ivec]->vec_sz=vec_sz;
        updt_vecs[ivec]->cpy_pst(list_,list_sz);
    }
    
    delete [] list;
    
    natms_ph=new_natms_ph;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
Update<_dimension_>::LoadUnLoadUpdateComm::
LoadUnLoadUpdateComm(Update* updt):
LoadUnLoadUpdate(),

snd_atms_lst(updt->snd_atms_lst),
snd_atms_lst_sz(updt->snd_atms_lst_sz),

rcv_atms_lst_sz(updt->rcv_atms_lst_sz),

snd_buff(updt->snd_buff),
snd_buff_sz(updt->snd_buff_sz),
snd_buff_cpcty(updt->snd_buff_cpcty),
snd_buff_grw(updt->snd_buff_grw),

rcv_buff(updt->rcv_buff),
rcv_buff_sz(updt->rcv_buff_sz),
rcv_buff_cpcty(updt->rcv_buff_cpcty),
rcv_buff_grw(updt->rcv_buff_grw),

updt_vecs(updt->updt_vecs),
nupdt_vecs(updt->nupdt_vecs),
tot_updt_vecs_sz(updt->tot_updt_vecs_sz)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Update<_dimension_>::LoadUnLoadUpdateComm::load
(int& iatm)
{
    if(snd_buff_cpcty<snd_buff_sz+tot_updt_vecs_sz)
    {
        byte* tmp_snd_buff=new byte[snd_buff_sz+tot_updt_vecs_sz+snd_buff_grw];
        memcpy(tmp_snd_buff,snd_buff,snd_buff_sz);
        delete [] snd_buff;
        snd_buff=tmp_snd_buff;
        snd_buff_cpcty=snd_buff_sz+tot_updt_vecs_sz+snd_buff_grw;
    }
    
    byte* tmp_snd_buff=snd_buff+snd_buff_sz;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        updt_vecs[ivec]->cpy(tmp_snd_buff,iatm);
    snd_buff_sz+=tot_updt_vecs_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Update<_dimension_>::LoadUnLoadUpdateComm::unload
(int& icomm,int& snd_p,int& rcv_p)
{
    
    MPI_Sendrecv(&snd_atms_lst_sz[icomm],1,MPI_INT,snd_p,0,
                 &rcv_atms_lst_sz[icomm],1,MPI_INT,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
    
    rcv_buff_sz=rcv_atms_lst_sz[icomm]*tot_updt_vecs_sz;
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=rcv_atms_lst_sz[icomm]*tot_updt_vecs_sz+rcv_buff_grw;
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
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        updt_vecs[ivec]->pst(tmp_rcv_buff0,tot_updt_vecs_sz,rcv_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Update<_dimension_>::LoadUnLoadUpdateComm::load_unload
(int& icomm,int& snd_p,int& rcv_p)
{
    snd_buff_sz=snd_atms_lst_sz[icomm]*tot_updt_vecs_sz;
    if(snd_buff_cpcty<snd_buff_sz)
    {
        delete [] snd_buff;
        snd_buff=new byte[snd_buff_sz+snd_buff_grw];
        snd_buff_cpcty=snd_buff_sz+snd_buff_grw;
    }
    
    byte* tmp_snd_buff=snd_buff;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        updt_vecs[ivec]->cpy(tmp_snd_buff,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
    
    MPI_Sendrecv(&snd_atms_lst_sz[icomm],1,MPI_INT,snd_p,0,
                 &rcv_atms_lst_sz[icomm],1,MPI_INT,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
    
    rcv_buff_sz=rcv_atms_lst_sz[icomm]*tot_updt_vecs_sz;
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=rcv_atms_lst_sz[icomm]*tot_updt_vecs_sz+rcv_buff_grw;
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
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        updt_vecs[ivec]->reserve(rcv_atms_lst_sz[icomm]);
        updt_vecs[ivec]->pst(tmp_rcv_buff,rcv_atms_lst_sz[icomm]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Update<_dimension_>::LoadUnLoadUpdateComm::update_mult
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
inline void Update<_dimension_>::LoadUnLoadUpdateComm::update_sing
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
inline void Update<_dimension_>::LoadUnLoadUpdateComm::xchng_buff
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
Update<_dimension_>::LoadUnLoadUpdateSelfComm::
LoadUnLoadUpdateSelfComm(Update* updt):
LoadUnLoadUpdate(),

snd_atms_lst(updt->snd_atms_lst),
snd_atms_lst_sz(updt->snd_atms_lst_sz),

rcv_atms_lst_sz(updt->rcv_atms_lst_sz),

updt_vecs(updt->updt_vecs),
nupdt_vecs(updt->nupdt_vecs)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Update<_dimension_>::LoadUnLoadUpdateSelfComm::
load(int& iatm)
{
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        updt_vecs[ivec]->cpy_pst(iatm);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Update<_dimension_>::LoadUnLoadUpdateSelfComm::unload
(int& icomm,int&,int&)
{
    rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Update<_dimension_>::LoadUnLoadUpdateSelfComm::load_unload
(int& icomm,int&,int&)
{
    rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm];
    
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        updt_vecs[ivec]->reserve(rcv_atms_lst_sz[icomm]);
        updt_vecs[ivec]->cpy_pst(snd_atms_lst[icomm],rcv_atms_lst_sz[icomm]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Update<_dimension_>::LoadUnLoadUpdateSelfComm::update_mult
(int& icomm,int&,int&,vec**& vecs
,int& nvecs,int&)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->cpy_pst(snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Update<_dimension_>::LoadUnLoadUpdateSelfComm::update_sing
(int& icomm,int&,int&,vec*& v)
{
    v->cpy_pst(snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Update<_dimension_>::LoadUnLoadUpdateSelfComm::xchng_buff
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
Exchange<_dimension_>::Exchange(Atoms* atoms
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
Exchange<_dimension_>::~Exchange()
{

    delete [] snd_buff[0];
    delete [] snd_buff[1];
    delete [] rcv_buff;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int _dimension_>
inline void Exchange<_dimension_>::load(int& iatm,int idir)
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
inline void Exchange<_dimension_>::load(byte*& buff,int& idir)
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
int Exchange<_dimension_>::xchng_buff(int idim,int idir)
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
void Exchange<_dimension_>::full_xchng()
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
Communincation<_dimension_>::Communincation(Atoms* atoms):
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
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<const int _dimension_>
Communincation<_dimension_>::~Communincation()
{
    delete [] *n_p_grid;
    delete [] n_p_grid;
    delete [] p_per_n;
}
/*--------------------------------------------
 automatically grid the domain
 --------------------------------------------*/
template<const int _dimension_>
void Communincation<_dimension_>::auto_grid(type0 (&H)[_dimension_][_dimension_])
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
void Communincation<_dimension_>::man_grid(int* n)
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




#endif
