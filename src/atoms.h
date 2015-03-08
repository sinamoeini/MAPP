#ifndef __MAPP__atoms__
#define __MAPP__atoms__

//swap list growth size
#define SWAPGROWTH 10

#include <iostream>
#include <exception>
#include <typeinfo>
#include <cmath>
#include <mpi.h>
#include "init.h"
#include "memory.h"
/*---------------------------------------------------------------------------
   _____   _          __      ___   _____   _       _____   _____
  /  ___/ | |        / /     /   | |  _  \ | |     /  ___/ |_   _|
  | |___  | |  __   / /     / /| | | |_| | | |     | |___    | |
  \___  \ | | /  | / /     / / | | |  ___/ | |     \___  \   | |
   ___| | | |/   |/ /     / /  | | | |     | |___   ___| |   | |
  /_____/ |___/|___/     /_/   |_| |_|     |_____| /_____/   |_|
 
 int** snd_list: snd_list[iswap][iatm]
 list of the atoms to send in each swap
 
 int* snd_list_capacity: snd_list_capacity[iswap]
 maximum size of the list in a swap
 
 int* snd_size: snd_size[iswap]
 number of atoms to send in each swap
 
 int* rcv_size: rcv[size]
 number of atoms to receive in each swap

 int* snd_p: snd_p[iswap]
 the processor id to send to in each swap
 
 int* rcv_p: rcv_p[iswap]
 the processor id to receive from in each swap
 
 int* pbc_correction: pbc_correction[iswap]
 the correction to x or s due to periodic boundary condition in each swap
 pbc_correction[iswap]=0:
 there is no cerrection 
 
 pbc_correction[iswap]=idim+1:
 s[iatm*x_dim+idim]++;
 for(int jdim=0;jdim<idim+1;jdim++)
    x[iatm*x_dim+jdim]+=H[idim][jdim]
 
 pbc_correction[iswap]=-(idim+1):
 s[iatm*x_dim+idim]--;
 for(int jdim=0;jdim<idim+1;jdim++)
    x[iatm*x_dim+jdim]-=H[idim][jdim]
 
 
 setup_ph(int,class VecLst*) assigns the following:
 snd_list[iswap]
 snd_size[iswap]
 snd_p[iswap]
 rcv_p[iswap]
 pbc_correction[iswap]
 
 xchng_ph(int,int,int*,int,class VecLst*) assigns the following:
 rcv_size[iswap]
 
 ---------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class SwapLst : protected InitPtrs
    {
    private:
        int* snd_list_capacity;        
        char* del_buff;
        int del_buff_capacity;
        
    protected:
    public:
        int** snd_list;
        
        int* snd_size;
        int* rcv_size;
        
        int* snd_p;
        int* rcv_p;
        
        int* pbc_correction;
        
        int grow_size;
        int no_swaps;
        
        SwapLst(MAPP*);
        SwapLst(MAPP*,int);
        ~SwapLst();
        
        void newlist(int);
        void add(int&,int&);
        void reset();
        void reset(int);
        
        int rectify(char*,int*);
        
    };
}
/*---------------------------------------------------------------------------
      ___   _     _   _____   _____
     /   | | |   / / | ____| /  ___|
    / /| | | |  / /  | |__   | |
   / / | | | | / /   |  __|  | |
  / /  | | | |/ /    | |___  | |___
 /_/   |_| |___/     |_____| \_____|
 
 
 pack(char*,int&,int):
 pack just 1 atom, not used
 
 pack(char*,int&,int*,int):
 pack multiple atoms, not used by
 atoms
 
 pack_n_mv(char*,int&,int):
 pack just 1 atom, and replace it
 with the last owned atom
 
 pack_n_mv(char*,int&,int*,int):
 pack multiple atoms, and replace the
 with the last owned atoms, not used
 by atoms
 
 unpack(char*,int&):
 unpack just 1 atom after last own atom
 
 unpack(char*,int&,int):
 unpack just multiple atoms after last
 own atom, not used by atoms
 
 unpack(char*,int&,int,int):
 unpack just multiple atoms after last
 own atom, but the buffer has stride.
 to clarify: values related to this
 vector have distances between (stride)
 them in the buffer
 
 unpack_ph(char*,int&):
 unpack just 1 atom after last phantom
 atom, not used  by atoms
 
 unpack_ph(char*,int&,int):
 unpack just multiple atoms after last
 phantom atom
 
 ---------------------------------------------------------------------------*/

namespace MAPP_NS {
    class Avec:protected InitPtrs
    {
    private:
    protected:
    public:
        int byte_size;
        int ph;
        int dim;
        char* name;
        char* print_format;
        
        Avec(MAPP* mapp):InitPtrs(mapp)
        {}
        
        template<typename T>
        void ret(T*& v)
        {v=static_cast<T*>(ret());}
        template<typename T>
        void ret_dump(T*& v)
        {v=static_cast<T*>(ret_dump());}
        
        
        virtual void* ret()=0;
        virtual void* ret(int)=0;
        virtual void* ret_dump()=0;
        
        virtual void move(int,int)=0;
        virtual void copy(int,int)=0;
        virtual void copy2bottom(int*,int)=0;
        
        virtual void gather_dump()=0;
        virtual void del_dump()=0;
        
        virtual void change_dimension(int)=0;

        virtual void pack(char*,int&,int)=0;
        virtual void pack(char*,int&,int*,int)=0;
        virtual void pack_n_mv(char*,int&,int)=0;
        virtual void pack_n_mv(char*,int&,int*,int)=0;
    
        virtual void unpack(char*,int&)=0;
        virtual void unpack(char*,int&,int)=0;
        virtual void unpack(char*,int&,int,int)=0;
        virtual void unpack_ph(char*,int&)=0;
        virtual void unpack_ph(char*,int&,int)=0;
        
        virtual void del(int*,int)=0;
        virtual void del(char*)=0;
        
        
        virtual void resize(int)=0;
        virtual void grow(int)=0;
        
        virtual void print_dump(FILE*,int)=0;
    };
}
/*---------------------------------------------------------------------------
      ___   _     _   _____   _____
     /   | | |   / / | ____| /  ___|
    / /| | | |  / /  | |__   | |
   / / | | | | / /   |  __|  | |
  / /  | | | |/ /    | |___  | |___   _____
 /_/   |_| |___/     |_____| \_____| |_____|
 ---------------------------------------------------------------------------*/
namespace MAPP_NS
{
    
    template<typename T>
    class AVec_: public Avec
    {
    private:
        T* vec;
        T* dump_vec;
    protected:
    public:
        AVec_(MAPP*,int,int,const char*);
        ~AVec_();
        
        void* ret();
        void* ret(int);
        void* ret_dump();
        
        void gather_dump();
        void del_dump();
        
        void change_dimension(int);
        
        void move(int,int);
        void copy(int,int);
        void copy2bottom(int*,int);
        
        void pack(char*,int&,int);
        void pack(char*,int&,int*,int);
        void pack_n_mv(char*,int&,int);
        void pack_n_mv(char*,int&,int*,int);
        
        void unpack(char*,int&);
        void unpack(char*,int&,int);
        void unpack(char*,int&,int,int);
        void unpack_ph(char*,int&);
        void unpack_ph(char*,int&,int);
        
        void resize(int);
        void grow(int);
        void print_dump(FILE*,int);
        
        void del(int*,int);
        void del(char*);
    };
}
/*---------------------------------------------------------------------------
  _     _   _____   _____   _       _____   _____
 | |   / / | ____| /  ___| | |     /  ___/ |_   _|
 | |  / /  | |__   | |     | |     | |___    | |
 | | / /   |  __|  | |     | |     \___  \   | |
 | |/ /    | |___  | |___  | |___   ___| |   | |
 |___/     |_____| \_____| |_____| /_____/   |_|
 ---------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class VecLst: protected InitPtrs
    {
    private:
    protected:
    public:
        int* vec_list;
        int no_vecs;
        int byte_size;
        
        int* ph_vec_list;
        int ph_no_vecs;
        int ph_byte_size;
        
        int* update_every_ph_vec_list;
        int update_every_ph_no_vecs;
        int update_every_ph_byte_size;
        
        VecLst(MAPP*,int,...);
        VecLst(MAPP*,int*,int);
        ~VecLst();
        void add_update(int);
        void del_update(int);
        
    };
    
}

/*---------------------------------------------------------------------------
      ___   _____   _____       ___  ___   _____
     /   | |_   _| /  _  \     /   |/   | /  ___/
    / /| |   | |   | | | |    / /|   /| | | |___
   / / | |   | |   | | | |   / / |__/ | | \___  \
  / /  | |   | |   | |_| |  / /       | |  ___| |
 /_/   |_|   |_|   \_____/ /_/        |_| /_____/
 
 ---------------------------------------------------------------------------*/
namespace MAPP_NS
{
    enum {COMM_MODE_0,COMM_MODE_1,COMM_MODE_2,COMM_MODE_3,COMM_MODE_4,COMM_MODE_5};
    class Atoms:protected InitPtrs
    {
    private:
        // snd and rcv buffers
        char* snd_buff_0;
        int snd_buff_0_capacity;
        char* snd_buff_1;
        int snd_buff_1_capacity;
        char* rcv_buff;
        int rcv_buff_capacity;
        
        char* snd_ph_buff;
        int snd_ph_buff_capacity;
        char* rcv_ph_buff;
        int rcv_ph_buff_capacity;
        
        int* old2new;
        int old2new_capacity;
        
        // to manage the size of above buffers
        template <typename TYPE>
        inline void buff_size_management(TYPE*&,int&,int);
  
        void grow(int,int,class VecLst*);
        
        int** comm_need;
        SwapLst* ph_lst;
        
        void xchng_cmplt(int,class VecLst*);
        void xchng_prtl(int,class VecLst*);
        
        void xchng_ph(int,class VecLst*);
        void xchng_ph(int,int* vec_list,int,int);
        
        // for xchng_cmplt(int,class VecLst*) & xchng_prtl(int,class VecLst*)
        inline void pack_cmplt_prtl(char*&,int&,int&,int,class VecLst*);
        // for xchng_cmplt(int,class VecLst*)
        inline void unpack_cmplt(char*&,int,class VecLst*);
        // for xchng_prtl(int,class VecLst*)
        inline void unpack_prtl(char*&,int,class VecLst*);
        
        // for xchng_ph(int*,int,int)
        inline void unpack_ph(char*&,class VecLst*,int);
        inline void pack_ph(char*&,class VecLst*,int*,int);
        
        // for update_ph(int*,int,int)
        inline void unpack_ph(char*&,int*,int,int);
        inline void pack_ph(char*&,int*,int,int*,int);

        
        inline void setup_comm_need();
               

        void setup_ph_n_neighbor(int,class VecLst*);
        
        void pre_pre_setup_ph(int);
        
        void pre_setup_ph();
        
        void setup_ph_0(class VecLst*);
        void setup_ph_1(class VecLst*);

        void post_setup_ph_0(class VecLst*);
        void post_setup_ph_1(class VecLst*);


        
    protected:
    public:
        int no_neigh_list_created;
        int comm_mode;
        /*-----------------------------------------------*/
        /* begining of atomic vectors related properties */
        // number of my processor owned atoms
        int natms;
        // number of my processor phantom atoms
        int natms_ph;
        //total number of owned atoms (sum of all processors natm)
        int tot_natms;

        // the maximim available size of non-phantom vectors (atms unit)
        // natms<=avec_size ALL THE TIME
        int avec_size;
        // the maximim available size of phantom vectors (atms unit)
        // natms+natms_ph<=avec_ph_size ALL THE TIME
        int avec_ph_size;
        /* end of atomic vectors related properties */
        
        // pointer to all the atomic vectors
        Avec** vectors;
        // number of atomic vectors
        int no_vecs;
        /* end of atomic vectors related properties */
        /*------------------------------------------*/
        
        
        /*----------------------------*/
        /* begining of box properties */
        //box dimension
        int dimension;
        // H matrix
        type0** H;
        // Hinv matrix
        type0** B;
        
        // these two define the domain under control of my processor
        // own atoms belong to this domain
        // lower bond of the local domain (s unit)
        type0* s_lo;
        // higher bond of the local domain (s unit)
        type0* s_hi;
        
        // these two define the domain under control of my processor + the domain that phanthom atoms of my processor should be inside
        // own atoms + phantom atoms belong to this domain
        // lower bond of the local domain + local phantom domain (s unit)
        type0* s_ph_lo;
        // higher bond of the local domain + local phantom domain (s unit)
        type0* s_ph_hi;
        
        
        int no_types;
        int type_n;
        type0* s_bound;
        type0** s_ph_lo_type;
        type0** s_ph_hi_type;
        type0** s_ph_lo_type_orig;
        type0** s_ph_hi_type_orig;

        // skin size which is the same as the skin size in neghborlist
        type0 skin;
        // thickness of extra layer to define the local phantom domain (length unit)
        // tot_cut_ph = skin + maximum cutoff
        type0 max_cut;
        // thickness of extra layer to define the local phantom domain (s unit)
        type0* max_cut_s;
        /* end of box box properties */
        /*---------------------------*/
        
        
        /*------------------------------------------*/
        /* begining of processor related properties */

        // total number of processors
        int tot_p;
        // my processor id
        int my_p_no;
        // the location of my processor in the grid (domain) for each dimension
        int* my_loc;
        // the number of my neighbor processors
        // 1st input: dimension
        // 2nd input: direction; 0 for behind me, 1 for infront of me
        int** neigh_p;
        // total number of processor in each dimension
        int* tot_p_grid;
        
        // total number of nodes
        int tot_n;
        // my node id
        int my_n_no;
        // number of processors per each node
        // 1st input: node number
        int* p_per_n;
        // the processor ids in each node
        // 1st input: node number
        // 2nd input: ith processor in that node
        // 0 <= i < p_per_n[node]
        int** n_p_grid;
        

        /* end of processor related properties */
        /*-------------------------------------*/
        
        Atoms(MAPP*);
        ~Atoms();
        void chng_dim(int);
        
        /*---------------------------------------------------------*/
        /* beginning of atomic vector related management functions */
        // add a new atomic vector
        // 1st input: 0 if the vector does not contain phantom atoms, 1 if it does
        // 2nd input: dimension of the vector
        // 3rd input: name of the vector
        template <typename TYPE>
        int add(int,int,const char*);
        // remove a vector by its number
        void del(int);
        // find an atomic vector by its name
        // if it cannot find it, it will produce an error and abort the execution
        int find(const char* name);
        // find an atomic vector by its name
        // if it cannot find it, it will return -1
        int find_exist(const char* name);
        /* end of atomic vector related management functions */
        /*---------------------------------------------------*/
        
       
        /*--------------------------------------*/
        /* beginning of grid related funcctions */
        void auto_grid_proc();
        // change the grid from command line
        void man_grid_proc(int*);
        /* end of grid related funcctions */
        /*--------------------------------*/
        
        void unpack_read(char*&,int,class VecLst*);
        
        void x2s(int);
        void s2x(int);

        // change the skin size from command line
        void chng_skin(type0);
        
        void restart();

        // init before a run
        void init(class VecLst*);

        /*
        update and make necessary changes
        after x is changed
         */
        void update(int,class VecLst*);
        
        void fin();
        
        /*
         exchange atoms untill everybody
         has their own atoms
         */
        void xchng_cmplt(class VecLst*);
        /*
         exchange atoms only once
         */
        void xchng_prtl(class VecLst*);
        /*
         update the phantom atoms of a 
         list of atomic vectors
         */
        void update_ph(int*,int,int);
        /*
         update the phantom atoms of an
         atomic vector
         */
        void update_ph(int);
    };
}

using namespace MAPP_NS;
/*---------------------------------------------------------------------------
      ___   _____   _____       ___  ___   _____
     /   | |_   _| /  _  \     /   |/   | /  ___/
    / /| |   | |   | | | |    / /|   /| | | |___
   / / | |   | |   | | | |   / / |__/ | | \___  \
  / /  | |   | |   | |_| |  / /       | |  ___| |
 /_/   |_|   |_|   \_____/ /_/        |_| /_____/

---------------------------------------------------------------------------*/
template<typename T>
int Atoms::add(int ph,int dim,const char* name)
{
    if(no_vecs==0)
    {
        if(typeid(T)!=typeid(type0))
        {
            error->abort("zeroth atomic vector shoud be of type %s"
                         ,typeid(type0).name());
        }
        if (ph!=1)
            error->abort("zeroth atomic vector shoud be phantom");
        if(dim<dimension)
            error->abort("dimension of zeroth atomic vector should "
                         "be more or equal to box dimension");
    }
    
    
    Avec** new_vectors=new Avec*[no_vecs+1];
    
    new_vectors[no_vecs]=new AVec_<T>(mapp,ph,dim,name);
    
    for(int ivec=0;ivec<no_vecs;ivec++)
        new_vectors[ivec]=vectors[ivec];
    if(no_vecs)
        delete [] vectors;
    vectors=new_vectors;
    no_vecs++;
    
    if(no_vecs==1)
    {
        int lngth=static_cast<int>(strlen(name))+1;
        char* x0name;
        CREATE1D(x0name,lngth+2);
        sprintf(x0name,"%s_0",name);
        add<T>(0,dim,x0name);
        delete [] x0name;
        return (no_vecs-2);
    }
    
    return no_vecs-1;
}
/*--------------------------------------------
 buffer size management
 --------------------------------------------*/
template<typename T>
inline void Atoms::buff_size_management(T*& buff
,int& curr_size,int size_needed)
{
    if(size_needed<curr_size)
        return;
    if(curr_size)
        delete [] buff;
    CREATE1D(buff,size_needed);
    curr_size=size_needed;
}
/*---------------------------------------------------------------------------
      ___   _     _   _____   _____
     /   | | |   / / | ____| /  ___|
    / /| | | |  / /  | |__   | |
   / / | | | | / /   |  __|  | |
  / /  | | | |/ /    | |___  | |___   _____
 /_/   |_| |___/     |_____| \_____| |_____|
 ---------------------------------------------------------------------------*/
/*--------------------------------------------
 constructor
 --------------------------------------------*/
template<typename T>
AVec_<T>::AVec_(MAPP* mapp,int phantom
,int dimension,const char* vec_name):Avec(mapp)
{
    int lngth=static_cast<int>(strlen(vec_name))+1;
    CREATE1D(name,lngth);
    memcpy(name,vec_name,lngth*sizeof(char));
    
    dim=dimension;
    ph=phantom;
    byte_size=sizeof(T)*dim;
    
    
    if(typeid(T)==typeid(char))
    {
        lngth=static_cast<int>(strlen((char*)"%c "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%c ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(unsigned char))
    {
        lngth=static_cast<int>(strlen((char*)"%c "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%c ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(short int))
    {
        lngth=static_cast<int>(strlen((char*)"%d "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%d ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(unsigned short int))
    {
        lngth=static_cast<int>(strlen((char*)"%d "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%d ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(int))
    {
        lngth=static_cast<int>(strlen((char*)"%d "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%d ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(unsigned int))
    {
        lngth=static_cast<int>(strlen((char*)"%d "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%d ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(long int))
    {
        lngth=static_cast<int>(strlen((char*)"%ld "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%ld ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(unsigned long int))
    {
        lngth=static_cast<int>(strlen((char*)"%ld "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%ld ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(long long int))
    {
        lngth=static_cast<int>(strlen((char*)"%lld "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%lld ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(unsigned long long int))
    {
        lngth=static_cast<int>(strlen((char*)"%lld "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%lld ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(float))
    {
        lngth=static_cast<int>(strlen((char*)"%f "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%f ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(double))
    {
        lngth=static_cast<int>(strlen((char*)"%lf "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%lf ",lngth*sizeof(char));
    }
    else if(typeid(T)==typeid(long double))
    {
        lngth=static_cast<int>(strlen((char*)"%Lf "))+1;
        CREATE1D(print_format,lngth);
        memcpy(print_format,(char*)"%Lf ",lngth*sizeof(char));
    }
    else
    {
        error->abort("this type is not supported for atomic vectors");
    }
    
    if(ph)
        CREATE1D(vec,dim*atoms->avec_ph_size);
    else
        CREATE1D(vec,dim*atoms->avec_size);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<typename T>
AVec_<T>::~AVec_()
{
    
    if(ph)
    {
        if(atoms->avec_ph_size)
            delete [] vec;
    }
    else
    {
        if(atoms->avec_size)
            delete [] vec;
    }
    
    delete [] print_format;
    delete [] name;
}
/*--------------------------------------------
 return the vector
 --------------------------------------------*/
template<typename T>
void* AVec_<T>::ret()
{
    return vec;
}
/*--------------------------------------------
 return the pointer to iatom
 --------------------------------------------*/
template<typename T>
void* AVec_<T>::ret(int iatm)
{
    return &vec[iatm*dim];
}
/*--------------------------------------------
 return the vector
 --------------------------------------------*/
template<typename T>
void* AVec_<T>::ret_dump()
{
    if(atoms->my_p_no==0)
        return dump_vec;
    else
        return NULL;
}
/*--------------------------------------------
 gather dump
 --------------------------------------------*/
template<typename T>
void AVec_<T>::gather_dump()
{
    int my_p_no=atoms->my_p_no;
    int tot_p=atoms->tot_p;
    int tot_natms=atoms->tot_natms;
    int natms=atoms->natms;
    int* rcv_size=NULL;
    
    if(my_p_no!=0)
        MPI_Send(&natms,1,MPI_INT,0,my_p_no,world);
    
    if(my_p_no==0)
    {
        CREATE1D(rcv_size,tot_p);
        MPI_Status status;
        for(int iproc=1;iproc<tot_p;iproc++)
            MPI_Recv(&rcv_size[iproc],1
            ,MPI_INT,iproc,iproc,world,&status);
    }
    
    if(my_p_no==0)
    {
        CREATE1D(dump_vec,dim*tot_natms);
        memcpy(dump_vec,vec,natms*byte_size);
    }
    
    if(my_p_no!=0)
        MPI_Send(vec,natms*byte_size,MPI_BYTE,0,my_p_no,world);
    
    if(my_p_no==0)
    {
        MPI_Status status;
        int tot_atoms=natms;
        for(int iproc=1;iproc<tot_p;iproc++)
        {
            MPI_Recv(&dump_vec[tot_atoms*dim]
            ,(byte_size*rcv_size[iproc]),MPI_BYTE
            ,iproc,iproc,world,&status);
            tot_atoms+=rcv_size[iproc];
        }
    }
    
    if(my_p_no==0)
    {
        delete [] rcv_size;
    }
    
    
}
/*--------------------------------------------
 delete dump
 --------------------------------------------*/
template<typename T>
void AVec_<T>::change_dimension(int d)
{
    if(d==dim)
        return;
    
    int d_min=MIN(d,dim);
    
    int tot;
    int max_size;
    if(ph==0)
    {
        tot=atoms->natms;
        max_size=atoms->avec_size;
    }
    else
    {
        tot=atoms->natms+atoms->natms_ph;
        max_size=atoms->avec_ph_size;
    }
    
    T* new_vec;
    CREATE1D(new_vec,max_size*tot);
    
    for(int i=0;i<tot;i++)
        for(int j=0;j<d_min;j++)
            new_vec[i*d+j]=vec[i*dim+j];
    
    if(max_size)
        delete [] vec;
    
    vec=new_vec;
    dim=d;
    byte_size=sizeof(T)*dim;
    
}
/*--------------------------------------------
 delete dump
 --------------------------------------------*/
template<typename T>
void AVec_<T>::del_dump()
{
    if(atoms->my_p_no==0 && atoms->tot_natms)
        delete [] dump_vec;
}
/*--------------------------------------------
 copy from jatm to iatm
 1st input: destination
 2nd input: source
 --------------------------------------------*/
template<typename T>
void AVec_<T>::copy(int iatm,int jatm)
{
    memcpy(&vec[iatm*dim],&vec[jatm*dim],byte_size);
}
/*--------------------------------------------
 copy from jatm to iatm
 1st input: destination
 2nd input: source
 --------------------------------------------*/
template<typename T>
void AVec_<T>::copy2bottom(int* lst,int no)
{
    
    no--;
    while(no>-1)
    {
        memcpy(&vec[dim*(atoms->natms+atoms->natms_ph+no)]
        ,&vec[dim*lst[no]],byte_size);
        no--;
    }
    
}
/*--------------------------------------------
 move from iatm to jatm
 --------------------------------------------*/
template<typename T>
void AVec_<T>::move(int iatm,int jatm)
{
    memcpy(&vec[jatm*dim],&vec[iatm*dim],byte_size);
}
/*--------------------------------------------
 pack a list without removing it from the
 atomic vector
 --------------------------------------------*/
template<typename T>
void AVec_<T>::pack(char* buff,int& buff_pos,int iatm)
{
    memcpy(&buff[buff_pos],&vec[iatm*dim],byte_size);
    buff_pos+=byte_size;
}
/*--------------------------------------------
 pack a list and remove them from the atomic
 vector
 --------------------------------------------*/
template<typename T>
void AVec_<T>::pack_n_mv(char* buff,int& buff_pos,int iatm)
{
    memcpy(&buff[buff_pos],&vec[iatm*dim],byte_size);
    memcpy(&vec[iatm*dim],&vec[(atoms->natms-1)*dim],byte_size);
    buff_pos+=byte_size;
}
/*--------------------------------------------
 unpack just one own atom
 --------------------------------------------*/
template<typename T>
void AVec_<T>::unpack(char* buff,int& buff_pos)
{
    memcpy(&vec[atoms->natms*dim],&buff[buff_pos],byte_size);
    buff_pos+=byte_size;
}
/*--------------------------------------------
 unpack just one phantom atoms
 --------------------------------------------*/
template<typename T>
void AVec_<T>::unpack_ph(char* buff,int& buff_pos)
{
    memcpy(&vec[(atoms->natms+atoms->natms_ph)*dim],&buff[buff_pos],byte_size);
    buff_pos+=byte_size;
}
/*----------------------------------------------------------------------------------------*/
/*--------------------------------------------
 pack a list without removing it from the
 atomic vector
 pack(buff,buff_pos,atm_list,atm_list_size)
 used by pack_ph pack_ph(char*&,class VecLst*
 ,int*,int)
 --------------------------------------------*/
template<typename T>
void AVec_<T>::pack(char* buff,int& buff_pos,int* lst,int no)
{
    for(int i=0;i<no;i++)
    {
        memcpy(&buff[buff_pos],&vec[lst[i]*dim],byte_size);
        buff_pos+=byte_size;
    }
}
/*--------------------------------------------
 pack a list and remove them from the atomic
 vector
 --------------------------------------------*/
template<typename T>
void AVec_<T>::pack_n_mv(char* buff,int& buff_pos,int* lst,int no)
{
    int my_natms=atoms->natms;
    for(int i=0;i<no;i++)
    {
        memcpy(&buff[buff_pos],&vec[lst[i]*dim],byte_size);
        memcpy(&vec[lst[i]*dim],&vec[(my_natms-1)*dim],byte_size);
        my_natms--;
        buff_pos+=byte_size;
    }
}
/*--------------------------------------------
 unpack the own atoms
 --------------------------------------------*/
template<typename T>
void AVec_<T>::unpack(char* buff,int& buff_pos,int xtra_natms)
{
    memcpy(&vec[atoms->natms*dim],&buff[buff_pos],byte_size*xtra_natms);
    buff_pos+=byte_size*xtra_natms;
}
/*--------------------------------------------
 unpack the phantom atoms
 used by unpack_ph(char*&,class VecLst*,int)
 --------------------------------------------*/
template<typename T>
void AVec_<T>::unpack_ph(char* buff,int& buff_pos,int xtra_natms)
{
    memcpy(&vec[(atoms->natms+atoms->natms_ph)*dim],&buff[buff_pos],byte_size*xtra_natms);
    buff_pos+=byte_size*xtra_natms;
}
/*--------------------------------------------
 resize the vector
 --------------------------------------------*/
template<typename T>
void AVec_<T>::resize(int xtra_natms)
{
    
    if(ph)
    {
        if(atoms->avec_ph_size)
            delete [] vec;
        xtra_natms+=atoms->avec_ph_size;
    }
    else
    {
        xtra_natms+=atoms->avec_size;
        if(atoms->avec_size)
            delete [] vec;
    }
    
    CREATE1D(vec,dim*xtra_natms);

}
/*--------------------------------------------
 grow the vector
 --------------------------------------------*/
template<typename T>
void AVec_<T>::grow(int xtra_natms)
{
    
    if(ph)
    {
        xtra_natms+=atoms->avec_ph_size;
        GROW(vec,dim*atoms->avec_ph_size,dim*xtra_natms);
    }
    else
    {
        xtra_natms+=atoms->avec_size;
        GROW(vec,dim*atoms->avec_size,dim*xtra_natms);
    }
    
}
/*--------------------------------------------
 print the dump values in a file
 --------------------------------------------*/
template<typename T>
void AVec_<T>::print_dump(FILE* fp,int iatm)
{
    for(int i=0;i<dim;i++)
        fprintf(fp,print_format,dump_vec[iatm*dim+i]);
}
/*--------------------------------------------
 unpack the own atoms with stride
 --------------------------------------------*/
template<typename T>
void AVec_<T>::unpack(char* buff,int& buff_pos,int stride,int xtra_natms)
{
    for(int i=0;i<xtra_natms;i++)
    {
        memcpy(&vec[(atoms->natms+i)*dim],&buff[buff_pos],byte_size);
        buff_pos+=stride;
    }
    buff_pos+=byte_size-stride*xtra_natms;
}
/*--------------------------------------------
 del a list of atoms
 --------------------------------------------*/
template<typename T>
void AVec_<T>::del(int* list,int list_size)
{
    if(list_size)
    {
        int start=list[0];
        int size;
        for(int i=0;i<list_size;i++)
        {
            if(i==list_size-1)
                size=atoms->natms+atoms->natms_ph-list[i]-1;
            else
                size=list[i+1]-list[i]-1;
            
            if(size)
            {
                memcpy(&vec[start*dim],&vec[(list[i]+1)*dim],size*byte_size);
            }
            start+=size;
        }

    }
}
/*--------------------------------------------
 del a list of atoms
 --------------------------------------------*/
template<typename T>
void AVec_<T>::del(char* mark)
{
    if(atoms->natms_ph)
    {
        int icurs,strt,istrt,isize;
        icurs=0;
        while (mark[icurs]=='1'&&icurs<atoms->natms_ph)
            icurs++;
        strt=icurs+atoms->natms;
        
        while(icurs<atoms->natms_ph)
        {
            while (mark[icurs]=='0'&&icurs<atoms->natms_ph)
                icurs++;
            istrt=icurs;
            
            while (mark[icurs]=='1'&&icurs<atoms->natms_ph)
                icurs++;
            isize=icurs-istrt;
            istrt+=atoms->natms;
            
            if(isize)
                memcpy(&vec[strt*dim],&vec[istrt*dim],isize*byte_size);
            strt+=isize;
        }
    }
}

#endif
