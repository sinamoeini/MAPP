#ifndef __MAPP__vec__
#define __MAPP__vec__
#include "atoms.h"
namespace MAPP_NS
{
    template<typename T>
    class VecTens
    {
    private:
        MPI_Datatype MPI_T;
        void assign_mpi_type();
    protected:
    public:
        Vec<T>* vec;
        T (*A)[__dim__];
        int* pnatms;
        int x_dim;
        int alloc_flag;
        
        bool box_chng;
        
        VecTens();
        VecTens(VecTens&&);
        ~VecTens();
        
        
        void init(Atoms*,bool);
        void init(Atoms*,Vec<T>*&,bool);
        void init(Atoms*,Vec<T>*& ,T(&)[__dim__][__dim__],bool);
        void fin();
        
    
        void add(T,const VecTens&);
        Vec<T>* operator()();
        T operator*(const VecTens&);
        VecTens& operator*= (T&);
        VecTens& operator=(const VecTens&);
        VecTens& operator=(VecTens&&);
        VecTens& operator-=(const VecTens&);
        VecTens& operator+=(const VecTens&);
        
    };
}

using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
VecTens<T>::VecTens():
alloc_flag(0)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
VecTens<T>::VecTens(VecTens&& other):
vec(other.vec),
box_chng(other.box_chng),
A(other.A),
alloc_flag(0)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
VecTens<T>::~VecTens()
{
    {
        if(alloc_flag==3)
        {
            delete vec;
            if(!box_chng) return;
            delete [] A;
            A=NULL;
        }
        else if(alloc_flag==2)
        {
            if(!box_chng) return;
            delete [] A;
            A=NULL;
        }
        else if(alloc_flag==1)
        {
            delete vec;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void VecTens<T>::assign_mpi_type()
{
    if(typeid(T)==typeid(type0))
        MPI_T=MPI_TYPE0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void VecTens<T>::init(Atoms* atoms,bool box_chng_)
{
    pnatms=&atoms->natms;
    x_dim=atoms->x->dim;
    box_chng=box_chng_;
    assign_mpi_type();
    alloc_flag=3;
    vec=new Vec<T>(atoms,x_dim);
    if(!box_chng) return;
    A=new T[__dim__][__dim__];
    for(int idim=0;idim<__dim__;idim++)
        for(int jdim=0;jdim<__dim__;jdim++)
            A[idim][jdim]=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void VecTens<T>::init(Atoms* atoms,Vec<T>*& v,bool box_chng_)
{
    pnatms=&atoms->natms;
    x_dim=atoms->x->dim;
    box_chng=box_chng_;
    assign_mpi_type();
    alloc_flag=2;
    vec=v;
    if(!box_chng) return;
    A=new T[__dim__][__dim__];
    for(int idim=0;idim<__dim__;idim++)
        for(int jdim=0;jdim<__dim__;jdim++)
            A[idim][jdim]=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void VecTens<T>::init(Atoms* atoms,Vec<T>*& v,T (&A_)[__dim__][__dim__],bool box_chng_)
{
    pnatms=&atoms->natms;
    x_dim=atoms->x->dim;
    box_chng=box_chng_;
    assign_mpi_type();
    
    alloc_flag=0;
    vec=v;
    if(!box_chng) return;
    A=A_;
}

/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void VecTens<T>::fin()
{
    if(alloc_flag==3)
    {
        delete vec;
        alloc_flag=0;
        if(!box_chng) return;
        delete [] A;
        A=NULL;
    }
    else if(alloc_flag==2)
    {
        alloc_flag=0;
        if(!box_chng) return;
        delete [] A;
        A=NULL;
    }
    else if(alloc_flag==1)
    {
        delete vec;
    }
    alloc_flag=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
void VecTens<T>::add(T a,const VecTens& other)
{
    T* vec0=this->vec->begin();
    T* vec1=other.vec->begin();
    for(int i=0;i<*pnatms*x_dim;i++)
        vec0[i]+=a*vec1[i];
    
    if(!box_chng) return;
    
    vec0=*this->A;
    vec1=*other.A;
    for(int i=0;i<__dim__*__dim__;i++)
        vec0[i]+=a*vec1[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
Vec<T>* VecTens<T>::operator()()
{
    return vec;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
T VecTens<T>::operator*(const VecTens<T> &rhs)
{
    T ans_lcl=0.0,ans;
    T* vec0=this->vec->begin();
    T* vec1=rhs.vec->begin();
    for(int i=0;i<*pnatms*x_dim;i++)
        ans_lcl+=vec0[i]*vec1[i];
    
    MPI_Allreduce(&ans_lcl,&ans,1,MPI_T,MPI_SUM,__world__);
    if(!box_chng) return ans;
    vec0=*this->A;
    vec1=*rhs.A;
    for(int i=0;i<__dim__*__dim__;i++)
        ans+=vec0[i]*vec1[i];
    return ans;
}

/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
VecTens<T>& VecTens<T>::operator*=(T& a)
{
    T* vec0=this->vec->begin();
    for(int i=0;i<*pnatms*x_dim;i++)
        vec0[i]*=a;
    
    if(!box_chng) return *this;
    
    vec0=*this->A;
    for(int i=0;i<__dim__*__dim__;i++)
        vec0[i]*=a;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
VecTens<T>& VecTens<T>::operator=(const VecTens& other)
{
    memcpy(this->vec->begin(),other.vec->begin(),*pnatms*x_dim*sizeof(T));
    if(!box_chng) return *this;
    memcpy(*this->A,*other.A,__dim__*__dim__*sizeof(T));
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
VecTens<T>& VecTens<T>::operator=(VecTens&& other)
{
    this->vec=other.vec;
    this->box_chng=other.box_chng;
    this->A=other.A;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
VecTens<T>& VecTens<T>::operator-=(const VecTens& rhs)
{
    T* vec0=this->vec->begin();
    T* vec1=rhs.vec->begin();
    for(int i=0;i<*pnatms*x_dim;i++)
        vec0[i]-=vec1[i];
    
    if(!box_chng) return *this;
    
    vec0=&(this->A[0][0]);
    vec1=&(rhs.A[0][0]);
    vec1=*rhs.A;
    for(int i=0;i<__dim__*__dim__;i++)
        vec0[i]-=vec1[i];
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T>
VecTens<T>& VecTens<T>::operator+=(const VecTens& rhs)
{
    T* vec0=this->vec->begin();
    T* vec1=rhs.vec->begin();
    for(int i=0;i<*pnatms*x_dim;i++)
        vec0[i]+=vec1[i];
    
    if(!box_chng) return *this;
    
    vec0=&(this->A[0][0]);
    vec1=&(rhs.A[0][0]);
    vec1=*rhs.A;
    for(int i=0;i<__dim__*__dim__;i++)
        vec0[i]+=vec1[i];
    
    return *this;
}



#endif
