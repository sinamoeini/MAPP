
#ifndef vec_h
#define vec_h
#include "atoms.h"
namespace MAPP_NS
{
    template<typename T>
    class VecTens
    {
    private:
    protected:
    public:
        Vec<T>* vec;
        T** A;
        int* pnatms;
        int x_dim;
        MPI_Comm* world;
        int dim;
        int alloc_flag;
        MPI_Datatype MPI_T;
        bool box_chng;
        
        VecTens(){alloc_flag=0;}
        
        ~VecTens()
        {
            if(alloc_flag==3)
            {
                delete vec;
                if(!box_chng) return;
                delete [] *A;
                delete [] A;
                A=NULL;
            }
            else if(alloc_flag==2)
            {
                if(!box_chng) return;
                delete [] *A;
                delete [] A;
                A=NULL;
            }
            else if(alloc_flag==1)
            {
                delete vec;
            }
        }
        
        void assign_mpi_type()
        {
            if(typeid(T)==typeid(type0))
                MPI_T=MPI_TYPE0;
        }
        
        void init(Atoms* atoms,bool bc)
        {
            pnatms=&atoms->natms;
            x_dim=atoms->x->dim;
            dim=atoms->dimension;
            world=&atoms->world;
            box_chng=bc;
            assign_mpi_type();
            alloc_flag=3;
            vec=new Vec<T>(atoms,x_dim);
            if(!box_chng) return;
            A=new T*[dim];
            *A=new T[dim*dim];
            for(int idim=1;idim<dim;idim++)
                A[idim]=A[idim-1]+dim;
        }
        void init(Atoms* atoms,Vec<T>*& v,bool bc)
        {
            pnatms=&atoms->natms;
            x_dim=atoms->x->dim;
            dim=atoms->dimension;
            world=&atoms->world;
            box_chng=bc;
            assign_mpi_type();
            alloc_flag=2;
            vec=v;
            if(!box_chng) return;
            A=new T*[dim];
            *A=new T[dim*dim];
            for(int idim=1;idim<dim;idim++)
                A[idim]=A[idim-1]+dim;
        }
        void init(Atoms* atoms,Vec<T>*& v,T**& A_,bool bc)
        {
            pnatms=&atoms->natms;
            x_dim=atoms->x->dim;
            dim=atoms->dimension;
            world=&atoms->world;
            box_chng=bc;
            assign_mpi_type();
            
            alloc_flag=0;
            vec=v;
            if(!box_chng) return;
            A=A_;
        }
        
        void fin()
        {
            if(alloc_flag==3)
            {
                delete vec;
                alloc_flag=0;
                if(!box_chng) return;
                delete [] *A;
                delete [] A;
                A=NULL;
            }
            else if(alloc_flag==2)
            {
                alloc_flag=0;
                if(!box_chng) return;
                delete [] *A;
                delete [] A;
                A=NULL;
            }
            else if(alloc_flag==1)
            {
                delete vec;
            }
            alloc_flag=0;
        }
        
        T operator*(const VecTens& rhs)
        {
            T ans_lcl=0.0,ans;
            T* vec0=this->vec->begin();
            T* vec1=rhs.vec->begin();
            for(int i=0;i<*pnatms*x_dim;i++)
                ans_lcl+=vec0[i]*vec1[i];
            
            MPI_Allreduce(&ans_lcl,&ans,1,MPI_T,MPI_SUM,*world);
            if(!box_chng) return ans;
            vec0=*this->A;
            vec1=*rhs.A;
            for(int i=0;i<dim*dim;i++)
                ans+=vec0[i]*vec1[i];
            return ans;
        }
        
        void add(T a,const VecTens& other)
        {
            T* vec0=this->vec->begin();
            T* vec1=other.vec->begin();
            for(int i=0;i<*pnatms*x_dim;i++)
                vec0[i]+=a*vec1[i];
            
            if(!box_chng) return;
            
            vec0=*this->A;
            vec1=*other.A;
            for(int i=0;i<dim*dim;i++)
                vec0[i]+=a*vec1[i];
        }

        /*********************************************************************************************************/
        
        Vec<T>* operator()()
        {return vec;}
        
        
        VecTens& operator*= (T& a)
        {
            T* vec0=this->vec->begin();
            for(int i=0;i<*pnatms*x_dim;i++)
                vec0[i]*=a;
            
            if(!box_chng) return *this;
            
            vec0=*this->A;
            for(int i=0;i<dim*dim;i++)
                vec0[i]*=a;
            return *this;
        }
        
        friend void operator*(T& a,VecTens& rhs)
        {
            T* vec0=rhs.vec->begin();
            for(int i=0;i<*rhs.pnatms*rhs.x_dim;i++)
                vec0[i]*=a;
            vec0=*rhs.A;
            for(int i=0;i<rhs.dim*rhs.dim;i++)
                vec0[i]*=a;
        }
        
        VecTens(const VecTens& other)
        {
            memcpy(vec->begin(),other.vec->begin(),*pnatms*x_dim*sizeof(T));
            if(!box_chng) return;
            memcpy(*A,*other.A,dim*dim*sizeof(T));
        }
        
        VecTens& operator=(const VecTens& other)
        {
            memcpy(this->vec->begin(),other.vec->begin(),*pnatms*x_dim*sizeof(T));
            if(!box_chng) return *this;
            memcpy(*this->A,*other.A,dim*dim*sizeof(T));
            return *this;
        }
        
        VecTens(VecTens&& other)
        {
            this->vec=other.vec;
            this->box_chng=other.box_chng;
            this->A=other.A;
            this->alloc_flag=0;
        }
        
        VecTens& operator =(VecTens&& other)
        {
            this->vec=other.vec;
            this->box_chng=other.box_chng;
            this->A=other.A;
            return *this;
        }
        
        VecTens& operator-(const VecTens& rhs)
        {
            T* vec0=this->vec->begin();
            T* vec1=rhs.vec->begin();
            T* my_vec=vec->begin();
            for(int i=0;i<*pnatms*x_dim;i++)
                my_vec[i]=vec0[i]-vec1[i];
            
            if(!box_chng) return *this;
            
            vec0=*this->A;
            vec1=*rhs.A;
            my_vec=*A;
            vec1=*rhs.A;
            for(int i=0;i<dim*dim;i++)
                my_vec[i]=vec0[i]-vec1[i];
            return *this;
        }
        
        void operator-=(const VecTens& rhs)
        {
            T* vec0=this->vec->begin();
            T* vec1=rhs.vec->begin();
            T* my_vec=vec->begin();
            for(int i=0;i<*pnatms*x_dim;i++)
                my_vec[i]=vec0[i]-vec1[i];
            
            if(!box_chng) return;
            
            vec0=*this->A;
            vec1=*rhs.A;
            my_vec=*A;
            vec1=*rhs.A;
            for(int i=0;i<dim*dim;i++)
                my_vec[i]=vec0[i]-vec1[i];

        }
        
        VecTens& operator+(const VecTens& rhs)
        {
            T* vec0=this->vec->begin();
            T* vec1=rhs.vec->begin();
            T* my_vec=vec->begin();
            for(int i=0;i<*pnatms*x_dim;i++)
                my_vec[i]=vec0[i]+vec1[i];
            
            if(!box_chng) return *this;
            
            vec0=*this->A;
            vec1=*rhs.A;
            my_vec=*A;
            vec1=*rhs.A;
            for(int i=0;i<dim*dim;i++)
                my_vec[i]=vec0[i]+vec1[i];
            return *this;
        }
        
        void operator+=(const VecTens& rhs)
        {
            T* vec0=this->vec->begin();
            T* vec1=rhs.vec->begin();
            T* my_vec=vec->begin();
            for(int i=0;i<*pnatms*x_dim;i++)
                my_vec[i]=vec0[i]+vec1[i];
            
            if(!box_chng) return;
            
            vec0=*this->A;
            vec1=*rhs.A;
            my_vec=*A;
            vec1=*rhs.A;
            for(int i=0;i<dim*dim;i++)
                my_vec[i]=vec0[i]+vec1[i];
        }
    };
}

#endif /* vec_h */
