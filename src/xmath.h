#ifndef __MAPP__xmath__
#define __MAPP__xmath__
#include <utility> 
#include "type_def.h"
namespace MAPP_NS
{
    namespace XMath
    {
        // recursive form of above functions that feeds into it
        void fac_list_rec(int,int,int,int*&,int*&,int&);
        template<typename,class>
        class XVec;
        // inverse square matrix
        void invert(type0**,type0**,int);
        // inverse lower triangle square matrix
        template<const int d>
        void invert_lower_triangle(type0 (&)[d][d],type0 (&)[d][d]);
        void invert_lower_triangle(type0**&,type0**&,int&);
        // return the list of all possible groups of integers that their products are equal to specific number
        void fac_list(int,int,int*&,int&);
        
        void square2lo_tri(type0(&)[3][3],type0(&)[3][3]);
        void quadrature_lg(int,type0*,type0*);
        void quadrature_hg(int,type0*,type0*);
        int M3sqroot(type0(&)[3][3],type0(&)[3][3]);
        
        template<typename T0,class COMP,class SWAP>
        void quicksort(T0,T0,COMP,SWAP);
        template<typename T0,class C0>
        void srch_lst_lst(T0*,int,C0*,T0*,int,C0*);
    };

    template<typename T0,class C0>
    class XMath::XVec
    {
    private:
        T0* data;
        int* key;
        int size;
        C0& act;
    protected:
    public:
        
        class iterator
        {
        private:
            XVec* xvec;
            int* ikey;
        protected:
        public:
            iterator(XVec* _xvec,int* _ikey)
            {
                xvec=_xvec;
                ikey=_ikey;
            }
            
            iterator(const iterator& other)
            {
                this->xvec=other.xvec;
                this->ikey=other.ikey;
            }
            
            iterator(iterator&& other)
            {
                this->xvec=other.xvec;
                this->ikey=other.ikey;
            }
            
            iterator& operator =(const iterator& other)
            {
                this->xvec=other.xvec;
                this->ikey=other.ikey;
                return *this;
            }
            
            iterator& operator =(iterator&& other)
            {
                this->xvec=other.xvec;
                this->ikey=other.ikey;
                return *this;
            }
            
            iterator& operator +=(int i)
            {
                ikey+=i;
                return *this;
            }
            
            iterator& operator ++()
            {
                ikey++;
                return *this;
            }
            
            iterator operator ++(int)
            {
                iterator old=*this;
                ikey++;
                return old;
            }
            
            iterator& operator --()
            {
                ikey--;
                return *this;
            }
            
            iterator operator --(int)
            {
                iterator old=*this;
                ikey--;
                return old;
            }
            
            T0& operator *()
            {
                return xvec->data[*ikey];
            }
            
            T0& operator [](int i)
            {
                return xvec->data[*(ikey+i)];
            }
            
            bool operator!=(const iterator& rhs)
            {
                return (this->ikey!=rhs.ikey);
            }
            
            bool operator==(const iterator& rhs)
            {
                return (this->ikey==rhs.ikey);
            }
            void neq()
            {
                xvec->act.neq(*ikey);
            }
            void eq()
            {
                xvec->act.eq(*ikey);
            }
        };
        
        XVec(T0* _data
             ,int _size,C0& _act):
        data(_data),
        size(_size),
        act(_act)
        {
            key=NULL;
            if(!size) return;
            
            key=new int[size];
            for(int i=0;i<size;i++) key[i]=i;
            quicksort(key,key+size,
            [this](int* ikey,int* jkey){return (data[*ikey]<data[*jkey]);},
            [](int* ikey,int* jkey){std::swap(*ikey,*jkey);});
        }
        
        XVec(XVec&& other):
        act(other.act)
        {
            this->data=other.data;
            this->key=other.key;
            this->size=other.size;
            other.key=NULL;
        }
        ~XVec()
        {
            delete [] key;
        }
        
        XVec& operator=(XVec&& other)
        {
            this->act=other.act;
            this->data=other.data;
            this->key=other.key;
            this->size=other.size;
            other.key=NULL;
            return *this;
        }
        
        iterator begin()
        {
            return iterator(this,key);
        }
        
        iterator end()
        {
            return iterator(this,key+size);
        }
        
        iterator rbegin()
        {
            return iterator(this,key+size-1);
        }
        
        iterator rend()
        {
            return iterator(this,key-1);
        }
    };
}
namespace MAPP_NS
{
    namespace XMatrixVector
    {
        template<const int dim>
        class _zero_
        {
        public:
            template<typename T>
            static inline void func(T* v)
            {
                *v=0;
                _zero_<dim-1>::func(v+1);
            }
        };
        template<>
        class _zero_<0>
        {
        public:
            template<typename T>
            static inline void func(T*)
            {
            }
        };
        
        template<const int i>
        class _rsq_
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return _rsq_<i-1>::func(V0+1,V1+1)+(*V0-*V1)*(*V0-*V1);
            }
        };
        
        template<>
        class _rsq_<1>
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return (*V0-*V1)*(*V0-*V1);
            }
        };
        
        template<const int i>
        class _V_V_
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return _V_V_<i-1>::func(V0+1,V1+1)+(*V0)*(*V1);
            }
        };
        
        template<>
        class _V_V_<1>
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return (*V0)*(*V1);
            }
        };
        template<>
        class _V_V_<0>
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return 0;
            }
        };
        
        template<const int stride,const int i>
        class _V_V_str_
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return _V_V_str_<stride,i-1>::func(V0+stride,V1+1)+(*V0)*(*V1);
            }
        };
        
        template<const int stride>
        class _V_V_str_<stride,1>
        {
        public:
            template<typename T>
            static inline T func(T* V0,T* V1)
            {
                return (*V0)*(*V1);
            }
        };
        
        
        template<const int dim,const int i>
        class _V_Mlt_
        {
        public:
            template<typename T>
            static inline void func(T* Mlt,T* V)
            {
                *V=_V_V_str_<dim,i>::func(Mlt,V);
                _V_Mlt_<dim,i-1>::func(Mlt+dim+1,V+1);
            }
        };
        
        template<const int dim>
        class _V_Mlt_<dim,1>
        {
        public:
            template<typename T>
            static inline void func(T* Mlt,T* V)
            {
                *V*=*Mlt;
            }
        };
        
        
        template<const int dim,const int i>
        class _x2s_
        {
        public:
            template<typename T>
            static inline void func(T* Mlt,T* V)
            {
                *V=_V_V_str_<dim,i>::func(Mlt,V);
                while(*V<0.0)
                    (*V)++;
                while(*V>=1.0)
                    (*V)--;
                _x2s_<dim,i-1>::func(Mlt+dim+1,V+1);
            }
        };
        
        template<const int dim>
        class _x2s_<dim,1>
        {
        public:
            template<typename T>
            static inline void func(T* Mlt,T* V)
            {
                *V*=*Mlt;
                while(*V<0.0)
                    (*V)++;
                while(*V>=1.0)
                    (*V)--;
            }
        };
        
        
        template<const int dim,const int i,const int j>
        class _Mlt_inv_
        {
        public:
            template<typename T>
            static inline void func(T* M,T* Minv)
            {
                *Minv=-_V_V_str_<dim,i-j>::func(Minv-(i-j)*dim,M)*(*(Minv+i-j));
                _Mlt_inv_<dim,i,j-1>::func(Minv-1,M-1);
            }
        };
        
        template<const int dim,const int i>
        class _Mlt_inv_<dim,i,i>
        {
        public:
            template<typename T>
            static inline void func(T* M,T* Minv)
            {
                *Minv=1.0/(*M);
                _Mlt_inv_<dim,i,i-1>::func(Minv-1,M-1);
            }
        };
        template<const int dim,const int i>
        class _Mlt_inv_<dim,i,0>
        {
        public:
            template<typename T>
            static inline void func(T* M,T* Minv)
            {
                *Minv=-_V_V_str_<dim,i>::func(Minv-i*dim,M)*(*(Minv+i));
                _Mlt_inv_<dim,i+1,i+1>::func(Minv+dim+i+1,M+dim+i+1);
            }
        };
        
        template<const int dim>
        class _Mlt_inv_<dim,dim,dim>
        {
        public:
            template<typename T>
            static inline void func(T* M,T* Minv)
            {}
        };
        
        template<const int dim,const int i,const int j>
        class _M_2_Mlt_
        {
        public:
            template<typename T>
            static inline void func(T* M_old,T* M_new)
            {
                
                M_new[i*dim+j]=(_V_V_<dim>::func(M_old+j*dim,M_old+i*dim)-
                _V_V_<j>::func(M_new+j*dim,M_new+i*dim))/M_new[j*(dim+1)];

                _M_2_Mlt_<dim,i,j+1>::func(M_old,M_new);
            }
        };
        
        template<const int dim,const int i>
        class _M_2_Mlt_<dim,i,i>
        {
        public:
            template<typename T>
            static inline void func(T* M_old,T* M_new)
            {
                M_new[i*dim+i]=sqrt(_V_V_<dim>::func(M_old+i*dim,M_old+i*dim)-
                    _V_V_<i>::func(M_new+i*dim,M_new+i*dim));
                
                _zero_<dim-i-1>::func(M_new+i*(dim+1)+1);
                _M_2_Mlt_<dim,i+1,0>::func(M_old,M_new);
            }
        };
        
        template<const int dim>
        class _M_2_Mlt_<dim,dim,0>
        {
        public:
            template<typename T>
            static inline void func(T*,T*)
            {}
        };
        
        template<const int dim,const int i,const int j>
        class _Mlt_Mlt_
        {
        public:
            template<typename T>
            static inline void func(T* M0,T* M1,T* B)
            {
                *B=_V_V_str_<dim,i-j+1>::func(M1-(i-j)*dim,M0);
                _Mlt_Mlt_<dim,i,j+1>::func(M0+1,M1+1,B+1);
            }
        };
        
        template<const int dim,const int i>
        class _Mlt_Mlt_<dim,i,i>
        {
        public:
            template<typename T>
            static inline void func(T* M0,T* M1,T* B)
            {
                *B=(*M0)*(*M1);
                _Mlt_Mlt_<dim,i+1,0>::func(M0+dim-i,M1+dim-i,B+dim-i);
            }
        };
        template<const int dim>
        class _Mlt_Mlt_<dim,dim,0>
        {
        public:
            template<typename T>
            static inline void func(T*,T*,T*)
            {}
        };
        
        
        /*-------------------------------------------------------------------*/
        template <const int i>
        class UnrolledLoop
        {
        public:
            template < typename FuncType >
            static inline void Do(FuncType func)
            {
                UnrolledLoop <i-1>::Do( func);
                func(i-1);
            }
        };
        template <>
        class UnrolledLoop<1>
        {
        public:
            template < typename FuncType >
            static inline void Do(FuncType func){func(0);}
        };
        
        template<const int dim,typename T>
        inline T rsq(T* V0,T* V1)
        {
            return _rsq_<dim>::func(V0,V1);
        }
        
        template<const int dim,typename T>
        inline T V_V(T* V0,T* V1)
        {
            return _V_V_<dim>(V0,V1);
        }
        
        template<const int dim,typename T>
        inline void V_Mlt(T* V,T** Mlt)
        {
            _V_Mlt_<dim,dim>::func(*Mlt,V);
        }
        template<const int dim,typename T>
        static inline void V_Mlt(T* V,T (&Mlt)[dim][dim])
        {
            _V_Mlt_<dim,dim>::func((T*)Mlt,V);
        }
        
        template<const int dim,typename T>
        inline void s2x(T* s,T (&H)[dim][dim])
        {
            _V_Mlt_<dim,dim>::func(&H[0][0],s);
        }
        
        template<const int dim,typename T>
        inline void s2x(T* s,T** H)
        {
            _V_Mlt_<dim,dim>::func(*H,s);
        }
        
        template<const int dim,typename T>
        inline void x2s(T* x,T (&B)[dim][dim])
        {
            _x2s_<dim,dim>::func(&B[0][0],x);
        }
        
        template<const int dim,typename T>
        inline void x2s(T* x,T** B)
        {
            _x2s_<dim,dim>::func(*B,x);
        }
        
        template<const int dim,typename T>
        inline void Mlt_inv(T** A,T** Ainv)
        {
            **Ainv=1.0/(**A);
            _Mlt_inv_<dim,1,1>::func(*A+dim+1,*Ainv+dim+1);
        }
        
        template<const int dim,typename T>
        inline void Mlt_inv(T (&A)[dim][dim],T (&Ainv)[dim][dim])
        {
            **Ainv=1.0/(**A);
            _Mlt_inv_<dim,1,1>::func((T*)A+dim+1,(T*)Ainv+dim+1);
        }
        
        template<const int dim,typename T>
        inline void Mlt_Mlt(T** M0,T** M1,T** B)
        {
            _Mlt_Mlt_<dim,0,0>::func(*M0,*M1,*B);
        }
        
        template<const int dim,typename T>
        inline void Mlt_Mlt(T (&M0)[dim][dim],T (&M1)[dim][dim],T (&B)[dim][dim])
        {
            _Mlt_Mlt_<dim,0,0>::func((T*)M0,(T*)M1,(T*)B);
        }
        
        template<const int dim,typename T>
        inline void M_2_Mlt(T** &M0,T** M1)
        {
            _M_2_Mlt_<dim,0,0>::func(*M0,*M1);
        }
        
        template<const int dim,typename T>
        inline void M_2_Mlt(T (&M0)[dim][dim],T (&M1)[dim][dim])
        {
            _M_2_Mlt_<dim,0,0>::func((T*)M0,(T*)M1);
        }
    }
}
/*--------------------------------------------
 quick sort algorithm
 
 example to test:
 
 int list_sz=1000;
 int* list=new int[list_sz];
 
 Random* rand=new Random(this,58245647);
 for(int i=0;i<list_sz;++i)
     list[i]=static_cast<int>(rand->uniform()*10000.0);
 
 XMath* xmath=new XMath();
 
 auto comp=
 [] (int* i,int* j)->bool
 {
     return (*i<*j);
 };
 auto swap=
 [] (int* i,int* j)->void
 {
     if(i==j)
         return;
     int k=*i;
     *i=*j;
     *j=k;
 };
 
 xmath->quicksort(list,list+list_sz,comp,swap);
 
 for(int i=1;i<list_sz;i++)
     if(list[i]<list[i-1])
         cout<<"error"<<endl;
 
 delete xmath;
 delete rand;
 delete [] list;
 --------------------------------------------*/
using namespace MAPP_NS;
template<typename T0,class COMP,class SWAP>
void XMath::quicksort(T0 start,T0 end,COMP comp,SWAP swap)
{
    if(start+1>=end)
        return;
    T0 pindex=start;
    for(T0 i=start;i!=end-1;i++)
        if(comp(i,end-1))
        {
            swap(i,pindex);
            pindex++;
        }
    swap(end-1,pindex);
    quicksort(start,pindex,comp,swap);
    quicksort(pindex+1,end,comp,swap);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0,class C0>
void XMath::srch_lst_lst(T0* ilst,int isize,C0* iact
,T0* jlst,int jsize,C0* jact)
{
    XVec<T0,C0> _ilst(ilst,isize,*iact);
    XVec<T0,C0> _jlst(jlst,jsize,*jact);
    
    if(isize==0 || jsize==0)
    {
        if(isize)
        {
            auto ipos=_ilst.begin();
            auto iend=_ilst.end();
            while (ipos!=iend)
            {
                ipos.neq();
                ++ipos;
            }
        }
        
        if(jsize)
        {
            auto jpos=_jlst.begin();
            auto jend=_jlst.end();
            while (jpos!=jend)
            {
                jpos.neq();
                ++jpos;
            }
        }
        
        return;
    }
    
    auto ipos=_ilst.begin();
    auto jpos=_jlst.begin();
    auto iend=_ilst.end();
    auto jend=_jlst.end();
    
    auto exit=
    [&]()->void
    {
        while (ipos!=iend)
        {
            ipos.neq();
            ++ipos;
        }
        
        while (jpos!=jend)
        {
            jpos.neq();
            ++jpos;
        }
        return;
    };
    
    while(1)
    {
        while (*ipos!=*jpos)
        {
            while(*ipos<*jpos)
            {
                ipos.neq();
                ++ipos;
                if(ipos==iend)
                    return exit();
            }
            std::swap(ipos,jpos);
            std::swap(iend,jend);
        }
        while(*ipos==*jpos)
        {
            ipos.eq();
            jpos.eq();
            ++ipos;
            ++jpos;
            if(ipos==iend || jpos==jend)
                return exit();
        }
    }    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<const int dim>
void XMath::invert_lower_triangle(type0 (&A)[dim][dim],type0 (&Ainv)[dim][dim])
{
    type0 ATA[dim][dim];
    type0 c[dim];
    type0 x[dim];
    type0 g[dim];
    type0 g0[dim];
    type0 h[dim];
    type0 a0,a1,alpha;
    type0 g0g0,gg,gg0,ratio;
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            Ainv[i][j]=ATA[i][j]=0.0;
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            for(int k=0;k<dim;k++)
                ATA[i][j]+=A[k][i]*A[k][j];
    
    for(int itry=0;itry<dim;itry++)
    {
        for(int i=0;i<dim;i++)
        {
            c[i]=A[itry][i];
            x[i]=c[i];
        }
        
        
        g0g0=0.0;
        for(int i=0;i<dim;i++)
        {
            h[i]=2.0*c[i];
            for(int j=0;j<dim;j++)
                h[i]-=2.0*ATA[i][j]*x[j];
            g[i]=h[i];
            g0g0+=h[i]*h[i];
        }
        
        int jtry=0;
        double error=1.0;
        while(jtry<dim+1 && error!=0.0)
        {
            
            if(g0g0==0.0)
            {
                error=0.0;
                continue;
            }
            
            
            a0=0.0;
            a1=0.0;
            for(int i=0;i<dim;i++)
            {
                a0+=h[i]*g[i];
                for(int j=0;j<dim;j++)
                    a1+=h[i]*ATA[i][j]*h[j];
            }
            if(a1==0.0)
            {
                error=0.0;
                continue;
            }
            alpha=0.5*a0/a1;
            
            for(int i=0;i<dim;i++)
                x[i]+=alpha*h[i];
            
            //cout << "chk 3" << endl;
            
            gg=0.0;
            gg0=0.0;
            for(int i=0;i<dim;i++)
            {
                g[i]=2.0*c[i];
                for(int j=0;j<dim;j++)
                    g[i]-=2.0*ATA[i][j]*x[j];
                gg+=g[i]*g[i];
                gg0+=g0[i]*g[i];
            }
            
            //cout << "chk 4" << endl;
            ratio=(gg-gg0)/g0g0;
            g0g0=gg;
            
            
            for(int i=0;i<dim;i++)
            {
                h[i]=ratio*h[i]+g[i];
                g0[i]=g[i];
            }
            
            
            error=0.0;
            for(int i=0;i<dim;i++)
            {
                for(int j=0;j<dim;j++)
                    error+=x[i]*ATA[i][j]*x[i];
                error-=2*c[i]*x[i];
            }
            error++;
            
            jtry++;
        }
        
        for(int i=0;i<dim;i++)
            Ainv[i][itry]=x[i];
    }
}


#endif
