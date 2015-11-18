/*--------------------------------------------
 Created by Sina on 2/5/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__xmath__
#define __MAPP__xmath__
#include "init.h"

namespace MAPP_NS
{
    class XMath
    {
    private:
        // recursive form of above functions that feeds into it
        void fac_rec(int,int,int,int*&,int**&,int&);
    protected:
    public:
        XMath();
        ~XMath();
        // inverse square matrix
        void invert(type0**,type0**,int);
        // inverse lower triangle square matrix
        void invert_lower_triangle(type0**&,type0**&,int&);
        // return the list of all possible groups of integers that their products are equal to specific number
        int fac(int,int,int**&);
        void square2lo_tri(type0**,type0**);
        void quadrature_lg(int,type0*,type0*);
        void quadrature_hg(int,type0*,type0*);
        int M3sqroot(type0**,type0**);
        
        template<typename T0,class COMP,class SWAP>
        void quicksort(T0,T0,COMP,SWAP);

    };
}
/*--------------------------------------------
 quick sort algorithm
 
    int list_sz=1000;
    int* list=new int[list_sz];
    
    Random* rand=new Random(this,58245647);
    for(int i=0;i<list_sz;++i)
    {
        list[i]=static_cast<int>(rand->uniform()*10000.0);
    }
    
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




namespace MAPP_NS
{
    class SPARSE : protected InitPtrs
    {
    private:
    protected:
    public:
        
        SPARSE(MAPP*,int,int);
        ~SPARSE();
        int no_elem;
        int h0,h1;
        int* idx0;
        int* idx1;
        double* mtrx;
        void add(int,int,double);
        void add(int,int*,int*,double*);
        void sort(int);
        void vectorize(int);
        
        int no_vec;
        int* idx_vec;
        int* init_pos;
        int* fin_pos;
        
    };
    
    // Parallel sparse decomposed by rows (h0)
    class SPARSE_P : protected InitPtrs
    {
    private:
    protected:
    public:
        SPARSE_P(MAPP*,int,int);
        ~SPARSE_P();
        
        int no_elem;
        int lcl_no_elem;
        int h0,h1;
        int lcl_h0;
        int lo_h0;
        int hi_h0;
        int totp;
        int myno;
        int* lcl_idx0;
        int* lcl_idx1;
        double* lcl_mtrx;
        void add(int,int,double);
        void addd(int,int,double);
        void add(int,int*,int*,double*);
        void sort(int);
        void vectorize(int);
        
        
        int lcl_no_vec;
        int* lcl_idx_vec;
        int* lcl_init_pos;
        int* lcl_fin_pos;
    };
    
    
    class SOLVEAXb : protected InitPtrs
    {
    private:
        SPARSE* A;
        int lsize;
        int tot_lsize;
        int totp;
        int myno;
        int lcl_lo,lcl_hi;
        
        int* lcl_idx0;
        int* lcl_idx1;
        double* lcl_mtrx;
        int lcl_no_elem;
        
        int* comm_snd_size;
        int* comm_rcv_size;
        int** comm_snd_lst;
        double** buff_snd;
        double* c;
        
        double** g;
        double* h;
        double d_sq;
    protected:
    public:
        SOLVEAXb(MAPP*,SPARSE*,double*,int);
        ~SOLVEAXb();
        void xchng(double*);
        void construct_ans();
        int solve(double);
        double* x;
        double* ans;
    };
    
    class COMB : protected InitPtrs
    {
        COMB(MAPP*);
        ~COMB();
        void clean();
        void comb(int*,int);
        void comb_rec(int,int,int,int,int*&,int*& ,int*,int,int);
        int** perm_list;
        int perm_list_size_0;
        int perm_list_size_1;
    };

}
#endif 
