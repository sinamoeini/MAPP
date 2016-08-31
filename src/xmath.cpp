/*--------------------------------------------
 Created by Sina on 2/5/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#include <cmath>
#include "xmath.h"
#include <limits>
#define TOLERANCE 1.0e-10
#include "error.h"
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
XMath::XMath()
{
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
XMath::~XMath()
{
}
/*--------------------------------------------
 return the list of all possible groups of 
 integers that their products are equal to a 
 specific number
 --------------------------------------------*/
void XMath::fac_list(int no,int dim,int*& list,int& list_sz)
{
    list_sz=0;
    list=NULL;
    int* tmp_list=new int[dim];
    fac_list_rec(no,dim,0,tmp_list,list,list_sz);
    delete [] tmp_list;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void XMath::fac_list_rec(int no,int idim,int ipos,
int*& tmp_list,int*& list,int& list_sz)
{
    
    if(idim>1)
    {
        for(int i=1;i<=no;i++)
            if(no%i==0)
            {
                tmp_list[ipos]=i;
                fac_list_rec(no/i,idim-1,ipos+1,tmp_list,list,list_sz);
            }
    }
    else
    {
        tmp_list[ipos]=no;
        int* list_=new int[(list_sz+1)*(ipos+1)];
        memcpy(list_,list,list_sz*(ipos+1)*sizeof(int));
        memcpy(list_+(ipos+1)*list_sz,tmp_list,(ipos+1)*sizeof(int));
        delete [] list;
        list=list_;
        list_sz++;
    }
}
/*--------------------------------------------
 autogrid the domain
 --------------------------------------------*/
void XMath::square2lo_tri(type0** H_old
,type0** H_new)
{
    type0** Q=new type0*[3];
    for(int idim=0;idim<3;idim++)
        Q[idim]=new type0[3];

    type0 H0H0;
    type0 H0H1;
    H0H0=H0H1=0.0;
    for(int i=0;i<3;i++)
    {
        H0H0+=H_old[0][i]*H_old[0][i];
        H0H1+=H_old[0][i]*H_old[1][i];
    }
    Q[0][0]=H_old[0][0];
    Q[0][1]=H_old[0][1];
    Q[0][2]=H_old[0][2];
    Q[1][0]=H0H0*H_old[1][0]-H0H1*H_old[0][0];
    Q[1][1]=H0H0*H_old[1][1]-H0H1*H_old[0][1];
    Q[1][2]=H0H0*H_old[1][2]-H0H1*H_old[0][2];
    Q[2][0]=H_old[0][1]*H_old[1][2]-H_old[0][2]*H_old[1][1];
    Q[2][1]=H_old[0][2]*H_old[1][0]-H_old[0][0]*H_old[1][2];
    Q[2][2]=H_old[0][0]*H_old[1][1]-H_old[0][1]*H_old[1][0];
    for(int i=0;i<3;i++)
    {
        H0H0=0.0;
        for(int j=0;j<3;j++)
            H0H0+=Q[i][j]*Q[i][j];
        
        H0H1=sqrt(H0H0);
        for(int j=0;j<3;j++)
            Q[i][j]/=H0H1;
    }
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            H_new[i][j]=0.0;
    
    for(int i=0;i<3;i++)
        for(int j=0;j<i+1;j++)
        {
            for(int k=0;k<3;k++)
                H_new[i][j]+=H_old[i][k]*Q[j][k];
        }
    
    for(int i=0;i<3;i++)
        delete [] Q[i];
    delete [] Q;
}
/*--------------------------------------------
 inversion funtion to calculate B whenever H
 is changed
 --------------------------------------------*/
void XMath::invert(type0** A,type0** Ainv,int dim)
{
    if(dim==0)
        return;
    
    
    type0** ATA=new type0*[dim];
    for(int idim=0;idim<dim;idim++)
        ATA[idim]=new type0[dim];
    type0* c=new type0[dim];
    type0* x=new type0[dim];
    type0* g=new type0[dim];
    type0* g0=new type0[dim];
    type0* h=new type0[dim];
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
    
    for(int i=0;i<dim;i++)
        delete [] ATA[i];
    delete [] ATA;
    delete [] c;
    delete [] x;
    delete [] g;
    delete [] g0;
    delete [] h;
}
/*--------------------------------------------
 inversion funtion to calculate B whenever H
 is changed
 test: 
     XMath* xmath=new XMath();
    
    int dim=4;
    type0** A;
    type0** A_inv;
    type0** I;
    CREATE_2D(A,dim,dim);
    CREATE_2D(A_inv,dim,dim);
    CREATE_2D(I,dim,dim);
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            A[i][j]=0.0;

    A[0][0]=1.25;
    A[1][0]=23.45; A[1][1]=0.87;
    A[2][0]=33.9;  A[2][1]=52.08; A[2][2]=7.85;
    A[3][0]=1.025; A[3][1]=75.9;  A[3][2]=9.06; A[3][3]=12.35;
    
    xmath->invert_lower_triangle(A,A_inv,dim);
    
    for(int i=0;i<dim;i++)
    {
        for(int j=0;j<dim;j++)
        {
            I[i][j]=0.0;
            for(int k=0;k<dim;k++)
                I[i][j]+=A[i][k]*A_inv[k][j];
            
            printf("%lf ",I[i][j]);
        }
        printf("\n");
    }
    
    DEL_2D(A);
    DEL_2D(A_inv);
    DEL_2D(I);
 --------------------------------------------*/
void XMath::invert_lower_triangle(type0**& A,type0**& A_inv,int& dim)
{
    for(int i=0;i<dim;i++)
    {
        A_inv[i][i]=1.0/A[i][i];
        for(int j=i+1;j<dim;j++)
            A_inv[i][j]=0.0;
        for(int j=i-1;j>-1;j--)
        {
            A_inv[i][j]=0.0;
            for(int k=j;k<i;k++)
                A_inv[i][j]-=A[i][k]*A_inv[k][j];
            A_inv[i][j]*=A_inv[i][i];
        }
    }
}
/*--------------------------------------------
 construct a legendre polynomial of degree n
 --------------------------------------------*/
void XMath::quadrature_lg(int n,type0* x,type0* w)
{
    int m=n/2+1;
    int iter,ord,icurs,jcurs;
    int max_iter=50;
    type0 a,u0,inv_u0,f,up,df,tmp0,tol;
    type0 ii,jj,del_u0;
    type0* p_coef=new type0[m];
    type0* dp_coef=new type0[m];
    
    up=0.0;
    for(int i=0;i<m;i++)
        p_coef[i]=dp_coef[i]=0.0;
    
    p_coef[0]=1.0;
    
    for(int i=1;i<n+1;i++)
    {
        ii=static_cast<type0>(i);
        if(i%2==0)
        {
            m=i/2+1;
            p_coef[m-1]=(2.0-1.0/ii)*p_coef[m-2];
            for(int j=m-2;j>0;j--)
            {
                jj=static_cast<type0>(j);
                p_coef[j]*=-(2.0*jj+1)/ii;
                p_coef[j]+=(1.0+(2*jj-1.0)/ii)*p_coef[j-1];
            }
            p_coef[0]*=-1.0/ii;
        }
        else
        {
            m=(i+1)/2;
            for(int j=0;j<m-1;j++)
            {
                jj=static_cast<type0>(j);
                p_coef[j]*=1.0+2.0*jj/ii;
                p_coef[j]-=2.0*(jj+1.0)*p_coef[j+1]/ii;
            }
            //2.0-1.0/ii=(1+2*(m-1)/i)
            p_coef[m-1]*=2.0-1.0/ii;
        }
    }
    
    m=n/2+1;
    
    for(int i=1;i<m;i++)
    {
        dp_coef[i-1]=p_coef[i]*static_cast<type0>(i);
    }
    
    tol=std::numeric_limits<type0>::epsilon();
    ord=m;
    icurs=n-1;
    a=p_coef[m-1];
    
    for(int i=0;i<m;i++)
        p_coef[i]/=a;
    for(int i=1;i<m;i++)
        dp_coef[i-1]/=a;
    
    while (ord>1)
    {
        u0=1.0;
        f=1.0;
        del_u0=0.0;
        iter=max_iter;
        while (fabs(f)>tol && iter)
        {
            u0+=del_u0;
            df=f=0.0;
            tmp0=1.0;
            for(int j=0;j<ord-1;j++)
            {
                f+=p_coef[j]*tmp0;
                df+=dp_coef[j]*tmp0;
                tmp0*=u0;
            }
            f+=p_coef[ord-1]*tmp0;
            del_u0=-f/df;
            iter--;
        }
        
        x[icurs]=sqrt(u0);
        
        inv_u0=1.0/u0;
        p_coef[0]*=-inv_u0;
        for(int i=1;i<ord-1;i++)
        {
            p_coef[i]*=-inv_u0;
            p_coef[i]+=inv_u0*p_coef[i-1];
            dp_coef[i-1]=p_coef[i]*static_cast<type0>(i);
        }
        
        ord--;
        icurs--;
    }
    delete [] p_coef;
    delete [] dp_coef;
    
    
    if(n%2==0)
    {
        icurs++;
        jcurs=icurs-1;
        for(int i=icurs;i<n;i++)
        {
            tmp0=a;
            for(int j=icurs;j<n;j++)
            {
                if(i!=j)
                    tmp0*=x[i]*x[i]-x[j]*x[j];
                else
                    tmp0*=x[i]+x[j];
            }
            w[i]=2.0/(tmp0*tmp0*(1.0-x[i]*x[i]));
            w[jcurs]=w[i];
            x[jcurs]=-x[i];
            jcurs--;
        }
    }
    else
    {
        x[icurs]=0.0;
        icurs++;
        tmp0=a;
        for(int i=icurs;i<n;i++)
            tmp0*=-x[i]*x[i];
        w[icurs-1]=2.0/(tmp0*tmp0);
        
        jcurs=icurs-2;
        for(int i=icurs;i<n;i++)
        {
            tmp0=a*x[i];
            for(int j=icurs;j<n;j++)
            {
                if(i!=j)
                    tmp0*=x[i]*x[i]-x[j]*x[j];
                else
                    tmp0*=x[i]+x[j];
            }
            w[i]=2.0/(tmp0*tmp0*(1.0-x[i]*x[i]));
            w[jcurs]=w[i];
            x[jcurs]=-x[i];
            jcurs--;
        }
    }
}
/*--------------------------------------------
 construct a legendre polynomial of degree n
 --------------------------------------------*/
void XMath::quadrature_hg(int n,type0* x,type0* w)
{
    int m=n/2+1;
    int iter,ord,icurs;
    int max_iter=50;
    type0 a,u0,inv_u0,f,up,df,tmp0,tol,tmp1;
    type0 ii,del_u0;

    
    type0* p_1=new type0[m];
    type0* p_2=new type0[m];
    type0* p_coef=new type0[m];
    type0* dp_coef=new type0[m];
    
    up=0.0;
    for(int i=0;i<m;i++)
        p_1[i]=p_2[i]=p_coef[i]=dp_coef[i]=0.0;
    
    p_1[0]=1.0;
    p_coef[0]=1.0;

    
    
    for(int i=1;i<n+1;i++)
    {
        ii=static_cast<type0>(i);
        m=i/2+1;
        if(i%2==0)
        {
            for(int j=0;j<m;j++)
                p_coef[j]=-2.0*(ii-1.0)*p_2[j];
            for(int j=1;j<m;j++)
                p_coef[j]+=2.0*p_1[j-1];
        }
        else
        {
            for(int j=0;j<m;j++)
                p_coef[j]=2.0*p_1[j]-2.0*(ii-1.0)*p_2[j];
        }
        
        for(int j=0;j<m;j++)
        {
            p_2[j]=p_1[j];
            p_1[j]=p_coef[j];
        }
    }
    
    m=n/2+1;

    
    for(int i=1;i<m;i++)
    {
        dp_coef[i-1]=p_coef[i]*static_cast<type0>(i);
    }
    
    tol=std::numeric_limits<type0>::epsilon();
    ord=m;
    icurs=n/2;
    if(n%2==1)
    {
        x[icurs]=0.0;
        icurs++;
    }
    
    a=p_coef[m-1];
    
    for(int i=0;i<m;i++)
        p_coef[i]/=a;
    for(int i=1;i<m;i++)
        dp_coef[i-1]/=a;
    
    

    while (ord>1)
    {
        u0=0.0;
        f=1.0;
        del_u0=0.0;
        iter=max_iter;
        while(fabs(f)>tol && iter)
        {
            u0+=del_u0;
            df=f=0.0;
            tmp0=1.0;
            for(int j=0;j<ord-1;j++)
            {
                f+=p_coef[j]*tmp0;
                df+=dp_coef[j]*tmp0;
                tmp0*=u0;
            }
            f+=p_coef[ord-1]*tmp0;
            del_u0=-f/df;
            iter--;
        }
        
        x[icurs]=sqrt(u0);
        
        tmp0=1.0;
        tmp1=0.0;
        for(int i=0;i<m;i++)
        {
            tmp1+=p_2[i]*tmp0;
            tmp0*=u0;
        }
        
        if(n%2==0)
            tmp1*=sqrt(u0);
        
        w[icurs]=1.0/(tmp1*tmp1);
        
        
        inv_u0=1.0/u0;
        p_coef[0]*=-inv_u0;
        for(int i=1;i<ord-1;i++)
        {
            p_coef[i]*=-inv_u0;
            p_coef[i]+=inv_u0*p_coef[i-1];
            dp_coef[i-1]=p_coef[i]*static_cast<type0>(i);
        }
        
        ord--;
        icurs++;
    }
    
    tmp0=sqrt(M_PI)/static_cast<type0>(n);
    for(int i=0;i<n-1;i++)
    {
        tmp0*=static_cast<type0>(2*(i+1));
    }
    if(n%2==1)
    {
        w[n/2]=tmp0/(p_2[0]*p_2[0]);
    }


    for(int i=0;i<n/2;i++)
    {
        icurs--;
        w[icurs]*=tmp0;
        w[i]=w[icurs];
        x[i]=-x[icurs];
    }
    
    delete [] p_1;
    delete [] p_2;
    delete [] p_coef;
    delete [] dp_coef;

}
/*--------------------------------------------
 calculates square root of 3x3 matrix
 ref: L. P. Franca
 An Algorithm to Compute The Square Root of
 a 3x3 Positive Definite Matrix
 Computers Math. Applic. Vol. 18, No. 5,
 pp. 459-466, 1989
 --------------------------------------------*/
int XMath::M3sqroot(type0** A,type0** Asq)
{
    type0 IA=0;
    for(int i=0;i<3;i++)
        IA+=A[i][i];
    type0 IIA=0;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            IIA-=A[i][j]*A[j][i];
    
    IIA+=IA*IA;
    IIA*=0.5;
    type0 IIIA=M3DET(A);
    type0 k=IA*IA-3*IIA;
    if(k<0.0)
        return 0;
    if(k<TOLERANCE)
    {
        if(IA<=0.0)
            return 0;
        M3ZERO(Asq);
        for(int i=0;i<3;i++)
            Asq[i][i]=sqrt(IA/3.0);
        return 1;
    }
    
    type0 l=IA*(IA*IA -4.5*IIA)+13.5*IIIA;
    
    type0 temp=l/(k*sqrt(k));
    if(temp>1.0||temp<-1.0)
        return 0;
    
    type0 phi=acos(temp);
    type0 lambda=sqrt((1.0/3.0)*(IA+2*sqrt(k)*cos(phi/3.0)));
    
    type0 IIIAsq=sqrt(IIIA);
    type0 y=-lambda*lambda+IA+2*IIIAsq/lambda;
    if(y<0.0)
        return 0;
    type0 IAsq=lambda+sqrt(y);
    type0 IIAsq=0.5*(IAsq*IAsq-IA);
    
    type0 coef0=IAsq*IIAsq-IIIAsq;
    if(coef0==0)
        return 0;
    coef0=1.0/coef0;
    
    M3ZERO(Asq);
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
                Asq[i][j]-=coef0*A[i][k]*A[k][j];
    
    type0 coef1=coef0*(IAsq*IAsq-IIAsq);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            Asq[i][j]+=coef1*A[i][j];
    
    type0 coef2=coef0*IAsq*IIIAsq;
    for(int i=0;i<3;i++)
        Asq[i][i]+=coef2;
    return 1;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
SPARSE::
SPARSE(MAPP* mapp,int n0,int n1): InitPtrs(mapp)
{
    h0=n0;
    h1=n1;
    no_vec=0;
    no_elem=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
SPARSE::~SPARSE()
{
    if(no_elem)
    {
        delete [] idx0;
        delete [] idx1;
        delete [] mtrx;
    }
    
    if(no_vec)
    {
        delete [] fin_pos;
        delete [] idx_vec;
        delete [] init_pos;
    }
    
}
/*--------------------------------------------
 add an element to the sparse
 --------------------------------------------*/
void SPARSE::add(int i0,int i1,double elem)
{
    if(i0>=h0 || i1>=h1 || i0<0 || i1<0)
        error->abort("in sparse the index "
                     "is higher than the number "
                     "columns/rows");
    
    if(elem==0.0)
        return;
    
    int no=0;
    while (no<no_elem && (idx0[no]!=i0
           || idx1[no]!=i1))
        no++;
    
    if(no<no_elem)
    {
        mtrx[no]=elem;
        return;
    }
    
    GROW(idx0,no_elem,no_elem+1);
    GROW(idx1,no_elem,no_elem+1);
    GROW(mtrx,no_elem,no_elem+1);
    
    idx0[no_elem]=i0;
    idx1[no_elem]=i1;
    mtrx[no_elem]=elem;
    no_elem++;
}
/*--------------------------------------------
 add elements
 --------------------------------------------*/
void SPARSE::add(int size,int* i0_lst,int* i1_lst,double* elem_lst)
{
    
    double elem;
    int i0,i1;
    for(int i=0; i<size;i++)
    {
        elem=elem_lst[i];
        i0=i0_lst[i];
        i1=i1_lst[i];
        
        if(elem==0.0)
            continue;

        
        int no=0;
        while (no<no_elem && (idx0[no]!=i0
                || idx1[no]!=i1))
            no++;
        
        if(no<no_elem)
        {
            mtrx[no]=elem;
            continue;
        }
        
        GROW(idx0,no_elem,no_elem+1);
        GROW(idx1,no_elem,no_elem+1);
        GROW(mtrx,no_elem,no_elem+1);
        
        idx0[no_elem]=i0;
        idx1[no_elem]=i1;
        mtrx[no_elem]=elem;
        no_elem++;
    }
}
/*--------------------------------------------
 sort the sparse matrix row/column wise
 --------------------------------------------*/
void SPARSE::sort(int r_c)
{
    int tmp0,tmp1;
    double tmpx;
    if(r_c==0)
    {
        for(int i=0;i<no_elem;i++)
        {
            for(int j=i+1;j<no_elem;j++)
            {
                if(idx0[j]<idx0[i])
                {
                    tmp0=idx0[j];
                    tmp1=idx1[j];
                    tmpx=mtrx[j];
                    
                    idx0[j]=idx0[i];
                    idx1[j]=idx1[i];
                    mtrx[j]=mtrx[i];
                    
                    idx0[i]=tmp0;
                    idx1[i]=tmp1;
                    mtrx[i]=tmpx;
                }
            }
        }
    }
    else if(r_c==1)
    {
        for(int i=0;i<no_elem;i++)
        {
            for(int j=i+1;j<no_elem;j++)
            {
                if(idx1[j]<idx1[i])
                {
                    tmp0=idx0[j];
                    tmp1=idx1[j];
                    tmpx=mtrx[j];
                    
                    idx0[j]=idx0[i];
                    idx1[j]=idx1[i];
                    mtrx[j]=mtrx[i];
                    
                    idx0[i]=tmp0;
                    idx1[i]=tmp1;
                    mtrx[i]=tmpx;
                }
            }
        }
    }
    
    int strt=0;
    int end;
    int pos=strt;
    
    if(r_c==0)
    {
        while(pos<no_elem)
        {
            while(idx0[pos]==idx0[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            end=pos;
            for(int j=strt;j<end;j++)
            {
                for(int k=j+1;k<end;k++)
                {
                    if(idx1[j]>idx1[k])
                    {
                        tmp1=idx1[k];
                        tmpx=mtrx[k];
                        
                        idx1[k]=idx1[j];
                        mtrx[k]=mtrx[j];
                        
                        idx1[j]=tmp1;
                        mtrx[j]=tmpx;
                    }
                }
            }
            strt=pos;
        }
    }
    else if(r_c==1)
    {
        while(pos<no_elem)
        {
            while(idx1[pos]==idx1[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            end=pos;
            for(int j=strt;j<end;j++)
            {
                for(int k=j+1;k<end;k++)
                {
                    if(idx0[j]>idx0[k])
                    {
                        tmp1=idx0[k];
                        tmpx=mtrx[k];
                        
                        idx0[k]=idx0[j];
                        mtrx[k]=mtrx[j];
                        
                        idx0[j]=tmp1;
                        mtrx[j]=tmpx;
                    }
                }
            }
            strt=pos;
        }
    }
}
/*--------------------------------------------
 vectorize the matrix row/column wise
 --------------------------------------------*/
void SPARSE::vectorize(int r_c)
{
    sort(r_c);
    
    if(no_vec)
    {
        delete [] fin_pos;
        delete [] idx_vec;
        delete [] init_pos;
    }
    
    
    no_vec=0;
    int strt=0;
    int pos=strt;
    if(r_c==0)
    {
        while(pos<no_elem)
        {
            while(idx0[pos]==idx0[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            GROW(init_pos,no_vec,no_vec+1);
            GROW(fin_pos,no_vec,no_vec+1);
            GROW(idx_vec,no_vec,no_vec+1);
            
            init_pos[no_vec]=strt;
            fin_pos[no_vec]=pos;
            idx_vec[no_vec]=idx0[strt];
            
            no_vec++;
            strt=pos;
        }
    }
    else if(r_c==1)
    {
        while(pos<no_elem)
        {
            while(idx1[pos]==idx1[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            
            GROW(init_pos,no_vec,no_vec+1);
            GROW(fin_pos,no_vec,no_vec+1);
            GROW(idx_vec,no_vec,no_vec+1);
            
            init_pos[no_vec]=strt;
            fin_pos[no_vec]=pos;
            idx_vec[no_vec]=idx1[strt];
            
            no_vec++;
            strt=pos;
        }
    }
}

/*--------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 --------------------------------------------*/
/*--------------------------------------------
 constructor: uses CG to solve Ax=b
 --------------------------------------------*/
SOLVEAXb::SOLVEAXb(MAPP* mapp,SPARSE* Mat,
double* b,int size): InitPtrs(mapp)
{
    A=Mat;

    if(A->h1!=size)
        error->abort("Dimensions of A and b does not match");

    myno=atoms->my_p;
    totp=atoms->tot_p;
    if(size<totp)
    error->abort("Size of vector cannot be "
        "less than number processor");
    
    
    lsize=size/totp;
    
    lcl_lo=myno*lsize;
    if(size%totp!=0)
        if(myno==totp-1)
            lsize+=size%totp;
    lcl_hi=lcl_lo+lsize;
    
    A->vectorize(1);
    
    int no_vec=A->no_vec;
    int* idx_vec=A->idx_vec;
    int* init_pos=A->init_pos;
    int* fin_pos=A->fin_pos;
    
    int pos=0;
    while (idx_vec[pos]<lcl_lo)
    {
        pos++;
        if(pos==no_vec)
            break;
    }
    int idx_pos_lo=pos;
    
    while (idx_vec[pos]<lcl_hi)
    {
        pos++;
        if(pos==no_vec)
            break;
    }
    int idx_pos_hi=pos;
    

    
    int i1,crsi0,crsi1;
    int comp0,comp1;
    double val;
    tot_lsize=lsize;
    int comm_idx;

    CREATE1D(comm_snd_size,totp-1);
    CREATE1D(comm_snd_lst,totp-1);
    CREATE1D(comm_rcv_size,totp-1);
    for(int i=0;i<totp-1;i++)
        comm_snd_size[i]=comm_rcv_size[i]=0;
    
    lcl_no_elem=0;
    for(int i0=idx_pos_lo;i0<idx_pos_hi;i0++)
    {
        for(int j=0;j<no_vec;j++)
        {
            i1=j+idx_pos_lo;
            if(i1>=no_vec)
                i1-=no_vec;
            val=0.0;
            
            crsi0=init_pos[i0];
            crsi1=init_pos[i1];
            while(crsi0<fin_pos[i0])
            {
                while(crsi1<fin_pos[i1] && A->idx0[crsi1]<A->idx0[crsi0])
                    crsi1++;
                if(crsi1<fin_pos[i1] && A->idx0[crsi1]==A->idx0[crsi0])
                    val+=A->mtrx[crsi0]*A->mtrx[crsi1];
                crsi0++;
            }
            if(val!=0.0)
            {
                comp0=idx_vec[i0]-lcl_lo;
                if(idx_vec[i1]>=lcl_lo && idx_vec[i1]<lcl_hi)
                {
                    comp1=idx_vec[i1]-lcl_lo;
                    
                    GROW(lcl_idx0,lcl_no_elem,lcl_no_elem+1);
                    GROW(lcl_idx1,lcl_no_elem,lcl_no_elem+1);
                    GROW(lcl_mtrx,lcl_no_elem,lcl_no_elem+1);
                    
                    lcl_idx0[lcl_no_elem]=comp0;
                    lcl_idx1[lcl_no_elem]=comp1;
                    lcl_mtrx[lcl_no_elem]=val;
                }
                else
                {
                    /*
                    comm_idx=0;
                    while (comm_idx<totp && idx_vec[i1]<comp_lim[comm_idx])
                        comm_idx++;
                    comm_idx--;*/
                    comm_idx=(idx_vec[i1]/(size/totp));
                    if(comm_idx>totp-1)
                        comm_idx=totp-1;
                    comm_idx-=myno+1;
                    if(comm_idx<0)
                        comm_idx+=totp;
                    
                    GROW(comm_snd_lst[comm_idx],comm_snd_size[comm_idx],comm_snd_size[comm_idx]+1);
                    comm_snd_lst[comm_idx][comm_snd_size[comm_idx]]=comp0;
                    comm_snd_size[comm_idx]++;
                    comm_rcv_size[totp-2-comm_idx]++;
                    comp1=tot_lsize;
                    GROW(lcl_idx0,lcl_no_elem,lcl_no_elem+1);
                    GROW(lcl_idx1,lcl_no_elem,lcl_no_elem+1);
                    GROW(lcl_mtrx,lcl_no_elem,lcl_no_elem+1);
                    
                    lcl_idx0[lcl_no_elem]=comp0;
                    lcl_idx1[lcl_no_elem]=comp1;
                    lcl_mtrx[lcl_no_elem]=val;
                    
                    tot_lsize++;
                }
                lcl_no_elem++;
            }
        }
    }
    /*
    int sndp,rcvp;
    MPI_Status status;
    
    for(int i=0;i<totp-1;i++)
    {
        rcvp=myno+i+1;
        sndp=myno-i-1;

        if(rcvp>=totp)
            rcvp-=totp;
        if(sndp<0)
            sndp+=totp;
        
        MPI_Sendrecv(&comm_snd_size[totp-2-i],1,
        MPI_INT,sndp,0,&comm_rcv_size[i],
                     1,MPI_INT,rcvp,0,world,
                     &status);
    }
    */

    CREATE1D(buff_snd,totp-1);
    for(int i=0;i<totp-1;i++)
    {
        CREATE1D(buff_snd[i],comm_snd_size[i]);
    }
 
    CREATE1D(c,lsize);
 
    for(int i=0;i<lsize;i++)
        c[i]=0.0;

    for(int i=idx_pos_lo;i<idx_pos_hi;i++)
        for(int j=init_pos[i];j<fin_pos[i];j++)
            c[idx_vec[i]-lcl_lo]+=b[A->idx0[j]]*(A->mtrx[j]);
    

    d_sq=0.0;
    for(int i=0;i<size;i++)
        d_sq+=b[i]*b[i];

    
    CREATE1D(x,tot_lsize);
    CREATE1D(g,2);
    CREATE1D(g[0],lsize);
    CREATE1D(g[1],lsize);
    CREATE1D(h,tot_lsize);
    
    for(int i=0;i<tot_lsize;i++)
        x[i]=0.0;
    for(int i=lcl_lo;i<lcl_hi;i++)
    {
        //x[i-lcl_lo]=b[i];
        x[i-lcl_lo]=0.0;
    }
    xchng(x);

    CREATE1D(ans,size);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
SOLVEAXb::~SOLVEAXb()
{
    if(totp-1)
    {
        for(int i=0;i<totp-1;i++)
        {
            if(comm_snd_size[i])
            {
                delete [] buff_snd[i];
                delete [] comm_snd_lst[i];
            }
        }
        
        delete [] buff_snd;
        delete [] comm_snd_lst;
        delete [] comm_rcv_size;
        delete [] comm_snd_size;
    }
    if(lcl_no_elem)
    {
        delete [] lcl_idx0;
        delete [] lcl_idx1;
        delete [] lcl_mtrx;
    }
    if(lsize)
    {
        delete [] x;
        delete [] h;
        delete [] c;
        delete [] g[0];
        delete [] g[1];
        delete [] g;
    }
    delete [] ans;
}
/*--------------------------------------------
packing
--------------------------------------------*/
void SOLVEAXb::xchng(double* vec)
{
    MPI_Request request;
    MPI_Status status;
    int rcvp,sndp;
    int curs=lsize;
    for(int i=0;i<totp-1;i++)
    {
        rcvp=myno-i-1;
        sndp=myno+i+1;
        
        if(sndp>=totp)
            sndp-=totp;
        if(rcvp<0)
            rcvp+=totp;
        for(int j=0;j<comm_snd_size[i];j++)
            buff_snd[i][j]=vec[comm_snd_lst[i][j]];
        if(comm_rcv_size[i])
        {
            MPI_Irecv(&vec[curs],comm_rcv_size[i],
                      MPI_DOUBLE,rcvp,0,world,
                      &request);
        }
        if(comm_snd_size[i])
            MPI_Send(&buff_snd[i][0],comm_snd_size[i],
                     MPI_DOUBLE,sndp,0,world);
        if(comm_rcv_size[i])
            MPI_Wait(&request,&status);
        curs+=comm_rcv_size[i];
    }
}
/*--------------------------------------------
 packing
 --------------------------------------------*/
int SOLVEAXb::solve(double tol)
{
    int niter=A->h1+1;
    int comp;
    double g1_sq_l,g01_sq_l,hBh_l,hg_l,xBx_l,cx_l;
    double g0_sq=0.0,g1_sq,g01_sq,hBh,hg,xBx,cx,coef;
    double error=0.0,alpha;
    
    for(int i=0;i<lsize;i++)
        h[i]=0.0;
    for(int iter=0;iter<niter;iter++)
    {
        
        if(iter)
        {
            comp=iter%2;
            for(int i=0;i<lsize;i++)
                g[comp][i]=0.0;
            
            for(int i=0;i<lcl_no_elem;i++)
                g[comp][lcl_idx0[i]]-=2.0*lcl_mtrx[i]*x[lcl_idx1[i]];
            
            
            g1_sq=g1_sq_l=0.0;
            g01_sq=g01_sq_l=0.0;
            for(int i=0;i<lsize;i++)
            {
                g[comp][i]+=2*c[i];
                g1_sq_l+=g[comp][i]*g[comp][i];
                g01_sq_l+=g[1-comp][i]*g[comp][i];
            }
            
            MPI_Allreduce(&g1_sq_l,&g1_sq,1,MPI_DOUBLE,MPI_SUM,world);
            MPI_Allreduce(&g01_sq_l,&g01_sq,1,MPI_DOUBLE,MPI_SUM,world);
            
            
            coef=(g1_sq-g01_sq)/g0_sq;
            for(int i=0;i<lsize;i++)
            {
                h[i]*=coef;
                h[i]+=g[comp][i];
            }
        }
        else
        {
            comp=iter%2;
            for(int i=0;i<lsize;i++)
                h[i]=0.0;
            
            for(int i=0;i<lcl_no_elem;i++)
                h[lcl_idx0[i]]-=2.0*lcl_mtrx[i]*x[lcl_idx1[i]];
            for(int i=0;i<lsize;i++)
                h[i]+=2*c[i];
            for(int i=0;i<lsize;i++)
                g[comp][i]=h[i];
            
            g1_sq=g1_sq_l=0.0;
            for(int i=0;i<lsize;i++)
                g1_sq_l+=h[i]*h[i];
            
            MPI_Allreduce(&g1_sq_l,&g1_sq,1,MPI_DOUBLE,MPI_SUM,world);
        }

        xchng(h);
        
        hBh=hBh_l=0.0;
        for(int i=0;i<lcl_no_elem;i++)
            hBh_l+=h[lcl_idx0[i]]*lcl_mtrx[i]*h[lcl_idx1[i]];
        MPI_Allreduce(&hBh_l,&hBh,1,MPI_DOUBLE,MPI_SUM,world);
        if(hBh<=0)
        {
            printf("The matrix is singular");
            return 1;
        }
        
        hg_l=hg=0.0;
        for(int i=0;i<lsize;i++)
            hg_l+=g[comp][i]*h[i];
        MPI_Allreduce(&hg_l,&hg,1,MPI_DOUBLE,MPI_SUM,world);
        
        alpha=0.5*hg/hBh;
        if(alpha<=0)
        {
            printf("The matrix is singular (alpha became negative)");
            return 1;
        }
        
        for(int i=0;i<lcl_no_elem;i++)
            x[i]+=alpha*h[i];
        
        xchng(x);
        
        xBx=xBx_l=0.0;
        for(int i=0;i<lcl_no_elem;i++)
            xBx_l+=x[lcl_idx0[i]]*lcl_mtrx[i]*x[lcl_idx1[i]];
        cx_l=cx=0.0;
        for(int i=0;i<lsize;i++)
            cx_l+=x[i]*c[i];
        
        MPI_Allreduce(&xBx_l,&xBx,1,MPI_DOUBLE,MPI_SUM,world);
        MPI_Allreduce(&cx_l,&cx,1,MPI_DOUBLE,MPI_SUM,world);
        
        error=abs(xBx-2.0*cx+d_sq);
        if(error<=tol)
        {
            if(atoms->my_p==0)
                printf("Converged in %d steps!\n",iter);
            construct_ans();
            return 0;
        }
        g0_sq=g1_sq;
    }
    printf("Did not converge error: %lf!\n",error);
    construct_ans();
    return 1;
}
/*--------------------------------------------
 packing
 --------------------------------------------*/
void SOLVEAXb::construct_ans()
{
    int size=A->h1;
    double* ans_l;
    CREATE1D(ans_l,size);
    for(int i=0;i<size;i++)
        ans_l[i]=ans[i]=0.0;

    for(int i=lcl_lo;i<lcl_hi;i++)
        ans_l[i]=x[i-lcl_lo];
    
    MPI_Allreduce(&ans_l[0],&ans[0],size,MPI_DOUBLE,MPI_SUM,world);
    delete [] ans_l;
}
/*--------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 --------------------------------------------*/
/*--------------------------------------------
 constructor
 --------------------------------------------*/
SPARSE_P::
SPARSE_P(MAPP* mapp,int n0,int n1): InitPtrs(mapp)
{
    h0=n0;
    h1=n1;
    no_elem=lcl_no_elem=0;
    myno=atoms->my_p;
    totp=atoms->tot_p;
    lcl_h0=static_cast<int>((1.0/totp)*h0+0.5);
    lo_h0=myno*lcl_h0;
    if(myno==totp-1)
        lcl_h0=h0-lcl_h0*(totp-1);
    hi_h0=lcl_h0+lcl_h0;
    
    lcl_no_vec=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
SPARSE_P::~SPARSE_P()
{

    if(lcl_no_elem)
    {
        delete [] lcl_idx0;
        delete [] lcl_idx1;
        delete [] lcl_mtrx;
    }
    if(lcl_no_vec)
    {
        delete [] lcl_idx_vec;
        delete [] lcl_init_pos;
        delete [] lcl_fin_pos;
    }
}
/*--------------------------------------------
 add elements
 --------------------------------------------*/
void SPARSE_P::add(int i0,int i1,double elem)
{
    if(i0>=h0 || i1>=h1 || i0<0 || i1<0)
        error->abort("in sparse the index "
                     "is higher than the number "
                     "columns/rows");
    
    if(elem==0.0)
        return;
    
    if(i0<lo_h0 || i0>=hi_h0)
        return;
    
    int lcl_i0=i0-lo_h0;
    
    int no=0;
    while (no<lcl_no_elem && (lcl_idx0[no]!=lcl_i0
                          || lcl_idx1[no]!=i1))
        no++;
    
    if(no<lcl_no_elem)
    {
        lcl_mtrx[no]=elem;
        return;
    }
    
    GROW(lcl_idx0,lcl_no_elem,no_elem+1);
    GROW(lcl_idx1,lcl_no_elem,no_elem+1);
    GROW(lcl_mtrx,lcl_no_elem,no_elem+1);
    
    lcl_idx0[no_elem]=lcl_i0;
    lcl_idx1[no_elem]=i1;
    lcl_mtrx[no_elem]=elem;
    lcl_no_elem++;
    no_elem++;
    
    MPI_Bcast(&no_elem,1,MPI_INT,0,world);
}
/*--------------------------------------------
 add elements
 --------------------------------------------*/
void SPARSE_P::addd(int i0,int i1,double elem)
{
    if(i0>=h0 || i1>=h1 || i0<0 || i1<0)
        error->abort("in sparse the index "
                     "is higher than the number "
                     "columns/rows");
    
    if(elem==0.0)
        return;
    
    if(i0<lo_h0 || i0>=hi_h0)
        return;
    
    int lcl_i0=i0-lo_h0;
    
    int no=0;
    while (no<lcl_no_elem && (lcl_idx0[no]!=lcl_i0
                              || lcl_idx1[no]!=i1))
        no++;
    
    if(no<lcl_no_elem)
    {
        lcl_mtrx[no]+=elem;
        return;
    }
    
    GROW(lcl_idx0,lcl_no_elem,no_elem+1);
    GROW(lcl_idx1,lcl_no_elem,no_elem+1);
    GROW(lcl_mtrx,lcl_no_elem,no_elem+1);
    
    lcl_idx0[no_elem]=lcl_i0;
    lcl_idx1[no_elem]=i1;
    lcl_mtrx[no_elem]=elem;
    lcl_no_elem++;
    no_elem++;
    
    MPI_Bcast(&no_elem,1,MPI_INT,0,world);
}
/*--------------------------------------------
 add elements
 --------------------------------------------*/
void SPARSE_P::add(int size,int* i0_lst,int* i1_lst,double* elem_lst)
{
    
    double elem;
    int i0,i1,lcl_i0;
    for(int i=0; i<size;i++)
    {
        elem=elem_lst[i];
        i0=i0_lst[i];
        i1=i1_lst[i];
        
        if(elem==0.0)
            continue;
        
        if(i0<lo_h0 || i0>=hi_h0)
            continue;
        
        lcl_i0=i0-lo_h0;
        
        int no=0;
        while (no<lcl_no_elem && (lcl_idx0[no]!=lcl_i0
                                  || lcl_idx1[no]!=i1))
            no++;
        
        if(no<lcl_no_elem)
        {
            lcl_mtrx[no]=elem;
            continue;
        }
        
        GROW(lcl_idx0,lcl_no_elem,no_elem+1);
        GROW(lcl_idx1,lcl_no_elem,no_elem+1);
        GROW(lcl_mtrx,lcl_no_elem,no_elem+1);
        
        lcl_idx0[no_elem]=lcl_i0;
        lcl_idx1[no_elem]=i1;
        lcl_mtrx[no_elem]=elem;
        lcl_no_elem++;
    }

    MPI_Allreduce(&lcl_no_elem,&no_elem,1,MPI_INT,MPI_SUM,world);
}
/*--------------------------------------------
 sort the sparse matrix row/column wise
 --------------------------------------------*/
void SPARSE_P::sort(int r_c)
{
    int tmp0,tmp1;
    double tmpx;
    if(r_c==0)
    {
        for(int i=0;i<lcl_no_elem;i++)
        {
            for(int j=i+1;j<lcl_no_elem;j++)
            {
                if(lcl_idx0[j]<lcl_idx0[i])
                {
                    tmp0=lcl_idx0[j];
                    tmp1=lcl_idx1[j];
                    tmpx=lcl_mtrx[j];
                    
                    lcl_idx0[j]=lcl_idx0[i];
                    lcl_idx1[j]=lcl_idx1[i];
                    lcl_mtrx[j]=lcl_mtrx[i];
                    
                    lcl_idx0[i]=tmp0;
                    lcl_idx1[i]=tmp1;
                    lcl_mtrx[i]=tmpx;
                }
            }
        }
    }
    else if(r_c==1)
    {
        for(int i=0;i<lcl_no_elem;i++)
        {
            for(int j=i+1;j<lcl_no_elem;j++)
            {
                if(lcl_idx1[j]<lcl_idx1[i])
                {
                    tmp0=lcl_idx0[j];
                    tmp1=lcl_idx1[j];
                    tmpx=lcl_mtrx[j];
                    
                    lcl_idx0[j]=lcl_idx0[i];
                    lcl_idx1[j]=lcl_idx1[i];
                    lcl_mtrx[j]=lcl_mtrx[i];
                    
                    lcl_idx0[i]=tmp0;
                    lcl_idx1[i]=tmp1;
                    lcl_mtrx[i]=tmpx;
                }
            }
        }
    }
    
    int strt=0;
    int end;
    int pos=strt;
    
    if(r_c==0)
    {
        while(pos<lcl_no_elem)
        {
            while(lcl_idx0[pos]==lcl_idx0[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            end=pos;
            for(int j=strt;j<end;j++)
            {
                for(int k=j+1;k<end;k++)
                {
                    if(lcl_idx1[j]>lcl_idx1[k])
                    {
                        tmp1=lcl_idx1[k];
                        tmpx=lcl_mtrx[k];
                        
                        lcl_idx1[k]=lcl_idx1[j];
                        lcl_mtrx[k]=lcl_mtrx[j];
                        
                        lcl_idx1[j]=tmp1;
                        lcl_mtrx[j]=tmpx;
                    }
                }
            }
            strt=pos;
        }
    }
    else if(r_c==1)
    {
        while(pos<no_elem)
        {
            while(lcl_idx1[pos]==lcl_idx1[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            end=pos;
            for(int j=strt;j<end;j++)
            {
                for(int k=j+1;k<end;k++)
                {
                    if(lcl_idx0[j]>lcl_idx0[k])
                    {
                        tmp1=lcl_idx0[k];
                        tmpx=lcl_mtrx[k];
                        
                        lcl_idx0[k]=lcl_idx0[j];
                        lcl_mtrx[k]=lcl_mtrx[j];
                        
                        lcl_idx0[j]=tmp1;
                        lcl_mtrx[j]=tmpx;
                    }
                }
            }
            strt=pos;
        }
    }
}
/*--------------------------------------------
 vectorize the matrix row/column wise
 --------------------------------------------*/
void SPARSE_P::vectorize(int r_c)
{
    sort(r_c);
    
    if(lcl_no_vec)
    {
        delete [] lcl_fin_pos;
        delete [] lcl_idx_vec;
        delete [] lcl_init_pos;
    }
    
    
    lcl_no_vec=0;
    int strt=0;
    int pos=strt;
    if(r_c==0)
    {
        while(pos<lcl_no_elem)
        {
            while(lcl_idx0[pos]==lcl_idx0[strt])
            {
                pos++;
                if(pos==lcl_no_elem)
                    break;
            }
            GROW(lcl_init_pos,lcl_no_vec,lcl_no_vec+1);
            GROW(lcl_fin_pos,lcl_no_vec,lcl_no_vec+1);
            GROW(lcl_idx_vec,lcl_no_vec,lcl_no_vec+1);
            
            lcl_init_pos[lcl_no_vec]=strt;
            lcl_fin_pos[lcl_no_vec]=pos;
            lcl_idx_vec[lcl_no_vec]=lcl_idx0[strt];
            
            lcl_no_vec++;
            strt=pos;
        }
    }
    else if(r_c==1)
    {
        while(pos<lcl_no_elem)
        {
            while(lcl_idx1[pos]==lcl_idx1[strt])
            {
                pos++;
                if(pos==lcl_no_elem)
                    break;
            }
            
            GROW(lcl_init_pos,lcl_no_vec,lcl_no_vec+1);
            GROW(lcl_fin_pos,lcl_no_vec,lcl_no_vec+1);
            GROW(lcl_idx_vec,lcl_no_vec,lcl_no_vec+1);
            
            lcl_init_pos[lcl_no_vec]=strt;
            lcl_fin_pos[lcl_no_vec]=pos;
            lcl_idx_vec[lcl_no_vec]=lcl_idx1[strt];
            
            lcl_no_vec++;
            strt=pos;
        }
    }
    
}
/*--------------------------------------------
 test field
 --------------------------------------------*/
COMB::COMB(MAPP* mapp):InitPtrs(mapp)
{
    perm_list_size_0=0;
    perm_list_size_1=0;
}
/*--------------------------------------------
 test field
 --------------------------------------------*/
COMB::~COMB()
{
    clean();
}
/*--------------------------------------------
 test field
 --------------------------------------------*/
void COMB::clean()
{
    if(perm_list_size_0)
    {
        for(int i=0;i<perm_list_size_1;i++)
            delete [] perm_list[i];
    }
    if(perm_list_size_1)
        delete [] perm_list;
    perm_list_size_0=0;
    perm_list_size_1=0;
}
/*--------------------------------------------
 test field
 --------------------------------------------*/
void COMB::comb(int* no_list, int list_size)
{
    int tot=0;
    for(int i=0;i<list_size;i++)
        tot+=no_list[i];

    clean();
    
    perm_list_size_0=tot;
    perm_list_size_1=0;
    
    int* tmp0;
    CREATE1D(tmp0,tot);
    int* tmp1;
    CREATE1D(tmp1,tot);
    
    comb_rec(no_list[0],tot,0,no_list[0]
    ,tmp0,tmp1,no_list,0,tot);
    delete [] tmp0;
    delete [] tmp1;
}
/*--------------------------------------------
 test field
 --------------------------------------------*/
void COMB::comb_rec(int no,int tot,int pos
,int max_no,int*& lvl,int*& tmp,int* y
,int y_pos,int y_tot)
{
    if(no!=0)
    {
        for(int i=pos;i<tot;i++)
        {
            lvl[no-1]=i;
            comb_rec(no-1,tot,i+1,max_no,lvl,tmp,y,y_pos,y_tot);
            
        }
    }
    else
    {
        if(y_pos==0)
            for(int i=0;i<tot;i++)
                tmp[i]=-1;
        int p=0;
        int pp=0;
        int ppp=max_no-1;
        while(ppp>=0)
        {
            if(tmp[p]==-1)
            {
                if(pp==lvl[ppp])
                {
                    tmp[p]=y_pos;
                    ppp--;
                }
                pp++;
            }
            p++;
        }
        
        if(y[y_pos]==tot)
        {

            GROW(perm_list,perm_list_size_1,perm_list_size_1+1);
            CREATE1D(perm_list[perm_list_size_1], perm_list_size_0);
            perm_list_size_1++;
            
            for(int i=0;i<y_tot;i++)
                if(tmp[i]==y_pos)
                    tmp[i]=-1;
        }
        else
        {
            int noo=y[y_pos+1];
            int tott=tot-y[y_pos];
            int* h=lvl+y[y_pos];
            comb_rec(noo,tott,0,noo,h,tmp,y,y_pos+1,y_tot);
            for(int i=0;i<y_tot;i++)
                if(tmp[i]==y_pos)
                    tmp[i]=-1;
        }
    }
}

