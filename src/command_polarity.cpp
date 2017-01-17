#include "command_polarity.h"
#include "atom_types.h"
#include "MAPP.h"
#include "neighbor_cna.h"
#include "dynamic.h"
#include "ff.h"
#include "atoms.h"
#include "comm.h"
#include "xmath.h"
#include "write_cfg.h"
#include "script_reader.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
Command_polarity::Command_polarity(int nargs,char** args)
{
    int nbins=atoi(args[1]);
    
    
    delete forcefield->neighbor;
    forcefield->neighbor=new NeighborCNA(atoms,mapp->type,forcefield->cut_sk_sq);
    
    Vec<type0>* polarity =new Vec<type0>(atoms,3);
    Vec<int>* polarity_xst =new Vec<int>(atoms,1,"BCC");
    
    Dynamic* dynamic=new Dynamic(atoms,comm,{mapp->type},{});
    dynamic->init(false);
    
    

    
    
    
    
    
    int** neighbor_list=forcefield->neighbor->neighbor_list;
    int* neighbor_list_size=forcefield->neighbor->neighbor_list_size;
    int natms=atoms->natms;
    type0* x=atoms->x->begin();
    type0* polar=polarity->begin();
    int* polar_xst=polarity_xst->begin();
    for(int i=0;i<natms;i++)
        polar_xst[i]=calc_polarity(x,i,neighbor_list[i],neighbor_list_size[i],polar+3*i);
    
    dynamic->update(polarity);
    dynamic->update(polarity_xst);
    
    Write_cfg* wr=NULL;
    ScriptReader::create(wr,"write cfg 1 dumps/dump BCC");
    wr->init();
    wr->write_file(0);
    wr->fin();
    delete wr;
    int jatm;
    int bin_no;
    int* count=new int[nbins];
    type0* polar_plot=new type0[nbins];
    type0 r,p;
    type0 dr=forcefield->max_cut()/static_cast<type0>(nbins);
    for(int i=0;i<nbins;i++)
    {
        count[i]=0;
        polar_plot[i]=0.0;
    }
    for(int iatm=0;iatm<natms;iatm++)
    {
        if(!polar_xst[iatm]) continue;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            
            r=sqrt(XMatrixVector::rsq<__dim__>(x+__dim__*iatm,x+__dim__*jatm));
            if(!polar_xst[jatm]) continue;
            bin_no=static_cast<int>(r/dr);
            count[bin_no]++;
            p=fabs(XMatrixVector::V_V<__dim__>(polar+__dim__*iatm,polar+__dim__*jatm));
            
            
            
            polar_plot[bin_no]+=p;
        }
    }
    
    delete polarity_xst;
    delete polarity;
    dynamic->fin();
    delete dynamic;
    
    int* count_glbl=new int[nbins];
    type0* polar_plot_glbl=new type0[nbins];
    
    MPI_Allreduce(polar_plot,polar_plot_glbl,nbins,MPI_TYPE0,MPI_SUM,comm->world);
    MPI_Allreduce(count,count_glbl,nbins,MPI_INT,MPI_SUM,comm->world);
    
    delete [] polar_plot;
    delete [] count;
    
    for(int i=0;i<nbins;i++)
    {
        if(count_glbl[i])
            polar_plot_glbl[i]/=static_cast<type0>(count_glbl[i]);
        else
            polar_plot_glbl[i]=0.0;
    }
    delete [] count_glbl;
    
    if(comm->my_p==0)
        for(int i=0;i<nbins;i++)
            printf("%lf %lf\n",dr*i,polar_plot_glbl[i]);
        
    
    delete [] polar_plot_glbl;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Command_polarity::~Command_polarity()
{
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline int Command_polarity::calc_polarity(type0* x,int iatm,int* lst,int sz,type0* polarity)
{
    
    if(sz<14) return 0;
    
    auto calc_rsq=
    [&x,&lst](int& i,int& j)->type0
    {
        return XMatrixVector::rsq<__dim__>(x+__dim__*lst[i],x+__dim__*lst[j]);
    };
    
    
    
    
    
    
    
    int lcl_neigh [14][6];
    type0 rsqs[14];
    type0 max_rsq[14];
    int key[14];
    
    
    for(int i=0;i<14;i++)
    {
        for(int j=0;j<14;j++)
        {
            key[j]=j;
            rsqs[j]=calc_rsq(i,j);
        }
        
        XMath::quicksort(key,key+14,
        [&rsqs](int* ikey,int* jkey){return (rsqs[*ikey] < rsqs[*jkey]); },
        [](int* ikey,int* jkey){std::swap(*ikey,*jkey);});
        
        memcpy(lcl_neigh[i],key+1,6*sizeof(int));
        max_rsq[i]=rsqs[key[6]];
    }
    
    for(int i=0;i<14;i++)
        key[i]=i;
    
    XMath::quicksort(key,key+14,
    [&max_rsq](int* ikey,int* jkey){return (max_rsq[*ikey] < max_rsq[*jkey]); },
    [](int* ikey,int* jkey){std::swap(*ikey,*jkey);});
    
    
    enum{H,S};
    int type[14];
    for(int i=0;i<8;i++)
        type[key[i]]=H;
    
    for(int i=8;i<14;i++)
        type[key[i]]=S;
    

    
    for(int i=0;i<14;i++)
    {
        if(type[i]==H)
        {
            int H_count=0;
            int S_count=0;
            for(int j=0;j<6;j++)
            {
                if(type[lcl_neigh[i][j]]==H)
                    H_count++;
                else
                    S_count++;
            }
            if(H_count!=3 || S_count!=3)
                return 0;
        }
        else
        {
            int H_count=0;
            for(int j=0;j<4;j++)
                if(type[lcl_neigh[i][j]]==H)
                    H_count++;
            
            if(H_count!=4)
                return 0;
        }
    }
    /*-------------------------------
     
     
         yz----------xyz
        /.            /|
       / .           / |
      /  .          /  |
     z------------zx   |
     |   .         |   |
     |   .         |   |
     |   y.........|..xy
     |  .          |  /
     | .           | /
     |.            |/
     o-------------x
     
     0:o    [-1,-1,-1]
     1:x    [+1,-1,-1]
     2:y    [-1,+1,-1]
     3:z    [-1,-1,+1]
     4:xyz  [+1,+1,+1]
     5:yz   [-1,+1,+1]
     6:zx   [+1,-1,+1]
     7:xy   [+1,+1,-1]
     -------------------------------*/
    
    
    
    int new_arrang[14];
    
    new_arrang[0]=key[0];
    int ic=1;
    for(int i=0,j;i<6;i++)
    {
        j=lcl_neigh[new_arrang[0]][i];
        if(type[j]==H)
            new_arrang[ic++]=j;
    }
    
    type0 HH[3][3];
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            HH[i][j]=x[lst[new_arrang[i+1]]*__dim__+j]-x[lst[new_arrang[0]]*__dim__+j];
        
    if(M3DET(HH)<0.0)
        std::swap(new_arrang[2],new_arrang[3]);
    
    ic=6;
    for(int i=0,j;i<6;i++)
    {
        j=lcl_neigh[new_arrang[1]][i];
        if(type[j]==H && j!=new_arrang[0])
            new_arrang[ic++]=j;
    }
    
    for(int i=0,j;i<6;i++)
    {
        j=lcl_neigh[new_arrang[2]][i];
        if(type[j]!=H)
            continue;
        
        if(new_arrang[6]==j)
            std::swap(new_arrang[6],new_arrang[7]);
        if(new_arrang[0]!=j && new_arrang[6]!=j && new_arrang[7]!=j)
            new_arrang[5]=j;
    }
    
    new_arrang[4]=-1;
    for(int i=0;i<14 && new_arrang[4]==-1;i++)
    {
        if(type[i]!=H) continue;
        if(i==new_arrang[0]) continue;
        if(i==new_arrang[1]) continue;
        if(i==new_arrang[2]) continue;
        if(i==new_arrang[3]) continue;
        if(i==new_arrang[5]) continue;
        if(i==new_arrang[6]) continue;
        if(i==new_arrang[7]) continue;
        new_arrang[4]=i;
    }
    
    
    
    type0 F[3][3];
    type0 a[4];
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<4;j++)
            a[j]=x[lst[new_arrang[j+4]]*__dim__+i]-x[lst[new_arrang[j]]*__dim__+i];
        
        F[i][0]=0.125*(a[0]-a[1]+a[2]+a[3]);
        F[i][1]=0.125*(a[0]+a[1]-a[2]+a[3]);
        F[i][2]=0.125*(a[0]+a[1]+a[2]-a[3]);
        
    }
    
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        {
            HH[i][j]=0.0;
            for(int k=0;k<3;k++)
                HH[i][j]+=F[i][k]*F[j][k];
        }
    type0 p[3];
    
    /*
    if(atoms->my_p==0 && iatm<10)
    {
        printf("{");
        for(int i=0;i<3;i++)
        {
            printf("{%lf",HH[i][0]);
            for(int j=1;j<3;j++)
                printf(",%lf",HH[i][j]);
            printf("}");
            if(i!=2)
                printf(",");
        }
        printf("}\n");
    }
    */
    test(HH,p);
    polarity[0]=p[0];
    polarity[1]=p[1];
    polarity[2]=p[2];
    
    /*
    if(atoms->my_p==0 && iatm<10)
    {
        printf("%lf | %lf %lf %lf \n",eig,p[0],p[1],p[2]);
        printf("----------------------------\n");

    }
    */
    /*
    
    */
    /*
    
    int ipair=0;
    int* Skey=key+8;
    type0 rsq_max_,rsq_;
    for(int i=0;i<6;i++)
    {
        rsq_max_=0.0;
        for(int j=i+1;j<6;j++)
        {
            rsq_=calc_rsq(Skey[i],Skey[j]);
            if(rsq_>rsq_max_)
            {
                rsq_max_=rsq_;
                ipair=j;
            }
        }
        
        rsqs[i/2]=rsq_max_;
        std::swap(Skey[i+1],Skey[ipair]);
        i++;
    }
    
    
    auto calc_polar=
    [&x,&lst](int& i,int& j,type0 rsq,type0*& pol)->void
    {
        type0 r=sqrt(rsq);
        for(int k=0;k<__dim__;k++)
        {
            pol[k]=(x[__dim__*lst[i]+k]-x[__dim__*lst[j]+k])/r;
        }
        
    };
    
   // printf("%lf %lf %lf\n",rsqs[0],rsqs[0],rsqs[0]);
    
    
    if(rsqs[0]>rsqs[1] && rsqs[0]>rsqs[2])
    {
        calc_polar(Skey[0],Skey[1],rsqs[0],polarity);
    }
    else if(rsqs[1]>rsqs[0] && rsqs[1]>rsqs[2])
    {
        calc_polar(Skey[2],Skey[3],rsqs[1],polarity);
    }
    else if(rsqs[2]>rsqs[1] && rsqs[2]>rsqs[0])
    {
        calc_polar(Skey[4],Skey[5],rsqs[2],polarity);
    }
    else
    {
        polarity[0]=polarity[1]=polarity[2]=0.0;
        return 0;
    }
    */
    return 1;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline type0 Command_polarity::test(type0 (&A)[3][3],type0 (&polar)[3])
{
    type0 p1=A[1][0]*A[1][0]+A[2][0]*A[2][0]+A[2][1]*A[2][1];
    
    
    if(p1==0.0)
    {
        polar[0]=polar[1]=polar[2]=0.0;
        if(A[0][0]>A[1][1] && A[0][0]>A[2][2])
        {
            polar[0]=1.0;
            return A[0][0];
        }
        else if(A[1][1]>A[0][0] && A[1][1]>A[2][2])
        {
            polar[1]=1.0;
            return A[1][1];
        }
        else if(A[2][2]>A[0][0] && A[2][2]>A[1][1])
        {
            polar[2]=1.0;
            return A[2][2];
        }
        
        
        
    }
    type0 eig[3];
    
    type0 q=(A[0][0]+A[1][1]+A[2][2])/3.0;
    type0 p2=(A[0][0]-q)*(A[0][0]-q)+(A[1][1]-q)*(A[1][1]-q)+(A[2][2]-q)*(A[2][2]-q)+2.0*p1;
    type0 p=sqrt(p2/6.0);
    type0 r=0.5*(
                 (A[0][0]-q)*(A[1][1]-q)*(A[2][2]-q)
                 +2.0*A[2][0]*A[1][0]*A[2][1]
                 -(A[0][0]-q)*A[2][1]*A[2][1]
                 -(A[1][1]-q)*A[2][0]*A[2][0]
                 -(A[2][2]-q)*A[1][0]*A[1][0])/(p*p*p);
    
    type0 phi;
    if(r<=-1.0)
        phi=M_PI/3.0;
    else if(r>=1.0)
        phi=0.0;
    else
        phi=acos(r)/3.0;
    
    eig[0]=q+2.0*p*cos(phi);
    eig[2]=q+2.0*p*cos(phi+2.0*M_PI/3.0);
    eig[1]=3.0*q-eig[0]-eig[2];
    
    
    auto inner=
    [] (type0(&x)[3],type0(&y)[3])->type0
    {
        
        return (x[0]*y[0]+x[1]*y[1]+x[2]*y[2]);
        
    };
    
    auto cross=
    [] (type0(&a)[3],type0(&b)[3],type0(&c)[3])->void
    {
        c[0]=a[1]*b[2]-a[2]*b[1];
        c[1]=a[2]*b[0]-a[0]*b[2];
        c[2]=a[0]*b[1]-a[1]*b[0];
    };
    
    auto normalize=
    [] (type0(&x)[3])->void
    {
        type0 r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        x[0]/=r;
        x[1]/=r;
        x[2]/=r;
    };
    
    
    
    

    A[0][0]-=eig[0];
    A[1][1]-=eig[0];
    A[2][2]-=eig[0];
    
    type0 prods[3];
    normalize(A[0]);
    normalize(A[1]);
    normalize(A[2]);
    prods[0]=fabs(inner(A[1],A[2]));
    prods[1]=fabs(inner(A[2],A[0]));
    prods[2]=fabs(inner(A[0],A[1]));
    
    if(prods[0]<prods[1] && prods[0]<prods[2])
    {
        cross(A[1],A[2],polar);
    }
    else if(prods[1]<prods[0] && prods[1]<prods[2])
    {
        cross(A[2],A[0],polar);
    }
    else
    {
        cross(A[0],A[1],polar);
    }
    
    normalize(polar);
    
    return eig[0];
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Command_polarity::QR(type0(&A)[3][3],type0(&Q)[3][3],type0(&R)[3][3])
{
    auto normalize=
    [] (type0(&x)[3])->void
    {
        type0 r=0.0;
        r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        x[0]/=r;
        x[1]/=r;
        x[2]/=r;
    };
    
    auto inner=
    [] (type0(&x)[3],type0(&y)[3])->type0
    {
        
        return (x[0]*y[0]+x[1]*y[1]+x[2]*y[2]);
        
    };

    
    memcpy(Q[0],A[0],3*sizeof(type0));
    normalize(Q[0]);
    
    
    type0 m,n=inner(A[1],Q[0]);
    for(int i=0;i<3;i++)
        Q[1][i]=A[1][i]-n*Q[0][i];
    normalize(Q[1]);
    
    n=inner(A[2],Q[0]);
    m=inner(A[2],Q[1]);
    for(int i=0;i<3;i++)
        Q[2][i]=A[2][i]-n*Q[0][i]-m*Q[1][i];
    normalize(Q[2]);
    
    for(int i=0;i<3;i++)
        for(int j=i;j<3;j++)
        {
            R[i][j]=0.0;
            for(int k=0;k<3;k++)
                R[i][j]+=Q[k][i]*A[k][j];
        }
    
    R[1][0]=R[2][0]=R[2][1]=0.0;
    
    
}
