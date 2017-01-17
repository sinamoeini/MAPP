//
//  comm.cpp
//  MAPP
//
//  Created by Sina on 11/21/16.
//  Copyright © 2016 MIT. All rights reserved.
//

#include "comm.h"
#include "atoms.h"
#include "exchange.h"
#include "xmath.h"
#include "error.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
Communication::Communication(MPI_Comm& __world):
world(__world),
my_p(get_rank(__world)),
tot_p(get_size(__world)),
skin(0.5),
xchng_id(0)
{
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
Communication::~Communication()
{
    delete [] *n_p_grid;
    delete [] n_p_grid;
    delete [] p_per_n;
}
/*--------------------------------------------
 automatically grid the domain
 --------------------------------------------*/
void Communication::auto_grid(type0(&H)[__dim__][__dim__])
{
    bool eq_p_per_n=true;
    for(int i=0;i<tot_n;i++)
        if(p_per_n[i]!=p_per_n[0])
            eq_p_per_n=false;
    
    int* fac_list;
    int fac_list_size;
    XMath::fac_list(tot_p,__dim__,fac_list,fac_list_size);
    
    type0 ratio=-1.0;
    type0 tmp_ratio;
    
    if(eq_p_per_n && tot_n>1)
    {
        int prin_dimension=0;
        type0 min_area=1.0/H[0][0];
        
        for(int i=1;i<__dim__;i++)
            if(1.0/H[i][i]<min_area)
            {
                min_area=1.0/H[i][i];
                prin_dimension=i;
            }
        
        int* ifac_list=fac_list;
        for(int ifac=0;ifac<fac_list_size;ifac++,ifac_list+=__dim__)
        {
            if(ifac_list[prin_dimension]!=tot_n) continue;
            
            tmp_ratio=0.0;
            for(int i=0;i<__dim__;i++)
                tmp_ratio+=static_cast<type0>(ifac_list[i])/H[i][i];
            
            if(tmp_ratio<ratio||ratio<0.0)
            {
                ratio=tmp_ratio;
                for(int i=0;i<__dim__;i++)
                    tot_p_grid[i]=ifac_list[i];
            }
            
        }
        
        delete [] fac_list;
        
        
        if(__dim__>1)
        {
            int* tmp_tot_p_grid=new int[__dim__-1];
            int* dim_indx=new int[__dim__-1];
            int* my_tmp_loc=new int[__dim__-1];
            int* my_tmp_denom=new int[__dim__-1];
            
            
            int pos=0;
            
            for(int i=0;i<__dim__;i++)
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
            
            for(int i=0;i<__dim__-1;i++)
            {
                int no=1;
                for(int j=0;j<i;j++)
                    no*=tmp_tot_p_grid[j];
                my_tmp_denom[i]=no;
            }
            
            
            int ttmp=my_p_in_my_node;
            for(int i=__dim__-2;i>-1;i--)
            {
                my_tmp_loc[i]=ttmp/my_tmp_denom[i];
                ttmp-=my_tmp_loc[i]*my_tmp_denom[i];
            }
            
            
            for(int i=0;i<__dim__-1;i++)
            {
                my_loc[dim_indx[i]]=my_tmp_loc[i];
            }
            my_loc[prin_dimension]=my_n;
            
            
            
            int nop;
            int nom;
            int tmp0,tmp1;
            for(int i=0;i<__dim__-1;i++)
            {
                
                tmp0=my_tmp_loc[i]-1;
                if(tmp0==-1)
                    tmp0=tmp_tot_p_grid[i]-1;
                
                tmp1=my_tmp_loc[i]+1;
                if(tmp1==tmp_tot_p_grid[i])
                    tmp1=0;
                
                nop=0;
                nom=0;
                
                for(int j=0;j<__dim__-1;j++)
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
        for(int ifac=0;ifac<fac_list_size;ifac++,ifac_list+=__dim__)
        {
            tmp_ratio=0.0;
            for(int i=0;i<__dim__;i++)
                tmp_ratio+=static_cast<type0>(ifac_list[i])/H[i][i];
            
            if(tmp_ratio<ratio||ratio<0.0)
            {
                ratio=tmp_ratio;
                for(int i=0;i<__dim__;i++)
                    tot_p_grid[i]=ifac_list[i];
            }
            
        }
        
        delete [] fac_list;
        
        int* list=new int[__dim__];
        
        MPI_Comm cartesian;
        for(int i=0;i<__dim__;i++)
            list[i]=1;
        MPI_Cart_create(world,__dim__,tot_p_grid,list,1,&cartesian);
        MPI_Cart_get(cartesian,__dim__,tot_p_grid,list,my_loc);
        for(int i=0;i<__dim__;i++)
            MPI_Cart_shift(cartesian,i,1,&neigh_p[i][0],&neigh_p[i][1]);
        
        MPI_Comm_free(&cartesian);
        
        delete [] list;
    }
    
    
    
    for(int i=0;i<__dim__;i++)
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
void Communication::__man_grid__(int(&n)[__dim__])
{
    for(int i=0;i<__dim__;i++)
        tot_p_grid[i]=n[i];
    
    int* list=new int[__dim__];
    for(int i=0;i<__dim__;i++)
        list[i]=1;
    
    MPI_Comm cartesian;
    MPI_Cart_create(world,__dim__,tot_p_grid,list,1,&cartesian);
    MPI_Cart_get(cartesian,__dim__,tot_p_grid,list,my_loc);
    for(int i=0;i<__dim__;i++)
        MPI_Cart_shift(cartesian,i,1,&neigh_p[i][0],&neigh_p[i][1]);
    
    MPI_Comm_free(&cartesian);
    
    delete [] list;
    
    for(int i=0;i<__dim__;i++)
    {
        s_lo[i]=
        static_cast<type0>(my_loc[i])
        /static_cast<type0>(tot_p_grid[i]);
        
        s_hi[i]=
        (static_cast<type0>(my_loc[i])+1.0)
        /static_cast<type0>(tot_p_grid[i]);
    }
}
/*--------------------------------------------
 manually grid the domain
 --------------------------------------------*/
void Communication::man_grid(int (&n)[__dim__])
{
     for(int i=0;i<__dim__;i++)
        tot_p_grid[i]=n[i];
    
    int* list=new int[__dim__];
    for(int i=0;i<__dim__;i++)
        list[i]=1;
    
    MPI_Comm cartesian;
    MPI_Cart_create(world,__dim__,tot_p_grid,list,1,&cartesian);
    MPI_Cart_get(cartesian,__dim__,tot_p_grid,list,my_loc);
    for(int i=0;i<__dim__;i++)
        MPI_Cart_shift(cartesian,i,1,&neigh_p[i][0],&neigh_p[i][1]);
    
    MPI_Comm_free(&cartesian);
    
    delete [] list;
    
    for(int i=0;i<__dim__;i++)
    {
        s_lo[i]=
        static_cast<type0>(my_loc[i])
        /static_cast<type0>(tot_p_grid[i]);
        
        s_hi[i]=
        (static_cast<type0>(my_loc[i])+1.0)
        /static_cast<type0>(tot_p_grid[i]);
    }
    
    reset();
    
    grid_established=true;
    
    if(my_p==0)
    {
        fprintf(output,"manual grid performed: ");
        
        for(int i=0;i<__dim__-1;i++)
            fprintf(output,"%dx",comm->tot_p_grid[i]);
        
        fprintf(output,"%d",comm->tot_p_grid[__dim__-1]);
        
        fprintf(output,"\n");
    }
}
/*--------------------------------------------
 automatically grid the domain
 --------------------------------------------*/
void Communication::auto_grid()
{
    type0 (&H)[__dim__][__dim__]=atoms->H;
    type0 det=1.0;
    for(int idim=0;idim<__dim__;idim++)
        det*=H[idim][idim]*H[idim][idim];
    
    if(det==0.0)
        Error::abort("cannot perform automatic "
        "domain griding, the domain has not been established");

    comm->auto_grid(H);
    grid_established=true;
    
    if(my_p==0)
    {
        fprintf(output,"automatic grid performed: ");
        
        for(int i=0;i<__dim__-1;i++)
            fprintf(output,"%dx",comm->tot_p_grid[i]);
        
        fprintf(output,"%d",comm->tot_p_grid[__dim__-1]);
        
        fprintf(output,"\n");
    }
    
}
/*--------------------------------------------
 manually grid the domain
 --------------------------------------------*/
void Communication::reset()
{
    Exchange xchng(atoms,this,atoms->nvecs);
    atoms->x2s_lcl();
    xchng.full_xchng();
    atoms->s2x_lcl();
}

