#include "line_search.h"
#include "ff.h"
#include "atom_types.h"
#include "thermo_dynamics.h"
using namespace std;
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearch::LineSearch(MAPP* mapp,VecLst*vec_list)
: InitPtrs(mapp)
{
    x_n=atoms->find("x");
    x_prev_n=atoms->find("x_prev");
    f_n=atoms->find("f");
    dim=atoms->dimension;
    x_dim=atoms->vectors[0]->dim;
    vecs_comm=vec_list;
    dof_n=atoms->find_exist("dof");
    
    CREATE2D(N,dim,dim);
    CREATE2D(M,dim,dim);
    xmath=new XMath(mapp);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
LineSearch::~LineSearch()
{
    delete xmath;
    for(int i=0;i<dim;i++)
    {
        delete [] N[i];
        delete [] M[i];
    }
    delete [] N;
    delete [] M;
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 LineSearch::inner_f_h()
{
    type0* f;
    type0* h;
    atoms->vectors[f_n]->ret(f);
    atoms->vectors[h_n]->ret(h);
    type0 inner=0.0;
    type0 inner_tot=0.0;

    for(int i=0;i<x_dim*atoms->natms;i++)
        inner+=f[i]*h[i];
    MPI_Allreduce(&inner,&inner_tot,1,MPI_TYPE0,MPI_SUM,world);

    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                inner_tot+=h_H[i][j]*f_H[i][j];
        return inner_tot;
    }
    else
    {
        return inner_tot;
    }
    
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
void LineSearch::normalize_h()
{
    type0* h;
    atoms->vectors[h_n]->ret(h);
    type0 inner=0.0;
    type0 inner_tot=0.0;
    
    for(int i=0;i<x_dim*atoms->natms;i++)
        inner+=h[i]*h[i];
    MPI_Allreduce(&inner,&inner_tot,1,MPI_TYPE0,MPI_SUM,world);
    
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                inner_tot+=h_H[i][j]*h_H[i][j];
        
        inner_tot=s_max/sqrt(inner_tot);
        
        for(int i=0;i<x_dim*atoms->natms;i++)
            h[i]*=inner_tot;
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                h_H[i][j]*=inner_tot;
    }
    else
    {
        
        inner_tot=d_max/sqrt(inner_tot);
        
        for(int i=0;i<x_dim*atoms->natms;i++)
            h[i]*=inner_tot;
    }
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 LineSearch::energy(type0 alpha)
{
    type0 energy;
    type0* x;
    type0* x_prev;
    type0* h;
    char* dof=NULL;
    
    atoms->vectors[x_n]->ret(x);
    atoms->vectors[x_prev_n]->ret(x_prev);
    atoms->vectors[h_n]->ret(h);
    if(dof_n>-1)
        atoms->vectors[dof_n]->ret(dof);

    if(chng_box)
    {
        type0** H=atoms->H;
        
        for(int i=0;i<dim;i++)
        {
            for(int j=0;j<dim;j++)
            {
                H[i][j]=H_prev[i][j]+alpha*h_H[i][j];
                N[i][j]=0.0;
                M[i][j]=0.0;
                for(int k=0;k<dim;k++)
                {
                    N[i][j]+=alpha*h_H[i][k]*B_prev[k][j];
                    M[i][j]+=alpha*H_prev[i][k]*(H_prev[k][j]+alpha*h_H[k][j]);
                }
            }
        }
        for(int i=0;i<dim;i++)
            N[i][i]++;

        
        if(dim==3)
            M3INV_TRI_LOWER(atoms->H,atoms->B);
        else
            xmath->invert_lower_triangle(atoms->H,atoms->B,dim);
       
        int icomp=0;
        for(int i=0;i<atoms->natms;i++)
        {
            for(int j=0;j<dim;j++)
            {
                x[icomp+j]=0.0;
                for(int k=j;k<dim;k++)
                    x[icomp+j]+=N[k][j]*x_prev[icomp+k]+alpha*H[k][j]*h[icomp+k];
            }
            for(int j=dim;j<x_dim;j++)
                x[icomp+j]=x_prev[icomp+j]+alpha*h[icomp+j];
                
            icomp+=x_dim;
        }
        
        if(dof_n>-1)
            for(int i=0;i<(atoms->natms)*x_dim;i++)
                if(dof[i]==1) x[i]=x_prev[i];

        thermo->start_comm_time();
        atoms->update(1,vecs_comm);
        thermo->stop_comm_time();
        
        thermo->start_force_time();
        energy=forcefield->energy_calc();
        thermo->stop_force_time();
        
        return energy;
    }
    else
    {
        
        for(int i=0;i<x_dim*atoms->natms;i++)
            x[i]=x_prev[i]+alpha*h[i];

        
        thermo->start_comm_time();
        atoms->update(0,vecs_comm);
        thermo->stop_comm_time();
        
        thermo->start_force_time();
        energy=forcefield->energy_calc();
        thermo->stop_force_time();
        
        return energy;
    }
}
/*--------------------------------------------
 find maximum h
 --------------------------------------------*/
type0 LineSearch::find_max_h()
{
    type0* h;
    type0 max_h=0.0;
    type0 max_h_tot=0.0;
    atoms->vectors[h_n]->ret(h);
    
    for(int i=0;i<atoms->natms*x_dim;i++)
        max_h=MAX(max_h,fabs(h[i]));
    MPI_Allreduce(&max_h,&max_h_tot,1,MPI_TYPE0,MPI_MAX,world);
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                max_h_tot=MAX(max_h_tot,fabs(h_H[i][j]));
    }
    
    return max_h_tot;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearch_BackTrack::LineSearch_BackTrack(MAPP* mapp,VecLst*vec_list)
:LineSearch(mapp,vec_list)
{
    c=0.4;
    rho=0.5;
    alpha_max=1.0;
    alpha_min=0.0;
    d_max=1.0;
    s_max=0.1;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
LineSearch_BackTrack::~LineSearch_BackTrack()
{
    
}
/*--------------------------------------------
 minimize line
 --------------------------------------------*/
int LineSearch_BackTrack::line_min(type0& nrgy,type0& alph)
{
    type0 inner;
    type0 alpha_m;
    type0 max_h;
    type0 min_h;
    type0 min_h_tot=0;
    type0 current_energy,ideal_energy;
    
    if(mapp->mode==DMD_mode)
    {
        type0* x;
        atoms->vectors[0]->ret(x);
        min_h=INFINITY;
        for(int i=0;i<atoms->natms;i++)
            for(int j=0;j<atom_types->no_types;j++)
                min_h=MIN(min_h,x[i*(3+atom_types->no_types)+3+j]);
        MPI_Allreduce(&min_h,&min_h_tot,1,MPI_TYPE0,MPI_MIN,world);
        min_h_tot*=0.99;
        
        alpha_min=MIN(alpha_min,min_h_tot);
    }
    
    alph=0.0;
    normalize_h();
    inner=inner_f_h();
    
    if(inner<=0)
        return LS_F_DOWNHILL;
    
    max_h=find_max_h();
    if(max_h==0)
        return LS_F_GRAD0;
    
    alpha_m=MIN(alpha_max,max_h);
    if(alpha_m<=alpha_min)
        return LS_F_ALPHAMIN;
    
    while (1)
    {
        ideal_energy=nrgy-alpha_m*c*inner;
        current_energy=energy(alpha_m);
        if(current_energy<=ideal_energy)
        {
            nrgy=current_energy;
            alph=alpha_m;
            return LS_S;
        }
        alpha_m*=rho;
        
        if(alpha_m<=alpha_min)
        {
            nrgy=energy(0);
            return LS_F_ALPHAMIN;
        }
    }
    
    
}

/*--------------------------------------------
 minimize line
 --------------------------------------------*/
/*
void LineSearch::ccc(type0** A,type0**TT,type0** AAA,type0 alpha)
{
    type0** Q;
    type0** H_tmp;
    type0 babs;
    type0** del_H;
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            H_tmp[i][j]=H_prev[i][j]+alpha*A[i][j];
    
    Q[0][0]=H_tmp[0][0];
    Q[0][1]=H_tmp[0][1];
    Q[0][2]=H_tmp[0][2];
    babs=sqrt(Q[0][0]*Q[0][0]+Q[0][1]*Q[0][1]+Q[0][2]*Q[0][2]);
    Q[0][0]/=babs;
    Q[0][1]/=babs;
    Q[0][2]/=babs;
    Q[2][0]=Q[0][1]*H_tmp[1][2]-Q[0][2]*H_tmp[1][1];
    Q[2][1]=Q[0][2]*H_tmp[1][0]-Q[0][0]*H_tmp[1][2];
    Q[2][2]=Q[0][0]*H_tmp[1][1]-Q[0][1]*H_tmp[1][0];
    babs=sqrt(Q[2][0]*Q[2][0]+Q[2][1]*Q[2][1]+Q[2][2]*Q[2][2]);
    Q[2][0]/=babs;
    Q[2][1]/=babs;
    Q[2][2]/=babs;
    Q[1][0]=Q[0][2]*Q[2][1]-Q[0][1]*Q[2][2];
    Q[1][1]=Q[0][0]*Q[2][2]-Q[0][2]*Q[2][0];
    Q[1][2]=Q[0][1]*Q[2][0]-Q[0][0]*Q[2][1];
    
    
    
    
    del_H[0][1]=del_H[0][2]=del_H[1][2]=0.0;
    del_H[0][0]=V3MULV3(Q[0],H_tmp[0])-H_prev[0][0];
    del_H[1][0]=V3MULV3(Q[0],H_tmp[1])-H_prev[1][0];
    del_H[1][1]=V3MULV3(Q[1],H_tmp[1])-H_prev[1][1];
    del_H[2][0]=V3MULV3(Q[0],H_tmp[2])-H_prev[2][0];
    del_H[2][1]=V3MULV3(Q[1],H_tmp[2])-H_prev[2][1];
    del_H[2][2]=V3MULV3(Q[2],H_tmp[2])-H_prev[2][2];
    
}
*/



