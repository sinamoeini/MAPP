/*--------------------------------------------
 Created by Sina on 2/1/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#include "min.h"
#include "ls_bt.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min::Min(MAPP* mapp):InitPtrs(mapp)
{
    ns_alloc=0;
    if(forcefield==NULL)
        error->abort("ff should be "
        "initiated before min");
    
    
    chng_box=0;
    dim=atoms->dimension;
    x_dim=atoms->vectors[0]->dim;
    err=LS_S;
    
    char** args;
    
    int narg;
    if(mapp->mode==MD_mode)
        narg=mapp->parse_line((char*)
        "PE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    else if (mapp->mode==DMD_mode)
        narg=mapp->parse_line((char*)
        "FE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    
    pe_idx=0;
    stress_idx=1;
    thermo=new ThermoDynamics(mapp,narg,args);
    for(int i=0;i<narg;i++)
        delete [] args[i];
    delete [] args;
    
    // default values 
    max_iter=10000;
    energy_tolerance=1.0e-8;
    affine=0;

    int dim=atoms->dimension;
    if(dim)
    {
        CREATE1D(nrgy_strss,dim*(dim+1)/2+1);
        ns_alloc=1;
    }
    CREATE2D(N,dim,dim);
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min::~Min()
{
    delete thermo;
    if(ns_alloc)
        delete [] nrgy_strss;
    for(int i=0;i<dim;i++)
        delete [] N[i];
    if(dim)
        delete [] N;
}
/*--------------------------------------------
 error messages
 --------------------------------------------*/
void Min::print_error()
{
    
    if(err==LS_F_DOWNHILL)
    {
        if(atoms->my_p_no==0)
            fprintf(output,"line search failed: not downhill direction\n");
    }
    else if(err==LS_F_GRAD0)
    {
        if(atoms->my_p_no==0)
            fprintf(output,"line search failed: gradient is zero\n");
    }
    else if(err==LS_MIN_ALPHA)
    {
        if(atoms->my_p_no==0)
            fprintf(output,"line search failed: minimum alpha reached\n");
    }
    else if (err==MIN_S_TOLERANCE)
    {
        if(atoms->my_p_no==0)
            fprintf(output,"minimization finished: energy tolerance reached\n");
    }
    else if (err==MIN_F_MAX_ITER)
    {
        if(atoms->my_p_no==0)
            fprintf(output,"minimization finished: maximum iteration reached\n");
    }
    else if (err==B_F_DOWNHILL)
    {
        if(atoms->my_p_no==0)
            fprintf(output,"bracketing failed: not downhill direction\n");
    }
    else if (err==B_F_MAX_ALPHA)
    {
        if(atoms->my_p_no==0)
            fprintf(output,"bracketing failed: maximum alpha reached\n");
    }

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::rectify(type0* a)
{
    if(dof_n==-1 && affine==0)
        return;
    else if(dof_n!=-1 && affine==0)
    {
        atoms->vectors[dof_n]->ret(dof);
        for(int i=0;i<atoms->natms*x_dim;i++)
            if(dof[i]==1) a[i]=0.0;
    }
    else if(dof_n==-1 && affine)
    {
        for(int i=0;i<atoms->natms;i++)
            for(int j=0;j<dim;j++)
                a[i*x_dim+j]=0.0;
    }
    else if(dof_n!=-1 && affine)
    {
        atoms->vectors[dof_n]->ret(dof);
        for(int i=0;i<atoms->natms;i++)
        {
            for(int j=0;j<dim;j++)
                a[i*x_dim+j]=0.0;
            for(int j=dim;j<x_dim;j++)
                if(dof[i*x_dim+j]==1) a[i*x_dim+j]=0.0;
        }

    }

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::reg_h_H(type0** T)
{
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            T[i][j]=0.0;
    
    if(H_dof[0][0])
        T[0][0]=nrgy_strss[1];
    if(H_dof[1][1])
        T[1][1]=nrgy_strss[2];
    if(H_dof[2][2])
        T[2][2]=nrgy_strss[3];
    if(H_dof[2][1])
        T[2][1]=nrgy_strss[4];
    if(H_dof[2][0])
        T[2][0]=nrgy_strss[5];
    if(H_dof[1][0])
        T[1][0]=nrgy_strss[6];

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::init_linesearch()
{
    if(ls==NULL)
    {
        char** arg=NULL;
        ls=new LineSearch_backtrack(mapp,0,arg);
    }
    ls->vecs_comm=vecs_comm;
    
    ls->dim=dim;
    ls->x_dim=x_dim;

    ls->x_prev_n=x_prev_n;
    ls->h_n=h_n;
    ls->f_n=f_n;
    ls->chng_box=chng_box;
    ls->nrgy_strss=nrgy_strss;
    ls->prev_val=-1.0;
    ls->affine=affine;

    if(chng_box)
    {
        ls->h_H=h_H;
        ls->f_H=f_H;
        ls->H_prev=H_prev;
        ls->B_prev=B_prev;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::prepare_affine_h(type0* x,type0* h)
{
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
        {
            N[i][j]=0.0;
            if(chng_box)
                for(int k=0;k<dim;k++)
                    N[i][j]+=B_prev[j][k]*h_H[k][i];
        }
    
    for(int iatm=0;iatm<atoms->natms;iatm++)
    {
        for(int j=0;j<dim;j++)
        {
            h[iatm*x_dim+j]=0.0;
            for(int k=0;k<dim;k++)
                h[iatm*x_dim+j]+=N[j][k]*x[iatm*x_dim+k];
        }
    }
}

