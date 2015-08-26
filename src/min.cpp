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
    
    int dim=atoms->dimension;
    if(dim)
    {
        CREATE1D(nrgy_strss,dim*(dim+1)/2+1);
        ns_alloc=1;
    }
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min::~Min()
{
    delete thermo;
    if(ns_alloc)
        delete [] nrgy_strss;
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
            fprintf(output,"line search failed: reached minimum alpha\n");
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
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::rectify_f(type0* f)
{
    if(dof_n==-1)
        return;
    atoms->vectors[dof_n]->ret(dof);
    int tot=(atoms->natms)*(atoms->vectors[0]->dim);
    for(int i=0;i<tot;i++) if(dof[i]==1) f[i]=0.0;
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
    //ls=new LineSearch_goldensection(mapp);
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

    if(chng_box)
    {
        ls->h_H=h_H;
        ls->f_H=f_H;
        ls->H_prev=H_prev;
        ls->B_prev=B_prev;
    }
}

