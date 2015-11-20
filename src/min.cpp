/*--------------------------------------------
 Created by Sina on 2/1/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#include "min.h"
#include "error.h"
#include "memory.h"
#include "timer.h"
#include "neighbor.h"
#include "thermo_dynamics.h"
#include "ls_styles.h"
#define INV_SQ_2 0.7071067811865475
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
    x_dim=mapp->x->dim;
    err=LS_S;
    
    char** args;
    int nargs;
    if(mapp->mode==MD_mode)
        nargs=mapp->parse_line(
        "PE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    else if (mapp->mode==DMD_mode)
        nargs=mapp->parse_line(
        "FE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    
    pe_idx=0;
    stress_idx=1;
    
    thermo=new ThermoDynamics(mapp,nargs,args);
    for(int i=0;i<nargs;i++)
        delete [] args[i];
    if(nargs)
        delete [] args;
    
    // default values 
    max_iter=10000;
    energy_tolerance=1.0e-8;
    affine=0;

    if(dim)
    {
        CREATE1D(nrgy_strss,dim*(dim+1)/2+1);
        ns_alloc=1;
    }
    
    N=new type0*[dim];
    *N=new type0[dim*dim];
    H_dof=new int*[dim];
    *H_dof=new int[dim*dim];

    for(int i=1;i<dim;i++)
    {
        N[i]=N[i-1]+dim;
        H_dof[i]=H_dof[i-1]+dim;
    }

    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            H_dof[i][j]=0;
    
    xmath=new XMath();
        
    if(mapp->mode==DMD_mode)
        forcefield_dmd=dynamic_cast<ForceFieldDMD*>(forcefield);
    max_dx=1.0;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min::~Min()
{
    delete thermo;
    if(ns_alloc)
        delete [] nrgy_strss;

    if(dim)
    {
        delete [] *H_dof;
        delete [] H_dof;
        delete [] *N;
        delete [] N;
    }
    
    delete xmath;
}
/*--------------------------------------------
 error messages
 --------------------------------------------*/
void Min::print_error()
{
    
    if(err==LS_F_DOWNHILL)
    {
        if(atoms->my_p==0)
            fprintf(output,"line search failed: not downhill direction\n");
    }
    else if(err==LS_F_GRAD0)
    {
        if(atoms->my_p==0)
            fprintf(output,"line search failed: gradient is zero\n");
    }
    else if(err==LS_MIN_ALPHA)
    {
        if(atoms->my_p==0)
            fprintf(output,"line search failed: minimum alpha reached\n");
    }
    else if (err==MIN_S_TOLERANCE)
    {
        if(atoms->my_p==0)
            fprintf(output,"minimization finished: energy tolerance reached\n");
    }
    else if (err==MIN_F_MAX_ITER)
    {
        if(atoms->my_p==0)
            fprintf(output,"minimization finished: maximum iteration reached\n");
    }
    else if (err==B_F_DOWNHILL)
    {
        if(atoms->my_p==0)
            fprintf(output,"bracketing failed: not downhill direction\n");
    }
    else if (err==B_F_MAX_ALPHA)
    {
        if(atoms->my_p==0)
            fprintf(output,"bracketing failed: maximum alpha reached\n");
    }

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::rectify(type0* a)
{
    if(mapp->dof==NULL && affine==0)
        return;
    else if(mapp->dof!=NULL && affine==0)
    {
        byte* dof=mapp->dof->begin();
        for(int i=0;i<atoms->natms*x_dim;i++)
            if(dof[i]==1) a[i]=0.0;
    }
    else if(mapp->dof==NULL && affine)
    {
        for(int i=0;i<atoms->natms;i++)
            for(int j=0;j<dim;j++)
                a[i*x_dim+j]=0.0;
    }
    else if(mapp->dof==NULL && affine)
    {
        byte* dof=mapp->dof->begin();
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
void Min::force_calc()
{
    type0* f=mapp->f->begin();
    if(chng_box && !affine)
        for(int i=0;i<x_dim*(atoms->natms+atoms->natms_ph);i++)
            f[i]=0.0;
    else
        for(int i=0;i<x_dim*atoms->natms;i++)
            f[i]=0.0;
    if(chng_box)
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                f_H[i][j]=0.0;

    forcefield->force_calc_timer(sts_flag,nrgy_strss);
    
    if(mapp->dof!=NULL && affine==0)
    {
        byte* dof=mapp->dof->begin();
        for(int i=0;i<atoms->natms*x_dim;i++)
            if(dof[i]) f[i]=0.0;
    }
    else if(mapp->dof==NULL && affine)
    {
        for(int i=0;i<atoms->natms;i++)
            for(int j=0;j<dim;j++)
                f[i*x_dim+j]=0.0;
    }
    else if(mapp->dof==NULL && affine)
    {
        byte* dof=mapp->dof->begin();
        for(int i=0;i<atoms->natms;i++)
        {
            for(int j=0;j<dim;j++)
                f[i*x_dim+j]=0.0;
            for(int j=dim;j<x_dim;j++)
                if(dof[i*x_dim+j]==1) f[i*x_dim+j]=0.0;
        }
    }
    
    if(!chng_box)
        return;
    
    if(affine)
    {
        type0** B=atoms->B;
        type0** H=atoms->H;
        if(H_dof[0][0])
            f_H[0][0]=-(nrgy_strss[1]*B[0][0]+nrgy_strss[6]*B[1][0]+nrgy_strss[5]*B[2][0])*H[0][0]*H[1][1]*H[2][2];
        if(H_dof[1][1])
            f_H[1][1]=-(nrgy_strss[2]*B[1][1]+nrgy_strss[4]*B[2][1])*H[0][0]*H[1][1]*H[2][2];
        if(H_dof[2][2])
            f_H[2][2]=-(nrgy_strss[3]*B[2][2])*H[0][0]*H[1][1]*H[2][2];
        if(H_dof[2][1])
            f_H[2][1]=-(nrgy_strss[4]*B[2][2])*H[0][0]*H[1][1]*H[2][2];
        if(H_dof[2][0])
            f_H[2][0]=-(nrgy_strss[6]*B[1][1]+nrgy_strss[5]*B[2][1])*H[0][0]*H[1][1]*H[2][2];
        if(H_dof[1][0])
            f_H[1][0]=-(nrgy_strss[5]*B[2][2])*H[0][0]*H[1][1]*H[2][2];
    }
    else
    {
        type0* f=mapp->f->begin();
        int tot_natms=atoms->natms+atoms->natms_ph;
        int dimension=atoms->dimension;
        type0* x=mapp->x->begin();
        int x_dim=mapp->x->dim;
        type0 xtmp[3];
        type0 image[3];
        type0 st_lcl[6];
        for(int i=0;i<6;i++)
            st_lcl[i]=0.0;
        type0 st[6];
        type0** H=atoms->H;
        
        
        for(int i=0;i<tot_natms;i++)
        {
            for(int j=0;j<dimension;j++)
                xtmp[j]=x[j];
            
            for(int j=dimension-1;j>-1;j--)
            {
                image[j]=0.0;
                
                while(xtmp[j]>=H[j][j])
                {
                    for(int k=0;k<j+1;k++)
                        xtmp[k]-=H[j][k];
                    image[j]++;
                }
                
                while(xtmp[j]<0.0)
                {
                    for(int k=0;k<j+1;k++)
                        xtmp[k]+=H[j][k];
                    image[j]--;
                }
            }
            
            
            st_lcl[0]+=f[0]*image[0];
            st_lcl[1]+=f[1]*image[1];
            st_lcl[2]+=f[2]*image[2];
            st_lcl[3]+=f[1]*image[2];
            st_lcl[4]+=f[0]*image[2];
            st_lcl[5]+=f[0]*image[1];
            
            x+=x_dim;
            f+=x_dim;
        }
        
        MPI_Allreduce(st_lcl,st,6,MPI_TYPE0,MPI_SUM,world);
        
        if(H_dof[0][0])
            f_H[0][0]=0.5*st[0];
        if(H_dof[1][1])
            f_H[1][1]=0.5*st[1];
        if(H_dof[2][2])
            f_H[2][2]=0.5*st[2];
        if(H_dof[2][1])
            f_H[2][1]=0.5*st[3];
        if(H_dof[2][0])
            f_H[2][0]=0.5*st[4];
        if(H_dof[1][0])
            f_H[1][0]=0.5*st[5];
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::prepare_affine_h(type0* x,type0* h)
{
    type0** B=atoms->B;
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
        {
            N[i][j]=0.0;
            if(chng_box)
                for(int k=0;k<dim;k++)
                    N[i][j]+=B[j][k]*h_H[k][i];
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
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::init()
{
    // determine if the there is chng_box
    chng_box=0;
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            if(H_dof[i][j])
                chng_box=1;
    
    sts_flag=1;
    if(chng_box && affine)
        sts_flag=1;
    else if(chng_box && !affine)
        sts_flag=2;
    
    if(chng_box)
    {
        h_H=new type0*[dim];
        f_H=new type0*[dim];
        H_prev=new type0*[dim];
        f_H_prev=new type0*[dim];
        
        *h_H=new type0[dim*dim];
        *f_H=new type0[dim*dim];
        *H_prev=new type0[dim*dim];
        *f_H_prev=new type0[dim*dim];
        
        for(int i=1;i<dim;i++)
        {
            h_H[i]=h_H[i-1]+dim;
            f_H[i]=f_H[i-1]+dim;
            H_prev[i]=H_prev[i-1]+dim;
            f_H_prev[i]=f_H_prev[i-1]+dim;
        }
    }
    
    if(mapp->f==NULL)
        mapp->f=new Vec<type0>(atoms,x_dim);
    
    x_prev_ptr=new Vec<type0>(atoms,x_dim);
    h_ptr=new Vec<type0>(atoms,x_dim);
    f_prev_ptr=new Vec<type0>(atoms,x_dim);

    
    vecs_comm=new VecLst(atoms);

    if(mapp->mode==MD_mode)
    {
        vecs_comm->add_updt(mapp->type);
        if(mapp->dof!=NULL)
            vecs_comm->add_xchng(mapp->dof);
    }
    else if(mapp->mode==DMD_mode)
    {
        vecs_comm->add_updt(mapp->c);
        vecs_comm->add_updt(mapp->ctype);
        if(mapp->dof!=NULL)
        {
            vecs_comm->add_xchng(mapp->dof);
            if(mapp->cdof!=NULL)
                vecs_comm->add_xchng(mapp->cdof);
        }
    }
    
    //vecs_comm->add_xchng(mapp->f);
    vecs_comm->add_xchng(h_ptr);    
    vecs_comm->add_xchng(x_prev_ptr);
    vecs_comm->add_xchng(f_prev_ptr);
    
    if(mapp->x_d!=NULL)
        vecs_comm->add_arch(mapp->x_d);
    
    if(ls==NULL)
    {
        char** args=NULL;
        ls=new LineSearch_backtrack<Min>(mapp,0,args);
    }
    
    ls->init(this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::fin()
{
    delete f_prev_ptr;
    delete h_ptr;
    delete x_prev_ptr;
    delete vecs_comm;
    timer->print_stats();
    neighbor->print_stats();
    
    if(chng_box)
    {
        delete [] *f_H;
        delete [] f_H;

        delete [] *h_H;
        delete [] h_H;

        delete [] *H_prev;
        delete [] H_prev;
        
        delete [] *f_H_prev;
        delete [] f_H_prev;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::zero_f()
{
    type0* f=mapp->f->begin();
    if(chng_box && ! affine)
        for(int i=0;i<x_dim*(atoms->natms+atoms->natms_ph);i++)
            f[i]=0.0;
    else
        for(int i=0;i<x_dim*atoms->natms;i++)
            f[i]=0.0;
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 Min::F(type0 alpha)
{
    type0* x=mapp->x->begin();
    type0* h=h_ptr->begin();
    type0* x_prev=x_prev_ptr->begin();
    
    for(int i=0;i<x_dim*atoms->natms;i++)
        x[i]=x_prev[i]+alpha*h[i];
    
    if(chng_box)
    {
        type0** H=atoms->H;
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
            {
                H[i][j]=H_prev[i][j]+alpha*h_H[i][j];
            }
        
        if(dim==3)
            M3INV_TRI_LOWER(atoms->H,atoms->B);
        else
            xmath->invert_lower_triangle(atoms->H,atoms->B,dim);
    }
    
    atoms->update(mapp->x);
    return forcefield->energy_calc_timer();
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 Min::dF(type0 alpha,type0& drev)
{
    type0* x=mapp->x->begin();
    type0* h=h_ptr->begin();
    type0* x_prev=x_prev_ptr->begin();
    type0* f=mapp->f->begin();
    
    type0 drev_lcl;
    
    for(int i=0;i<x_dim*atoms->natms;i++)
        x[i]=x_prev[i]+alpha*h[i];
    
    if(chng_box)
    {
        type0** H=atoms->H;
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
            {
                H[i][j]=H_prev[i][j]+alpha*h_H[i][j];
            }
        
        if(dim==3)
            M3INV_TRI_LOWER(atoms->H,atoms->B);
        else
            xmath->invert_lower_triangle(atoms->H,atoms->B,dim);
    }
    
    atoms->update(mapp->x);
    
    force_calc();
    
    drev_lcl=0.0;
    if(affine==0)
    {
        for(int i=0;i<atoms->natms*x_dim;i++)
            drev_lcl+=h[i]*f[i];
    }
    else
    {
        for(int i=dim;i<atoms->natms*x_dim;i+=dim)
        {
            for(int j=dim;j<x_dim;j++)
            {
                drev_lcl+=h[i]*f[i];
                i++;
            }
        }
    }
    MPI_Allreduce(&drev_lcl,&drev,1,MPI_TYPE0,MPI_SUM,world);
    
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
            {
                drev+=h_H[i][j]*f_H[i][j];
            }
    }
    drev*=-1.0;
    
    return nrgy_strss[0];
}
/*--------------------------------------------
 find maximum h
 --------------------------------------------*/
void Min::ls_prep(type0& dfa,type0& h_norm,type0& max_a)
{
    type0* h=h_ptr->begin();
    type0* f=mapp->f->begin();
    type0 max_a_lcl=max_dx;
    type0 h_norm_lcl,dfa_lcl;
    
    h_norm_lcl=0.0;
    dfa_lcl=0.0;
    if(affine==0)
    {
        for(int i=0;i<atoms->natms*x_dim;i++)
        {
            h_norm_lcl+=h[i]*h[i];
            dfa_lcl+=h[i]*f[i];
        }
    }
    else
    {
        for(int i=dim;i<atoms->natms*x_dim;i+=dim)
        {
            for(int j=dim;j<x_dim;j++)
            {
                h_norm_lcl+=h[i]*h[i];
                dfa_lcl+=h[i]*f[i];
                i++;
            }
        }
    }

    MPI_Allreduce(&h_norm_lcl,&h_norm,1,MPI_TYPE0,MPI_SUM,world);
    MPI_Allreduce(&dfa_lcl,&dfa,1,MPI_TYPE0,MPI_SUM,world);
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
            {
                h_norm+=h_H[i][j]*h_H[i][j];
                dfa+=h_H[i][j]*f_H[i][j];
            }
    }
    
    dfa*=-1.0;
    
    if(h_norm==0.0)
    {
        max_a=0.0;
        dfa=0.0;
        return;
    }
    
    if(dfa>=0.0)
    {
        max_a=0.0;
        dfa=1.0;
        return;
    }
    
    h_norm=sqrt(h_norm);
    
    
    if(mapp->mode==MD_mode)
    {
        for(int i=0;i<atoms->natms*x_dim;i++)
            max_a_lcl=MIN(max_a_lcl,0.999*fabs(max_dx/h[i]));
    }
    else if(mapp->mode==DMD_mode)
    {
        int icmp=0;
        type0 a_max=0.99*(forcefield_dmd->alpha_max)/sqrt(2.0);
        type0 a_min=forcefield_dmd->alpha_min;

        max_a_lcl=a_max;
        type0* x=mapp->x->begin();
        for(int i=0;i<atoms->natms;i++)
        {
            for(int j=0;j<dim;j++)
            {
                max_a_lcl=MIN(max_a_lcl,0.999*fabs(max_dx/h[icmp]));
                icmp++;
            }
            for(int j=dim;j<x_dim;j++)
            {
                if(h[icmp]>0.0)
                    max_a_lcl=MIN(max_a_lcl,0.999*(a_max-x[icmp])/h[icmp]);
                if(h[icmp]<0.0)
                    max_a_lcl=MIN(max_a_lcl,0.999*(a_min-x[icmp])/h[icmp]);
                icmp++;
            }
        }
    }

    MPI_Allreduce(&max_a_lcl,&max_a,1,MPI_TYPE0,MPI_MIN,world);
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                if(h_H[i][j]!=0.0)
                    max_a=MIN(max_a,0.999*fabs(max_dx/h_H[i][j]));
    }
    
    
}
/*--------------------------------------------
 reset to initial position
 --------------------------------------------*/
void Min::F_reset()
{
    type0* x=mapp->x->begin();
    type0* x_prev=x_prev_ptr->begin();
    
    memcpy(x,x_prev,atoms->natms*x_dim*sizeof(type0));
    
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                atoms->H[i][j]=H_prev[i][j];
        
        if(dim==3)
            M3INV_TRI_LOWER(atoms->H,atoms->B);
        else
            xmath->invert_lower_triangle(atoms->H,atoms->B,dim);
    }
    
    
    
    atoms->update(mapp->x);
}



