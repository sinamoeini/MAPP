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
    
    H_dof=new int*[dim];
    *H_dof=new int[dim*dim];

    for(int i=1;i<dim;i++)
        H_dof[i]=H_dof[i-1]+dim;

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
        bool* dof=mapp->dof->begin();
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
        bool* dof=mapp->dof->begin();
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
    type0* fvec=f()->begin();
    for(int i=0;i<x_dim*atoms->natms;i++)
        fvec[i]=0.0;

    forcefield->force_calc_timer(sts_flag,nrgy_strss);
    
    if(!chng_box)
    {
        if(mapp->dof!=NULL)
        {
            bool* dof=mapp->dof->begin();
            for(int i=0;i<atoms->natms*x_dim;i++)
                if(dof[i]) fvec[i]=0.0;
        }
        return;
    }
    
    if(affine)
    {
        type0** B=atoms->B;
        type0** H=atoms->H;
        type0 vol=H[0][0]*H[1][1]*H[2][2];
        if(H_dof[0][0])
            f.A[0][0]=-(nrgy_strss[1]*B[0][0]+nrgy_strss[6]*B[1][0]+nrgy_strss[5]*B[2][0])*vol;
        if(H_dof[1][1])
            f.A[1][1]=-(nrgy_strss[2]*B[1][1]+nrgy_strss[4]*B[2][1])*vol;
        if(H_dof[2][2])
            f.A[2][2]=-(nrgy_strss[3]*B[2][2])*vol;
        if(H_dof[2][1])
            f.A[2][1]=-(nrgy_strss[4]*B[2][2])*vol;
        if(H_dof[2][0])
            f.A[2][0]=-(nrgy_strss[5]*B[2][2])*vol;
        if(H_dof[1][0])
            f.A[1][0]=-(nrgy_strss[6]*B[1][1]+nrgy_strss[5]*B[2][1])*vol;

        if(mapp->dof==NULL)
        {
            fvec=f()->begin();
            
            if(x_dim!=dim)
            {
                for(int iatm=0;iatm<atoms->natms;iatm++,fvec+=x_dim)
                    for(int i=0;i<dim;i++)
                        fvec[i]=0.0;
            }
            else
                for(int i=0;i<atoms->natms*dim;i++)
                    fvec[i]=0.0;
            return;
        }

        type0 st_lcl[6];
        type0 st[6];
        for(int i=0;i<6;i++)
            st_lcl[i]=0.0;
        
        type0* xvec=x()->begin();
        fvec=f()->begin();
        bool* dof=mapp->dof->begin();
        for(int i=0;i<atoms->natms;i++)
        {
            if(dof[0])
            {
                st_lcl[0]+=fvec[0]*xvec[0];
                st_lcl[4]+=fvec[0]*xvec[2];
                st_lcl[5]+=fvec[0]*xvec[1];
            }
            
            if(dof[1])
            {
                st_lcl[1]+=fvec[1]*xvec[1];
                st_lcl[3]+=fvec[1]*xvec[2];
            }
            
            if(dof[2])
            {
                st_lcl[2]+=fvec[2]*xvec[2];
            }
            
            for(int j=0;j<dim;j++)
                fvec[j]=0.0;
            
            for(int j=dim;j<x_dim;j++)
                if(dof[j])
                    fvec[j]=0.0;
            
            fvec+=x_dim;
            xvec+=x_dim;
            dof+=x_dim;
        }
        
        MPI_Allreduce(st_lcl,st,6,MPI_TYPE0,MPI_SUM,world);
        
        st_lcl[0]=-(st[0]*B[0][0]+st[5]*B[1][0]+st[4]*B[2][0]);
        st_lcl[1]=-(st[1]*B[1][1]+st[3]*B[2][1]);
        st_lcl[2]=-(st[2]*B[2][2]);
        st_lcl[3]=-(st[3]*B[2][2]);
        st_lcl[4]=-(st[4]*B[2][2]);
        st_lcl[5]=-(st[5]*B[1][1]+st[4]*B[2][1]);
        
        if(H_dof[0][0])
            f.A[0][0]-=st_lcl[0];
        if(H_dof[1][1])
            f.A[1][1]-=st_lcl[1];
        if(H_dof[2][2])
            f.A[2][2]-=st_lcl[2];
        if(H_dof[2][1])
            f.A[2][1]-=st_lcl[3];
        if(H_dof[2][0])
            f.A[2][0]-=st_lcl[4];
        if(H_dof[1][0])
            f.A[1][0]-=st_lcl[5];
        
    }
    else
    {
        type0** B=atoms->B;
        type0** H=atoms->H;
        
        type0 st_lcl[6];
        type0 st[6];
        for(int i=0;i<6;i++)
            st_lcl[i]=0.0;
        
        fvec=f()->begin();
        type0* xvec=x()->begin();
        for(int iatm=0;iatm<atoms->natms;iatm++)
        {
            st_lcl[0]+=fvec[0]*xvec[0];
            st_lcl[1]+=fvec[1]*xvec[1];
            st_lcl[2]+=fvec[2]*xvec[2];
            st_lcl[3]+=fvec[1]*xvec[2];
            st_lcl[4]+=fvec[0]*xvec[2];
            st_lcl[5]+=fvec[0]*xvec[1];
            
            fvec+=x_dim;
            xvec+=x_dim;
        }
        

        MPI_Allreduce(st_lcl,st,6,MPI_TYPE0,MPI_SUM,world);
        
        type0 vol=H[0][0]*H[1][1]*H[2][2];
        for(int i=0;i<6;i++)
            st_lcl[i]=st[i]+nrgy_strss[i+1]*vol;
        
        st[0]=-(st_lcl[0]*B[0][0]+st_lcl[5]*B[1][0]+st_lcl[4]*B[2][0]);
        st[1]=-(st_lcl[1]*B[1][1]+st_lcl[3]*B[2][1]);
        st[2]=-(st_lcl[2]*B[2][2]);
        st[3]=-(st_lcl[3]*B[2][2]);
        st[4]=-(st_lcl[4]*B[2][2]);
        st[5]=-(st_lcl[5]*B[1][1]+st_lcl[4]*B[2][1]);

        if(H_dof[0][0])
            f.A[0][0]=st[0];
        if(H_dof[1][1])
            f.A[1][1]=st[1];
        if(H_dof[2][2])
            f.A[2][2]=st[2];
        if(H_dof[2][1])
            f.A[2][1]=st[3];
        if(H_dof[2][0])
            f.A[2][0]=st[4];
        if(H_dof[1][0])
            f.A[1][0]=st[5];

        if(mapp->dof==NULL)
            return;
        
        bool* dof=mapp->dof->begin();
        fvec=f()->begin();
        for(int i=0;i<atoms->natms*x_dim;i++)
            if(dof[i]) fvec[i]=0.0;
    }

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::prepare_affine_h()
{
    
    if(chng_box)
    {
        type0** B=atoms->B;
        type0 N[3][3];
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
            {
                N[i][j]=0.0;
                for(int k=0;k<dim;k++)
                    N[i][j]+=B[j][k]*h.A[k][i];
            }
        
        type0* xvec=x0()->begin();
        type0* hvec=h()->begin();
        if(dim==3)
        {
            for(int iatm=0;iatm<atoms->natms;iatm++)
            {
                M3V_TRI_UPPER(N,xvec,hvec);
                xvec+=x_dim;
                hvec+=x_dim;
            }
        }
        else
        {
            for(int iatm=0;iatm<atoms->natms;iatm++)
            {
                for(int j=0;j<dim;j++)
                {
                    hvec[j]=0.0;
                    for(int k=j;k<dim;k++)
                        hvec[j]+=N[j][k]*xvec[k];
                }
                xvec+=x_dim;
                hvec+=x_dim;
            }
        }
    }
    else
    {
        type0* hvec=h()->begin();
        for(int iatm=0;iatm<atoms->natms;iatm++,hvec+=x_dim)
            for(int j=0;j<dim;j++)
                hvec[j]=0.0;
    }
    

    if(mapp->dof==NULL)
        return;
    
    bool* dof=mapp->dof->begin();
    type0*hvec=h()->begin();
    for(int i=0;i<atoms->natms*x_dim;i++)
        if(dof[i]) hvec[i]=0.0;
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
        sts_flag=2;
    else if(chng_box && !affine)
        sts_flag=2;
    
    if(chng_box)
    {
        H_prev=new type0*[dim];
        f_H_prev=new type0*[dim];
        
        *H_prev=new type0[dim*dim];
        *f_H_prev=new type0[dim*dim];
        
        for(int i=1;i<dim;i++)
        {
            H_prev[i]=H_prev[i-1]+dim;
            f_H_prev[i]=f_H_prev[i-1]+dim;
        }
    }
    
    if(mapp->f==NULL)
        mapp->f=new Vec<type0>(atoms,x_dim);
    
    x_prev_ptr=new Vec<type0>(atoms,x_dim);
    h_ptr=new Vec<type0>(atoms,x_dim);
    f_prev_ptr=new Vec<type0>(atoms,x_dim);
    
    x.init(atoms,atoms->x,atoms->H,chng_box);
    x0.init(atoms,x_prev_ptr,H_prev,chng_box);
    f.init(atoms,mapp->f,chng_box);
    f0.init(atoms,f_prev_ptr,f_H_prev,chng_box);
    h.init(atoms,h_ptr,chng_box);

    
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
    x.fin();
    x0.fin();
    f.fin();
    f0.fin();
    h.fin();
    delete f_prev_ptr;
    delete h_ptr;
    delete x_prev_ptr;
    delete vecs_comm;
    timer->print_stats();
    neighbor->print_stats();
    
    if(chng_box)
    {
        delete [] *H_prev;
        delete [] H_prev;
        
        delete [] *f_H_prev;
        delete [] f_H_prev;
    }
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 Min::F(type0 alpha)
{
    x=x0;
    x.add(alpha,h);

    if(chng_box)
    {
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
    
    x=x0;
    x.add(alpha,h);
    
    if(chng_box)
    {
        if(dim==3)
            M3INV_TRI_LOWER(atoms->H,atoms->B);
        else
            xmath->invert_lower_triangle(atoms->H,atoms->B,dim);
    }
    
    atoms->update(mapp->x);
    force_calc();

    if(affine)
    {
        drev=0.0;
        if(x_dim!=dim)
        {
            
            type0* hvec=h()->begin()+dim;
            type0* fvec=f()->begin()+dim;
            type0 dfa_lcl=0.0;
            
            for(int iatm=0;iatm<atoms->natms;iatm++)
            {
                for(int i=0;i<x_dim-dim;i++)
                    dfa_lcl+=hvec[i]*fvec[i];
                
                fvec+=x_dim;
                hvec+=x_dim;
            }
            
            MPI_Allreduce(&dfa_lcl,&drev,1,MPI_TYPE0,MPI_SUM,world);
        }
        
        if(chng_box)
        {
            for(int i=0;i<dim;i++)
                for(int j=0;j<dim;j++)
                    drev+=h.A[i][j]*f.A[i][j];
        }
        
        drev*=-1.0;
    }
    else
    {
        drev=-(f*h);
    }
    
    return nrgy_strss[0];
}
/*--------------------------------------------
 find maximum h
 --------------------------------------------*/
void Min::ls_prep(type0& dfa,type0& h_norm,type0& max_a)
{
    type0 max_a_lcl=max_dx;

    if(affine)
    {
        h_norm=0.0;
        
        if(x_dim!=dim)
        {
            type0* hvec=h()->begin()+dim;
            type0 h_norm_lcl=0.0;
            for(int iatm=0;iatm<atoms->natms;iatm++)
            {
                for(int i=0;i<x_dim-dim;i++)
                    h_norm_lcl+=hvec[i]*hvec[i];
                hvec+=x_dim;
            }
            MPI_Allreduce(&h_norm_lcl,&h_norm,1,MPI_TYPE0,MPI_SUM,world);
        }
        
        if(chng_box)
        {
            for(int i=0;i<dim;i++)
                for(int j=0;j<dim;j++)
                    h_norm+=h.A[i][j]*h.A[i][j];
        }
    }
    else
    {
        h_norm=h*h;
    }

    dfa=-f_h;
    
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
        type0* hvec=h()->begin();
        for(int i=0;i<atoms->natms*x_dim;i++)
            max_a_lcl=MIN(max_a_lcl,0.999*fabs(max_dx/hvec[i]));
    }
    else if(mapp->mode==DMD_mode)
    {
        type0* hvec=h()->begin();
        int icmp=0;
        type0 a_max=0.99*(forcefield_dmd->alpha_max)/sqrt(2.0);
        type0 a_min=forcefield_dmd->alpha_min;

        max_a_lcl=a_max;
        type0* x=mapp->x->begin();
        for(int i=0;i<atoms->natms;i++)
        {
            for(int j=0;j<dim;j++)
            {
                max_a_lcl=MIN(max_a_lcl,0.999*fabs(max_dx/hvec[icmp]));
                icmp++;
            }
            for(int j=dim;j<x_dim;j++)
            {
                if(hvec[icmp]>0.0)
                    max_a_lcl=MIN(max_a_lcl,0.999*(a_max-x[icmp])/hvec[icmp]);
                if(hvec[icmp]<0.0)
                    max_a_lcl=MIN(max_a_lcl,0.999*(a_min-x[icmp])/hvec[icmp]);
                icmp++;
            }
        }
    }

    MPI_Allreduce(&max_a_lcl,&max_a,1,MPI_TYPE0,MPI_MIN,world);
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                if(h.A[i][j]!=0.0)
                    max_a=MIN(max_a,0.999*fabs(max_dx/h.A[i][j]));
    }
}
/*--------------------------------------------
 reset to initial position
 --------------------------------------------*/
void Min::F_reset()
{
    x=x0;
    atoms->update(mapp->x);
}



