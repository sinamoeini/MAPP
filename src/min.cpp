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
#include "MAPP.h"
#include "script_reader.h"
#define INV_SQ_2 0.7071067811865475
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min::Min():
nrgy_strss(forcefield->nrgy_strss)
{
    output_flag=true;
    if(forcefield==NULL)
        error->abort("ff should be "
        "initiated before min");
    
    
    chng_box=false;
    x_dim=atoms->x->dim;
    err=LS_S;
    
    char** args;
    int nargs;
    if(mode==MD_mode)
        nargs=ScriptReader::parse_line(
        "PE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    else if(mode==DMD_mode)
        nargs=ScriptReader::parse_line(
        "FE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    
    pe_idx=0;
    stress_idx=1;
    
    thermo=new ThermoDynamics(nargs,args);
    for(int i=0;i<nargs;i++)
        delete [] args[i];
    if(nargs)
        delete [] args;
    
    // default values 
    max_iter=10000;
    energy_tolerance=1.0e-8;
    affine=0;

    
    H_dof=new bool*[__dim__];
    *H_dof=new bool[__dim__*__dim__];

    for(int i=1;i<__dim__;i++)
        H_dof[i]=H_dof[i-1]+__dim__;

    for(int i=0;i<__dim__;i++)
        for(int j=0;j<__dim__;j++)
            H_dof[i][j]=false;
    
    
    if(mode==DMD_mode)
        forcefield_dmd=dynamic_cast<ForceFieldDMD*>(forcefield);
    max_dx=1.0;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min::~Min()
{
    delete thermo;

    if(__dim__)
    {
        delete [] *H_dof;
        delete [] H_dof;
    }
    
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
    else if(err==MIN_S_TOLERANCE)
    {
        if(atoms->my_p==0)
            fprintf(output,"minimization finished: energy tolerance reached\n");
    }
    else if(err==MIN_F_MAX_ITER)
    {
        if(atoms->my_p==0)
            fprintf(output,"minimization finished: maximum iteration reached\n");
    }
    else if(err==B_F_DOWNHILL)
    {
        if(atoms->my_p==0)
            fprintf(output,"bracketing failed: not downhill direction\n");
    }
    else if(err==B_F_MAX_ALPHA)
    {
        if(atoms->my_p==0)
            fprintf(output,"bracketing failed: maximum alpha reached\n");
    }

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::force_calc()
{
    type0* fvec=f()->begin();
    for(int i=0;i<x_dim*atoms->natms;i++)
        fvec[i]=0.0;

    forcefield->force_calc_timer(sts_flag);
    
    if(output_flag)
    {
        thermo->update(pe_idx,nrgy_strss[0]);
        thermo->update(stress_idx,6,&nrgy_strss[1]);
    }

    
    if(!chng_box)
    {
        if(mapp->x_dof)
        {
            bool* dof=mapp->x_dof->begin();
            for(int i=0;i<atoms->natms*x_dim;i++)
                fvec[i]*=dof[i];
        }
        return;
    }
    
    if(affine)
    {
        type0 (&H)[__dim__][__dim__]=atoms->H;
        type0 (&B)[__dim__][__dim__]=atoms->B;
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

        if(mapp->x_dof==NULL)
        {
            fvec=f()->begin();
            
            if(x_dim!=__dim__)
            {
                for(int iatm=0;iatm<atoms->natms;iatm++,fvec+=x_dim)
                    for(int i=0;i<__dim__;i++)
                        fvec[i]=0.0;
            }
            else
                for(int i=0;i<atoms->natms*__dim__;i++)
                    fvec[i]=0.0;
            return;
        }

        type0 st_lcl[6];
        type0 st[6];
        for(int i=0;i<6;i++)
            st_lcl[i]=0.0;
        
        type0* xvec=x()->begin();
        fvec=f()->begin();
        bool* dof=mapp->x_dof->begin();
        for(int i=0;i<atoms->natms;i++)
        {
            if(!dof[0])
            {
                st_lcl[0]+=fvec[0]*xvec[0];
                st_lcl[4]+=fvec[0]*xvec[2];
                st_lcl[5]+=fvec[0]*xvec[1];
            }
            
            if(!dof[1])
            {
                st_lcl[1]+=fvec[1]*xvec[1];
                st_lcl[3]+=fvec[1]*xvec[2];
            }
            
            if(!dof[2])
            {
                st_lcl[2]+=fvec[2]*xvec[2];
            }
            
            for(int j=0;j<__dim__;j++)
                fvec[j]=0.0;
            
            for(int j=__dim__;j<x_dim;j++)
                fvec[j]*=dof[j];
            
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
        type0 (&B)[__dim__][__dim__]=atoms->B;
        type0 (&H)[__dim__][__dim__]=atoms->H;
        
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

        if(mapp->x_dof==NULL)
            return;
        
        bool* dof=mapp->x_dof->begin();
        fvec=f()->begin();
        for(int i=0;i<atoms->natms*x_dim;i++)
            fvec[i]*=dof[i];
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::prepare_affine_h()
{
    
    if(chng_box)
    {
        type0 (&B)[__dim__][__dim__]=atoms->B;
        type0 N[3][3];
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
            {
                N[i][j]=0.0;
                for(int k=0;k<__dim__;k++)
                    N[i][j]+=B[j][k]*h.A[k][i];
            }
        
        type0* xvec=x0()->begin();
        type0* hvec=h()->begin();
        if(__dim__==3)
        {
            for(int iatm=0;iatm<atoms->natms;iatm++)
            {
                M3V_UT(N,xvec,hvec);
                xvec+=x_dim;
                hvec+=x_dim;
            }
        }
        else
        {
            for(int iatm=0;iatm<atoms->natms;iatm++)
            {
                for(int j=0;j<__dim__;j++)
                {
                    hvec[j]=0.0;
                    for(int k=j;k<__dim__;k++)
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
            for(int j=0;j<__dim__;j++)
                hvec[j]=0.0;
    }
    

    if(mapp->x_dof==NULL)
        return;
    
    bool* dof=mapp->x_dof->begin();
    type0* hvec=h()->begin();
    for(int i=0;i<atoms->natms*x_dim;i++)
        hvec[i]*=dof[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 Min::calc_ave_f_norm()
{
    force_calc();
    return sqrt(f*f)/ndofs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 Min::calc_ndofs()
{
    type0 ndofs_lcl=0.0;
    if(mode==DMD_mode)
    {
        type0* c=mapp->c->begin();
        int c_dim=mapp->c->dim;
        if(mapp->x_dof==NULL)
        {
            for(int iatm=0;iatm<atoms->natms;iatm++)
            {
                for(int i=0;i<c_dim;i++)
                    if(c[i]>=0.0)
                        ndofs_lcl++;
                
                c+=c_dim;
            }
        }
        else
        {
            bool* dof=mapp->x_dof->begin()+__dim__;
            for(int iatm=0;iatm<atoms->natms;iatm++)
            {
                for(int i=0;i<c_dim;i++)
                    if(c[i]>=0.0 && dof[i])
                        ndofs_lcl++;
                
                c+=c_dim;
                dof+=x_dim;
            }
            
        }
    }
    
    if(!affine)
    {
        if(mapp->x_dof)
        {
            bool* dof=mapp->x_dof->begin();
            for(int iatm=0;iatm<atoms->natms;iatm++)
            {
                for(int i=0;i<__dim__;i++)
                    if(dof[i])
                        ndofs_lcl++;
                dof+=x_dim;
            }
        }
        else
        {
            ndofs_lcl+=static_cast<type0>(atoms->natms*__dim__);
        }
    }
    type0 ndofs;
    MPI_Allreduce(&ndofs_lcl,&ndofs,1,MPI_TYPE0,MPI_SUM,world);
    for(int i=0;i<__dim__;i++)
        for(int j=0;j<__dim__;j++)
            if(H_dof[i][j])
                ndofs++;
    
    return ndofs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::init()
{
    // determine if the there is chng_box
    chng_box=false;
    for(int i=0;i<__dim__;i++)
        for(int j=0;j<__dim__;j++)
            if(H_dof[i][j])
                chng_box=true;
    
    ndofs=calc_ndofs();
    if(ndofs==0)
        error->abort("for minimization to "
        "perform the degrees of freedom "
        "must be greater than 0.0");
    
    sts_flag=true;

    
    x.init(atoms,atoms->x,atoms->H,chng_box);
    x0.init(atoms,chng_box);
    f.init(atoms,forcefield->f,chng_box);
    f0.init(atoms,chng_box);
    h.init(atoms,chng_box);

    
    vecs_comm=new VecLst(atoms);

    if(mode==MD_mode)
    {
        vecs_comm->add_updt(mapp->type);
        if(mapp->x_dof)
            vecs_comm->add_xchng(mapp->x_dof);
    }
    else if(mode==DMD_mode)
    {
        vecs_comm->add_updt(mapp->c);
        vecs_comm->add_updt(mapp->ctype);
        if(mapp->x_dof)
        {
            vecs_comm->add_xchng(mapp->x_dof);
        }
        if(mapp->c_dof)
            vecs_comm->add_xchng(mapp->c_dof);
    }
    
    vecs_comm->add_xchng(h());
    vecs_comm->add_xchng(x0());
    vecs_comm->add_xchng(f0());
    
    if(mapp->x_d)
        vecs_comm->add_arch(mapp->x_d);
    
    if(ls==NULL)
    {
        LineSearch_backtrack<Min>* ls_=NULL;
        ScriptReader::create(ls_,"ls bt");
        ls=ls_;
        
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
    delete vecs_comm;
    if(output_flag)
    {
        print_error();
        timer->print_stats();
        neighbor->print_stats();
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
        if(__dim__==3)
            M3INV_LT(atoms->H,atoms->B);
        else
            XMatrixVector::Mlt_inv<__dim__>(atoms->H,atoms->B);
        
    }
    
    atoms->update(atoms->x);
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
        if(__dim__==3)
            M3INV_LT(atoms->H,atoms->B);
        else
            XMatrixVector::Mlt_inv<__dim__>(atoms->H,atoms->B);
    }
    
    atoms->update(atoms->x);
    force_calc();

    if(affine)
    {
        drev=0.0;
        if(x_dim!=__dim__)
        {
            
            type0* hvec=h()->begin()+__dim__;
            type0* fvec=f()->begin()+__dim__;
            type0 dfa_lcl=0.0;
            
            for(int iatm=0;iatm<atoms->natms;iatm++)
            {
                for(int i=0;i<x_dim-__dim__;i++)
                    dfa_lcl+=hvec[i]*fvec[i];
                
                fvec+=x_dim;
                hvec+=x_dim;
            }
            
            MPI_Allreduce(&dfa_lcl,&drev,1,MPI_TYPE0,MPI_SUM,world);
        }
        
        if(chng_box)
        {
            for(int i=0;i<__dim__;i++)
                for(int j=0;j<__dim__;j++)
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
        
        if(x_dim!=__dim__)
        {
            type0* hvec=h()->begin()+__dim__;
            type0 h_norm_lcl=0.0;
            for(int iatm=0;iatm<atoms->natms;iatm++)
            {
                for(int i=0;i<x_dim-__dim__;i++)
                    h_norm_lcl+=hvec[i]*hvec[i];
                hvec+=x_dim;
            }
            MPI_Allreduce(&h_norm_lcl,&h_norm,1,MPI_TYPE0,MPI_SUM,world);
        }
        
        if(chng_box)
        {
            for(int i=0;i<__dim__;i++)
                for(int j=0;j<__dim__;j++)
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
    
    
    if(mode==MD_mode)
    {
        type0* hvec=h()->begin();
        for(int i=0;i<atoms->natms*x_dim;i++)
            max_a_lcl=MIN(max_a_lcl,0.999*fabs(max_dx/hvec[i]));
    }
    else if(mode==DMD_mode)
    {
        type0* hvec=h()->begin();
        int icmp=0;
        type0 a_max=0.99*(forcefield_dmd->alpha_max)/sqrt(2.0);
        type0 a_min=forcefield_dmd->alpha_min;

        max_a_lcl=a_max;
        type0* x=atoms->x->begin();
        for(int i=0;i<atoms->natms;i++)
        {
            for(int j=0;j<__dim__;j++)
            {
                max_a_lcl=MIN(max_a_lcl,0.999*fabs(max_dx/hvec[icmp]));
                icmp++;
            }
            for(int j=__dim__;j<x_dim;j++)
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
        for(int i=0;i<__dim__;i++)
            for(int j=0;j<__dim__;j++)
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
    atoms->update(atoms->x);
    if(chng_box)
    {
        if(__dim__==3)
            M3INV_LT(atoms->H,atoms->B);
        else
            XMatrixVector::Mlt_inv<__dim__>(atoms->H,atoms->B);
    }
}



