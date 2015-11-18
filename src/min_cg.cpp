#include <stdlib.h>
#include "memory.h"
#include "error.h"
#include "min_cg.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min_cg::Min_cg(MAPP* mapp,int nargs,char** args):Min(mapp)
{
    int icmp;
    int jcmp;
    int iarg=2;
    while(iarg<nargs)
    {
        if(!strcmp(args[iarg],"max_iter"))
        {
            iarg++;
            if(iarg==nargs)
                error->abort("max_iter in min cg should at least have 1argument");
            max_iter=atoi(args[iarg]);
            iarg++;
        }
        else if(!strcmp(args[iarg],"e_tol"))
        {
            iarg++;
            if(iarg==nargs)
                error->abort("e_tol in min cg should at least have 1argument");
            energy_tolerance=atof(args[iarg]);
            iarg++;
        }
        else if(sscanf(args[iarg],"H[%d][%d]",&icmp,&jcmp)==2)
        {
            if(icmp<0 || icmp>=atoms->dimension)
                error->abort("wrong component in min cg for H[%i][%i]",icmp,jcmp);
            if(jcmp<0 || jcmp>=atoms->dimension)
                error->abort("wrong component in min cg for H[%i][%i]",icmp,jcmp);
            
            if(icmp<=jcmp)
                H_dof[jcmp][icmp]=1;
            else
                H_dof[icmp][jcmp]=1;
            iarg++;
        }
        else if(!strcmp(args[iarg],"affine"))
        {
            affine=1;
            iarg++;
        }
        else
            error->abort("unknown keyword for min cg: %s",args[iarg]);
    }
    
    if(max_iter<0)
        error->abort("max_iter in min cg should be greater than 0");
    if(energy_tolerance<0.0)
         error->abort("e_tol in min cg should be greater than 0.0");
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Min_cg::~Min_cg()
{
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void Min_cg::init()
{
    Min::init();
    atoms->init(vecs_comm,chng_box);
    force_calc();
    curr_energy=nrgy_strss[0];
    thermo->update(pe_idx,nrgy_strss[0]);
    thermo->update(stress_idx,6,&nrgy_strss[1]);
    
    thermo->init();
    
    if(write!=NULL)
        write->init();
    
}
/*--------------------------------------------
 min
 --------------------------------------------*/
void Min_cg::run()
{
    if(max_iter==0) return;
    
    type0* x=mapp->x->begin();
    type0* f=mapp->f->begin();
    type0* x_0=x_prev_ptr->begin();
    type0* f_0=f_prev_ptr->begin();
    type0* h=h_ptr->begin();
    type0 prev_energy;
    type0 alpha;
    int size;
    int istp=0;
    type0 f0_f0;
    type0 f_f;
    type0 f_h;
    type0 f_f0;
    type0 ratio;
    type0 inner;
    err=LS_S;
    
    memcpy(h,f,x_dim*atoms->natms*sizeof(type0));
    
    inner=0.0;
    for(int i=0;i<atoms->natms*x_dim;i++)
        inner+=f[i]*f[i];
    f0_f0=0.0;
    MPI_Allreduce(&inner,&f0_f0,1,MPI_TYPE0,MPI_SUM,world);
    
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                h_H[i][j]=f_H[i][j];
        
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                f0_f0+=f_H[i][j]*f_H[i][j];
    }
    
    while(err==LS_S)
    {
        
        if(f0_f0==0.0)
        {
            err=LS_F_GRAD0;
            continue;
        }
        
        x=mapp->x->begin();
        f=mapp->f->begin();
        x_0=x_prev_ptr->begin();
        f_0=f_prev_ptr->begin();
        size=atoms->natms*x_dim*sizeof(type0);
        
        
        memcpy(x_0,x,size);
        memcpy(f_0,f,size);
        
        if(chng_box)
        {
            for(int i=0;i<dim;i++)
                for(int j=0;j<dim;j++)
                {
                    H_prev[i][j]=atoms->H[i][j];
                    B_prev[i][j]=atoms->B[i][j];
                    f_H_prev[i][j]=f_H[i][j];
                }
        }
        
        prev_energy=curr_energy;
        
        
        if(write!=NULL)
            write->write();
        
        thermo->thermo_print();
        
        if(affine) prepare_affine_h(x_0,h);
        err=ls->line_min(curr_energy,alpha,1);
        h=h_ptr->begin();
        if(affine) rectify(h);
        
        if(err!=LS_S)
        {
            thermo->update(pe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            continue;
        }
        
        
        if(prev_energy-curr_energy<energy_tolerance)
            err=MIN_S_TOLERANCE;
        
        if(istp+1==max_iter)
            err=MIN_F_MAX_ITER;
        
        
        force_calc();
        
        
        
        if(thermo->test_prev_step() || err)
        {
            thermo->update(pe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
        }
        
        istp++;
        step_no++;
        
        if(err)
            continue;
        
        f=mapp->f->begin();
        h=h_ptr->begin();
        f_0=f_prev_ptr->begin();
        inner=0.0;
        for(int i=0;i<x_dim*atoms->natms;i++)
            inner+=f[i]*f[i];
        
        f_f=0.0;
        MPI_Allreduce(&inner,&f_f,1,MPI_TYPE0,MPI_SUM,world);
        if(chng_box)
        {
            for(int i=0;i<dim;i++)
                for(int j=0;j<dim;j++)
                    f_f+=f_H[i][j]*f_H[i][j];
        }
        
        inner=0.0;
        for(int i=0;i<x_dim*atoms->natms;i++)
            inner+=f[i]*f_0[i];
        f_f0=0.0;
        MPI_Allreduce(&inner,&f_f0,1,MPI_TYPE0,MPI_SUM,world);
        if(chng_box)
        {
            for(int i=0;i<dim;i++)
                for(int j=0;j<dim;j++)
                    f_f0+=f_H_prev[i][j]*f_H[i][j];
        }
        
        ratio=(f_f-f_f0)/(f0_f0);
        
        inner=0.0;
        
        for(int i=0;i<x_dim*atoms->natms;i++)
        {
            h[i]*=ratio;
            h[i]+=f[i];
            inner+=h[i]*f[i];
        }
        
        MPI_Allreduce(&inner,&f_h,1,MPI_TYPE0,MPI_SUM,world);
        
        if(chng_box)
        {
            for(int i=0;i<dim;i++)
                for(int j=0;j<dim;j++)
                {
                    h_H[i][j]*=ratio;
                    h_H[i][j]+=f_H[i][j];
                    f_h+=h_H[i][j]*f_H[i][j];
                }
        }
        
        
        
        if(f_h<0.0)
        {
            size=atoms->natms*x_dim*sizeof(type0);
            memcpy(h,f,size);
            
            if(chng_box)
            {
                for(int i=0;i<dim;i++)
                    for(int j=0;j<dim;j++)
                    {
                        h_H[i][j]=f_H[i][j];
                    }
            }
        }
        
        f0_f0=f_f;
        
    }

}
/*--------------------------------------------
 finishing minimization
 --------------------------------------------*/
void Min_cg::fin()
{
    
    if(write!=NULL)
        write->fin();
     
    thermo->fin();
    print_error();
    atoms->fin();
    
    Min::fin();
    
}

