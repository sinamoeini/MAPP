/*--------------------------------------------
 Created by Sina on 07/02/14.
 Copyright (c) 2013 MIT. All rights reserved.
 
 L-BFGS minimization is written based on
 Numerical Optimization written by Nocedal & 
 Wright, second edition, pages 177-179, 
 Algorithm 7.4 & 7.5 Equation (7.20)
 
 with respect to notations in Nocedal:
 new_y_i=y_i
 new_rho_i=rho_i
 new_alpha_i=alpha_i
 new_beta=beta
 --------------------------------------------*/
#include <stdlib.h>
#include "min_l-bfgs.h"
#include "ff.h"
#include "memory.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min_lbfgs::Min_lbfgs(MAPP* mapp,int nargs,char** args):Min(mapp)
{
    m_it=2;
    int icmp;
    int jcmp;
    int iarg=2;
    while(iarg<nargs)
    {
        if(!strcmp(args[iarg],"max_iter"))
        {
            iarg++;
            if(iarg==nargs)
                error->abort("max_iter in min l-bfgs should at least have 1argument");
            max_iter=atoi(args[iarg]);
            iarg++;
        }
        else if(!strcmp(args[iarg],"e_tol"))
        {
            iarg++;
            if(iarg==nargs)
                error->abort("e_tol in min l-bfgs should at least have 1argument");
            energy_tolerance=atof(args[iarg]);
            iarg++;
        }
        else if(!strcmp(args[iarg],"m"))
        {
            iarg++;
            if(iarg==nargs)
                error->abort("m in min l-bfgs should at least have 1argument");
            m_it=atoi(args[iarg]);
            iarg++;
        }
        else if(sscanf(args[iarg],"H[%d][%d]",&icmp,&jcmp)==2)
        {
            if(icmp<0 || icmp>=dim)
                error->abort("wrong component in min l-bfgs for H[%i][%i]",icmp,jcmp);
            if(jcmp<0 || jcmp>=dim)
                error->abort("wrong component in min l-bfgs for H[%i][%i]",icmp,jcmp);
            
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
            error->abort("unknown keyword for min l-bfgs: %s",args[iarg]);
    }
    
    if(max_iter<0)
        error->abort("max_iter in min l-bfgs should be greater than 0");
    if(energy_tolerance<=0.0)
        error->abort("e_tol in min l-bfgs should be greater than 0.0");

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Min_lbfgs::~Min_lbfgs()
{
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void Min_lbfgs::init()
{
     Min::init();
    
    if(chng_box)
    {
        CREATE1D(y_H,m_it);
        CREATE1D(s_H,m_it);
        for(int i=0;i<m_it;i++)
        {
            CREATE2D(y_H[i],dim,dim);
            CREATE2D(s_H[i],dim,dim);
        }

    }
    
    CREATE1D(rho,m_it);
    CREATE1D(alpha,m_it);
    
    
    /* begining of chekcking if the primary atomic vectors exist */
    /* if not, add them */

    s=new Vec<type0>*[m_it];
    y=new Vec<type0>*[m_it];
    for(int i=0;i<m_it;i++)
    {
        s[i]=new Vec<type0>(atoms,x_dim);
        y[i]=new Vec<type0>(atoms,x_dim);
    }
    
    for(int i=0;i<m_it;i++)
        vecs_comm->add_xchng(s[i]);
    for(int i=0;i<m_it;i++)
        vecs_comm->add_xchng(y[i]);

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
 run
 --------------------------------------------*/
void Min_lbfgs::run()
{
    if(max_iter==0) return;
    
    type0* f;
    type0* x;
    type0* x_0;
    type0* f_0;
    type0* h;
    type0* tmp_vec0;
    type0* tmp_vec1;
    type0 prev_energy;
    type0 alpha_m,beta,gamma;
    type0 inner0_lcl,inner1_lcl,inner0,inner1;
    type0** H_y_tmp=NULL;
    type0** H_s_tmp=NULL;
    
    
    int k_it=0;
    int size;
    int istp=0;
    
    gamma=1.0;
    err=LS_S;
    
    while(err==LS_S)
    {
        /* beginning of storing the parameters */
        x=mapp->x->begin();
        f=mapp->f->begin();
        x_0=x_prev_ptr->begin();
        f_0=f_prev_ptr->begin();
        h=h_ptr->begin();
        size=atoms->natms*x_dim*sizeof(type0);
        
        memcpy(h,f,size);
        memcpy(x_0,x,size);
        memcpy(f_0,f,size);
        
        if(chng_box)
        {
            for(int i=0;i<dim;i++)
            {
                for(int j=0;j<dim;j++)
                {
                    h_H[i][j]=f_H[i][j];
                    f_H_prev[i][j]=f_H[i][j];
                    H_prev[i][j]=atoms->H[i][j];
                    B_prev[i][j]=atoms->B[i][j];
                }
            }
        }
        
        /* end of storing the parameters */
        
        
        /* beginning of determining new h & h_H */
        
        /* level 0 */
        for(int i=0;i<k_it;i++)
        {
            tmp_vec0=s[i]->begin();
            inner0_lcl=0.0;
            for(int j=0;j<atoms->natms*x_dim;j++)
                inner0_lcl+=h[j]*tmp_vec0[j];
            
            MPI_Allreduce(&inner0_lcl,&inner0,1,MPI_TYPE0,MPI_SUM,world);
            
            if(chng_box)
            {
                for(int j=0;j<dim;j++)
                    for(int k=0;k<dim;k++)
                        inner0+=h_H[j][k]*s_H[i][j][k];
            }
            
            
            alpha[i]=-rho[i]*inner0;
            
            tmp_vec0=y[i]->begin();
            for(int j=0;j<atoms->natms*x_dim;j++)
                h[j]+=alpha[i]*tmp_vec0[j];
            
            if(chng_box)
                for(int j=0;j<dim;j++)
                    for(int k=0;k<dim;k++)
                        h_H[j][k]+=alpha[i]*y_H[i][j][k];
            
        }
        
        /* level 1 */
        for(int j=0;j<atoms->natms*x_dim;j++)
            h[j]*=gamma;
        
        if(chng_box)
            for(int j=0;j<dim;j++)
                for(int k=0;k<dim;k++)
                    if(H_dof[j][k])
                        h_H[j][k]*=gamma;
        
        /* level 2 */
        for(int i=k_it-1;i>-1;i--)
        {
            
            tmp_vec0=y[i]->begin();
            inner0_lcl=0.0;
            for(int j=0;j<atoms->natms*x_dim;j++)
                inner0_lcl+=h[j]*tmp_vec0[j];
            
            MPI_Allreduce(&inner0_lcl,&inner0,1,MPI_TYPE0,MPI_SUM,world);
            
            if(chng_box)
                for(int j=0;j<dim;j++)
                    for(int k=0;k<dim;k++)
                        inner0+=h_H[j][k]*y_H[i][j][k];
            
            beta=-rho[i]*inner0;
            
            tmp_vec0=s[i]->begin();
            for(int j=0;j<atoms->natms*x_dim;j++)
                h[j]-=(alpha[i]-beta)*tmp_vec0[j];
            
            if(chng_box)
                for(int j=0;j<dim;j++)
                    for(int k=0;k<dim;k++)
                        h_H[j][k]-=(alpha[i]-beta)*s_H[i][j][k];
            
        }
        
        /* end of determining new h & h_H */
        
        
        /* beginning of book-keeping & linesearch */
        prev_energy=curr_energy;
        
    
        if(write!=NULL)
            write->write();

        
        thermo->thermo_print();
        
        if(affine) prepare_affine_h(x_0,h);
        err=ls->line_min(curr_energy,alpha_m,0);
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
        
        /* end of book-keeping & linesearch */
        
        
        /* beginning of calculating new f & f_H */
        force_calc();
        /* end of calculating new f & f_H */
        
        /* beginning of updating thermodynamics variables */
        if(thermo->test_prev_step() || err)
        {
            thermo->update(pe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
        }
        /* end of updating thermodynamics variables */

        istp++;
        step_no++;
        
        if(err)
            continue;
    
        
        /* beginning of storing vectors & calculating ratio for the next step */
        if(m_it)
        {
            Vec<type0>* s_tmp=s[m_it-1];
            Vec<type0>* y_tmp=y[m_it-1];

            if(chng_box)
            {
                H_s_tmp=s_H[m_it-1];
                H_y_tmp=y_H[m_it-1];
            }
            
            for(int i=m_it-1;i>0;i--)
            {
                s[i]=s[i-1];
                y[i]=y[i-1];
                if(chng_box)
                {
                    s_H[i]=s_H[i-1];
                    y_H[i]=y_H[i-1];
                }
                rho[i]=rho[i-1];
            }
            s[0]=s_tmp;
            y[0]=y_tmp;
            if(chng_box)
            {
                s_H[0]=H_s_tmp;
                y_H[0]=H_y_tmp;
            }
            
            if(k_it!=m_it)
                k_it++;
            
            x=mapp->x->begin();
            x_0=x_prev_ptr->begin();
            f=mapp->f->begin();
            f_0=f_prev_ptr->begin();
            tmp_vec0=s[0]->begin();
            tmp_vec1=y[0]->begin();
            
            inner0_lcl=0.0;
            inner1_lcl=0.0;
            for(int i=0;i<x_dim*atoms->natms;i++)
            {
                tmp_vec0[i]=x[i]-x_0[i];
                tmp_vec1[i]=f_0[i]-f[i];
                inner0_lcl+=tmp_vec0[i]*tmp_vec1[i];
                inner1_lcl+=tmp_vec1[i]*tmp_vec1[i];
            }
            MPI_Allreduce(&inner0_lcl,&inner0,1,MPI_TYPE0,MPI_SUM,world);
            MPI_Allreduce(&inner1_lcl,&inner1,1,MPI_TYPE0,MPI_SUM,world);
            
            if(chng_box)
                for(int i=0;i<dim;i++)
                    for(int j=0;j<dim;j++)
                    {
                        s_H[0][i][j]=atoms->H[i][j]-H_prev[i][j];
                        y_H[0][i][j]=f_H_prev[i][j]-f_H[i][j];
                        inner0+=s_H[0][i][j]*y_H[0][i][j];
                        inner1+=s_H[0][i][j]*s_H[0][i][j];
                    }
            gamma=inner0/inner1;
            rho[0]=1.0/inner0;
        }
        else
        {
            x=mapp->x->begin();
            x_0=x_prev_ptr->begin();
            f=mapp->f->begin();
            f_0=f_prev_ptr->begin();
            inner0_lcl=0.0;
            inner1_lcl=0.0;
            for(int i=0;i<x_dim*atoms->natms;i++)
            {
                inner0_lcl+=(x[i]-x_0[i])*(f_0[i]-f[i]);
                inner1_lcl+=(f_0[i]-f[i])*(f_0[i]-f[i]);
            }
            MPI_Allreduce(&inner0_lcl,&inner0,1,MPI_TYPE0,MPI_SUM,world);
            MPI_Allreduce(&inner1_lcl,&inner1,1,MPI_TYPE0,MPI_SUM,world);
            
            if(chng_box)
                for(int i=0;i<dim;i++)
                    for(int j=0;j<dim;j++)
                    {
                        inner0+=(atoms->H[i][j]-H_prev[i][j])*(f_H_prev[i][j]-f_H[i][j]);
                        inner1+=(f_H_prev[i][j]-f_H[i][j])*(f_H_prev[i][j]-f_H[i][j]);
                    }
            
            gamma=inner0/inner1;
        }
        /* end of storing vectors & calculating ration for the next step */
        
    }

    
}
/*--------------------------------------------
 fin after a run
 --------------------------------------------*/
void Min_lbfgs::fin()
{
    
    if(m_it)
    {
        delete [] rho;
        delete [] alpha;
    }
    
    if(chng_box)
    {
        for(int i=0;i<m_it;i++)
        {
            for(int j=0;j<dim;j++)
            {
                delete [] y_H[i][j];
                delete [] s_H[i][j];
            }
            
            delete [] y_H[i];
            delete [] s_H[i];
        }
        delete [] y_H;
        delete [] s_H;
    }
    
    if(write!=NULL)
        write->fin();
    
    thermo->fin();
    print_error();
    atoms->fin();
    
    for(int i=0;i<m_it;i++)
    {
        delete s[i];
        delete y[i];
    }
    if(m_it)
    {
        delete [] s;
        delete [] y;
    }
    
    Min::fin();

}

