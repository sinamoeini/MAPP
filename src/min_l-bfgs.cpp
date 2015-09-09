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
#include "write.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min_lbfgs::Min_lbfgs(MAPP* mapp,int narg,char** arg):Min(mapp)
{

    m_it=2;
    
    CREATE2D(H_dof,dim,dim);
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            H_dof[i][j]=0;
    
    
    int icmp;
    int jcmp;
    int iarg=2;
    while(iarg<narg)
    {
        if(!strcmp(arg[iarg],"max_iter"))
        {
            iarg++;
            if(iarg==narg)
                error->abort("max_iter in min l-bfgs should at least have 1 arguement");
            max_iter=atoi(arg[iarg]);
            iarg++;
        }
        else if(!strcmp(arg[iarg],"e_tol"))
        {
            iarg++;
            if(iarg==narg)
                error->abort("e_tol in min l-bfgs should at least have 1 arguement");
            energy_tolerance=atof(arg[iarg]);
            iarg++;
        }
        else if(!strcmp(arg[iarg],"m"))
        {
            iarg++;
            if(iarg==narg)
                error->abort("m in min l-bfgs should at least have 1 arguement");
            m_it=atoi(arg[iarg]);
            iarg++;
        }
        else if(sscanf(arg[iarg],"H[%d][%d]",&icmp,&jcmp)==2)
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
        else if(!strcmp(arg[iarg],"affine"))
        {
            affine=1;
            iarg++;
        }
        else
            error->abort("unknown keyword for min l-bfgs: %s",arg[iarg]);
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
    for(int i=0;i<dim;i++)
        delete [] H_dof[i];
    delete [] H_dof;
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void Min_lbfgs::init()
{
    
    // determine if the there is chng_box
    chng_box=0;
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            if(H_dof[i][j])
                chng_box=1;
    
    /*
     add the tensors for box if necessary
     */
    
    if(chng_box)
    {
        CREATE2D(h_H,dim,dim);
        CREATE2D(f_H,dim,dim);
        CREATE2D(f_H_prev,dim,dim);
        CREATE2D(H_prev,dim,dim);
        CREATE2D(B_prev,dim,dim);
        CREATE1D(H_y,m_it);
        CREATE1D(H_s,m_it);
        for(int i=0;i<m_it;i++)
        {
            CREATE2D(H_y[i],dim,dim);
            CREATE2D(H_s[i],dim,dim);
        }

    }
    
    CREATE1D(rho,m_it);
    CREATE1D(alpha,m_it);
    CREATE1D(s_list,m_it);
    CREATE1D(y_list,m_it);
    
    
    /* begining of chekcking if the primary atomic vectors exist */
    /* if not, add them */
    f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<type0>(0,x_dim,"f");
    
    if(mapp->mode==DMD_mode)
        c_type_n=atoms->find("c");
    else
        c_type_n=atoms->find("type");

    id_n=atoms->find("id");
    dof_n=atoms->find_exist("dof");
    /* end of chekcking if the primary atomic vectors exist */
    
    /* beginning of adding of the new atomic vectors for this min scheme */
    char* vec_name;
    
    for (int i=0;i<m_it;i++)
    {
        CREATE1D(vec_name,20);
        sprintf(vec_name,"s_%d",i);
        s_list[i]=atoms->add<type0>(0,x_dim,vec_name);
        delete [] vec_name;
        
        CREATE1D(vec_name,20);
        sprintf(vec_name,"y_%d",i);
        y_list[i]=atoms->add<type0>(0,x_dim,vec_name);
        delete [] vec_name;
    }
    
    x_prev_n=atoms->add<type0>(0,x_dim,"x_prev");
    f_prev_n=atoms->add<type0>(0,x_dim,"f_prev");
    h_n=atoms->add<type0>(0,x_dim,"h");
    /* end of adding of the new atomic vectors for this min scheme */
    
    /* begining of creating the new VecLst */
    int* tmp_lst;
    int icurs=0;
    
    if(dof_n<0)
        CREATE1D(tmp_lst,2*m_it+7);
    else
    {
        if(mapp->mode==MD_mode)
        {
            CREATE1D(tmp_lst,2*m_it+8);
            tmp_lst[icurs++]=dof_n;
            
        }
        else if (mapp->mode==DMD_mode)
        {
            CREATE1D(tmp_lst,2*m_it+9);
            tmp_lst[icurs++]=dof_n;
            cdof_n=atoms->find("cdof");
            tmp_lst[icurs++]=cdof_n;
        }
    }
    
    tmp_lst[icurs++]=0;
    tmp_lst[icurs++]=c_type_n;
    tmp_lst[icurs++]=id_n;
    tmp_lst[icurs++]=f_n;
    for (int i=0;i<m_it;i++)
    {
        tmp_lst[icurs++]=s_list[i];
        tmp_lst[icurs++]=y_list[i];
    }
    tmp_lst[icurs++]=x_prev_n;
    tmp_lst[icurs++]=f_prev_n;
    tmp_lst[icurs++]=h_n;
    
    vecs_comm=new VecLst(mapp,tmp_lst,icurs);
    delete [] tmp_lst;
    /* end of creating the new VecLst */
    
    /*
     add the update to vec list
     */
    vecs_comm->add_update(0);
    
    /*
     initiate the run
     */
    atoms->init(vecs_comm);
    
    
    /*
     do the first force calculatation for
     thermo and write
     */
    type0* f;
    atoms->vectors[f_n]->ret(f);
    for(int i=0;i<x_dim*atoms->natms;i++)
        f[i]=0.0;
    forcefield->force_calc_timer(2+affine,nrgy_strss);
    rectify(f);
    if(chng_box)
        reg_h_H(f_H);
    
    curr_energy=nrgy_strss[0];
    thermo->update(pe_idx,nrgy_strss[0]);
    thermo->update(stress_idx,6,&nrgy_strss[1]);
    
    thermo->init();
    
    if(write!=NULL)
        write->init();
    

    /*
     create the linesearch
     */

    init_linesearch();
    
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
    int s_list_tmp,y_list_tmp;
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
        atoms->vectors[h_n]->ret(h);
        atoms->vectors[0]->ret(x);
        atoms->vectors[f_n]->ret(f);
        atoms->vectors[x_prev_n]->ret(x_0);
        atoms->vectors[f_prev_n]->ret(f_0);
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
            
            atoms->vectors[s_list[i]]->ret(tmp_vec0);
            inner0_lcl=0.0;
            for(int j=0;j<atoms->natms*x_dim;j++)
                inner0_lcl+=h[j]*tmp_vec0[j];
            
            MPI_Allreduce(&inner0_lcl,&inner0,1,MPI_TYPE0,MPI_SUM,world);
            
            if(chng_box)
            {
                for(int j=0;j<dim;j++)
                    for(int k=0;k<dim;k++)
                        inner0+=h_H[j][k]*H_s[i][j][k];
            }
            
            
            alpha[i]=-rho[i]*inner0;
            
            atoms->vectors[y_list[i]]->ret(tmp_vec0);
            for(int j=0;j<atoms->natms*x_dim;j++)
                h[j]+=alpha[i]*tmp_vec0[j];
            
            if(chng_box)
                for(int j=0;j<dim;j++)
                    for(int k=0;k<dim;k++)
                        h_H[j][k]+=alpha[i]*H_y[i][j][k];
            
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
            
            atoms->vectors[y_list[i]]->ret(tmp_vec0);
            inner0_lcl=0.0;
            for(int j=0;j<atoms->natms*x_dim;j++)
                inner0_lcl+=h[j]*tmp_vec0[j];
            
            MPI_Allreduce(&inner0_lcl,&inner0,1,MPI_TYPE0,MPI_SUM,world);
            
            if(chng_box)
                for(int j=0;j<dim;j++)
                    for(int k=0;k<dim;k++)
                        inner0+=h_H[j][k]*H_y[i][j][k];
            
            beta=-rho[i]*inner0;
            
            atoms->vectors[s_list[i]]->ret(tmp_vec0);
            for(int j=0;j<atoms->natms*x_dim;j++)
                h[j]-=(alpha[i]-beta)*tmp_vec0[j];
            
            if(chng_box)
                for(int j=0;j<dim;j++)
                    for(int k=0;k<dim;k++)
                        h_H[j][k]-=(alpha[i]-beta)*H_s[i][j][k];
            
        }
        
        /* end of determining new h & h_H */
        
        
        /* beginning of book-keeping & linesearch */
        prev_energy=curr_energy;
        
        if(write!=NULL)
            write->write();
        
        thermo->thermo_print();
        
        if(affine) prepare_affine_h(x_0,h);
        err=ls->line_min(curr_energy,alpha_m,0);
        atoms->vectors[h_n]->ret(h);
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
        atoms->vectors[f_n]->ret(f);
        for(int i=0;i<x_dim*atoms->natms;i++)
            f[i]=0.0;
        forcefield->force_calc_timer(2+affine,nrgy_strss);
        rectify(f);
        if(chng_box)
            reg_h_H(f_H);
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
    
        
        /* beginning of storing vectors & calculating ration for the next step */
        if(m_it)
        {
            
            s_list_tmp=s_list[m_it-1];
            y_list_tmp=y_list[m_it-1];
            if(chng_box)
            {
                H_s_tmp=H_s[m_it-1];
                H_y_tmp=H_y[m_it-1];
            }
            
            for(int i=m_it-1;i>0;i--)
            {
                s_list[i]=s_list[i-1];
                y_list[i]=y_list[i-1];
                if(chng_box)
                {
                    H_s[i]=H_s[i-1];
                    H_y[i]=H_y[i-1];
                }
                rho[i]=rho[i-1];
            }
            s_list[0]=s_list_tmp;
            y_list[0]=y_list_tmp;
            if(chng_box)
            {
                H_s[0]=H_s_tmp;
                H_y[0]=H_y_tmp;
            }
            
            if(k_it!=m_it)
                k_it++;
            
            atoms->vectors[0]->ret(x);
            atoms->vectors[x_prev_n]->ret(x_0);
            atoms->vectors[f_n]->ret(f);
            atoms->vectors[f_prev_n]->ret(f_0);
            atoms->vectors[s_list[0]]->ret(tmp_vec0);
            atoms->vectors[y_list[0]]->ret(tmp_vec1);
            
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
                        H_s[0][i][j]=atoms->H[i][j]-H_prev[i][j];
                        H_y[0][i][j]=f_H_prev[i][j]-f_H[i][j];
                        inner0+=H_s[0][i][j]*H_y[0][i][j];
                        inner1+=H_s[0][i][j]*H_s[0][i][j];
                    }
            gamma=inner0/inner1;
            rho[0]=1.0/inner0;
        }
        else
        {
            atoms->vectors[0]->ret(x);
            atoms->vectors[x_prev_n]->ret(x_0);
            atoms->vectors[f_n]->ret(f);
            atoms->vectors[f_prev_n]->ret(f_0);
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
                delete [] H_y[i][j];
                delete [] H_s[i][j];
            }
            
            delete [] H_y[i];
            delete [] H_s[i];
        }
        delete [] H_y;
        delete [] H_s;
        
    }
    
    
    if(write!=NULL)
        write->fin();
    thermo->fin();
    print_error();
    atoms->fin();
    
    delete vecs_comm;
    
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            delete [] f_H_prev[i];
        delete [] f_H_prev;
        
        for(int i=0;i<dim;i++)
            delete [] f_H[i];
        delete [] f_H;
        
        for(int i=0;i<dim;i++)
            delete [] h_H[i];
        delete [] h_H;
        
        for(int i=0;i<dim;i++)
            delete [] H_prev[i];
        delete [] H_prev;
        
        for(int i=0;i<dim;i++)
            delete [] B_prev[i];
        delete [] B_prev;
        
    }
    
    
    atoms->del(h_n);
    atoms->del(f_prev_n);
    atoms->del(x_prev_n);
    
    int tmp;
    
    for(int i=0;i<m_it;i++)
    {
        for(int j=i+1;j<m_it;j++)
        {
            if(s_list[j]<s_list[i])
            {
                tmp=s_list[j];
                s_list[j]=s_list[i];
                s_list[i]=tmp;
            }
            
            if(y_list[j]<y_list[i])
            {
                tmp=y_list[j];
                y_list[j]=y_list[i];
                y_list[i]=tmp;
            }
        }
    }
    
    
    if(m_it)
    {
        for(int i=m_it-1;i>-1;i--)
        {
            atoms->del(y_list[i]);
            atoms->del(s_list[i]);
        }
    }

    if(m_it)
    {
        delete [] s_list;
        delete [] y_list;
    }

    

}

