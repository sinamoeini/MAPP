/*--------------------------------------------
 Created by Sina on 06/20/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <stdlib.h>
#include "md_nh.h"
#include "atoms.h"
#include "ff.h"
#include "random.h"
#include "atom_types.h"
#include "error.h"
#include "memory.h"
#include "write.h"
#include "timer.h"
#include "neighbor.h"
#include "thermo_dynamics.h"
#include "cmd.h"
#include "pgcmc.h"
#include "MAPP.h"


using namespace MAPP_NS;
enum {NONE,X,Y,XY,Z,ZX,YZ,XYZ,TAU};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MD_nh::MD_nh(int nargs,char** args):MD()
{
    //the defaults
    gcmc=NULL;
    seed=0;
    xchng_seed=0;
    drag=0.0;
    nreset=1;
    no_it_eta=1;
    no_ch_eta=3;
    no_it_peta=1;
    no_ch_peta=3;
    crt_vel=false;
    dof_adj[0]=dof_adj[1]=dof_adj[2]=true;
    t_tar=0.0;
    stress_mode=NONE;
    for(int i=0;i<6;i++)
    {
        H_dof[i]=false;
        tau_tar[i]=tau_freq[i]=0.0;
    }

    cmd(nargs,args);
    if(xchng_seed)
    {
        gcmc=new PGCMC(1,gas_type,mu,t_tar,xchng_seed);
        count_idx=thermo->add("count");
    }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
MD_nh::~MD_nh()
{
    delete gcmc;
}
/*--------------------------------------------
 setup before start of run
 --------------------------------------------*/
void MD_nh::init()
{
    
    if(boltz==0.0)
        error->abort("boltzmann should be set after md nh and before run");
    if(gcmc && hplanck==0.0)
        error->abort("plank should be set after md nh and before run");
    if(dt==0.0)
        error->abort("time_step should be set after md nh and before run");
    
    tdrag=1.0-dt*t_freq*drag/static_cast<type0>(no_ch_eta);
    
    if(stress_mode)
    {
        CREATE1D(peta_d,no_ch_peta);
        CREATE1D(peta_dd,no_ch_peta);
        CREATE1D(peta_m,no_ch_peta);
        for(int i=0;i<no_ch_peta;i++)
            peta_d[i]=0.0;
        tau_freq_m=0.0;
        for(int i=0;i<6;i++)
            if(H_dof[i])
                tau_freq_m=MAX(tau_freq_m,tau_freq[i]);
        
        pdrag=1.0-dt*tau_freq_m*drag/static_cast<type0>(no_ch_peta);
    }
    
    CREATE1D(eta_d,no_ch_eta);
    CREATE1D(eta_dd,no_ch_eta);
    CREATE1D(eta_m,no_ch_eta);

    for(int i=0;i<no_ch_eta;i++)
        ke_vec[i]=eta_d[i]=0.0;
    
    

    dt2=0.5*dt;
    dt4=0.25*dt;
    dt8=0.125*dt;

    if(mapp->x_d==NULL)
    {
        mapp->x_d=new Vec<type0>(atoms,3);
        memset(mapp->x_d->begin(),0,atoms->natms*3*sizeof(type0));
    }


    dof_xst=false;
    if(mapp->x_dof!=NULL)
    {
        dof_xst=true;
        int ndofs_lcl[3];
        ndofs_lcl[0]=ndofs_lcl[1]=ndofs_lcl[2]=0;
        bool* dof=mapp->x_dof->begin();

        for(int i=0;i<atoms->natms;i++)
        {
            if(dof[0]) ndofs_lcl[0]++;
            if(dof[1]) ndofs_lcl[1]++;
            if(dof[2]) ndofs_lcl[2]++;
            
            dof+=3;
        }
        int ndofs[3];
        MPI_Allreduce(ndofs_lcl,ndofs,3,MPI_INT,MPI_SUM,world);
        
        for(int i=0;i<3;i++)
            if(ndofs[i])
                dof_adj[i]=false;
        
        no_dof=static_cast<type0>(ndofs[0]+ndofs[1]+ndofs[2])-3.0;
        
        
        omega_denom=0.0;
        for(int i=0;i<3;i++)
            if(H_dof[i])
                omega_denom+=static_cast<type0>(ndofs[i]);
    }
    else
    {
        no_dof=static_cast<type0>(atoms->tot_natms*3-3);
        if(stress_mode)
        {
            omega_denom=0.0;
            for(int i=0;i<3;i++)
                if(H_dof[i])
                    omega_denom+=static_cast<type0>(atoms->tot_natms);
        }
    }

    ke_tar=t_tar*boltz*no_dof;
    
    if(no_dof==0.0)
        error->abort("degrees of freedom shoud be greater than 0 for md nh");
    
    vecs_comm=new VecLst(atoms);
    vecs_comm->add_updt(mapp->type);
    vecs_comm->add_xchng(mapp->x_d);
    
    if(mapp->x_dof!=NULL)
        vecs_comm->add_xchng(mapp->x_dof);
    
    if(stress_mode)
        atoms->init(vecs_comm,true);
    else
        atoms->init(vecs_comm,false);

    if(crt_vel)
    {
        create_vel(seed,t_tar);
        crt_vel=false;
    }
    
    vol=1.0;
    for(int idim=0;idim<3;idim++)
        vol*=atoms->H[idim][idim];
    
    memset(forcefield->f->begin(),0,atoms->natms*3*sizeof(type0));
    forcefield->force_calc_timer(true);

    if(stress_mode)
    {
        for(int i=0;i<6;i++)
            virial_pe[i]=-nrgy_strss[i+1];
        if(dof_xst) modify_vrial();
    }

    
    for(int j=0;j<6;j++)
    {
        nrgy_strss[1+j]*=-1.0;
        nrgy_strss[1+j]+=ke_vec[j]/vol;
    }
    thermo->update(stress_idx,6,&nrgy_strss[1]);
    thermo->update(pe_idx,nrgy_strss[0]);
    thermo->update(ke_idx,ke_cur);
    thermo->update(temp_idx,t_cur);
    if(gcmc)
    {
        gcmc->init();
        thermo->update(count_idx,static_cast<type0>(gcmc->tot_ngas));
    }
    
    thermo->init();
    
    if(write!=NULL)
        write->init();
    
    for(int i=0;i<no_ch_eta;i++)
        eta_m[i]=boltz*t_tar/(t_freq*t_freq);
    eta_m[0]*=no_dof;

    for(int i=0;i<no_ch_eta;i++)
        eta_d[i]=0.0;
    
    for(int i=1;i<no_ch_eta;i++)
        eta_dd[i] =(eta_m[i-1]*eta_d[i-1]*eta_d[i-1]-boltz*t_tar)/eta_m[i];
    
    if(stress_mode)
    {
        for(int i=0;i<no_ch_peta;i++)
            peta_m[i]=boltz*t_tar/(tau_freq_m*tau_freq_m);
        
        for(int i=0;i<no_ch_peta;i++)
            peta_d[i]=0.0;
        
        for(int i=1;i<no_ch_peta;i++)
            peta_dd[i]=(peta_m[i-1]*peta_d[i-1]*peta_d[i-1]-boltz*t_tar)/peta_m[i];
        for(int i=0;i<6;i++)
            omega_d[i]=omega_m[i]=0.0;
    }
    
    
}
/*--------------------------------------------
 finalize after the run is complete
 --------------------------------------------*/
void MD_nh::fin()
{
    if(gcmc)
        gcmc->fin();
    if(write!=NULL)
        write->fin();
    thermo->fin();
    atoms->fin();
    timer->print_stats();
    neighbor->print_stats();
    
    delete vecs_comm;
    
    delete [] eta_d;
    delete [] eta_dd;
    delete [] eta_m;
    
    if(stress_mode)
    {
        delete [] peta_d;
        delete [] peta_dd;
        delete [] peta_m;
    }
}
/*--------------------------------------------
 MDrun
 --------------------------------------------*/
void MD_nh::run(int no_stps)
{
    
    if(stress_mode)
    {
        for(int i=0;i<no_stps;i++)
        {
            
            update_NH_tau(dt2);
            update_NH_T(dt2);
            calc_sigma(i);
            update_omega_d(dt2);
            update_x_d_xpnd(dt2);
            update_x_d(dt2);
            
            update_x(dt);
            atoms->update(atoms->x);
            
            thermo->thermo_print();
            if(write!=NULL) write->write();

            memset(forcefield->f->begin(),0,atoms->natms*3*sizeof(type0));
            forcefield->force_calc_timer(true);
            
            vol=1.0;
            for(int idim=0;idim<3;idim++)
                vol*=atoms->H[idim][idim];
            
            for(int j=0;j<6;j++)
                virial_pe[j]=-nrgy_strss[1+j];
            if(dof_xst) modify_vrial();
            
            update_x_d(dt2);
            update_x_d_xpnd(dt2);
            update_omega_d(dt2);
            update_NH_T(dt2);
            update_NH_tau(dt2);
            
            if(thermo->test_prev_step()|| i==no_stps-1)
            {
                for(int j=0;j<6;j++)
                {
                    nrgy_strss[1+j]*=-1.0;
                    nrgy_strss[1+j]+=ke_vec[j]/vol;
                }
                
                thermo->update(stress_idx,6,&nrgy_strss[1]);
                thermo->update(pe_idx,nrgy_strss[0]);
                thermo->update(ke_idx,ke_cur);
                thermo->update(temp_idx,t_cur);
                
            }
            step_no++;
        }
    }
    else
    {

        for(int i=0;i<no_stps;i++)
        {
            update_NH_T(dt2);
            update_x_d(dt2);
            update_x(dt);
            
            if(gcmc && i%xchng_every==0)
            {
                gcmc->xchng(false,nxchng_attmpts);
                no_dof+=static_cast<type0>(gcmc->dof_diff);
                ke_tar=t_tar*boltz*no_dof;
                thermo->update(count_idx,static_cast<type0>(gcmc->tot_ngas));
                if(no_dof<=0.0)
                    error->abort("degrees of freedom shoud be greater than 0 for md nh");
            }
            else
               atoms->update(atoms->x);
            
            memset(forcefield->f->begin(),0,atoms->natms*3*sizeof(type0));
            thermo->thermo_print();
            
            if(write!=NULL)
                write->write();
            
            if(thermo->test_prev_step()|| i==no_stps-1)
                 forcefield->force_calc_timer(true);
            else
                forcefield->force_calc_timer(false);

            update_x_d(dt2);
            update_NH_T(dt2);
            
            if(thermo->test_prev_step()|| i==no_stps-1)
            {
                for(int j=0;j<6;j++)
                {
                    nrgy_strss[1+j]*=-1.0;
                    nrgy_strss[1+j]+=ke_vec[j]/vol;
                }
                
                thermo->update(stress_idx,6,&nrgy_strss[1]);
                thermo->update(pe_idx,nrgy_strss[0]);
                thermo->update(ke_idx,ke_cur);
                thermo->update(temp_idx,t_cur);
            }
            
            step_no++;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_nh::update_H(type0 dlt)
{
    type0 dlt2,dlt4,dlt8,exfac;
    type0 H0[3][3];
    type0 B0[3][3];
    dlt2=0.5*dlt;
    dlt4=0.25*dlt;
    dlt8=0.125*dlt;
    
    
    M3EQV(atoms->H,H0);
    for(int i=0;i<2;i++)
    {
        //M3INV_LT(H0,B0);
        
        if(H_dof[4])
        {
            exfac=exp(dlt8*omega_d[0]);
            H0[2][0]*=exfac;
            H0[2][0]+=dlt4*(omega_d[5]*H0[2][1]
            +omega_d[4]*H0[2][2]);
            H0[2][0]*=exfac;
        }
        
        if(H_dof[3])
        {
            exfac=exp(dlt4*omega_d[1]);
            H0[2][1]*=exfac;
            H0[2][1]+=dlt2*omega_d[3]*H0[2][2];
            H0[2][1]*=exfac;
        }
        
        if(H_dof[5])
        {
            exfac=exp(dlt4*omega_d[0]);
            H0[1][0]*=exfac;
            H0[1][0]+=dlt2*omega_d[5]*H0[1][1];
            H0[1][0]*=exfac;
        }
        
        if(H_dof[4])
        {
            exfac=exp(dlt8*omega_d[0]);
            H0[2][0]*=exfac;
            H0[2][0]+=dlt4*(omega_d[5]*H0[2][1]+omega_d[4]*H0[2][2]);
            H0[2][0]*=exfac;
        }
        
        if(H_dof[0])
        {
            exfac=exp(dlt*omega_d[0]);
            H0[0][0]*=exfac;
        }
        
        if(H_dof[1])
        {
            exfac=exp(dlt*omega_d[1]);
            H0[1][1]*=exfac;
            H0[1][0]*=exfac;
        }
        
        if(H_dof[2])
        {
            exfac=exp(dlt*omega_d[2]);
            H0[2][2]*=exfac;
            H0[2][1]*=exfac;
            H0[2][0]*=exfac;
        }
        
        if(H_dof[4])
        {
            exfac=exp(dlt8*omega_d[0]);
            H0[2][0]*=exfac;
            H0[2][0]+=dlt4*(omega_d[5]*H0[2][1]+omega_d[4]*H0[2][2]);
            H0[2][0]*=exfac;
        }
        
        if(H_dof[5])
        {
            exfac=exp(dlt4*omega_d[0]);
            H0[1][0]*=exfac;
            H0[1][0]+=dlt2*omega_d[5]*H0[1][1];
            H0[1][0]*=exfac;
        }
        
        if(H_dof[3])
        {
            exfac=exp(dlt4*omega_d[1]);
            H0[2][1]*=exfac;
            H0[2][1]+=dlt2*omega_d[3]*H0[2][2];
            H0[2][1]*=exfac;
        }
        
        if(H_dof[4])
        {
            exfac=exp(dlt8*omega_d[0]);
            H0[2][0]*=exfac;
            H0[2][0]+=dlt4*(omega_d[5]*H0[2][1]+omega_d[4]*H0[2][2]);
            H0[2][0]*=exfac;
        }
        if(i==0)
            M3INV_LT(H0,B0);
    }
    
    M3MUL_LT(atoms->B,H0,M[0]);
    M3MUL_LT(B0,H0,M[1]);
    
    M3EQV(H0,atoms->H);
    M3INV_LT(atoms->H,atoms->B);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_nh::update_x(type0 dlt)
{
    type0* x=atoms->x->begin();
    type0* x_d=mapp->x_d->begin();
    int natms=atoms->natms;

    for(int i=0;i<3;i++)
        dx_ave_lcl[i]=0.0;
    
    if(stress_mode)
    {
        type0 x_0[3];
        type0 x_1[3];
        update_H(0.5*dlt);
        
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                M[1][i][j]*=dlt;
        
        if(dof_xst)
        {
            bool* dof=mapp->x_dof->begin();
            
            for(int i=0;i<natms;i++)
            {
                for(int j=0;j<3;j++)
                    dx_ave_lcl[j]-=x[j];
                
                V3M_LT(x,M[0],x_0);
                V3M_LT(x_d,M[1],x_1);
                
                for(int j=0;j<3;j++)
                    if(dof[j])
                        x[j]=x_0[j]+x_1[j];
                
                for(int j=0;j<3;j++)
                    dx_ave_lcl[j]+=x[j];
                
                x+=3;
                x_d+=3;
                dof+=3;
            }
        }
        else
        {
            for(int i=0;i<natms;i++)
            {
                for(int j=0;j<3;j++)
                    dx_ave_lcl[j]-=x[j];
                
                V3M_LT(x,M[0],x_0);
                V3M_LT(x_d,M[1],x_1);
                
                for(int j=0;j<3;j++)
                    x[j]=x_0[j]+x_1[j];
                
                for(int j=0;j<3;j++)
                    dx_ave_lcl[j]+=x[j];
                
                x+=3;
                x_d+=3;
            }
        }

    }
    else
    {
        if(dof_xst)
        {
            bool* dof=mapp->x_dof->begin();
            
            for(int i=0;i<natms;i++)
            {
                for(int j=0;j<3;j++)
                    if(dof[j])
                    {
                        x[j]+=x_d[j]*dlt;
                        dx_ave_lcl[j]+=x_d[j]*dlt;
                    }
                
                x+=3;
                x_d+=3;
                dof+=3;
            }
        }
        else
        {
            for(int i=0;i<natms;i++)
            {
                for(int j=0;j<3;j++)
                {
                    x[j]+=x_d[j]*dlt;
                    dx_ave_lcl[j]+=x_d[j]*dlt;
                }
                
                x+=3;
                x_d+=3;
            }
        }
    }

    
    MPI_Allreduce(dx_ave_lcl,dx_ave,3,MPI_TYPE0,MPI_SUM,world);
    dx_ave[0]*=1.0/(atoms->tot_natms);
    dx_ave[1]*=1.0/(atoms->tot_natms);
    dx_ave[2]*=1.0/(atoms->tot_natms);
    

    if(dof_adj[0] || dof_adj[1] || dof_adj[2])
    {
        for(int i=0;i<3;i++)
            if(!dof_adj[i])
                dx_ave[i]=0.0;
            
        x=atoms->x->begin();
        for(int i=0;i<natms;i++)
        {
            x[0]-=dx_ave[0];
            x[1]-=dx_ave[1];
            x[2]-=dx_ave[2];
            x+=3;
        }
    }

}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_nh::update_x_d(type0 dlt)
{
    
    type0* x_d=mapp->x_d->begin();
    type0* f=forcefield->f->begin();
    md_type* type=mapp->type->begin();
    type0* mass=atom_types->mass;
    type0 m;
    int natms=atoms->natms;
    
    for(int i=0;i<6;i++) ke_vec_lcl[i]=0.0;

    if(dof_xst)
    {
        bool* dof=mapp->x_dof->begin();
        
        for(int i=0;i<natms;i++)
        {
            m=mass[type[i]];
            
            if(dof[0])
            {
                x_d[0]+=f[0]*dlt/m;
                ke_vec_lcl[0]+=m*x_d[0]*x_d[0];
            }
            if(dof[1])
            {
                x_d[1]+=f[1]*dlt/m;
                ke_vec_lcl[1]+=m*x_d[1]*x_d[1];
            }
            if(dof[2])
            {
                x_d[2]+=f[2]*dlt/m;
                ke_vec_lcl[2]+=m*x_d[2]*x_d[2];
            }
            if(dof[1] && dof[2]) ke_vec_lcl[3]+=m*x_d[1]*x_d[2];
            if(dof[0] && dof[2]) ke_vec_lcl[4]+=m*x_d[0]*x_d[2];
            if(dof[0] && dof[1]) ke_vec_lcl[5]+=m*x_d[0]*x_d[1];
            
            x_d+=3;
            f+=3;
            dof+=3;
        }
    }
    else
    {
        for(int i=0;i<natms;i++)
        {
            
            m=mass[type[i]];
            
            x_d[0]+=f[0]*dlt/m;
            x_d[1]+=f[1]*dlt/m;
            x_d[2]+=f[2]*dlt/m;

            ke_vec_lcl[0]+=m*x_d[0]*x_d[0];
            ke_vec_lcl[1]+=m*x_d[1]*x_d[1];
            ke_vec_lcl[2]+=m*x_d[2]*x_d[2];
            ke_vec_lcl[3]+=m*x_d[1]*x_d[2];
            ke_vec_lcl[4]+=m*x_d[0]*x_d[2];
            ke_vec_lcl[5]+=m*x_d[0]*x_d[1];
            
            x_d+=3;
            f+=3;
        }
    }

    MPI_Allreduce(ke_vec_lcl,ke_vec,6,MPI_TYPE0,MPI_SUM,world);
    ke_cur=(ke_vec[0]+ke_vec[1]+ke_vec[2]);

    t_cur=ke_cur/(boltz*no_dof);
}
/*--------------------------------------------
 Nosé–Hoover thermostat chains 
 --------------------------------------------*/
void MD_nh::update_NH_T(type0 dlt)
{
    type0 dltm,dltm2,dltm4,exfac,velfac;
    
    for(int i=0;i<no_ch_eta;i++)
        eta_m[i]=boltz*t_tar/(t_freq*t_freq);
    eta_m[0]*=no_dof;
    
    dltm=dlt/static_cast<type0>(no_it_eta);
    dltm2=0.5*dltm;
    dltm4=0.25*dltm;
    velfac=1.0;
    
    eta_dd[0]=(ke_cur-ke_tar)/eta_m[0];
    for(int it=0;it<no_it_eta;it++)
    {
        exfac=1.0;
        for(int ich=no_ch_eta-1;ich>-1;ich--)
        {
            eta_d[ich]*=exfac;
            eta_d[ich]+=eta_dd[ich]*dltm2;
            eta_d[ich]*=tdrag;
            eta_d[ich]*=exfac;
            exfac=exp(-dltm4*eta_d[ich]);
        }

        //rescale x_d dlt & claculate the new temperature
        exfac=exp(-dltm*eta_d[0]);
        velfac*=exfac;
        t_cur*=exfac*exfac;
        ke_cur*=exfac*exfac;
        ke_vec[0]*=exfac*exfac;
        ke_vec[1]*=exfac*exfac;
        ke_vec[2]*=exfac*exfac;
        ke_vec[3]*=exfac*exfac;
        ke_vec[4]*=exfac*exfac;
        ke_vec[5]*=exfac*exfac;
        
        
        for(int ich=0;ich<no_ch_eta;ich++)
        {
            if(ich==no_ch_eta-1) exfac=1.0;
            else exfac=exp(-dltm4*eta_d[ich+1]);
            
            eta_d[ich]*=exfac;
            if(ich) eta_dd[ich]=(eta_m[ich-1]*eta_d[ich-1]*eta_d[ich-1]-boltz*t_tar)/eta_m[ich];
            else eta_dd[ich]=(ke_cur-ke_tar)/eta_m[ich];
            eta_d[ich]+=eta_dd[ich]*dltm2;
            eta_d[ich]*=exfac;
        }
    }

    int natms=atoms->natms;
    type0* x_d=mapp->x_d->begin();
    if(dof_xst)
    {
        bool* dof=mapp->x_dof->begin();
        for(int i=0;i<natms*3;i++)
            if(dof[i])
                x_d[i]*=velfac;
    }
    else
    {
        for(int i=0;i<natms*3;i++)
            x_d[i]*=velfac;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_nh::update_NH_tau(type0 dlt)
{
    type0 dltm,dltm2,dltm4,exfac,kec;
    type0 dof=0;
    
    for(int i=0;i<6;i++)
        if(H_dof[i]) dof++;
    
    for(int i=0;i<6;i++)
        if(H_dof[i])
            omega_m[i]=boltz*t_tar/(tau_freq[i]*tau_freq[i]);
    
    for(int i=0;i<no_ch_peta;i++)
        peta_m[i]=boltz*t_tar/(tau_freq_m*tau_freq_m);
    
    for(int i=1;i<no_ch_peta;i++)
        peta_dd[i]=(peta_m[i-1]*peta_d[i-1]*peta_d[i-1]-boltz*t_tar)/peta_m[i];
    
    kec=0.0;
    for(int i=0;i<6;i++)
        if(H_dof[i])
            kec+=omega_m[i]*omega_d[i]*omega_d[i];
    
    peta_dd[0]=(kec-dof*boltz*t_tar)/peta_m[0];
    
    dltm=dlt/static_cast<type0>(no_it_peta);
    dltm2=0.5*dltm;
    dltm4=0.25*dltm;
    
    for(int it=0;it<no_it_peta;it++)
    {
        exfac=1.0;
        for(int ich=no_ch_peta-1;ich>-1;ich--)
        {
            peta_d[ich]*=exfac;
            peta_d[ich]+=peta_dd[ich]*dltm2;
            peta_d[ich]*=pdrag;
            peta_d[ich]*=exfac;
            exfac=exp(-dltm4*peta_d[ich]);
        }
        
        exfac=exp(-dltm*peta_d[0]);
        for(int i=0;i<6;i++)
            if(H_dof[i])
                omega_d[i]*=exfac;
        
        kec*=exfac*exfac;
        
        for(int ich=0;ich<no_ch_peta;ich++)
        {
            if(ich==no_ch_peta-1) exfac=1.0;
            else exfac=exp(-dltm4*peta_d[ich+1]);
            
            peta_d[ich]*=exfac;
            if(ich) peta_dd[ich]=(peta_m[ich-1]*peta_d[ich-1]*peta_d[ich-1]-boltz*t_tar)/peta_m[ich];
            else peta_dd[ich]=(kec-dof*boltz*t_tar)/peta_m[ich];
            peta_d[ich]+=peta_dd[ich]*dltm2;
            peta_d[ich]*=exfac;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_nh::update_omega_d(type0 dlt)
{
    MTK_1=0.0;
    if(omega_denom)
    {
        for(int i=0;i<3;i++)
            if(H_dof[i])
                MTK_1+=ke_vec[i];
        MTK_1/=omega_denom;
    }

    if(stress_mode==XYZ)
    {
        type0 tmp=(virial_pe[0]+virial_pe[1]+virial_pe[2])/3.0;
        virial_pe[0]=virial_pe[1]=virial_pe[2]=tmp;
        tmp=(ke_vec[0]+ke_vec[1]+ke_vec[2])/3.0;
        ke_vec[0]=ke_vec[1]=ke_vec[2]=tmp;
    }
    else if(stress_mode==YZ)
    {
        type0 tmp=(virial_pe[1]+virial_pe[2])/2.0;
        virial_pe[1]=virial_pe[2]=tmp;
        tmp=(ke_vec[1]+ke_vec[2])/2.0;
        ke_vec[1]=ke_vec[2]=tmp;
    }
    else if(stress_mode==ZX)
    {
        type0 tmp=(virial_pe[2]+virial_pe[0])/2.0;
        virial_pe[2]=virial_pe[0]=tmp;
        tmp=(ke_vec[2]+ke_vec[0])/2.0;
        ke_vec[2]=ke_vec[0]=tmp;
    }
    else if(stress_mode==XY)
    {
        type0 tmp=(virial_pe[0]+virial_pe[1])/2.0;
        virial_pe[0]=virial_pe[1]=tmp;
        tmp=(ke_vec[0]+ke_vec[1])/2.0;
        ke_vec[0]=ke_vec[1]=tmp;
    }

    for(int i=0;i<3;i++)
        if(H_dof[i])
            omega_d[i]+=((virial_pe[i]-tau_hydro)*vol+ke_vec[i]
            -sigma_tar[i]*vol0-MTK_1)*dlt/omega_m[i];
    
    for(int i=3;i<6;i++)
        if(H_dof[i])
            omega_d[i]+=(virial_pe[i]*vol+ke_vec[i]
            -sigma_tar[i]*vol0-MTK_1)*dlt/omega_m[i];

    MTK_2=0.0;
    if(omega_denom)
    {
        for(int i=0;i<3;i++)
            if(H_dof[i])
                MTK_2+=omega_d[i];
        MTK_2/=omega_denom;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_nh::update_x_d_xpnd(type0 dlt)
{
    type0 fac[3];
    int natms=atoms->natms;
    
    fac[0]=exp(0.5*dlt*(omega_d[0]+MTK_2));
    fac[1]=exp(0.5*dlt*(omega_d[1]+MTK_2));
    fac[2]=exp(0.5*dlt*(omega_d[2]+MTK_2));
    
    for(int i=0;i<6;i++) ke_vec_lcl[i]=0.0;
    
    type0* x_d=mapp->x_d->begin();
    md_type* type=mapp->type->begin();
    type0* mass=atom_types->mass;
    type0 m;
    
    if(dof_xst)
    {
        bool* dof=mapp->x_dof->begin();
        type0 x_d_[3];
        
        for(int i=0;i<natms;i++)
        {
            if(!dof[0]) x_d_[0]=x_d[0];
            if(!dof[1]) x_d_[1]=x_d[1];
            if(!dof[2]) x_d_[2]=x_d[2];
            
            x_d[0]*=fac[0];
            x_d[1]*=fac[1];
            x_d[2]*=fac[2];
            
            x_d[0]-=dlt*(x_d[1]*omega_d[5]
            +x_d[2]*omega_d[4]);
            x_d[1]-=dlt*x_d[2]*omega_d[3];
            
            x_d[0]*=fac[0];
            x_d[1]*=fac[1];
            x_d[2]*=fac[2];
            
            if(!dof[0]) x_d[0]=x_d_[0];
            if(!dof[1]) x_d[1]=x_d_[1];
            if(!dof[2]) x_d[2]=x_d_[2];
            
            m=mass[type[i]];
            if(dof[0]) ke_vec_lcl[0]+=m*x_d[0]*x_d[0];
            if(dof[1]) ke_vec_lcl[1]+=m*x_d[1]*x_d[1];
            if(dof[2]) ke_vec_lcl[2]+=m*x_d[2]*x_d[2];
            if(dof[1] && dof[2]) ke_vec_lcl[3]+=m*x_d[1]*x_d[2];
            if(dof[0] && dof[2]) ke_vec_lcl[4]+=m*x_d[0]*x_d[2];
            if(dof[0] && dof[1]) ke_vec_lcl[5]+=m*x_d[0]*x_d[1];
            
            x_d+=3;
            dof+=3;
        }
    }
    else
    {
        for(int i=0;i<natms;i++)
        {
            x_d[0]*=fac[0];
            x_d[1]*=fac[1];
            x_d[2]*=fac[2];
            
            x_d[0]-=dlt*(x_d[1]*omega_d[5]
            +x_d[2]*omega_d[4]);
            x_d[1]-=dlt*x_d[2]*omega_d[3];
            
            x_d[0]*=fac[0];
            x_d[1]*=fac[1];
            x_d[2]*=fac[2];
            
            m=mass[type[i]];
            
            ke_vec_lcl[0]+=m*x_d[0]*x_d[0];
            ke_vec_lcl[1]+=m*x_d[1]*x_d[1];
            ke_vec_lcl[2]+=m*x_d[2]*x_d[2];
            ke_vec_lcl[3]+=m*x_d[1]*x_d[2];
            ke_vec_lcl[4]+=m*x_d[0]*x_d[2];
            ke_vec_lcl[5]+=m*x_d[0]*x_d[1];
            
            x_d+=3;
        }
    }

    MPI_Allreduce(ke_vec_lcl,ke_vec,6,MPI_TYPE0,MPI_SUM,world);
    ke_cur=(ke_vec[0]+ke_vec[1]+ke_vec[2]);
    t_cur=ke_cur/(boltz*no_dof);
}
/*--------------------------------------------
 create initial velocity 
 --------------------------------------------*/
void MD_nh::create_vel(int seed,type0 temperature)
{
    type0* x_d=mapp->x_d->begin();
    md_type* type=mapp->type->begin();
    type0* mass=atom_types->mass;
    type0 m;
    
    int natms=atoms->natms;
    for(int i=0;i<6;i++)
        ke_vec_lcl[i]=0.0;
    
    class Random* random=new Random(seed);

    if(dof_xst)
    {
        bool* dof=mapp->x_dof->begin();
        for(int i=0;i<natms;i++)
        {
            for(int j=0;j<3;j++)
                if(dof[j])
                    x_d[j]=random->gaussian()/(sqrt(mass[type[i]]));
            
            m=mass[type[i]];
            if(dof[0]) ke_vec_lcl[0]+=m*x_d[0]*x_d[0];
            if(dof[1]) ke_vec_lcl[1]+=m*x_d[1]*x_d[1];
            if(dof[2]) ke_vec_lcl[2]+=m*x_d[2]*x_d[2];
            if(dof[1] && dof[2]) ke_vec_lcl[3]+=m*x_d[1]*x_d[2];
            if(dof[0] && dof[2]) ke_vec_lcl[4]+=m*x_d[0]*x_d[2];
            if(dof[0] && dof[1]) ke_vec_lcl[5]+=m*x_d[0]*x_d[1];
            
            x_d+=3;
            dof+=3;
        }
    }
    else
    {
        for(int i=0;i<natms;i++)
        {
            for(int j=0;j<3;j++)
                x_d[j]=random->gaussian()/(sqrt(mass[type[i]]));
            m=mass[type[i]];
            ke_vec_lcl[0]+=m*x_d[0]*x_d[0];
            ke_vec_lcl[1]+=m*x_d[1]*x_d[1];
            ke_vec_lcl[2]+=m*x_d[2]*x_d[2];
            ke_vec_lcl[3]+=m*x_d[1]*x_d[2];
            ke_vec_lcl[4]+=m*x_d[0]*x_d[2];
            ke_vec_lcl[5]+=m*x_d[0]*x_d[1];
            
            x_d+=3;
        }
    }

    delete random;

    MPI_Allreduce(ke_vec_lcl,ke_vec,6,MPI_TYPE0,MPI_SUM,world);
    
    ke_cur=(ke_vec[0]+ke_vec[1]+ke_vec[2]);
    t_cur=ke_cur/(boltz*no_dof);
    
    type0 ke_des=(boltz*no_dof)*temperature;
    type0 factor=sqrt(ke_des/ke_cur);
    type0 facsq=ke_des/ke_cur;

    if(dof_xst)
    {
        bool* dof=mapp->x_dof->begin();
        x_d=mapp->x_d->begin();
        for(int i=0;i<natms*3;i++)
            if(dof[i]) x_d[i]*=factor;
    }
    else
    {
        x_d=mapp->x_d->begin();
        for(int i=0;i<natms*3;i++)
            x_d[i]*=factor;
    }

    for(int i=0;i<6;i++)
        ke_vec[i]*=facsq;

    ke_cur*=facsq;
    t_cur*=facsq;

}
/*--------------------------------------------
 create initial velocity
 --------------------------------------------*/
void MD_nh::init_vel(type0 temperature)
{

    type0* x_d=mapp->x_d->begin();
    md_type* type=mapp->type->begin();
    type0* mass=atom_types->mass;
    type0 m;
    
    int natms=atoms->natms;
    
    for(int i=0;i<6;i++)
        ke_vec_lcl[i]=0.0;
    
    if(dof_xst)
    {
        bool* dof=mapp->x_dof->begin();
        for(int i=0;i<natms;i++)
        {
            m=mass[type[i]];
            if(dof[0]) ke_vec_lcl[0]+=m*x_d[0]*x_d[0];
            if(dof[1]) ke_vec_lcl[1]+=m*x_d[1]*x_d[1];
            if(dof[2]) ke_vec_lcl[2]+=m*x_d[2]*x_d[2];
            if(dof[1] && dof[2]) ke_vec_lcl[3]+=m*x_d[1]*x_d[2];
            if(dof[0] && dof[2]) ke_vec_lcl[4]+=m*x_d[0]*x_d[2];
            if(dof[0] && dof[1]) ke_vec_lcl[5]+=m*x_d[0]*x_d[1];
            
            x_d+=3;
            dof+=3;
        }
    }
    else
    {
        for(int i=0;i<natms;i++)
        {
            m=mass[type[i]];
            ke_vec_lcl[0]+=m*x_d[0]*x_d[0];
            ke_vec_lcl[1]+=m*x_d[1]*x_d[1];
            ke_vec_lcl[2]+=m*x_d[2]*x_d[2];
            ke_vec_lcl[3]+=m*x_d[1]*x_d[2];
            ke_vec_lcl[4]+=m*x_d[0]*x_d[2];
            ke_vec_lcl[5]+=m*x_d[0]*x_d[1];
            
            x_d+=3;
        }
    }

    MPI_Allreduce(ke_vec_lcl,ke_vec,6,MPI_TYPE0,MPI_SUM,world);
    
    ke_cur=(ke_vec[0]+ke_vec[1]+ke_vec[2]);
    t_cur=ke_cur/(boltz*no_dof);
    
    if(ke_cur==0.0)
        error->abort("kinetic energy of the system should be "
        "greater than 0.0 for md nh, please assign velocities "
        "or use create_vel keyword");
    
}
/*--------------------------------------------
 calculate the target sigma
 --------------------------------------------*/
void MD_nh::calc_sigma(int istep)
{
    type0** H=atoms->H;

    if((nreset && istep%nreset==0) || istep==0)
    {
        M3EQV(atoms->B,B_ref);
        vol0=H[0][0]*H[1][1]*H[2][2];
    }
    
    type0 m[6];
    M3MUL_LT_2_voigt(B_ref,H,m);
    
    tau_hydro=(tau_tar[0]+tau_tar[1]+tau_tar[2])/3.0;
    
    sigma_tar[0]
    =m[0]*m[0]*(tau_tar[0]-tau_hydro)
    +m[5]*m[5]*(tau_tar[1]-tau_hydro)
    +m[4]*m[4]*(tau_tar[2]-tau_hydro)
    +m[4]*m[5]*tau_tar[3]*2.0
    +m[0]*m[4]*tau_tar[4]*2.0
    +m[0]*m[5]*tau_tar[5]*2.0;
    
    sigma_tar[1]
    =m[1]*m[1]*(tau_tar[1]-tau_hydro)
    +m[3]*m[3]*(tau_tar[2]-tau_hydro)
    +m[1]*m[3]*tau_tar[3]*2.0;
    
    sigma_tar[2]
    =m[2]*m[2]*(tau_tar[2]-tau_hydro);
    
    sigma_tar[3]
    =m[2]*m[3]*(tau_tar[2]-tau_hydro)
    +m[1]*m[2]*tau_tar[3];
    
    sigma_tar[4]
    =m[2]*m[4]*(tau_tar[2]-tau_hydro)
    +m[2]*m[5]*tau_tar[3]
    +m[2]*m[0]*tau_tar[4];
    
    sigma_tar[5]
    =m[1]*m[5]*(tau_tar[1]-tau_hydro)
    +m[3]*m[4]*(tau_tar[2]-tau_hydro)
    +(m[1]*m[4]+m[3]*m[5])*tau_tar[3]
    +m[0]*m[3]*tau_tar[4]
    +m[0]*m[1]*tau_tar[5];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_nh::modify_vrial()
{
    type0* xvec=atoms->x->begin();
    type0* fvec=forcefield->f->begin();
    bool* dof=mapp->x_dof->begin();
    type0 st_lcl[6];
    type0 st[6];
    for(int i=0;i<6;i++)
        st_lcl[i]=0.0;
    
    for(int i=0;i<atoms->natms;i++)
    {
        if(!dof[0])
            st_lcl[0]+=fvec[0]*xvec[0];
            st_lcl[4]+=fvec[0]*xvec[2];
            st_lcl[5]+=fvec[0]*xvec[1];
        
        if(!dof[1])
            st_lcl[1]+=fvec[1]*xvec[1];
            st_lcl[3]+=fvec[1]*xvec[2];
        
        if(!dof[2])
            st_lcl[2]+=fvec[2]*xvec[2];
        
        fvec+=3;
        xvec+=3;
        dof+=3;
    }
    
    MPI_Allreduce(st_lcl,st,6,MPI_TYPE0,MPI_SUM,world);
    
    for(int i=0;i<6;i++)
        virial_pe[i]-=st[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_nh::cmd(int nargs,char** args)
{
    type0 tau_tar_ave,tau_freq_ave;
    char* ensemble=NULL;
    char* couple=NULL;
    char* gcmc_type=NULL;
    
    Pattern cmd;
    
    /*----------------------------*/
    cmd.cmd("nh");
    cmd.add_var(ensemble,"ensemble");
    cmd.add_vdesc(0,"defines the ensemble of the md simulation");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=vlogic("eq","nvt")
    +vlogic("eq","npt")+vlogic("eq","ntaut");
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("couple");
    cmd.add_var(couple,"form");
    cmd.add_vdesc(0,"defines coupling form of npt ensemble");
    /*--------------------------------------------------------*/
    cmd.add_clog()=logic(ensemble,"!eq","npt")/vlogic("!set");
    cmd.add_vlog(0)=
    vlogic("eq","xyz")+vlogic("eq","xzy")
    +vlogic("eq","yzx")+vlogic("eq","yxz")
    +vlogic("eq","zxy")+vlogic("eq","zxy")
    +vlogic("eq","xy")+vlogic("eq","yx")
    +vlogic("eq","xz")+vlogic("eq","zx")
    +vlogic("eq","yz")+vlogic("eq","zy")
    +vlogic("eq","x")+vlogic("eq","y")
    +vlogic("eq","z");
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("temp");
    cmd.add_var(t_tar,"t_tar",t_freq,"t_per");
    cmd.add_vdesc(0,"defines the target temperature");
    cmd.add_vdesc(1,"defines the period of nose hoover bath");
    /*--------------------------------------------------------*/
    cmd.add_clog()=vlogic("set");
    cmd.add_vlog(0)=vlogic("gt",0.0);
    cmd.add_vlog(1)=vlogic("gt",0.0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("stress");
    cmd.add_var(tau_tar_ave,"tau_tar_ave",tau_freq_ave,"tau_per");
    cmd.add_vdesc(0,"defines the target average stress");
    cmd.add_vdesc(1,"defines the period of nose hoover bath for average stress");
    /*--------------------------------------------------------*/
    cmd.add_clog()=logic(ensemble,"eq","npt")-vlogic("set");
    cmd.add_vlog(1)=vlogic("gt",0.0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd_voigt("tau",3);
    cmd.add_var(tau_tar,"tau_tar",tau_freq,"tau_per");
    cmd.add_vdesc(0,"defines the target stress");
    cmd.add_vdesc(1,"defines the period of nose hoover bath for said stress");
    /*--------------------------------------------------------*/
    cmd.add_clog()=logic(ensemble,"eq","ntaut")-vlogic("set");
    cmd.add_vlog(1)=vlogic("gt",0.0);
    /*------------------------------------------------------------------------------------*/
    
    
    /*----------------------------*/
    cmd.cmd("eta_iter");
    cmd.add_var(no_it_eta,"eta_iters");
    cmd.add_vdesc(0,"defines no. of iterations in nose hoover bath");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=vlogic("gt",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("eta_chains");
    cmd.add_var(no_ch_eta,"neta_chains");
    cmd.add_vdesc(0,"defines no. of masses in nose hoover chains bath");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=vlogic("gt",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("peta_iter");
    cmd.add_var(no_it_peta,"eta_iters");
    cmd.add_vdesc(0,"defines no. of iterations in nose hoover bath for stress/pressure");
    /*--------------------------------------------------------*/
    cmd.add_clog()=(logic(ensemble,"!eq","npt")+logic(ensemble,"!eq","ntaut"))/vlogic("!set");
    cmd.add_vlog(0)=vlogic("gt",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("peta_chains");
    cmd.add_var(no_ch_peta,"neta_chains");
    cmd.add_vdesc(0,"defines no. of masses in nose hoover chains bath for stress/pressure");
    /*--------------------------------------------------------*/
    cmd.add_clog()=(logic(ensemble,"!eq","npt")+logic(ensemble,"!eq","ntaut"))/vlogic("!set");
    cmd.add_vlog(0)=vlogic("gt",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("create_vel");
    cmd.add_var(seed,"random_seed");
    cmd.add_vdesc(0,"defines the random seed for gaussian distribution");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=vlogic("gt",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("nreset");
    cmd.add_var(nreset,"n");
    cmd.add_vdesc(0,"defines the number of steps to restart");
    /*--------------------------------------------------------*/
    cmd.add_clog()=(logic(ensemble,"!eq","npt")+logic(ensemble,"!eq","ntaut"))/vlogic("!set");
    cmd.add_vlog(0)=vlogic("gt",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("drag");
    cmd.add_var(drag,"fac");
    cmd.add_vdesc(0,"defines the drag factor for damped md");
    /*--------------------------------------------------------*/
    cmd.add_clog()=(logic(ensemble,"!eq","npt")+logic(ensemble,"!eq","ntaut"))/vlogic("!set");
    cmd.add_vlog(0)=vlogic("gt",0.0);
    /*------------------------------------------------------------------------------------*/
    
    
    
    /*----------------------------*/
    cmd.cmd("mu");
    cmd.add_var(mu,"chemical_pot",gcmc_type,"exchange_element",xchng_seed,"exchange_seed",xchng_every,"N",
    nxchng_attmpts,"nxchng");
    cmd.add_vdesc(0,"defines the chemical potential");
    cmd.add_vdesc(1,"defines element to be exchanged");
    cmd.add_vdesc(2,"defines the random seed for exchange");
    cmd.add_vdesc(3,"defines every number steps to do exchange");
    cmd.add_vdesc(4,"defines number of exchange attempts");
    /*--------------------------------------------------------*/
    cmd.add_clog()=logic(ensemble,"!eq","nvt")/vlogic("!set");
    //cmd.add_vlog(1)=vlogic("eq",atom_types->atom_names[0]);
    //for(int i=1;i<atom_types->no_types;i++) cmd.add_vlog(1)+=vlogic("eq",atom_types->atom_names[i]);
    cmd.add_vlog(2)=vlogic("gt",0);
    cmd.add_vlog(3)=vlogic("gt",0);
    /*------------------------------------------------------------------------------------*/
    
    
    args++;
    nargs--;
    cmd.scan(args,nargs);
    
    t_freq=1.0/t_freq;
    tau_freq_ave=1.0/tau_freq_ave;
    for(int i=0;i<6;i++) tau_freq[i]=1.0/tau_freq[i];
    
    if(!strcmp(ensemble,"npt") && couple!=NULL)
    {
        
        stress_mode=0;
        for(int i=0;i<strlen(couple);i++)
        {
            if(couple[i]=='x')
            {
                stress_mode+=1;
                H_dof[0]=true;
                tau_tar[0]=tau_tar_ave;
                tau_freq[0]=tau_freq_ave;
            }
            else if(couple[i]=='y')
            {
                stress_mode+=2;
                H_dof[1]=true;
                tau_tar[1]=tau_tar_ave;
                tau_freq[1]=1.0/tau_freq_ave;
            }
            else if(couple[i]=='z')
            {
                stress_mode+=4;
                H_dof[2]=true;
                tau_tar[2]=tau_tar_ave;
                tau_freq[2]=tau_freq_ave;
            }
        }
    }
    else
    {
        if(!strcmp(ensemble,"ntaut"))
            stress_mode=TAU;
        else if(!strcmp(ensemble,"npt"))
        {
            stress_mode=XYZ;
            H_dof[0]=H_dof[1]=H_dof[2]=true;
            tau_tar[0]=tau_tar[1]=tau_tar[2]=tau_tar_ave;
            tau_freq[0]=tau_freq[1]=tau_freq[2]=tau_freq_ave;
        }
        else
            stress_mode=NONE;
        
    }
    
    if(seed)
        crt_vel=true;
    if(gcmc_type!=NULL)
        gas_type=atom_types->find_type(gcmc_type);
}

