/*--------------------------------------------
 Created by Sina on 06/20/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <stdlib.h>
#include "md_nh.h"
#include "ff.h"
#include "random.h"
#include "atom_types.h"
#include "error.h"
#include "memory.h"
#include "write.h"
#include "timer.h"
#include "neighbor.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
enum {NONE,TAU,XYZ,YZ,ZX,XY};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MD_nh::MD_nh(MAPP* mapp,int nargs,char** args)
: MD(mapp)
{
    //the defaults
    nreset=10;
    no_it_eta=1;
    no_ch_eta=3;
    no_it_peta=1;
    no_ch_peta=3;
    crt_vel=false;
    t_tar=boltz=dt=0.0;
    stress_mode=NONE;
    for(int i=0;i<6;i++)
    {
        H_dof[i]=false;
        tau_tar[i]=tau_freq[i]=0.0;
    }
    
    if(atoms->dimension!=3)
        error->abort("dimension of the box should be 3 for md nh");
    
    int iarg=2;
    if(!strcmp(args[iarg],"nvt"))
    {
        stress_mode=NONE;
        iarg++;
    }
    else if(!strcmp(args[iarg],"ntaut"))
    {
        stress_mode=TAU;
        iarg++;
    }
    else if(!strcmp(args[iarg],"npt")&&((iarg+2)<=nargs))
    {
        
        iarg++;
        if(!strcmp(args[iarg],"xyz")
           || !strcmp(args[iarg],"xzy")
           || !strcmp(args[iarg],"yzx")
           || !strcmp(args[iarg],"yxz")
           || !strcmp(args[iarg],"zxy")
           || !strcmp(args[iarg],"zyx"))
        {
            stress_mode=XYZ;
            iarg++;
        }
        else if(!strcmp(args[iarg],"yz")
                || !strcmp(args[iarg],"zy"))
        {
            stress_mode=YZ;
            iarg++;
        }
        else if(!strcmp(args[iarg],"zx")
                || !strcmp(args[iarg],"xz"))
        {
            stress_mode=ZX;
            iarg++;
        }
        else if(!strcmp(args[iarg],"xy")
                || !strcmp(args[iarg],"yx"))
        {
            stress_mode=XY;
            iarg++;
        }
        else
            error->abort("unknown coupling style for npt in md nh: %s",args[iarg]);
        
    }
    else error->abort("unknown ensemble for md nh: %s",args[iarg]);


    while(iarg<nargs)
        keywords(nargs,args,iarg);


}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
MD_nh::~MD_nh()
{
}
/*--------------------------------------------
 setup before start of run
 --------------------------------------------*/
void MD_nh::init()
{
    
    if(t_tar==0)
        error->abort("temp was not set by md nh");
    if(stress_mode!=NONE)
    {
        int tmp=1;
        if(stress_mode==TAU)
        {
            tmp=0;
            for(int i=0;i<6;i++)
                if(H_dof[i])
                    tmp++;
        }
        else if(stress_mode==XYZ)
        {
            if(!H_dof[0])
                tmp=0;
            if(!H_dof[1])
                tmp=0;
            if(!H_dof[2])
                tmp=0;
        }
        else if(stress_mode==YZ)
        {
            if(!H_dof[1])
                tmp=0;
            if(!H_dof[2])
                tmp=0;
            
        }
        else if(stress_mode==ZX)
        {
            if(!H_dof[0])
                tmp=0;
            if(!H_dof[2])
                tmp=0;
            
        }
        else if(stress_mode==XY)
        {
            if(!H_dof[0])
                tmp=0;
            if(!H_dof[1])
                tmp=0;
            
        }
        
        if(tmp==0)
        {
            if(stress_mode==TAU)
                error->abort("no stress was not set by md nh");
            else
                error->abort("average stress was not set by md nh");
        }
    }
    
    if(boltz==0.0)
        error->abort("boltzmann should be set after md nh and before run");
    if(dt==0.0)
        error->abort("time_step should be set after md nh and before run");
    
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
    }
    
    CREATE1D(eta_d,no_ch_eta);
    CREATE1D(eta_dd,no_ch_eta);
    CREATE1D(eta_m,no_ch_eta);

    for(int i=0;i<no_ch_eta;i++)
        ke_vec[i]=eta_d[i]=0.0;
    
    
    no_dof=static_cast<type0>(atoms->tot_natms*3-3);
    ke_tar=t_tar*boltz*no_dof;
    dt2=0.5*dt;
    dt4=0.25*dt;
    dt8=0.125*dt;

    if(mapp->x_d==NULL)
        mapp->x_d=new Vec<type0>(atoms,3);


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
        
        no_dof=static_cast<type0>(ndofs[0]+ndofs[1]+ndofs[2])-3.0;
        
        
        omega_denom=0.0;
        for(int i=0;i<3;i++)
            if(H_dof[i])
                omega_denom+=static_cast<type0>(ndofs[i]);
    }
    else
    {
        if(stress_mode)
        {
            omega_denom=0.0;
            for(int i=0;i<3;i++)
                if(H_dof[i])
                    omega_denom+=static_cast<type0>(atoms->tot_natms);
        }
    }
    
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
        create_vel(seed,t_tar);
    else
        init_vel(t_tar);
    
    vol=1.0;
    for(int idim=0;idim<3;idim++)
        vol*=atoms->H[idim][idim];
    
    memset(forcefield->f->begin(),0,atoms->natms*3*sizeof(type0));
    forcefield->force_calc_timer(true);
    for(int i=0;i<6;i++)
        virial_pe[i]=-nrgy_strss[i+1];
    if(dof_xst) modify_vrial();    
    
    for(int j=0;j<6;j++)
    {
        nrgy_strss[1+j]*=-1.0;
        nrgy_strss[1+j]+=ke_vec[j]/vol;
    }
    thermo->update(stress_idx,6,&nrgy_strss[1]);
    thermo->update(pe_idx,nrgy_strss[0]);
    thermo->update(ke_idx,ke_cur);
    thermo->update(temp_idx,t_cur);
    
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
            peta_m[i] = boltz * t_tar/(tau_freq_m*tau_freq_m);
        
        for(int i=0;i<no_ch_peta;i++)
            peta_d[i]=0.0;
        
        for(int i=1;i<no_ch_peta;i++)
            peta_dd[i] =(peta_m[i-1]*peta_d[i-1]*peta_d[i-1]-boltz*t_tar)/peta_m[i];
        for(int i=0;i<6;i++)
            omega_d[i]=omega_m[i]=0.0;
    }
    

}
/*--------------------------------------------
 finalize after the run is complete
 --------------------------------------------*/
void MD_nh::fin()
{
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
            atoms->update(mapp->x);
            
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
            
            atoms->update(mapp->x);
            
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
    type0* x=mapp->x->begin();
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
    
    x=mapp->x->begin();
    for(int i=0;i<natms;i++)
    {
        x[0]-=dx_ave[0];
        x[1]-=dx_ave[1];
        x[2]-=dx_ave[2];
        x+=3;
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
            
            if(!dof[0]) x_d[0]+=f[0]*dlt/m;
            if(!dof[1]) x_d[1]+=f[1]*dlt/m;
            if(!dof[2]) x_d[2]+=f[2]*dlt/m;
            
            if(dof[0]) ke_vec_lcl[0]+=m*x_d[0]*x_d[0];
            if(dof[1]) ke_vec_lcl[1]+=m*x_d[1]*x_d[1];
            if(dof[2]) ke_vec_lcl[2]+=m*x_d[2]*x_d[2];
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
        
        exfac=exp(-dltm4*eta_d[1]);
        eta_d[0]*=exfac;
        eta_dd[0]=(ke_cur-ke_tar)/eta_m[0];
        eta_d[0]+=eta_dd[0]*dltm2;
        eta_d[0]*=exfac;
        
        for(int ich=1;ich<no_ch_eta;ich++)
        {
            if(ich==no_ch_eta-1) exfac=1.0;
            else exfac=exp(-dltm4*eta_d[ich+1]);
            
            eta_d[ich]*=exfac;
            eta_dd[ich]=(eta_m[ich-1]*eta_d[ich-1]*eta_d[ich-1]-boltz*t_tar)/eta_m[ich];
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
            peta_d[ich]*=exfac;
            exfac=exp(-dltm4*peta_d[ich]);
        }
        
        exfac=exp(-dltm*peta_d[0]);
        for(int i=0;i<6;i++)
            if(H_dof[i])
                omega_d[i]*=exfac;
        
        kec*=exfac*exfac;
        
        exfac=exp(-dltm4*peta_d[1]);
        peta_d[0]*=exfac;
        peta_dd[0]=(kec-dof*boltz*t_tar)/peta_m[0];
        peta_d[0]+=peta_dd[0]*dltm2;
        peta_d[0]*=exfac;
        
        for(int ich=1;ich<no_ch_peta;ich++)
        {
            if(ich==no_ch_peta-1) exfac=1.0;
            else exfac=exp(-dltm4*peta_d[ich+1]);
            
            peta_d[ich]*=exfac;
            peta_dd[ich]=(peta_m[ich-1]*peta_d[ich-1]*peta_d[ich-1]-boltz*t_tar)/peta_m[ich];
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
    
    class Random* random=new Random(mapp,seed);

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
void MD_nh::keywords(int nargs,char** args,int& iarg)
{
    int icmp,jcmp;
    if(!strcmp(args[iarg],"temp")&&((iarg+2)<=nargs))
    {
        iarg++;
        if(nargs-iarg<2)
            error->abort("temp in md nh should at least have 2 arguments");
        
        t_tar=atof(args[iarg]);
        if(t_tar<=0)
            error->abort("temp in md nh should be greater than 0.0");
        iarg++;
        if(atof(args[iarg])<=0)
            error->abort("temp frequency (arguement 2 after temp) in md nh should be greater than 0.0");
        t_freq=1.0/atof(args[iarg]);
        iarg++;
    }
    else if(!strcmp(args[iarg],"stress")&&((iarg+7)<=nargs))
    {
        iarg++;
        if(stress_mode==NONE)
            error->abort("stress in md nh is valid for ntaut or npt ensemble");
        if(stress_mode==TAU)
        {
            if(nargs-iarg<12)
                error->abort("stress in md nh should at least have 12 arguments");
            for(int i=0;i<6;i++)
            {
                tau_tar[i]=atof(args[iarg]);
                iarg++;
                
                if(atof(args[iarg])<=0)
                    error->abort("stress frequency (arguement %d after stress) in md nh should be greater than 0.0",(i+1)*2);
                tau_freq[i]=1.0/atof(args[iarg]);
                
                H_dof[i]=true;
                iarg++;
            }
        }
        else
        {
            if(nargs-iarg<2)
                error->abort("stress in md nh should at least have 2 arguments");
            if(stress_mode==XYZ)
            {
                tau_tar[0]=tau_tar[1]=tau_tar[2]=atof(args[iarg]);
                H_dof[0]=H_dof[1]=H_dof[2]=true;
                iarg++;
                if(atof(args[iarg])<=0)
                    error->abort("stress frequency (arguement 2 after ave) in md nh should be greater than 0.0");
                tau_freq[0]=tau_freq[1]=tau_freq[2]=1.0/atof(args[iarg]);
                iarg++;
            }
            else if(stress_mode==YZ)
            {
                tau_tar[1]=tau_tar[2]=atof(args[iarg]);
                H_dof[1]=H_dof[2]=true;
                iarg++;
                if(atof(args[iarg])<=0)
                    error->abort("stress frequency (arguement 2 after stress) in md nh should be greater than 0.0");
                tau_freq[1]=tau_freq[2]=1.0/atof(args[iarg]);
                iarg++;
            }
            else if(stress_mode==ZX)
            {
                tau_tar[0]=tau_tar[2]=atof(args[iarg]);
                H_dof[0]=H_dof[2]=true;
                iarg++;
                if(atof(args[iarg])<=0)
                    error->abort("stress frequency (arguement 2 after stress) in md nh should be greater than 0.0");
                tau_freq[0]=tau_freq[2]=1.0/atof(args[iarg]);
                iarg++;
            }
            else if(stress_mode==XY)
            {
                tau_tar[0]=tau_tar[1]=atof(args[iarg]);
                H_dof[0]=H_dof[1]=true;
                iarg++;
                if(atof(args[iarg])<=0)
                    error->abort("stress frequency (arguement 2 after stress) in md nh should be greater than 0.0");
                tau_freq[0]=tau_freq[1]=1.0/atof(args[iarg]);
                iarg++;
            }
        }
        
    }
    else if(!strcmp(args[iarg],"eta_iter")&&((iarg+2)<=nargs))
    {
        iarg++;
        no_it_eta=atoi(args[iarg]);
        if(no_it_eta<1)
            error->abort("eta_iter in md nh should be greater than 0");
        iarg++;
    }
    else if(!strcmp(args[iarg],"peta_iter")&&((iarg+2)<=nargs))
    {
        if(stress_mode==NONE)
            error->abort("peta_iter in md nh is valid for npt ntaut ensemble");
        iarg++;
        no_it_peta=atoi(args[iarg]);
        if(no_it_peta<1)
            error->abort("peta_iter in md nh should be greater than 0");
        iarg++;
    }
    else if(!strcmp(args[iarg],"eta_chains")&&((iarg+2)<=nargs))
    {
        iarg++;
        no_ch_eta=atoi(args[iarg]);
        if(no_ch_eta<3)
            error->abort("eta_chains in md nh should be greater than 2");
        iarg++;
    }
    else if(!strcmp(args[iarg],"peta_chains")&&((iarg+2)<=nargs))
    {
        if(stress_mode==NONE)
            error->abort("peta_chains in md nh is valid for npt ntaut ensemble");
        iarg++;
        no_ch_peta=atoi(args[iarg]);
        if(no_ch_peta<3)
            error->abort("peta_chains in md nh should be greater than 2");
        iarg++;
    }
    else if(!strcmp(args[iarg],"create_vel")&&((iarg+2)<=nargs))
    {
        crt_vel=true;
        iarg++;
        seed=atoi(args[iarg]);
        iarg++;
    }
    else if(sscanf(args[iarg],"T[%d][%d]",&icmp,&jcmp)==2)
    {
        if(stress_mode!=TAU)
            error->abort("T[%d][%d] in md nh must be used with ntaut ensemble",icmp,jcmp);
        if(icmp<0 || icmp>=3)
            error->abort("wrong component in md nh for T[%d][%d]",icmp,jcmp);
        if(jcmp<0 || jcmp>=3)
            error->abort("wrong component in md nh for T[%d][%d]",icmp,jcmp);
        iarg++;
        if(nargs-iarg<2)
            error->abort("T[%d][%d] in md nh must have 2 arguments",icmp,jcmp);
        
        int voigt[3][3]={{0,5,4},{5,1,3},{4,3,2}};
        int idx=voigt[icmp][jcmp];
        
        tau_tar[idx]=atof(args[iarg]);
        iarg++;
        if(atof(args[iarg])<=0)
            error->abort("stress frequency (arguement 2 after T[%d][%d]) in md nh should be greater than 0.0",icmp,jcmp);
        tau_freq[idx]=1.0/atof(args[iarg]);
        iarg++;
        H_dof[idx]=true;
        
    }
    else if(!strcmp(args[iarg],"nreset"))
    {
        iarg++;
        if(nargs-iarg<1)
            error->abort("nreset in md nh must have 1 rguments");
        if(atoi(args[iarg])<0)
            error->abort("nreset in md nh should be greater than 0");
        nreset=atoi(args[iarg]);
        iarg++;
    }
    else error->abort("unknown keyword for md nh: %s",args[iarg]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_nh::modify_vrial()
{
    type0* xvec=mapp->x->begin();
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

