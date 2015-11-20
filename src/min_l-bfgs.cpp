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
    
    x.init(atoms,atoms->x,atoms->H,chng_box);
    x0.init(atoms,x_prev_ptr,H_prev,chng_box);
    f.init(atoms,mapp->f,f_H,chng_box);
    f0.init(atoms,f_prev_ptr,f_H_prev,chng_box);
    h.init(atoms,h_ptr,h_H,chng_box);
    
    s=new vvec<type0>[m_it];
    y=new vvec<type0>[m_it];
    for(int i=0;i<m_it;i++)
    {
        s[i].init(atoms,chng_box);
        y[i].init(atoms,chng_box);
    }
    
    CREATE1D(rho,m_it);
    CREATE1D(alpha,m_it);
    
    for(int i=0;i<m_it;i++)
        vecs_comm->add_xchng(s[i]());
    for(int i=0;i<m_it;i++)
        vecs_comm->add_xchng(y[i]());
    
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
    
    type0 prev_energy;
    type0 alpha_m,gamma;
    type0 inner0,inner1;
    
    
    
    int k_it=0;
    int istp=0;
    
    gamma=1.0;
    err=LS_S;
    
    while(err==LS_S)
    {
        x0=x;
        h=f0=f;
        
        for(int i=0;i<k_it;i++)
        {
            alpha[i]=-rho[i]*(s[i]*h);
            h.add(alpha[i],y[i]);
        }
        
        h*=gamma;
        
        for(int i=k_it-1;i>-1;i--)
            h.add(-alpha[i]-rho[i]*(y[i]*h),s[i]);
        
        prev_energy=curr_energy;
        
        if(write!=NULL)
            write->write();
        
        thermo->thermo_print();
        
        if(affine) prepare_affine_h(x_prev_ptr->begin(),h_ptr->begin());
        err=ls->line_min(curr_energy,alpha_m,0);
        if(affine) rectify(h_ptr->begin());
        
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
        
        if(err) continue;
        
        if(m_it)
        {
            vvec<type0> s_tmp=std::move(s[m_it-1]);
            vvec<type0> y_tmp=std::move(y[m_it-1]);
            
            for(int i=m_it-1;i>0;i--)
            {
                s[i]=std::move(s[i-1]);
                y[i]=std::move(y[i-1]);
                rho[i]=rho[i-1];
            }
            s[0]=std::move(s_tmp);
            y[0]=std::move(y_tmp);
            
            if(k_it!=m_it)
                k_it++;
            
            s[0]=x;
            s[0]-=x0;
            y[0]=f0;
            y[0]-=f;
            
            inner0=s[0]*y[0];
            inner1=y[0]*y[0];
            
            gamma=inner0/inner1;
            rho[0]=1.0/inner0;
        }
        else
        {
            gamma=(x*f0-x*f-x0*f0+x0*f)/(f*f+f0*f0-2.0*(f*f0));
        }
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
    
    if(m_it)
    {
        delete [] s;
        delete [] y;
    }
    
    x.fin();
    x0.fin();
    f.fin();
    f0.fin();
    h.fin();
    
    
    if(write!=NULL)
        write->fin();
    
    thermo->fin();
    print_error();
    atoms->fin();
    Min::fin();

}

