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
#include "cmd.h"
#include "ls.h"
#include "dynamic.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min_lbfgs::Min_lbfgs(int nargs,char** args):Min()
{
    m_it=2;
    char* min_style;
    
    Pattern cmd;
    
    /*----------------------------*/
    cmd.cmd("min");
    cmd.add_var(min_style,"style");
    cmd.add_vdesc(0,"defines the style of minimization");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=VLogics("eq","l-bfgs");
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("e_tol");
    cmd.add_var(energy_tolerance,"tol");
    cmd.add_vdesc(0,"defines energy tolerance to stop the minimization");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=VLogics("ge",0.0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("max_iter");
    cmd.add_var(max_iter,"niter");
    cmd.add_vdesc(0,"defines maximum number of iterations");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=VLogics("ge",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("m");
    cmd.add_var(m_it,"n");
    cmd.add_vdesc(0,"defines maximum number of vectors in memory");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=VLogics("ge",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd(affine,"affine");
    cmd.add_cdesc("when used, atoms displacements will be affine");
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd_2d_lt<__dim__>(H_dof,"H");
    cmd.add_cdesc("when used, this degree of the freedom of the box would be considered in addition");
    /*------------------------------------------------------------------------------------*/
    
    cmd.scan(args,nargs);

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
    
    s=new VecTens<type0>[m_it];
    y=new VecTens<type0>[m_it];
    for(int i=0;i<m_it;i++)
    {
        s[i].init(atoms,chng_box);
        y[i].init(atoms,chng_box);
    }
    
    CREATE1D(rho,m_it);
    CREATE1D(alpha,m_it);
    
    for(int i=0;i<m_it;i++)
        dynamic->add_xchng(s[i]());
    for(int i=0;i<m_it;i++)
        dynamic->add_xchng(y[i]());
    
    dynamic->init(chng_box);
    force_calc();
    curr_energy=nrgy_strss[0];

    if(output_flag)
    {
        thermo->init();
        if(__write__)
            __write__->init();
    }
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
        
        if(output_flag)
        {
            thermo->thermo_print();
            if(__write__)
                __write__->write();
        }
        
        
        f_h=f*h;
        if(f_h<0.0)
        {
            h=f;
            k_it=0;
            f_h=f*f;
        }
        if(affine) prepare_affine_h();
        err=ls->line_min(curr_energy,alpha_m,0);
        
        if(err!=LS_S)
            continue;
        
        if(prev_energy-curr_energy<energy_tolerance)
            err=MIN_S_TOLERANCE;
        
        if(istp+1==max_iter)
            err=MIN_F_MAX_ITER;
        
        force_calc();

        
        istp++;
        if(output_flag) step_no++;
        
        if(err) continue;
        
        if(m_it)
        {
            VecTens<type0> s_tmp=std::move(s[m_it-1]);
            VecTens<type0> y_tmp=std::move(y[m_it-1]);
            
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
        
    if(output_flag)
    {
        if(__write__)
            __write__->fin();
        thermo->fin();
    }
    dynamic->fin();
    Min::fin();
}

