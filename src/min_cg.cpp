#include <stdlib.h>
#include "memory.h"
#include "error.h"
#include "min_cg.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
#include "cmd.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min_cg::Min_cg(MAPP* mapp,int nargs,char** args):Min(mapp)
{
    char* min_style=NULL;
    
    Pattern cmd(error);
    
    /*----------------------------*/
    cmd.cmd("min");
    cmd.add_var(min_style,"style");
    cmd.add_vdesc(0,"defines the style of minimization");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=vlogic("eq","cg");
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("e_tol");
    cmd.add_var(energy_tolerance,"tol");
    cmd.add_vdesc(0,"defines energy tolerance to stop the minimization");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=vlogic("ge",0.0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd("max_iter");
    cmd.add_var(max_iter,"niter");
    cmd.add_vdesc(0,"defines maximum number of iterations");
    /*--------------------------------------------------------*/
    cmd.add_vlog(0)=vlogic("ge",0);
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd(affine,"affine");
    cmd.add_cdesc("when used, atoms displacements will be affine");
    /*------------------------------------------------------------------------------------*/
    
    /*----------------------------*/
    cmd.cmd_2d_lt(H_dof,"H",3);
    cmd.add_cdesc("when used, this degree of the freedom of the box would be considered in addition");
    /*------------------------------------------------------------------------------------*/
    
    cmd.scan(args,nargs);
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
    if(output_flag)
    {
        thermo->init();
        if(write!=NULL)
            write->init();
    }
}
/*--------------------------------------------
 min
 --------------------------------------------*/
void Min_cg::run()
{
    if(max_iter==0) return;
    
    curr_energy=nrgy_strss[0];
    type0 prev_energy;
    type0 alpha;
    int istp=0;
    type0 f0_f0,f_f,f_f0;
    type0 ratio;
    err=LS_S;

    h=f;
    f_h=f0_f0=f*f;
    
    while(err==LS_S)
    {
        if(f0_f0==0.0)
        {
            err=LS_F_GRAD0;
            continue;
        }
        
        x0=x;
        f0=f;
        
        prev_energy=curr_energy;
        
        if(output_flag)
        {
            thermo->thermo_print();
            if(write!=NULL)
                write->write();
        }
        
        f_h=f*h;
        if(f_h<0.0)
        {
            h=f;
            f_h=f0_f0;
        }
        if(affine) prepare_affine_h();
        err=ls->line_min(curr_energy,alpha,1);
       
        if(err!=LS_S)
            continue;
        
        if(prev_energy-curr_energy<energy_tolerance)
            err=MIN_S_TOLERANCE;
        
        if(istp+1==max_iter)
            err=MIN_F_MAX_ITER;
        
        force_calc();
        
        istp++;
        if(output_flag) step_no++;
        
        if(err)
            continue;

        f_f=f*f;
        f_f0=f*f0;
        
        ratio=(f_f-f_f0)/(f0_f0);
        
        h*=ratio;
        h+=f;

        f0_f0=f_f;
    }
}
/*--------------------------------------------
 finishing minimization
 --------------------------------------------*/
void Min_cg::fin()
{
    if(output_flag)
    {
        if(write!=NULL)
            write->fin();
        thermo->fin();
    }
        
    atoms->fin();
    Min::fin();
}

