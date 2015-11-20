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
    x.init(atoms,atoms->x,atoms->H,chng_box);
    x0.init(atoms,x_prev_ptr,H_prev,chng_box);
    f.init(atoms,mapp->f,f_H,chng_box);
    f0.init(atoms,f_prev_ptr,f_H_prev,chng_box);
    h.init(atoms,h_ptr,h_H,chng_box);
    
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
    
    type0 prev_energy;
    type0 alpha;
    int istp=0;
    type0 f0_f0;
    type0 f_f;
    type0 f_h;
    type0 f_f0;
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
        
        if(write!=NULL)
            write->write();
        
        thermo->thermo_print();
        
        if(affine) prepare_affine_h(x0()->begin(),h()->begin());
        err=ls->line_min(curr_energy,alpha,1);
        if(affine) rectify(h()->begin());
        
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
        
        f_f=f*f;
        f_f0=f*f0;
        
        ratio=(f_f-f_f0)/(f0_f0);
        
        h*=ratio;
        h+=f;
        f_h*=ratio;
        f_h+=f_f;
        
        if(f_h<0.0)
        {
            h=f;
            f_h=f_f;
        }
        f0_f0=f_f;
    }
}
/*--------------------------------------------
 finishing minimization
 --------------------------------------------*/
void Min_cg::fin()
{
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

