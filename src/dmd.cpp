#include "dmd.h"
#include "thermo_dynamics.h"
#include "atom_types.h"
#include "write.h"
#include "error.h"
#include "memory.h"
#include "timer.h"
#include "min_styles.h"
#include "neighbor_dmd.h"
#include "gmres.h"
#include "MAPP.h"
#include "script_reader.h"
#include <limits>
//#define DMD_DEBUG
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMD::DMD():
nrgy_strss(forcefield->nrgy_strss)
{
    if(forcefield==NULL)
        error->abort("ff should be "
        "initiated before dmd");
    
    if(mode!=DMD_mode)
        error->abort("dmd works only "
        "for md mode");
    
    char** args;
    int nargs=ScriptReader::parse_line(
    "Time FE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    
    time_idx=0;
    //msd_idx=1;
    fe_idx=1;
    stress_idx=2;
    
    thermo=new ThermoDynamics(nargs,args);
    for(int i=0;i<nargs;i++)
        delete [] args[i];
    if(nargs)
        delete [] args;

    
    
    c_dim=mapp->c->dim;
    forcefield_dmd=dynamic_cast<ForceFieldDMD*>(forcefield);
    neighbor_dmd=dynamic_cast<Neighbor_dmd*>(neighbor);
    
    a_tol=sqrt(2.0*std::numeric_limits<type0>::epsilon());
    dt_min=std::numeric_limits<type0>::epsilon();
    eps=std::numeric_limits<type0>::epsilon();
    eps_sqr=sqrt(2.0*std::numeric_limits<type0>::epsilon());
    inf=std::numeric_limits<type0>::infinity();
    t_fin=1.0e7;
    f_tol=1.0e-5;
    cd_tol=1.0e-4;
    min_flag=NO_FLAG;
    min=NULL;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
DMD::~DMD()
{
    delete thermo;
    delete min;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD::coef(int nargs,char** args)
{
    if(!strcmp(args[1],"f_tol"))
    {
        f_tol=atof(args[2]);
        min_flag=F_FLAG;
    }
    else if(!strcmp(args[1],"cd_tol"))
    {
        cd_tol=atof(args[2]);
        min_flag=CD_FLAG;
    }
    else if(!strcmp(args[1],"step_tol"))
    {
        step_tol=atoi(args[2]);
        min_flag=STEP_FLAG;
    }
    else
        error->abort("unknown keyword");
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void DMD::dmd_min(int nargs,char** args)
{
    char** args_=new char*[nargs];
    size_t len=strlen("min")+1;
    args_[0]=new char[len];
    memcpy(args_[0],"min",len);
    for(int i=1;i<nargs;i++)
    {
        len=strlen(args[i])+1;
        args_[i]=new char[len];
        memcpy(args_[i],args[i],len);
    }
    
    #define Min_Style
    #define MinStyle(class_name,style_name)   \
    else if(strcmp(args[1],#style_name)==0)   \
    {if(min)delete min;                 \
    min= new class_name(nargs,args_);}

    if(0){}
    #include "min_styles.h"
    else
        error->abort("wrong style of minimization"
        ": %s",args[1]);
    #undef Min_Style
    #undef MinStyle
    min->output_flag=false;
    
    for(int i=0;i<nargs;i++)
        delete [] args_[i];
    delete [] args_;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD::run(type0 t)
{
    if(t<=0.0)
        error->abort("t_fin in dmd should be greater than 0.0");
    if(dt_min>t)
        error->abort("t_fin in dmd should be greater than dt_min");
    t_fin=t;
    run();
}
/*--------------------------------------------

 --------------------------------------------*/
void DMD::do_min()
{
    forcefield_dmd->dynamic_flag=true;
    min->force_calc();
    min->run();
    nmin++;
    if(min_flag==F_FLAG)
        f_norm0=min->calc_ave_f_norm();
}
/*--------------------------------------------
 rectify
 --------------------------------------------*/
bool DMD::decide_min(int& istep,type0& del_t)
{
    max_succ_dt=MAX(max_succ_dt,del_t);
    t_cur+=del_t;
    if(write)
        write->write();
    thermo->thermo_print();
    
    bool thermo_flag=(thermo->test_prev_step()
    || istep==max_step-1
    || t_cur>=t_fin);
    istep++;
    step_no++;
    if(min_flag==F_FLAG)
    {
        type0 f_norm=min->calc_ave_f_norm();
        if(thermo_flag)
        {
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            thermo->update(time_idx,t_cur);
        }

        if(f_norm-f_norm0>=f_tol)
        {
            return true;
        }
    }
    else
    {
        if(thermo_flag)
        {
            forcefield_dmd->force_calc_timer(true);
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            thermo->update(time_idx,t_cur);
        }
        if(min_flag==CD_FLAG)
        {
            if(c_d_norm<cd_tol)
            {
                return true;
            }
        }
        else if(min_flag==STEP_FLAG)
        {
            if(dstep+1==step_tol)
            {
                dstep=0;
                return true;
            }
            dstep++;
        }
    }
    return false;
}
/*--------------------------------------------
 rectify
 --------------------------------------------*/
void DMD::rectify(type0* f)
{
    if(mapp->c_dof==NULL)
        return;
    
    bool* cdof=mapp->c_dof->begin();
    for(int i=0;i<atoms->natms*c_dim;i++)
        f[i]*=cdof[i];
}
/*--------------------------------------------
 clac ndofs
 --------------------------------------------*/
type0 DMD::calc_nc_dofs()
{
    type0* c=mapp->c->begin();
    type0 ndofs_lcl=0.0;
    if(mapp->c_dof==NULL)
    {
        for(int idof=0;idof<atoms->natms*c_dim;idof++)
            if(c[idof]>=0.0)
                ndofs_lcl++;
    }
    else
    {
        bool* c_dof=mapp->c_dof->begin();
        for(int idof=0;idof<atoms->natms*c_dim;idof++)
            if(c[idof]>=0.0 && c_dof[idof])
                ndofs_lcl++;
    }
    
    type0 ndofs=0.0;
    MPI_Allreduce(&ndofs_lcl,&ndofs,1,MPI_TYPE0,MPI_SUM,world);

    
    return ndofs;
}
/*--------------------------------------------
 default init
 --------------------------------------------*/
void DMD::init()
{
    dstep=-1;
    nmin=0;
    t_cur=0.0;
    nc_dofs=calc_nc_dofs();
    if(nc_dofs==0.0)
        error->abort("cannot start a dmd with no degrees of freedom");

    if(mapp->c_d==NULL)
    {
        mapp->c_d=new Vec<type0>(atoms,mapp->c->orig_dim,"c_d");
        memset(mapp->c_d->begin(),0,atoms->natms*mapp->c->orig_dim*sizeof(type0));
        if(c_dim!=mapp->c->orig_dim)
            mapp->c_d->change_dimension(0.0,mapp->c->orig_dim,mapp->c->orig_dim-c_dim);
        
        delete [] mapp->c_d->print_format;
        mapp->c_d->print_format=new char[strlen("%e ")+1];
        memcpy(mapp->c_d->print_format,"%e ",strlen("%e ")+1);
    }
    
    if(min==NULL)
    {
        old_skin=atoms->get_skin();
        atoms->set_skin(0.001);
        vecs_comm=new VecLst(atoms);
        vecs_comm->add_updt(mapp->c);
        vecs_comm->add_updt(mapp->ctype);
        atoms->init(vecs_comm,false);
    }
    else
    {
        min->init();
        min->run();
        if(min_flag==F_FLAG)
            f_norm0=min->calc_ave_f_norm();
        nmin++;
    }
    forcefield_dmd->force_calc_timer(true);
    thermo->update(fe_idx,nrgy_strss[0]);
    thermo->update(stress_idx,6,&nrgy_strss[1]);
    thermo->update(time_idx,0.0);
    thermo->init();
    if(write)
        write->init();

}
/*--------------------------------------------
 default fin
 --------------------------------------------*/
void DMD::fin()
{
    if(write)
        write->fin();
    thermo->fin();
    
    if(min==NULL)
    {
        atoms->fin();
        delete vecs_comm;
        atoms->set_skin(old_skin);
    }
    else
    {
        min->fin();
    }
    
    print_stats();
    timer->print_stats();
    neighbor->print_stats();
    forcefield_dmd->dynamic_flag=true;
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void DMD::print_stats()
{
    if(atoms->my_p==0)
    {
        fprintf(output,"dmd stats:\n");
        fprintf(output,"efficiancy fac: %e\n",t_cur/timer->tot_time);
        fprintf(output,"no. minimizations performed: %d\n",nmin);
    }
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
type0 DMD::vac_msd()
{
    
    type0* c=mapp->c->begin();
    type0* x=atoms->x->begin();
    int x_dim=atoms->x->dim;
    type0* sum_lcl;
    type0* sum;
    CREATE1D(sum_lcl,__dim__+2);
    CREATE1D(sum,__dim__+2);
    
    for(int i=0;i<__dim__+2;i++)
        sum_lcl[i]=sum[i]=0.0;
    
    type0 c_vac,rsq;
    
    for(int i=0;i<atoms->natms;i++)
    {
        c_vac=1.0;
        for(int j=0;j<c_dim;j++)
            if(c[j]>=0.0)
                c_vac-=c[j];
        
        rsq=0.0;
        for(int j=0;j<__dim__;j++)
        {
            sum_lcl[j]+=c_vac*x[j];
            rsq+=x[j]*x[j];
        }
        
        sum_lcl[__dim__]+=c_vac*rsq;
        sum_lcl[__dim__+1]+=c_vac;
        
        c+=c_dim;
        x+=x_dim;
    }
    
    MPI_Allreduce(sum_lcl,sum,__dim__+2,MPI_TYPE0,MPI_SUM,world);
    type0 ans=0.0;
    for(int i=0;i<__dim__+1;i++)
        sum[i]/=sum[__dim__+1];
        
    for(int i=0;i<__dim__;i++)
        ans-=sum[i]*sum[i];
    ans+=sum[__dim__];
    delete [] sum_lcl;
    delete [] sum;
    return ans;
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void DMD::reset()
{
    ncs=atoms->natms*c_dim;
    forcefield_dmd->dynamic_flag=false;
    neighbor_dmd->create_2nd_list();
    forcefield_dmd->init_static();
    forcefield_dmd->dc_timer();
    rectify(mapp->c_d->begin());
    c_d_norm=1.0;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMDImplicit::DMDImplicit():DMD()
{
    max_iter=20;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
DMDImplicit::~DMDImplicit()
{
}
/*--------------------------------------------
 solve the implicit equation
 --------------------------------------------*/
bool DMDImplicit::solve_non_lin()
{
    type0 res_tol=0.005*a_tol*sqrt(nc_dofs)/err_fac;
    type0 cost,cost_p;
    type0 denom=1.0*a_tol*sqrt(nc_dofs)/err_fac;
    type0 r,r_lcl,norm=1.0,del=0.0,delp=0.0;
    int iter=0,solver_iter;
    type0* c=mapp->c->begin();
    

    
    if(max_c_d*dt>0.0)
    {
        memcpy(c,y_0,ncs*sizeof(type0));
        atoms->update(mapp->c);
    }
    
    type0* c_d=mapp->c_d->begin();
    cost=cost_p=forcefield_dmd->update_J(beta,a,F)/res_tol;
    
    //printf("dt %e cost %e max_c_d %e\n",dt,cost,max_c_d);
    while(cost>=1.0 && iter<max_iter)
    {
        for(int i=0;i<ncs;i++) del_c[i]=0.0;
        
        gmres->solve(res_tol,F,del_c,solver_iter,norm);
        
        rectify(del_c);
        
        
        type0 tmp;
        r_lcl=1.0;
        for(int i=0;i<ncs;i++)
        {
            if(c[i]>=0.0)
            {
                --++del_c[i];
                tmp=c[i]+r_lcl*del_c[i];
                if(tmp>1.0)
                {
                    r_lcl=(1.0-c[i])/del_c[i];
                    while(c[i]+r_lcl*del_c[i]>1.0)
                        r_lcl=nextafter(r_lcl,0.0);
                }
                if(tmp<0.0)
                {
                    r_lcl=-c[i]/del_c[i];
                    while(c[i]+r_lcl*del_c[i]<0.0)
                        r_lcl=nextafter(r_lcl,0.0);
                }
            }
            else
                del_c[i]=0.0;
        }
        
        
        MPI_Allreduce(&r_lcl,&r,1,MPI_TYPE0,MPI_MIN,world);
        
        for(int i=0;i<ncs;i++)
        {
                c[i]+=r*del_c[i];
                --++c[i];
        }
        
        del=fabs(r*norm/denom);

        atoms->update(mapp->c);
        cost_p=forcefield_dmd->update_J(beta,a,F)/res_tol;
        cost=MIN(cost_p,del*err_fac*10.0);
        //printf("iter %d max_c_d %e cost_p %e cost %e norm %e r %e| %e\n",iter,max_c_d,cost_p,cost,norm,r,dt);
        delp=del;
        iter++;
    }
    
    
    rectify(c_d);

    if(cost<1.0)
    {
        if(iter)
            solve_acc++;
        return true;
    }
    
    
    return false;
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void DMDImplicit::intg_fail()
{
    intg_rej++;
    if(t_fin-t_cur<=2.0*dt_min)
    {
        if(q>1)
            q--;
        else
            error->abort("reached minimum order & del_t (%e) at line %d",dt,__LINE__);
    }
    else
    {
        if(dt==dt_min)
        {
            if(q>1)
                q--;
            else
                error->abort("reached minimum order & del_t (%e) at line %d",dt,__LINE__);
        }
        else
        {
            type0 r=pow(0.5/err,1.0/static_cast<type0>(q+1));
            
            if(r*dt<dt_min)
                dt=dt_min;
            else if(r*dt>t_fin-t_cur-dt_min)
                dt=t_fin-t_cur-dt_min;
            else
                dt*=r;
        }
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void DMDImplicit::nonl_fail()
{
    solve_rej++;
    if(t_fin-t_cur<=2.0*dt_min)
    {
        if(q>1)
            q--;
        else
            error->abort("reached minimum order & del_t (%e) at line %d",dt,__LINE__);
    }
    else
    {
        if(dt==dt_min)
        {
            if(q>1)
                q--;
            else
                error->abort("reached minimum order & del_t (%e) at line %d",dt,__LINE__);
        }
        else
        {
            type0 r=0.25;
            if(r*dt<dt_min)
                dt=dt_min;
            else if(r*dt>t_fin-t_cur-dt_min)
                dt=t_fin-t_cur-dt_min;
            else
                dt*=r;
        }
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void DMDImplicit::run()
{
    if(max_step==0)
        return;

    int istep;
    bool min_run=false;
    
    istep=0;
    while (istep<max_step && t_cur<t_fin)
    {
        
        restart();
        start();
        q_p=-1;
        dt_p=0.0;
        const_q=const_dt=0;
        intp_failure=0;

        
        while (istep<max_step && t_cur<t_fin && !min_run)
        {
            err=1.0;
            for(;;)
            {
                if(!interpolate())
                {
                    intp_failure=1;
                    intp_rej++;
                    interpolate_fail();
                }
                else
                    intp_failure--;
                
                err_fac_calc();
                
                if(!solve_non_lin())
                {
                    nonl_fail();
                    continue;
                }
                
                err_calc();
                
                if(err>=1.0)
                {
                    intg_fail();
                    continue;
                }
                break;
            }
            
            max_succ_q=MAX(max_succ_q,q);
            min_run=decide_min(istep,dt);
            if(min_run) continue;
            ord_dt();
            q_p=q;
            dt_p=dt;
            update_for_next();
        }

        if(!min_run) continue;
        do_min();
        min_run=false;
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
inline void DMDImplicit::ord_dt()
{
    if(dt==dt_p)
        const_dt++;
    else
        const_dt=0;

    if(q==q_p)
        const_q++;
    else
        const_q=0;
    
    type0 r=1.0;
    dq=0;

    if(intp_failure<=-q)
        ord_dt(r);

    if(10.0<r)
        r=10.0;
    else if(1.0<=r && r<1.5)
        r=1.0;
    else if(r<=0.5)
        r=0.5;
        
    if(r*dt>t_fin-t_cur)
        dt_new=t_fin-t_cur;
    else
    {
        if(t_fin-t_cur<=2.0*dt_min)
        {
            dt_new=2.0*dt_min;
        }
        else
        {
            if(r*dt<dt_min)
                dt_new=dt_min;
            else if(r*dt>=t_fin-t_cur-dt_min)
                dt_new=t_fin-t_cur-dt_min;
            else
            {
                dt_new=r*dt;
            }
        }
    }
   
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMDImplicit::reset()
{
    DMD::reset();
    a=vecs_0[0]->begin();
    y_0=vecs_0[1]->begin();
    F=vecs_0[2]->begin();
    del_c=vecs_0[3]->begin();
    gmres->refresh();
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void DMDImplicit::print_stats()
{
    DMD::print_stats();
    if(atoms->my_p==0)
    {
        fprintf(output,"final   average velocity: %e\n",c_d_norm);
        fprintf(output,"max timestep: %e\n",max_succ_dt);
        fprintf(output,"max    order: %d\n",max_succ_q);
        fprintf(output,"rejected integration   attempts: %d\n",intg_rej);
        fprintf(output,"rejected interpolation attempts: %d\n",intp_rej);
        fprintf(output,"total implicit equations: %d = accepted(%d) + rejected(%d)\n",solve_rej+solve_acc,solve_acc,solve_rej);
    }
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void DMDImplicit::init()
{
    gmres=new GMRES<type0,ForceFieldDMD>(max_iter,c_dim,*forcefield_dmd);
    DMD::init();
    max_succ_q=1;
    max_succ_dt=0.0;
    solve_rej=0;
    solve_acc=0;
    intg_rej=0;
    intp_rej=0;
    
    vecs_0=new Vec<type0>*[4];
    for(int ivec=0;ivec<4;ivec++)
        vecs_0[ivec]=new Vec<type0>(atoms,c_dim);
    allocate();
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void DMDImplicit::fin()
{
    deallocate();
    for(int ivec=0;ivec<4;ivec++)
        delete vecs_0[ivec];
    delete [] vecs_0;
    DMD::fin();
    delete gmres;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMDExplicit::DMDExplicit():DMD()
{
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
DMDExplicit::~DMDExplicit()
{
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMDExplicit::reset()
{
    DMD::reset();
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void DMDExplicit::init()
{
    DMD::init();
    max_succ_dt=0.0;
    intg_rej=0;
    intp_rej=0;
    allocate();
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void DMDExplicit::fin()
{
    deallocate();
    DMD::fin();
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void DMDExplicit::print_stats()
{
    DMD::print_stats();
    if(atoms->my_p==0)
    {
        fprintf(output,"max timestep: %e\n",max_succ_dt);
        fprintf(output,"rejected integration   attempts: %d\n",intg_rej);
        fprintf(output,"rejected interpolation attempts: %d\n",intp_rej);
    }
}
/*--------------------------------------------
 step addjustment after failure
 --------------------------------------------*/
inline void DMDExplicit::fail_stp_adj(type0 err,type0& del_t)
{
    if(t_fin-t_cur<=2.0*dt_min)
    {
        error->abort("reached minimum order & del_t (%e) at line %d",del_t,__LINE__);
    }
    else
    {
        if(del_t==dt_min)
        {
            
            error->abort("reached minimum order & del_t (%e) at line %d",del_t,__LINE__);
        }
        else
        {
            type0 r=0.9/err;
            r=MAX(r,0.5);
            
            if(r*del_t<dt_min)
                del_t=dt_min;
            else if(r*del_t>t_fin-t_cur-dt_min)
                del_t=t_fin-t_cur-dt_min;
            else
                del_t*=r;
        }
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void DMDExplicit::run()
{
    if(max_step==0)
        return;
    
    type0 del_t,del_t_tmp;
    type0 err=0.0;
    int q;
    int istep;
    bool min_run=false;
    istep=0;
    while (istep<max_step && t_cur<t_fin)
    {
        restart(del_t,q);
        
        while (istep<max_step && t_cur<t_fin && !min_run)
        {
            err=1.0;
            while(err>=1.0)
            {
                interpolate_n_err(err,del_t);
                if(err<1.0)
                    continue;
                
                fail_stp_adj(err,del_t);
            }
            
            min_run=decide_min(istep,del_t);
            if(min_run) continue;
            
            del_t_tmp=del_t;
            ord_dt(del_t,err);
            store_vecs(del_t_tmp);
            
        }
        
        if(!min_run) continue;
        do_min();
        min_run=false;
    }
}
