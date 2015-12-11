#include "dmd.h"
#include "thermo_dynamics.h"
#include "atom_types.h"
#include "write.h"
#include "error.h"
#include "memory.h"
#include "timer.h"
#include "min_styles.h"
#include <limits>
//#define DMD_DEBUG
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMD::DMD(MAPP* mapp):InitPtrs(mapp),
nrgy_strss(forcefield->nrgy_strss)
{
    if(forcefield==NULL)
        error->abort("ff should be "
        "initiated before dmd");
    
    if(mapp->mode!=DMD_mode)
        error->abort("dmd works only "
        "for md mode");
    
    char** args;
    int nargs=mapp->parse_line(
    "Time FE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    
    fe_idx=1;
    stress_idx=2;
    time_idx=0;
    thermo=new ThermoDynamics(mapp,nargs,args);
    for(int i=0;i<nargs;i++)
        delete [] args[i];
    if(nargs)
        delete [] args;

    
    
    c_dim=mapp->c->dim;
    forcefield_dmd=dynamic_cast<ForceFieldDMD*>(forcefield);
    neighbor_dmd=dynamic_cast<Neighbor_dmd*>(neighbor);
    
    a_tol=sqrt(2.0*numeric_limits<type0>::epsilon());
    min_del_t=std::numeric_limits<type0>::epsilon();
    max_t=1.0e7;
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
    #define Min_Style
    #define MinStyle(class_name,style_name)   \
    else if(strcmp(args[1],#style_name)==0)   \
    {if(min!=NULL)delete min;                 \
    min= new class_name(mapp,nargs,args);}

    if(0){}
    #include "min_styles.h"
    else
        error->abort("wrong style of minimization"
        ": %s",args[1]);
    #undef Min_Style
    #undef MinStyle
    min->output_flag=false;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMD::run(type0 t)
{
    if(t<=0.0)
        error->abort("max_t in dmd should be greater than 0.0");
    if(min_del_t>t)
        error->abort("max_t in dmd should be greater than min_del_t");
    max_t=t;
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
    tot_t+=del_t;
    if(write!=NULL)
        write->write();
    thermo->thermo_print();
    
    bool thermo_flag=(thermo->test_prev_step()
    || istep==max_step-1
    || tot_t>=max_t);
    istep++;
    step_no++;
    if(min_flag==F_FLAG)
    {
        type0 f_norm=min->calc_ave_f_norm();
        if(thermo_flag)
        {
            thermo->update(fe_idx,nrgy_strss[0]);
            thermo->update(stress_idx,6,&nrgy_strss[1]);
            thermo->update(time_idx,tot_t);
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
            thermo->update(time_idx,tot_t);
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
    tot_t=0.0;
    nc_dofs=calc_nc_dofs();
    if(nc_dofs==0.0)
        error->abort("cannot start a dmd with no degrees of freedom");

    if(mapp->c_d==NULL)
    {
        mapp->c_d=new Vec<type0>(atoms,c_dim,"c_d");
        for(int i=0;i<atoms->natms*c_d_norm;i++)
            mapp->c_d->begin()[i]=0.0;
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
    if(write!=NULL)
        write->init();
}
/*--------------------------------------------
 default fin
 --------------------------------------------*/
void DMD::fin()
{
    if(write!=NULL)
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
        fprintf(output,"efficiancy fac: %e\n",tot_t/timer->tot_time);
        fprintf(output,"no. minimizations performed: %d\n",nmin);
    }
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void DMD::reset()
{
    ncs=atoms->natms*c_dim;
    forcefield_dmd->dynamic_flag=false;
    neighbor_dmd->create_2nd_list();
    forcefield_dmd->dc_timer();
    rectify(mapp->c_d->begin());
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMDImplicit::DMDImplicit(MAPP* mapp):DMD(mapp)
{
    iter_dcr_thrsh=3;
    max_iter=20;
    m_tol=sqrt(numeric_limits<type0>::epsilon());
    
    
    char** args;
    int nargs=mapp->parse_line("ls brent max_iter 10",args);
    
    ls_dmd=new LineSearch_brent<DMDImplicit>(mapp,nargs,args);
    
    for(int i=0;i<nargs;i++)
        delete [] args[i];
    if(nargs)
        delete [] args;
    
    ls_dmd->init(this);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
DMDImplicit::~DMDImplicit()
{
    delete ls_dmd;
}
/*--------------------------------------------
 solve the implicit equation
 --------------------------------------------*/
int DMDImplicit::solve_n_err(type0& cost,type0& err)
{
    type0* c=mapp->c->begin();
    type0* c_d=mapp->c_d->begin();
    
    type0 gamma;
    type0 inner[2],inner_lcl[2],tmp0;
    type0 err_lcl;
    type0 ratio;
    type0 g0_g0=-1.0,g_g,g_g0;
    type0 curr_cost;
    type0 tol=m_tol;
    type0 mm_tol=tol/sqrt(nc_dofs);
    int ls_succ=LS_S,iter=0;

    
    memcpy(c,y_0,ncs*sizeof(type0));
    atoms->update(mapp->c);
    curr_cost=forcefield_dmd->imp_cost_grad_timer(true,mm_tol,beta,a,g);
    rectify(g);
    inner_lcl[0]=0.0;
    for(int i=0;i<ncs;i++) inner_lcl[0]+=g[i]*g[i];
    MPI_Allreduce(&inner_lcl,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
    
    memcpy(h,g,ncs*sizeof(type0));
    h_h=g_h=g0_g0;
    
    while(g0_g0!=0.0 && iter<max_iter && ls_succ==LS_S)
    {
        memcpy(g0,g,ncs*sizeof(type0));
        memcpy(c0,c,ncs*sizeof(type0));
        
        
        ls_succ=ls_dmd->line_min(curr_cost,gamma,1);
        iter++;
        
        if(ls_succ!=LS_S)
            continue;
        curr_cost=forcefield_dmd->imp_cost_grad_timer(true,mm_tol,beta,a,g);
        rectify(g);
        
        inner_lcl[0]=inner_lcl[1]=0.0;
        for(int i=0;i<ncs;i++)
        {
            inner_lcl[0]+=g[i]*g0[i];
            inner_lcl[1]+=g[i]*g[i];
        }
        
        MPI_Allreduce(inner_lcl,inner,2,MPI_TYPE0,MPI_SUM,world);
        g_g0=inner[0];
        g_g=inner[1];
        ratio=(g_g-g_g0)/g0_g0;
        
        if(g_g==0.0)
        {
            g0_g0=0.0;
            continue;
        }
        
        int bound_chk_lcl=1,bound_chk;
        for(int i=0;i<ncs;i++)
        {
            h[i]*=ratio;
            h[i]+=g[i];
            if((c[i]==0.0 && h[i]<0.0) || (c[i]==1.0 && h[i]>0.0))
                bound_chk_lcl=0;
            inner_lcl[0]+=h[i]*g[i];
            inner_lcl[1]+=h[i]*h[i];
        }
        
        MPI_Allreduce(&bound_chk_lcl,&bound_chk,1,MPI_INT,MPI_MIN,world);
        if(bound_chk)
        {
            MPI_Allreduce(inner_lcl,inner,2,MPI_TYPE0,MPI_SUM,world);
            g_h=inner[0];
            h_h=inner[1];
        }
        
        if(g_h<0.0  || bound_chk==0)
        {
            memcpy(h,g,ncs*sizeof(type0));
            h_h=g_h=g_g;
        }
        
        g0_g0=g_g;
    }
    
    rectify(c_d);
    
    err_lcl=0.0;
    type0 c_d_norm_lcl=0.0;
    type0 c_d_max_lcl=0.0;
    type0 c_d_max;
    type0 err_max_lcl=0.0;
    type0 err_max;
    type0 cost_max_lcl=0.0;
    type0 cost_max;
    
    for(int i=0;i<ncs;i++)
    {
        if(c[i]>=0.0)
        {
            tmp0=(y_0[i]-c[i]);
            err_lcl+=tmp0*tmp0;
            c_d_norm_lcl+=c_d[i]*c_d[i];
            c_d_max_lcl=MAX(c_d_max_lcl,fabs(c_d[i]));
            err_max_lcl=MAX(err_max_lcl,fabs(tmp0));
            cost_max_lcl=MAX(cost_max_lcl,fabs(c[i]+a[i]-beta*c_d[i]));
        }
    }
    
    err=0.0;
    MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
    MPI_Allreduce(&c_d_max_lcl,&c_d_max,1,MPI_TYPE0,MPI_MAX,world);
    MPI_Allreduce(&err_max_lcl,&err_max,1,MPI_TYPE0,MPI_MAX,world);
    MPI_Allreduce(&cost_max_lcl,&cost_max,1,MPI_TYPE0,MPI_MAX,world);
    MPI_Allreduce(&c_d_norm_lcl,&c_d_norm,1,MPI_TYPE0,MPI_SUM,world);
    err=sqrt(err/nc_dofs)/a_tol;
    c_d_norm=sqrt(c_d_norm/nc_dofs);
    err*=err_prefac;

    cost=cost_max/m_tol;
    
    if(iter)
    {
        if(cost<1.0)
            solve_acc++;
        else
            solve_rej++;
    }
    if(err>=1.0)
        intg_rej++;
    
    
#ifdef DMD_DEBUG
    if(atoms->my_p==0)
        printf("%lf %lf %d %e %d %lf %e | Max %e %e\n",err,cost,iter,c_d_norm,ls_succ,g0_g0,max_a,c_d_max,err_prefac*err_max/a_tol);
#endif
    

    return iter;
}
/*--------------------------------------------
 find the the cost function given gamma
 --------------------------------------------*/
void DMDImplicit::ls_prep(type0& dfa,type0& h_norm,type0& max_a_)
{
    type0 max_a_lcl;
    max_a_lcl=numeric_limits<type0>::infinity();

    for(int i=0;i<ncs;i++)
    {
        if(h[i]>0.0)
        {
            c1[i]=(1.0-c0[i])/h[i];
            max_a_lcl=MIN(c1[i],max_a_lcl);
        }
        else if(h[i]<0.0)
        {
            c1[i]=(0.0-c0[i])/h[i];
            max_a_lcl=MIN(c1[i],max_a_lcl);
        }

    }
    MPI_Allreduce(&max_a_lcl,&max_a_,1,MPI_TYPE0,MPI_MIN,world);
    h_norm=sqrt(h_h);
    dfa=-g_h;
    if(max_a_==numeric_limits<type0>::infinity())
        error->abort("h_norm is 0.0");
    max_a=max_a_;
    
    for(int i=0;i<ncs;i++)
    {
        if(c0[i]>=0.0)
        {
            if(h[i]>0.0)
            {
                if(c1[i]>max_a)
                    c1[i]=c0[i]+max_a*h[i];
                else
                    c1[i]=1.0;
            }
            else if(h[i]<0.0)
            {
                if(c1[i]>max_a)
                    c1[i]=c0[i]+max_a*h[i];
                else
                    c1[i]=0.0;
            }
            else
                c1[i]=c0[i];
        }
    }
}
/*--------------------------------------------
 find the the cost function given gamma
 --------------------------------------------*/
type0 DMDImplicit::F(type0 alpha)
{
    type0* c=mapp->c->begin();
    
    if(alpha<max_a)
    {
        for(int i=0;i<ncs;i++)
            if(c0[i]>=0.0)
                c[i]=c0[i]+alpha*h[i];
    }
    else
    {
        for(int i=0;i<ncs;i++)
            if(c0[i]>=0.0)
                c[i]=c1[i];
    }
    
    atoms->update(mapp->c);
    
    return forcefield_dmd->imp_cost_grad_timer(false,0.0,beta,a,g);
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 DMDImplicit::dF(type0 alpha,type0& drev)
{
    type0* c=mapp->c->begin();
    
    if(alpha<max_a)
    {
        for(int i=0;i<ncs;i++)
            if(c0[i]>=0.0)
                c[i]=c0[i]+alpha*h[i];
    }
    else
    {
        for(int i=0;i<ncs;i++)
            if(c0[i]>=0.0)
                c[i]=c1[i];
    }
    
    atoms->update(mapp->c);
    
    type0 cost=forcefield_dmd->imp_cost_grad_timer(true,0.0,beta,a,g);
    type0 inner0=0.0;
    for(int i=0;i<ncs;i++)
        if(c[i]>=0.0)
            inner0+=g[i]*h[i];
    
    MPI_Allreduce(&inner0,&drev,1,MPI_TYPE0,MPI_SUM,world);
    drev*=-1.0;
    return cost;
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
void DMDImplicit::F_reset()
{
    type0* c=mapp->c->begin();
    memcpy(c,c0,ncs*sizeof(type0));
    atoms->update(mapp->c);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void DMDImplicit::reset()
{
    DMD::reset();
    a=vecs_0[0]->begin();
    y_0=vecs_0[1]->begin();
    g=vecs_0[2]->begin();
    h=vecs_0[3]->begin();
    g0=vecs_0[4]->begin();
    c0=vecs_0[5]->begin();
    c1=vecs_0[6]->begin();
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
    DMD::init();
    max_succ_q=1;
    max_succ_dt=0.0;
    solve_rej=0;
    solve_acc=0;
    intg_rej=0;
    intp_rej=0;
    
    vecs_0=new Vec<type0>*[7];
    for(int ivec=0;ivec<7;ivec++)
        vecs_0[ivec]=new Vec<type0>(atoms,c_dim);
    allocate();
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void DMDImplicit::fin()
{
    deallocate();
    for(int ivec=0;ivec<7;ivec++)
        delete vecs_0[ivec];
    delete [] vecs_0;
    DMD::fin();
}
/*--------------------------------------------
 step addjustment after failure
 --------------------------------------------*/
inline void DMDImplicit::fail_stp_adj(type0 err,type0 m_err,type0& del_t,int& q)
{
    const_stps=0;
    
    if(max_t-tot_t<=2.0*min_del_t)
    {
        if(q>1)
            q--;
        else
            error->abort("reached minimum order & del_t (%e)",del_t);
    }
    else
    {
        if(del_t==min_del_t)
        {
            if(q>1)
                q--;
            else
                error->abort("reached minimum order & del_t (%e)",del_t);
        }
        else
        {
            type0 r=pow(0.5/MAX(err,m_err),1.0/static_cast<type0>(q+1));
            
            r=MIN(r,0.9);
            
            if(r*del_t<min_del_t)
                del_t=min_del_t;
            else if(r*del_t>max_t-tot_t-min_del_t)
                del_t=max_t-tot_t-min_del_t;
            else
                del_t*=r;
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

    type0 del_t,del_t_tmp;
    type0 cost,err;
    int q;
    int istep;
    bool min_run=false;
    
    istep=0;
    while (istep<max_step && tot_t<max_t)
    {
        
        restart(del_t,q);
        const_stps=0;
        int iter=-1,iter_0=-1;

        while (istep<max_step && tot_t<max_t && !min_run)
        {
            err=1.0;
            cost=1.0;
            while (MAX(cost,err)>=1.0)
            {
                interpolate(del_t,q);
#ifdef DMD_DEBUG
                if(atoms->my_p==0) printf("%e %d: ",del_t,q);
#endif
                iter=solve_n_err(cost,err);
                if(MAX(cost,err)<1.0)
                    continue;
                if(cost>=1.0)
                    iter_0=-1;
                
                fail_stp_adj(err,cost,del_t,q);
            }
            
            if(iter_0!=-1 && (iter<iter_0 || iter==0))
                iter_dcr_cntr++;
            else
                iter_dcr_cntr=0;
            
            iter_0=iter;
            
            max_succ_q=MAX(max_succ_q,q);
            
            min_run=decide_min(istep,del_t);
            if(min_run) continue;

            del_t_tmp=del_t;
            ord_dt(err,del_t,q);
            store_vecs(del_t_tmp);

        }

        if(!min_run) continue;
        do_min();
        min_run=false;
    }
}
/*--------------------------------------------
 run
 --------------------------------------------*/
inline void DMDImplicit::ord_dt(type0 err,type0& del_t,int& q)
{
    type0 r=1.0;
    int del_q=0;
   
    if(iter_dcr_cntr>iter_dcr_thrsh && const_stps>q+1)
        ord_dt(err,del_t,q,r,del_q);

    
    if(r>=2.0)
        r=2.0;
    else if(r<=0.5)
        r=0.5;
    else
        r=1.0;
    
    if(r>1.0)
        iter_dcr_cntr=0;
    
    bool const_stp_chk=false;
    if(r*del_t>max_t-tot_t)
        del_t=max_t-tot_t;
    else
    {
        if(max_t-tot_t<=2.0*min_del_t)
        {
            del_t=2.0*min_del_t;
        }
        else
        {
            if(r*del_t<min_del_t)
                del_t=min_del_t;
            else if(r*del_t>=max_t-tot_t-min_del_t)
                del_t=max_t-tot_t-min_del_t;
            else
            {
                del_t*=r;
                if(r==1.0)
                    const_stp_chk=true;
            }
        }
    }
    
    if(const_stp_chk && del_q==0)
        const_stps++;
    else
        const_stps=0;
    
    q+=del_q;
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
int DMDImplicit::test()
{
    type0 dfa,h_norm;
    ls_prep(dfa,h_norm,max_a);
    dfa=-g_h;
    int no=100;
    type0 frac=1.0e-2*max_a;
    type0 fu,u=0.0;
    type0 fa=F(0.0);
    
    if(atoms->my_p==0) printf("--------------------- %d\n",step_no);
    if(atoms->my_p==0) printf("u fu f_x u*dfa\n");
    
    for(int i=0;i<no;i++)
    {
        fu=F(u);
        if(atoms->my_p==0)printf("%22.20lf %22.20lf %22.20lf \n",u,fu-fa,u*dfa);
        u+=frac;
    }
    F_reset();
    if(atoms->my_p==0) printf("---------------------n");
    error->abort("");
    return 0;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
DMDExplicit::DMDExplicit(MAPP* mapp):
DMD(mapp)
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
}
/*--------------------------------------------
 given the direction h do the line search
 --------------------------------------------*/
void DMDExplicit::fin()
{
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
