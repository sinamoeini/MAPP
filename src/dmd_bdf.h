#ifdef DMD_Style
    DMDStyle(DMD_bdf,bdf)
#else
#ifndef __MAPP__dmd_bdf__
#define __MAPP__dmd_bdf__
#include "dmd.h"
namespace MAPP_NS
{
    class DMD_bdf :public DMDImplicit
    {
    private:
    protected:
        // stuff for book keeping
        type0* t;
        type0** y;
        type0* dy;
        Vec<type0>** vecs_1;
        void reset();
        
        // stuff for coefficients
        type0 alpha_dy;
        type0* alpha_y;
        type0 dalpha_dy;
        type0* dalpha_y;
        
        //stuff for error calculation
        type0* e_n;
        
        void interpolate(type0&,int&);
        void ord_dt(type0&,int&,type0);
        void ratio_calc(int,type0,type0&,type0&);
        
        void allocate();
        void deallocate();
        inline void fail_stp_adj(type0,type0&,int&);
        inline void init_stp_adj(type0&);
        type0 est_dt();
    public:
        DMD_bdf(MAPP *,int,char**);
        ~DMD_bdf();
        void run();
        void run0();
        void init();
        void fin();
    };
}

#endif 
#endif 



/*
 del_t_tmp=del_t;
 if(min==NULL)
 {
 if(thermo->test_prev_step()|| istep==max_step-1 || tot_t>=max_t)
 {
 forcefield_dmd->force_calc_timer(true);
 thermo->update(fe_idx,nrgy_strss[0]);
 thermo->update(stress_idx,6,&nrgy_strss[1]);
 thermo->update(time_idx,tot_t);
 }
 ord_dt(del_t,q,err);
 }
 else
 {
 type0 f_norm=min->calc_ave_f_norm();
 if(f_norm-init_f_norm>=f_tol)
 {
 forcefield_dmd->dynamic_flag=true;
 min->run();
 forcefield_dmd->dynamic_flag=false;
 init_f_norm=min->calc_ave_f_norm();
 reset();
 type0 sum=forcefield_dmd->ddc_norm_timer()/sqrt(static_cast<type0>(dof_tot));
 type0 del_t_=MIN(sqrt(2.0*a_tol/sum),1.0e-3*max_t);
 init_stp_adj(del_t_);
 }
 else
 {
 
 ord_dt(del_t,q,err);
 }
 }
*/