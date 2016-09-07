#ifdef MD_Style
    MDStyle(MD_nh,nh)
#else
#ifndef __MAPP__md_nh__
#define __MAPP__md_nh__
#include "md.h"
namespace MAPP_NS
{
    class MD_nh : public MD
    {
    private:
        int x_d_n,f_n,type_n,id_n;
        bool dof_xst;
        type0 MTK_1,MTK_2;
        type0 ke_cur,ke_tar,t_tar,t_cur;
        type0 t_freq,tau_freq_m;
        type0 dt2,dt4,dt8;
        type0 no_dof;
        int stress_mode;
        int seed;
        bool crt_vel;

        type0 M[2][3][3];
        type0 B_ref[3][3];
        type0 vol0;
        
        type0 ke_vec[6];
        type0 ke_vec_lcl[6];
        bool H_dof[6];
        bool dof_adj[3];
        type0 virial_pe[6];
        type0 tau_freq[6];
        type0 tau_tar[6];
        type0 sigma_tar[6];
        type0 omega_m[6];
        type0 omega_d[6];

        
        int no_it_eta,no_ch_eta;
        type0* eta_m;
        type0* eta_d;
        type0* eta_dd;

        int no_it_peta,no_ch_peta;
        type0* peta_m;
        type0* peta_d;
        type0* peta_dd;
        
        void update_H(type0);
        void update_x(type0);
        void update_x_d(type0);
        void update_NH_T(type0);
        void update_NH_tau(type0);
        void update_omega_d(type0);
        void update_x_d_xpnd(type0);
                
        void create_vel(int,type0);
        void init_vel(type0);
        
        type0 dx_ave_lcl[3];
        type0 dx_ave[3];
        
        type0 omega_denom;
        
        void calc_sigma(int);
        int nreset;
        type0 tau_hydro;
        type0 vol;
        type0 drag,pdrag,tdrag;
        
        void cmd(int,char**);
        void modify_vrial();
        
        // gcmc related values
        type0 mu;
        md_type gas_type;
        int xchng_every;
        int nxchng_attmpts;
        int xchng_seed;
        int count_idx;
        class GCMC* gcmc;
    protected:
    public:
        MD_nh(int,char**);
        ~MD_nh();
        void init();
        void fin();
        void run(int);
    };
}
#endif
#endif
