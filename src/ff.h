/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__ff__
#define __MAPP__ff__
#include "type_def.h"
namespace MAPP_NS
{
    template<typename T>
    class Vec;
    class ForceField {
    private:
        int cut_sz;
        int ns_alloc;

        
    protected:
        type0* nrgy_strss_lcl;
        void cut_off_alloc();
        void cut_off_dealloc();
        virtual void force_calc(bool)=0;
        virtual type0 energy_calc()=0;
    public:
        ForceField();
        virtual ~ForceField();

        
        virtual void init()=0;
        virtual void fin()=0;
        virtual void coef(int,char**)=0;
        type0* rsq_crd;
        type0** cut;
        type0** cut_sq;
        type0** cut_sk_sq;
        type0 max_cut();

        void force_calc_timer(bool);
        type0 energy_calc_timer();
        Vec<type0>* f;
        type0* nrgy_strss;
        
        int gcmc_n_vars;
        int gcmc_n_cutoff;
        bool gcmc_tag_enabled;
        
        
    };
}

namespace MAPP_NS
{
    class ForceFieldMD : public ForceField
    {
    private:
    protected:
        virtual void pre_xchng_energy(class GCMC*)=0;
        virtual type0 xchng_energy(class GCMC*)=0;
        virtual void post_xchng_energy(class GCMC*)=0;
    public:
        ForceFieldMD():ForceField(){}
        virtual ~ForceFieldMD(){}
        virtual void init_xchng()=0;
        virtual void fin_xchng()=0;
        
        void pre_xchng_energy_timer(class GCMC*);
        type0 xchng_energy_timer(class GCMC*);
        void post_xchng_energy_timer(class GCMC*);
    };
}

namespace MAPP_NS
{
    class ForceFieldDMD : public ForceField
    {
    private:
    protected:
        virtual void force_calc(bool)=0;
        virtual type0 energy_calc()=0;
        virtual void dc()=0;
        virtual type0 ddc_norm()=0;
        virtual void ddc(type0*)=0;
    public:
        ForceFieldDMD():ForceField()
        {dynamic_flag=true;}
        virtual ~ForceFieldDMD(){}
        virtual void init()=0;
        virtual void fin()=0;
        virtual void coef(int,char**)=0;
        virtual void set_temp(type0)=0;
        
        void dc_timer();
        type0 ddc_norm_timer();
        void ddc_timer(type0*);
        
        type0 alpha_min,alpha_max;
        bool dynamic_flag;
        
        
        
        
        virtual void init_static()=0;
        virtual type0 update_J(type0,type0*,type0*)=0;
        virtual void operator()(Vec<type0>*,Vec<type0>*)=0;
    };
}
#endif

