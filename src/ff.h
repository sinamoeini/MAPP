/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__ff__
#define __MAPP__ff__
#include "init.h"
namespace MAPP_NS {
    class ForceField : protected InitPtrs{
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
        ForceField(MAPP *);
        virtual ~ForceField();

        
        virtual void init()=0;
        virtual void fin()=0;
        virtual void coef(int,char**)=0;
        type0* rsq_crd;
        type0* cut;
        type0* cut_sq;
        type0* cut_sk_sq;
        type0 max_cut();

        void force_calc_timer(bool);
        type0 energy_calc_timer();
        
        Vec<type0>* f;
        type0* nrgy_strss;
    };
}

namespace MAPP_NS
{
    class ForceFieldMD : public ForceField
    {
    private:
    protected:
    public:
        ForceFieldMD(MAPP* mapp):ForceField(mapp){}
        virtual ~ForceFieldMD(){}
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
        virtual type0 dc_en_proj(bool,type0*,type0&)=0;
        virtual type0 ddc_norm()=0;
        virtual type0 imp_cost_grad(bool,type0,type0,type0,type0*,type0*)=0;
        virtual type0 dc_norm_grad(bool,type0*,type0*)=0;
        virtual type0 en_grad(bool,type0*,type0*)=0;
    public:
        ForceFieldDMD(MAPP* mapp):ForceField(mapp)
        {dynamic_flag=true;}
        virtual ~ForceFieldDMD(){}
        virtual void init()=0;
        virtual void fin()=0;
        virtual void coef(int,char**)=0;
        
        type0 imp_cost_grad_timer(bool,type0,type0,type0,type0*,type0*);
        type0 dc_norm_grad_timer(bool,type0*,type0*);
        type0 en_grad_timer(bool,type0*,type0*);
        void dc_timer();
        type0 dc_en_proj_timer(bool,type0*,type0&);
        type0 ddc_norm_timer();
    
        
        type0 alpha_min,alpha_max;
        bool dynamic_flag;
    };
}
#endif

