/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__ff__
#define __MAPP__ff__
#include "init.h"
#include "gcmc.h"
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
        virtual void xchng_energy(GCMC*)=0;
    public:
        ForceField(MAPP *);
        virtual ~ForceField();

        
        virtual void init()=0;
        virtual void fin()=0;
        virtual void init_xchng()=0;
        virtual void fin_xchng()=0;
        virtual void coef(int,char**)=0;
        type0* rsq_crd;
        type0** cut;
        type0** cut_sq;
        type0** cut_sk_sq;
        type0 max_cut();

        void force_calc_timer(bool);
        type0 energy_calc_timer();
        void xchng_energy_timer(GCMC*);
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
        virtual type0 ddc_norm()=0;
        virtual void ddc(type0*)=0;
    public:
        ForceFieldDMD(MAPP* mapp):ForceField(mapp)
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

