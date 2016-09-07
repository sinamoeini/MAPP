#ifndef __MAPP__thermo_dynamics__
#define __MAPP__thermo_dynamics__
#include "type_def.h"
namespace MAPP_NS
{
    
    class ThermoQuantity
    {
    private:
    protected:
    public:
        void init(const char*);
        ThermoQuantity();
        ThermoQuantity(ThermoQuantity&);
        ThermoQuantity(ThermoQuantity&&);
        ~ThermoQuantity();
        ThermoQuantity& operator =(ThermoQuantity&);
        ThermoQuantity& operator =(ThermoQuantity&&);
        void mod(int);
        char* hdr_name;
        int hdr_name_lngth;
        char* mod_hdr_name;
        int mod_hdr_name_lngth;
        type0 value;
    };
    
    
    class ThermoDynamics 
    {
    private:
    protected:
        char* step_name;
        int step_name_lngth;
        int init_step;
        void hdr_print();
        void tail_print();
        void val_print();
        
        int mod_lngth;
        ThermoQuantity* quantities;
        int no_quantities;
        char* qform;
        char* sform;
    public:
        ThermoDynamics(int,char**);
        ~ThermoDynamics();
        int add(const char*);
        void update(int,type0);
        void update(int,int,type0*);
        
        
        void init();
        void thermo_print();
        void fin();
        
        int test_prev_step();
        int print_step;
    };
}

#endif 
