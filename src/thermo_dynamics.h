#ifndef __MAPP__thermo_dynamics__
#define __MAPP__thermo_dynamics__

#include "init.h"
#include "atoms.h"


namespace MAPP_NS {
    
    class ThermoQuantity
    {
    private:
    protected:
    public:
        void init(const char*);
        ~ThermoQuantity();
        void mod(int);
        char* hdr_name;
        int hdr_name_lngth;
        char* mod_hdr_name;
        int mod_hdr_name_lngth;
        type0 value;
    };
    
    
    class ThermoDynamics : protected InitPtrs
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
        /*
        type0 run_time;
        type0 force_time;
        type0 comm_time;
         */
        char* qform0;
        char* qform1;
        char* sform;
    public:
        ThermoDynamics(MAPP*,int,char**);
        ~ThermoDynamics(); 
        /*
        void start_force_time();
        void stop_force_time();
        void start_comm_time();
        void stop_comm_time();
        */
        void update(int,type0);
        void update(int,int,type0*);
        
        
        void init();
        void thermo_print();
        void fin();
        
        int test_prev_step();
        int print_step;
    };

    
}

#endif /* defined(__MAPP__thermo_dynamics__) */
