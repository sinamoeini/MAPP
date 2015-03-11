#ifndef __MAPP__timer__
#define __MAPP__timer__
#include "init.h"
namespace MAPP_NS {
    enum{COMM_TIME_mode,NEIGH_TIME_mode,FORCE_TIME_mode,WRITE_TIME_mode,OTHER_TIME_mode};
    
    class Timer:protected InitPtrs
    {
    private:
        int curr_time_mode;
        int prev_time_mode;
    protected:
    public:
        int* level2mode;
        int level2mode_capacity;
        int curr_level;
        
        int no_of_modes;
        
        type0 tot_time;
        type0* time;
        
        void start(int);
        void stop(int);
        type0 mode_time(int);
        void init();
        void fin();
        void print_time_stats();
        
        Timer(MAPP*);
        ~Timer();

    };
}
#endif
