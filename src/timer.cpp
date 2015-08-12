#include "timer.h"
#include "atoms.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor of timer
 --------------------------------------------*/
Timer::Timer(MAPP* mapp):InitPtrs(mapp)
{
    no_of_modes=5;
    curr_level=-1;
    level2mode_capacity=0;
    
    CREATE1D(time,no_of_modes);
    
    tot_time=0.0;
    
    for(int i=0;i<no_of_modes;i++)
        time[i]=0.0;
}
/*--------------------------------------------
 destructor of timer
 --------------------------------------------*/
Timer::~Timer()
{
    if(no_of_modes)
        delete [] time;
    
    if (level2mode_capacity)
        delete [] level2mode;
}
/*--------------------------------------------
 start communication time;
 --------------------------------------------*/
void Timer::start(int mode)
{
    type0 t=MPI_Wtime();
    if(curr_level!=-1)
        time[level2mode[curr_level]]+=t;
    
    curr_level++;
    if(curr_level==level2mode_capacity)
    {
        GROW(level2mode,level2mode_capacity,level2mode_capacity+1);
        level2mode_capacity++;
    }
    level2mode[curr_level]=mode;
    
    time[mode]-=t;
    
}
/*--------------------------------------------
 start communication time;
 --------------------------------------------*/
void Timer::stop(int mode)
{
    if(mode!=level2mode[curr_level])
        error->abort("timer error");
        
        
    type0 t=MPI_Wtime();
    time[mode]+=t;
    curr_level--;
    if(curr_level!=-1)
        time[level2mode[curr_level]]-=t;
}
/*--------------------------------------------
 start communication time;
 --------------------------------------------*/
type0 Timer::mode_time(int mode)
{
    return time[mode];
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Timer::init()
{
    type0 t=MPI_Wtime();
    tot_time=-t;
    for(int i=0;i<no_of_modes;i++)
        time[i]=0.0;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Timer::fin()
{
    type0 t=MPI_Wtime();
    tot_time+=t;

    t=tot_time;
    for(int i=0;i<no_of_modes;i++)
        t-=time[i];
    
    time[OTHER_TIME_mode]+=t;
}

/*--------------------------------------------
 print
 --------------------------------------------*/
void Timer::print_time_stats()
{
    if(atoms->my_p_no==0)
    {
        fprintf(output,"total time: %lf secs\n"
        ,tot_time);
        
        fprintf(output,"force time: %lf secs (%05.2lf%%)\n"
        ,time[FORCE_TIME_mode],time[FORCE_TIME_mode]*100.0/tot_time);
        
        fprintf(output,"neigh time: %lf secs (%05.2lf%%)\n"
        ,time[NEIGH_TIME_mode],time[NEIGH_TIME_mode]*100.0/tot_time);
        
        fprintf(output,"comm  time: %lf secs (%05.2lf%%)\n"
        ,time[COMM_TIME_mode],time[COMM_TIME_mode]*100.0/tot_time);
        
        fprintf(output,"write time: %lf secs (%05.2lf%%)\n"
        ,time[WRITE_TIME_mode],time[WRITE_TIME_mode]*100.0/tot_time);
        
        fprintf(output,"other time: %lf secs (%05.2lf%%)\n"
        ,time[OTHER_TIME_mode],time[OTHER_TIME_mode]*100.0/tot_time);
    }

    
}

/*--------------------------------------------
 print
 --------------------------------------------*/
type0 Timer::tst_start()
{
    tst_time=-MPI_Wtime();
    return 0.0;
}
/*--------------------------------------------
 print
 --------------------------------------------*/
type0 Timer::tst_stop()
{
    tst_time+=MPI_Wtime();
    return tst_time;
}
