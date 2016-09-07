#ifdef Command_Style
    CommandStyle(Command_time_step,time_step)
#else
#ifndef __MAPP__command_time_step__
#define __MAPP__command_time_step__
namespace MAPP_NS
{
    class Command_time_step
    {
    private:
    protected:
    public:
        Command_time_step(int,char**);
        ~Command_time_step();
    };
}
#endif
#endif

