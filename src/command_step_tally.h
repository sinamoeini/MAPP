#ifdef Command_Style
    CommandStyle(Command_step_tally,step_tally)
#else
#ifndef __MAPP__command_step_tally__
#define __MAPP__command_step_tally__
namespace MAPP_NS
{
    class Command_step_tally
    {
    private:
    protected:
    public:
        Command_step_tally(int,char**);
        ~Command_step_tally();
    };
}
#endif
#endif

