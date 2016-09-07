#ifdef Command_Style
    CommandStyle(Command_run,run)
#else
#ifndef __MAPP__command_run__
#define __MAPP__command_run__
namespace MAPP_NS
{
    class Command_run
    {
    private:
    protected:
    public:
        Command_run(int,char**);
        ~Command_run();
    };
}
#endif
#endif