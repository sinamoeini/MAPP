#ifdef Command_Style
    CommandStyle(Command_reset,reset)
#else
#ifndef __MAPP__command_reset__
#define __MAPP__command_reset__
namespace MAPP_NS
{
    class Command_reset
    {
    private:
    protected:
    public:
        Command_reset(int,char**);
        ~Command_reset();
    };
}
#endif
#endif
