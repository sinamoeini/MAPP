#ifdef Command_Style
    CommandStyle(Command_mode,mode)
#else
#ifndef __MAPP__command_mode__
#define __MAPP__command_mode__
namespace MAPP_NS
{
    class Command_mode
    {
    private:
    protected:
    public:
        Command_mode(int,char**);
        ~Command_mode();
    };
}
#endif
#endif
