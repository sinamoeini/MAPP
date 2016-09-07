#ifdef Command_Style
    CommandStyle(Command_dmd,dmd)
#else
#ifndef __MAPP__command_dmd__
#define __MAPP__command_dmd__
namespace MAPP_NS
{
    class Command_dmd
    {
    private:
    protected:
    public:
        Command_dmd(int,char**);
        ~Command_dmd();
    };
}
#endif
#endif

