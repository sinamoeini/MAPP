#ifdef Command_Style
    CommandStyle(Command_fix,fix)
#else
#ifndef __MAPP__command_fix__
#define __MAPP__command_fix__
namespace MAPP_NS
{
    class Command_fix
    {
    private:
    protected:
    public:
        Command_fix(int,char**);
        ~Command_fix();
    };
}
#endif 
#endif
