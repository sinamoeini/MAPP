#ifdef Command_Style
    CommandStyle(Command_ff_coef,ff_coef)
#else
#ifndef __MAPP__command_ff_coef__
#define __MAPP__command_ff_coef__
namespace MAPP_NS
{
    class Command_ff_coef
    {
    private:
    protected:
    public:
        Command_ff_coef(int,char**);
        ~Command_ff_coef();
    };
}
#endif
#endif

