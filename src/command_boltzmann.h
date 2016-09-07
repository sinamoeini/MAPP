#ifdef Command_Style
    CommandStyle(Command_boltzmann,boltzmann)
#else
#ifndef __MAPP__command_boltzmann__
#define __MAPP__command_boltzmann__
namespace MAPP_NS
{
    class Command_boltzmann
    {
    private:
    protected:
    public:
        Command_boltzmann(int,char**);
        ~Command_boltzmann();
    };
}
#endif
#endif

