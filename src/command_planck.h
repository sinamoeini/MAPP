#ifdef Command_Style
    CommandStyle(Command_planck,planck)
#else
#ifndef __MAPP__command_planck__
#define __MAPP__command_planck__
namespace MAPP_NS
{
    class Command_planck
    {
    private:
    protected:
    public:
        Command_planck(int,char**);
        ~Command_planck();
    };
}
#endif
#endif

