#ifdef Command_Style
    CommandStyle(Command_temp,temp)
#else
#ifndef __MAPP__command_temp__
#define __MAPP__command_temp__
namespace MAPP_NS
{
    class Command_temp
    {
    private:
    protected:
    public:
        Command_temp(int,char**);
        ~Command_temp();
    };
}
#endif
#endif