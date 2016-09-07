#ifdef Command_Style
    CommandStyle(Command_group,group)
#else
#ifndef __MAPP__command_group__
#define __MAPP__command_group__
namespace MAPP_NS
{
    class Command_group
    {
    private:
    protected:
    public:
        Command_group(int,char**);
        ~Command_group();
    };
}
#endif 
#endif
