#ifdef Command_Style
    CommandStyle(Command_del_write,del_write)
#else
#ifndef __MAPP__command_del_write__
#define __MAPP__command_del_write__
namespace MAPP_NS
{
    class Command_del_write
    {
    private:
    protected:
    public:
        Command_del_write(int,char**);
        ~Command_del_write();
    };
}
#endif
#endif

