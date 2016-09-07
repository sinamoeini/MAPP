#ifdef Command_Style
    CommandStyle(Command_write,write)
#else
#ifndef __MAPP__command_write__
#define __MAPP__command_write__
namespace MAPP_NS
{
    class Command_write
    {
    private:
    protected:
    public:
        Command_write(int,char**);
        ~Command_write();
    };
}
#endif
#endif

