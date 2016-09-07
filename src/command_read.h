#ifdef Command_Style
    CommandStyle(Command_read,read)
#else
#ifndef __MAPP__command_read__
#define __MAPP__command_read__
namespace MAPP_NS
{
    class Command_read
    {
    private:
    protected:
    public:
        Command_read(int,char**);
        ~Command_read();
    };
}
#endif
#endif

