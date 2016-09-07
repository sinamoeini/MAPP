#ifdef Command_Style
    CommandStyle(Command_precision,precision)
#else
#ifndef __MAPP__command_precision__
#define __MAPP__command_precision__
namespace MAPP_NS
{
    class Command_precision
    {
    private:
    protected:
    public:
        Command_precision(int,char**);
        ~Command_precision();
    };
}
#endif
#endif
