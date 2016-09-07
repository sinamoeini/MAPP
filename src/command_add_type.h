#ifdef Command_Style
    CommandStyle(Command_add_type,add_type)
#else
#ifndef __MAPP__command_add_type__
#define __MAPP__command_add_type__
namespace MAPP_NS
{
    class Command_add_type
    {
    private:
    protected:
    public:
        Command_add_type(int,char**);
        ~Command_add_type();
    };
}
#endif
#endif

