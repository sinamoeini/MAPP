#ifdef Command_Style
    CommandStyle(Command_skin,skin)
#else
#ifndef __MAPP__command_skin__
#define __MAPP__command_skin__
namespace MAPP_NS
{
    class Command_skin
    {
    private:
    protected:
    public:
        Command_skin(int,char**);
        ~Command_skin();
    };
}
#endif
#endif

