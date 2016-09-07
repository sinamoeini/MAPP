#ifdef Command_Style
    CommandStyle(Command_change_box,change_box)
#else
#ifndef __MAPP__command_change_box__
#define __MAPP__command_change_box__
namespace MAPP_NS
{
    class Command_change_box
    {
    private:
    protected:
    public:
        Command_change_box(int,char**);
        ~Command_change_box();
    };
}
#endif
#endif
