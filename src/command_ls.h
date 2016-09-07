#ifdef Command_Style
    CommandStyle(Command_ls,ls)
#else
#ifndef __MAPP__command_ls__
#define __MAPP__command_ls__
namespace MAPP_NS
{
    class Command_ls
    {
    private:
    protected:
    public:
        Command_ls(int,char**);
        ~Command_ls();
    };
}
#endif
#endif

