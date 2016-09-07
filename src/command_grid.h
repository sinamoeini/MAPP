#ifdef Command_Style
    CommandStyle(Command_grid,grid)
#else
#ifndef __MAPP__command_grid__
#define __MAPP__command_grid__
namespace MAPP_NS
{
    class Command_grid
    {
    private:
    protected:
    public:
        Command_grid(int,char**);
        ~Command_grid();
    };
}
#endif
#endif

