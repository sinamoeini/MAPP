#ifdef Command_Style
    CommandStyle(Command_displace,displace)
#else
#ifndef __MAPP__command_displace__
#define __MAPP__command_displace__
namespace MAPP_NS
{
    class Command_displace
    {
    private:
    protected:
    public:
        Command_displace(int,char**);
        ~Command_displace();
    };
}
#endif 
#endif
