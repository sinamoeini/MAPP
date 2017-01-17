#ifdef Command_Style
    CommandStyle(Command_polarity,polarity)
#else

#ifndef __MAPP__command_polarity__
#define __MAPP__command_polarity__
#include "global.h"
namespace MAPP_NS
{
    class Command_polarity
    {
    private:
        int calc_polarity(type0*,int,int*,int,type0*);
        type0 test(type0 (&)[3][3],type0(&)[3]);
        void QR(type0(&)[3][3],type0(&)[3][3],type0(&)[3][3]);
    protected:
    public:
        Command_polarity(int,char**);
        ~Command_polarity();
        
    };
}

#endif
#endif
