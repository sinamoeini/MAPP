#ifdef LS_Style
    LSStyle(LineSearch_backtrack,backtrack)
#else
#ifndef __MAPP__ls_bt__
#define __MAPP__ls_bt__

#include "ls.h"
namespace MAPP_NS
{
    class LineSearch_backtrack : public LineSearch
    {
    private:
    protected:
        type0 c,rho,min_alpha;
    public:
        LineSearch_backtrack(MAPP*,int,char**);
        ~LineSearch_backtrack();
        int line_min(type0&,type0&,int);
    };
}
#endif
#endif