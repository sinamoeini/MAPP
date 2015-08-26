#ifdef LS_Style
    LSStyle(LineSearch_brent,brent)
#else
#ifndef __MAPP__ls_brent__
#define __MAPP__ls_brent__

#include "ls.h"
namespace MAPP_NS
{
    class LineSearch_brent : public LineSearch
    {
    private:
    protected:
        type0 tol,zeps;
        int max_iter;
    public:
        LineSearch_brent(MAPP*,int,char**);
        ~LineSearch_brent();
        int line_min(type0&,type0&,int);
    };
    
}
#endif
#endif
