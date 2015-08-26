#ifdef LS_Style
    LSStyle(LineSearch_goldensection,goldensection)
#else
#ifndef __MAPP__ls_gs__
#define __MAPP__ls_gs__

#include "ls.h"
namespace MAPP_NS
{
    class LineSearch_goldensection : public LineSearch
    {
    private:
    protected:
        type0 tol;
    public:
        LineSearch_goldensection(MAPP*,int,char**);
        ~LineSearch_goldensection();
        int line_min(type0&,type0&,int);
    };
}
#endif
#endif
