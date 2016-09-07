/*--------------------------------------------
 Created by Sina on 05/11/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__error__
#define __MAPP__error__
#include <string.h>
namespace MAPP_NS
{
    class Error
    {
    private:
    protected:
        int my_no;
    public:
        Error();
        ~Error();
        void abort(const char*,...);
        void abort_sing(const char*,...);
        void warning(const char*,...);
        void print(const char*,...);
        void abort();        
    };
}
#endif
