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
    public:
        static int my_p;

        Error();
        ~Error();
        static void abort(const char*,...);
        static void abort_sing(const char*,...);
        static void warning(const char*,...);
        static void print(const char*,...);
        static void abort();        
    };
}
#endif
