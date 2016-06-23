/*--------------------------------------------
 Created by Sina on 05/11/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__error__
#define __MAPP__error__
#include "init.h"
namespace MAPP_NS {
    class Error : protected InitPtrs {
    private:
    protected:
        int my_no;
    public:
        Error(MAPP *);
        ~Error();
        void abort(const char*,...);
        void warning(const char*,...);
        void print(const char*,...);
        void abort();
        
        template <typename... T0>
        void append(char*& buff,const char* format,T0... vals)
        {
            size_t len=snprintf(NULL,0,format,vals...)+1;
            if(buff!=NULL)
            {
                size_t old_len=strlen(buff);
                char* buff_=new char[old_len+len];
                memcpy(buff_,buff,old_len);
                sprintf(buff_+old_len,format,vals...);
                delete [] buff;
                buff=buff_;
            }
            else
            {
                buff=new char[len];
                sprintf(buff,format,vals...);
            }
        }
        
        void append(char*& buff,const char* format)
        {
            size_t len=strlen(format)+1;
            if(buff!=NULL)
            {
                size_t old_len=strlen(buff);
                char* buff_=new char[old_len+len];
                memcpy(buff_,buff,old_len);
                memcpy(buff_+old_len,format,len);
                delete [] buff;
                buff=buff_;
            }
            else
            {
                buff=new char[len];
                memcpy(buff,format,len);
            }
        }
        
    };
}
#endif
