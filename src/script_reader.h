#ifndef __MAPP__script_reader__
#define __MAPP__script_reader__
#include <stdio.h>
namespace MAPP_NS
{
    class ScriptReader
    {
    private:
        char* line;
        int line_len;
        int sz_inc;
        FILE*& fp;
        bool get_line_proc0(char*&,int&,int&);
        bool is_char(char*);
    protected:
    public:
        ScriptReader(FILE*&);
        ~ScriptReader();
        int operator()(char**&,int&);
        //extra functions that can be used by other files
        static int read_line(FILE*,char*&);
        static int read_line(FILE*,char*&,int&,int&);
        static int parse_line(const char*,char**&);
        static int parse_line(char*,char**&,int&);
        static int hash_remover(char*,char*&);
        static int hash_remover(char*&);
        static int concatenate(int,char**,char*&);
        static void open_file(FILE*&,const char*,const char*);
        
        
        template<class C>
        static void create(C*& ptr,const char*line)
        {
            if(ptr!=NULL)
                delete ptr;
            ptr=NULL;
            
            char** args;
            int nargs=parse_line(line,args);
            ptr=new C(nargs,args);
            for(int i=0;i<nargs;i++)
                delete [] args[i];
            if(nargs)
                delete [] args;
        }
    };
}
#endif
