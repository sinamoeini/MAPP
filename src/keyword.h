#ifndef __MAPP__keyword__
#define __MAPP__keyword__
#include "var.h"
#include "logics.h"
namespace MAPP_NS
{
    class KeyWord
    {
    private:
    protected:
    public:
        int nvars;
        Var** vars;
        VarManager* var_mngr;
        Logics* cmd_logic;
        char* name;
        bool st_set;
        
        KeyWord();
        KeyWord(VarManager&);
        ~KeyWord();

        
        void init(Logics*,VarManager&,const char*);
        template <typename... T0>
        void init(Logics*,VarManager&,const char*,T0...);
        void init(Logics* clog,VarManager&,bool&,const char*);
        template <typename... T0>
        void init(Logics*,VarManager&,bool&,const char*,T0...);
        
        
        template<typename T0>
        void add_var(T0&,const char*);
        
        char* scan(char**& args,int& nargs);
        char* finalize();
        Var& operator [] (int i);
        
    };
    
    
    
    class FileKeyWord
    {
    private:
    protected:
    public:
        Logics* cmd_logic;
        Var* keyw;
        Var* var;
        VarManager* var_mngr;
        char* name;
        bool st_set;
        
        FileKeyWord();
        ~FileKeyWord();
        
        template <typename T0>
        void init(Logics*,VarManager&,T0&,const char*);
        template <typename T0,typename... T1>
        void init(Logics* clog,VarManager& vm,T0& v,const char* format,T1... vals);
        
        template <typename T0,typename... T1>
        void init(Logics*,VarManager&,T0&,bool&,const char*,T1...);
        template <typename T0>
        void init(Logics*,VarManager&,T0&,bool&, const char*);
        
        
        char* scan(char**,int);
        char* finalize();
        
        
    };
    
};
using namespace MAPP_NS;
/*----------------------------------------------------------------
 _   _    _____  __    __  _          __  _____   _____    _____  
| | / /  | ____| \ \  / / | |        / / /  _  \ |  _  \  |  _  \ 
| |/ /   | |__    \ \/ /  | |  __   / /  | | | | | |_| |  | | | | 
| |\ \   |  __|    \  /   | | /  | / /   | | | | |  _  /  | | | | 
| | \ \  | |___    / /    | |/   |/ /    | |_| | | | \ \  | |_| | 
|_|  \_\ |_____|  /_/     |___/|___/     \_____/ |_|  \_\ |_____/
 ----------------------------------------------------------------*/
/*--------------------------------------------
 
 --------------------------------------------*/
template <typename... T0>
void KeyWord::init(Logics* clog,VarManager& vm,const char* format,T0... vals)
{
    cmd_logic=clog;
    var_mngr=&vm;
    size_t len=snprintf(NULL,0,format,vals...)+1;
    char* name_=new char[len];
    sprintf(name_,format,vals...);
    
    Var** vars_=new Var*[nvars+1];
    memcpy(vars_,vars,nvars*sizeof(Var*));
    delete [] vars;
    vars=vars_;
    
    vars[nvars]=var_mngr->adddd_var(st_set,name_);
    name=vars[0]->name;
    nvars++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <typename... T0>
void KeyWord::init(Logics* clog,VarManager& vm,bool& set,const char* format,T0... vals)
{
    cmd_logic=clog;
    var_mngr=&vm;
    size_t len=snprintf(NULL,0,format,vals...)+1;
    char* name_=new char[len];
    sprintf(name_,format,vals...);
    
    Var** vars_=new Var*[nvars+1];
    memcpy(vars_,vars,nvars*sizeof(Var*));
    delete [] vars;
    vars=vars_;
    
    vars[nvars]=var_mngr->adddd_var(set,name_);
    name=vars[0]->name;
    nvars++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void KeyWord::add_var(T0& v,const char* var_name)
{
    Var** vars_=new Var*[nvars+1];
    memcpy(vars_,vars,nvars*sizeof(Var*));
    delete [] vars;
    vars=vars_;
    
    char* v_name=NULL;
    if(var_name==NULL)
        Var::append(v_name,"entry #%d after %s",nvars+1,vars[0]->name);
    else
        Var::append(v_name,"%s (entry #%d after %s)",var_name,nvars+1,vars[0]->name);
    
    vars[nvars]=var_mngr->adddd_var(v,v_name);
    delete [] v_name;
    nvars++;
}

/*-------------------------------------------------------------------------------------------
 _____   _   _       _____   _   _    _____  __    __  _          __  _____   _____    _____  
|  ___| | | | |     | ____| | | / /  | ____| \ \  / / | |        / / /  _  \ |  _  \  |  _  \ 
| |__   | | | |     | |__   | |/ /   | |__    \ \/ /  | |  __   / /  | | | | | |_| |  | | | | 
|  __|  | | | |     |  __|  | |\ \   |  __|    \  /   | | /  | / /   | | | | |  _  /  | | | | 
| |     | | | |___  | |___  | | \ \  | |___    / /    | |/   |/ /    | |_| | | | \ \  | |_| | 
|_|     |_| |_____| |_____| |_|  \_\ |_____|  /_/     |___/|___/     \_____/ |_|  \_\ |_____/
 -------------------------------------------------------------------------------------------*/
/*--------------------------------------------
 
 --------------------------------------------*/
template <typename T0,typename... T1>
void FileKeyWord::init(Logics* clog,VarManager& vm,T0& v,const char* format,T1... vals)
{
    cmd_logic=clog;
    var_mngr=&vm;
    size_t len=snprintf(NULL,0,format,vals...)+1;
    char* name_=new char[len];
    sprintf(name_,format,vals...);
    
    keyw=var_mngr->adddd_var(st_set,name_);
    name=keyw->name;
    var=var_mngr->adddd_var(v,name);
}

/*--------------------------------------------
 
 --------------------------------------------*/
template <typename T0,typename... T1>
void FileKeyWord::init(Logics* clog,VarManager& vm,T0& v,bool& set,const char* format,T1... vals)
{
    cmd_logic=clog;
    var_mngr=&vm;
    size_t len=snprintf(NULL,0,format,vals...)+1;
    char* name_=new char[len];
    sprintf(name_,format,vals...);
    
    keyw=var_mngr->adddd_var(set,name_);
    name=keyw->name;
    var=var_mngr->adddd_var(v,name);
}

/*--------------------------------------------
 
 --------------------------------------------*/
template <typename T0>
void FileKeyWord::init(Logics* clog,VarManager& vm,T0& v,const char* name_)
{
    cmd_logic=clog;
    var_mngr=&vm;
    keyw=var_mngr->adddd_var(st_set,name_);
    name=keyw->name;
    var=var_mngr->adddd_var(v,name);
}

/*--------------------------------------------
 
 --------------------------------------------*/
template <typename T0>
void FileKeyWord::init(Logics* clog,VarManager& vm,T0& v,bool& set, const char* name_)
{
    cmd_logic=clog;
    keyw=var_mngr->adddd_var(set,name_);
    name=keyw->name;
    var=var_mngr->adddd_var(v,name);
}









#endif

