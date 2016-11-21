#include "keyword.h"
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
KeyWord::KeyWord()
{
    nvars=0;
    vars=NULL;
    st_set=false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
KeyWord::KeyWord(VarManager& vm):
var_mngr(&vm)
{
    nvars=0;
    vars=NULL;
    st_set=false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
KeyWord::~KeyWord()
{
    delete [] vars;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void KeyWord::init(Logics* clog,VarManager& vm,const char* name_)
{
    cmd_logic=clog;
    var_mngr=&vm;
    Var** vars_=new Var*[nvars+1];
    memcpy(vars_,vars,nvars*sizeof(Var*));
    delete [] vars;
    vars=vars_;
    
    vars[nvars]=var_mngr->adddd_var(st_set,name_);
    
    nvars++;
    
    name=vars[0]->name;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void KeyWord::init(Logics* clog,VarManager& vm,bool& set, const char* name_)
{
    cmd_logic=clog;
    var_mngr=&vm;
    Var** vars_=new Var*[nvars+1];
    memcpy(vars_,vars,nvars*sizeof(Var*));
    delete [] vars;
    vars=vars_;
    
    vars[nvars]=var_mngr->adddd_var(set,name_);
    
    nvars++;
    
    name=vars[0]->name;
}
/*--------------------------------------------
 
 --------------------------------------------*/
char* KeyWord::scan(char**& args,int& nargs)
{
    if(nargs<nvars)
    {
        char* err_msg=NULL;
        Var::append(err_msg,"%d arguments are expected for %s",nvars-1,name);
        return err_msg;
    }
    
    for(int i=0;i<nvars;i++)
    {
        char* err_msg=vars[i]->scan(*args);
        if(err_msg)
            return err_msg;
        args++;
        nargs--;
    }
    
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
char* KeyWord::finalize()
{
    return (*cmd_logic)(vars[0]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Var& KeyWord::operator [] (int i)
{
    if(i+1>nvars)
    {
        printf("exceeding dimension of the keyword %d %d\n",i+1,nvars);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    return *vars[i+1];
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
FileKeyWord::FileKeyWord()
{
    var=NULL;
    keyw=NULL;
    st_set=false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
FileKeyWord::~FileKeyWord()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
char* FileKeyWord::scan(char** args,int nargs)
{
    char* err_msg=keyw->scan(args[0]);
    if(err_msg)
        return err_msg;
    err_msg=var->scan(args[1]);
    if(err_msg)
        return err_msg;
    
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
char* FileKeyWord::finalize()
{
    return (*cmd_logic)(keyw);
}
