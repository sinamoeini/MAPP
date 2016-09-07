#include "var.h"
#include "atom_types.h"

using namespace MAPP_NS;
const char* verb_aux ="should";
const char* negate="not";
const char* verb="be";
const char* s_verb="is";
const char* p_verb="are";


bool FindReplace::is_same(char* src,char* key)
{
    while(1)
    {
        if(*key=='\0')
            return true;
        
        if(*src!=*key)
            return false;
        src++;
        key++;
    }
    
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int FindReplace::look_up(char* s)
{
    for(int i=0;i<n;i++)
        if(is_same(s,srch[i]))
            return i;
    
    return -1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
FindReplace::FindReplace()
{
    srch=rplc=NULL;
    srch_len=rplc_len=NULL;
    n=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
FindReplace::~FindReplace()
{
    for(int i=0;i<n;i++)
    {
        delete [] srch[i];
        delete [] rplc[i];
    }
    delete [] srch;
    delete [] rplc;
    delete [] srch_len;
    delete [] rplc_len;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void FindReplace::reverse()
{
    std::swap(srch,rplc);
    std::swap(srch_len,rplc_len);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void FindReplace::operator()(char*& buff)
{
    if(buff==NULL)
        return;
        
    char* s=buff;
    char* str=NULL;
    size_t str_len=0;
    
    int ilook;
    while(*s!='\0')
    {
        ilook=look_up(s);
        if(ilook==-1)
        {
            attach(str,str_len,*s);
            s++;
            continue;
        }
        attach(str,str_len,ilook);
        s+=srch_len[ilook];
    }
    
    attach(str,str_len,'\0');
    delete [] buff;
    buff=str;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void FindReplace::operator()(const char* s,const char* r)
{
    char** srch_=new char*[n+1];
    memcpy(srch_,srch,n*sizeof(char*));
    delete [] srch;
    srch=srch_;
    
    char** rplc_=new char*[n+1];
    memcpy(rplc_,rplc,n*sizeof(char*));
    delete [] rplc;
    rplc=rplc_;
    
    size_t* srch_len_=new size_t[n+1];
    memcpy(srch_len_,srch_len,n*sizeof(size_t));
    delete [] srch_len;
    srch_len=srch_len_;
    
    size_t* rplc_len_=new size_t[n+1];
    memcpy(rplc_len_,rplc_len,n*sizeof(size_t));
    delete [] rplc_len;
    rplc_len=rplc_len_;
    
    srch_len[n]=strlen(s);
    rplc_len[n]=strlen(r);
    srch[n]=new char[srch_len[n]+1];
    memcpy(srch[n],s,srch_len[n]+1);
    rplc[n]=new char[rplc_len[n]+1];
    memcpy(rplc[n],r,rplc_len[n]+1);
    
    n++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void FindReplace::attach(char*& str,size_t& str_len,int i)
{
    
    char* str_=new char[str_len+rplc_len[i]];
    memcpy(str_,str,str_len);
    memcpy(str_+str_len,rplc[i],rplc_len[i]);
    delete [] str;
    str=str_;
    str_len+=rplc_len[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void FindReplace::attach(char*& str,size_t& str_len,char c)
{
    
    char* str_=new char[str_len+1];
    memcpy(str_,str,str_len);
    str_[str_len]=c;
    delete [] str;
    
    str=str_;
    str_len++;
}

/*-------------------------------------------------------------------------------------------
 _     _   ___   _____        ___  ___       ___   __   _       ___   _____   _____   _____   
| |   / / /   | |  _  \      /   |/   |     /   | |  \ | |     /   | /  ___| | ____| |  _  \  
| |  / / / /| | | |_| |     / /|   /| |    / /| | |   \| |    / /| | | |     | |__   | |_| |  
| | / / / / | | |  _  /    / / |__/ | |   / / | | | |\   |   / / | | | |  _  |  __|  |  _  /  
| |/ / / /  | | | | \ \   / /       | |  / /  | | | | \  |  / /  | | | |_| | | |___  | | \ \  
|___/ /_/   |_| |_|  \_\ /_/        |_| /_/   |_| |_|  \_| /_/   |_| \_____/ |_____| |_|  \_\
 -------------------------------------------------------------------------------------------*/
VarManager::VarManager()
{
    MAPP_NS::g_vm=this;
    
    old_vm=MAPP_NS::g_vm;
    st_vars=dy_vars=NULL;
    nst_vars=ndy_vars=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
VarManager::~VarManager()
{   
    for(int i=0;i<ndy_vars;i++)
        delete dy_vars[i];
    delete [] dy_vars;
    
    for(int i=0;i<nst_vars;i++)
        delete st_vars[i];
    delete [] st_vars;

    MAPP_NS::g_vm=old_vm;
}
/*----------------------
 _     _   ___   _____
| |   / / /   | |  _  \
| |  / / / /| | | |_| |
| | / / / / | | |  _  /
| |/ / / /  | | | | \ \  
|___/ /_/   |_| |_|  \_\
 -----------------------*/
Var::Var(const char* name_,void* addr_):
addr(addr_)
{
    set=false;
    size_t len=strlen(name_)+1;
    name=new char[len];
    memcpy(name,name_,len);
    loc=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Var::Var():
addr(NULL)
{
    set=true;    
    name=NULL;
    loc=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Var::~Var()
{
    delete [] name;
    delete [] loc;
}
/*--------------------------------------------
 
 --------------------------------------------*/
char* Var::scan(char*)
{
    set=true;
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::operator < (Var&)
{
    printf("this is a dummy function %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::operator > (Var&)
{
    printf("this is a dummy function %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::operator <= (Var&)
{
    printf("this is a dummy function %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::operator >= (Var&)
{
    printf("this is a dummy function %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Var::operator == (Var&)
{
    printf("this is a dummy function %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Var* Var::clone()
{
    printf("this is a dummy function %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Var::det_loc()
{
    if(name==NULL)
        return;
    
    char* p=strchr(name,'(');
    if(p==NULL)
    {
        p=strchr(name,'#');
        if(p==NULL)
            return;
        p=name;
        name=loc;
        loc=p;
        return;
    }
    
    
    
    p--;
    
    size_t len=strlen(p)+1;
    delete [] loc;
    loc=new char[len];
    memcpy(loc,p,len);
    len=p-name+1;
    char* name_=new char[len];
    memcpy(name_,name,len);
    name_[len-1]='\0';
    delete [] name;
    name=name_;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Var::att_loc()
{
    if(loc==NULL)
        return;
    if(strlen(name)==0)
    {
        delete [] name;
        name=loc;
        loc=NULL;
        return;
    }
    
    char* name_=NULL;
    append(name_,"%s%s",name,loc);
    delete [] name;
    delete [] loc;
    name=name_;
    loc=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
const char* Var::get_type_name()
{
    printf("this is a dummy function %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return NULL;
}
/*----------------------------------------------------------------
  _     _   ___   ____         __  _   __   _   _____  __
 | |   / / /   | |  _  \      / / | | |  \ | | |_   _| \ \
 | |  / / / /| | | |_| |     / /  | | |   \| |   | |    \ \
 | | / / / / | | |  _  /    { {   | | | |\   |   | |     } }
 | |/ / / /  | | | | \ \     \ \  | | | | \  |   | |    / /
 |___/ /_/   |_| |_|  \_\     \_\ |_| |_|  \_|   |_|   /_/
 ----------------------------------------------------------------*/
template<>
var<int>::var(const int val_):
Var(),
val(st_val)
{
    size_t len=snprintf(NULL,0,"%d",val_)+1;
    name=new char[len];
    sprintf(name,"%d",val_);
    val=val_;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
char* var<int>::sub_scan(char* inp)
{
    size_t len=strlen(inp);
    if(len==0)
    {
        char* err_msg=NULL;
        append(err_msg,"%s cannot be empty",name);
        return err_msg;
    }
    bool sign=true;
    bool E=true;
    size_t ipos=0;
    while (ipos!=len)
    {
        if(sign)
        {
            if(inp[ipos]=='+' || inp[ipos]=='-')
            {
                ipos++;
                continue;
            }
            
            sign=false;
        }
        
        
        if(E && (inp[ipos]=='e' || inp[ipos]=='E'))
        {
            E=false;
            sign=true;
            ipos++;
            continue;
        }
        
        if(isdigit(inp[ipos]))
        {
            ipos++;
            continue;
        }
        
        char* err_msg=NULL;
        append(err_msg,"illegal value: %s\n%s should be int",inp,name);
        return err_msg;
        
    }
    val=atoi(inp);
    
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
const char* var<int>::get_type_name()
{
    return "int";
}
/*------------------------------------------------------------------------------------------
  _     _   ___   ____         __  _____   _____   _   _   _____   _       _____  __ 
 | |   / / /   | |  _  \      / / |  _  \ /  _  \ | | | | |  _  \ | |     | ____| \ \ 
 | |  / / / /| | | |_| |     / /  | | | | | | | | | | | | | |_| | | |     | |__    \ \ 
 | | / / / / | | |  _  /    { {   | | | | | | | | | | | | |  _  { | |     |  __|    } }
 | |/ / / /  | | | | \ \     \ \  | |_| | | |_| | | |_| | | |_| | | |___  | |___   / / 
 |___/ /_/   |_| |_|  \_\     \_\ |_____/ \_____/ \_____/ |_____/ |_____| |_____| /_/ 
 ------------------------------------------------------------------------------------------*/
template<>
var<double>::var(const double val_):
Var(),
val(st_val)
{
    size_t len=snprintf(NULL,0,"%g",val_)+1;
    name=new char[len];
    sprintf(name,"%g",val_);
    val=val_;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
char* var<double>::sub_scan(char* inp)
{
    size_t len=strlen(inp);
    if(len==0)
    {
        char* err_msg=NULL;
        append(err_msg,"%s cannot be empty",name);
        return err_msg;
    }
    
    bool sign=true;
    bool E=true;
    bool dot=true;
    size_t ipos=0;
    while (ipos!=len)
    {
        if(sign)
        {
            if(inp[ipos]=='+' || inp[ipos]=='-')
            {
                ipos++;
                continue;
            }
            
            sign=false;
        }
        
        if(dot && inp[ipos]=='.')
        {
            dot=false;
            ipos++;
            continue;
        }
        
        if(E && (inp[ipos]=='e' || inp[ipos]=='E'))
        {
            E=false;
            sign=true;
            ipos++;
            continue;
        }
        
        if(isdigit(inp[ipos]))
        {
            ipos++;
            continue;
        }
        
        char* err_msg=NULL;
        append(err_msg,"illegal value: %s\n%s should be double",inp,name);
        return err_msg;
    }
    
    if(dot)
    {
        char* err_msg=NULL;
        append(err_msg,"illegal value: %s\n%s should be double",inp,name);
        return err_msg;
    }
    
    val=atof(inp);
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
const char* var<double>::get_type_name()
{
    return "double";
}
/*------------------------------------------------------------------------------------
  _     _   ___   ____         __  _____   _   _       ___   _____   __   _   __ __     
 | |   / / /   | |  _  \      / / /  ___| | | | |     /   | |  _  \  \ \ | | / / \ \  
 | |  / / / /| | | |_| |     / /  | |     | |_| |    / /| | | |_| |   \ \| |/ /   \ \ 
 | | / / / / | | |  _  /    { {   | |     |  _  |   / / | | |  _  /    }     {     } }
 | |/ / / /  | | | | \ \     \ \  | |___  | | | |  / /  | | | | \ \   / /| |\ \   / / 
 |___/ /_/   |_| |_|  \_\     \_\ \_____| |_| |_| /_/   |_| |_|  \_\ /_/ |_| \_\ /_/  
 -------------------------------------------------------------------------------------*/
template<>
var<char*>::var(const char* v):
Var(),
val(st_val)
{
    size=static_cast<int>(strlen(v))+1;
    st_val=new char[size];
    name=new char[size];
    memcpy(name,v,size);
    memcpy(st_val,v,size);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
char* var<char*>::sub_scan(char* inp)
{
    size_t len=strlen(inp);
    if(len==0)
    {
        char* err_msg=NULL;
        append(err_msg,"%s cannot be empty",name);
        return err_msg;
    }
    len++;
    size=static_cast<int>(len);
    delete [] st_val;
    delete [] val;
    st_val=new char[len];
    memcpy(st_val,inp,len);
    val=st_val;
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
bool var<char*>::operator < (Var& right)
{
    int other_size;
    char* other_val=get_val(right,other_size);
    if(this->size!=other_size)
        return false;
    for(int i=0;i<this->size;i++)
        if(val[i]>=other_val[i])
            return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
bool var<char*>::operator > (Var& right)
{
    int other_size;
    char* other_val=get_val(right,other_size);
    if(this->size!=other_size)
        return false;
    for(int i=0;i<this->size;i++)
        if(val[i]<=other_val[i])
            return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
bool var<char*>::operator <= (Var& right)
{
    int other_size;
    char* other_val=get_val(right,other_size);
    if(this->size!=other_size)
        return false;
    for(int i=0;i<this->size;i++)
        if(val[i]>other_val[i])
            return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
bool var<char*>::operator >= (Var& right)
{
    int other_size;
    char* other_val=get_val(right,other_size);
    if(this->size!=other_size)
        return false;
    for(int i=0;i<this->size;i++)
        if(val[i]<other_val[i])
            return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
bool var<char*>::operator == (Var& right)
{
    int other_size;
    char* other_val=get_val(right,other_size);
    if(this->size!=other_size)
        return false;
    for(int i=0;i<this->size;i++)
        if(val[i]!=other_val[i])
            return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<>
const char* var<char*>::get_type_name()
{
    return "string";
}
/*----------------------------------------------------------------------
 _       _____   _____   _   _____   _____        _       _____   _____
| |     /  _  \ /  ___| | | /  ___| /  ___/      | |     /  _  \ /  ___|
| |     | | | | | |     | | | |     | |___       | |     | | | | | |
| |     | | | | | |  _  | | | |     \___  \      | |     | | | | | |  _
| |___  | |_| | | |_| | | | | |___   ___| |      | |___  | |_| | | |_| |
|_____| \_____/ \_____/ |_| \_____| /_____/      |_____| \_____/ \_____/
 ----------------------------------------------------------------------*/
Logics::Log::Log(const char* op_name_,const bool _is):
_is_(_is)
{
    size_t len=strlen(op_name_)+1;
    op_name=new char[len];
    memcpy(op_name,op_name_,len);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Log::~Log()
{
    delete [] op_name;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Logics::Log::eq(const char* st0,const char* st1)
{
    if(st0==NULL && st1==NULL)
        return 1;
    if(st0==NULL)
        return 0;
    if(st1==NULL)
        return 0;
    if(strcmp(st0,st1))
        return 0;
    
    return 1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::Log::finish_sntnc(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    if(L==NULL)
        return;
    if(R!=NULL)
        Var::append(err_msg,"%s %s %s %s",L,_verb_,op,R);
    else
        Var::append(err_msg,"%s %s %s",L,_verb_,op);
    
    
    delete [] L;
    delete [] _verb_;
    delete [] op;
    delete [] R;
    
    L=NULL;
    _verb_=NULL;
    op=NULL;
    R=NULL;
    p_flag=false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::Log::start_sntnc(char*& L,char*& _verb_,char*& op,char*& R,int v_flag,char* l_name,char* r_name)
{
    
    
    delete [] L;
    delete [] _verb_;
    delete [] op;
    delete [] R;
    
    L=NULL;
    _verb_=NULL;
    op=NULL;
    R=NULL;
    
    Var::append(L,"%s",l_name);
    
    if(v_flag==AUX_VERB && _is_)
        Var::append(_verb_,"%s %s",verb_aux,verb);
    else if(v_flag==AUX_VERB && !_is_)
        Var::append(_verb_,"%s %s %s",verb_aux,negate,verb);
    else if(v_flag==VERB && _is_)
        Var::append(_verb_,"%s",s_verb);
    else if(v_flag==VERB && !_is_)
        Var::append(_verb_,"%s %s",s_verb,negate);
    
    Var::append(op,"%s",op_name);
    if(r_name!=NULL)
        Var::append(R,"%s",r_name);
    else
        r_name=NULL;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::Log::process_sntnc(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,char* l_name,char* r_name)
{
    if(L==NULL)
    {
        start_sntnc(L,_verb_,op,R,v_flag,l_name,r_name);
        __is__=_is_;
        p_flag=false;
        
        return;
    }
    
    int case_no=eq(R,r_name)+2*eq(op,op_name)+4*eq(L,l_name);
    if(_is_!=__is__ || case_no<3 || case_no==7)
    {
        finish_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
        Var::append(err_msg,", %s ",del);
        print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
        return;
    }
    
    if(case_no==3)
    {
        if(!p_flag)
        {
            if(v_flag==VERB && _is_)
            {
                delete [] _verb_;
                _verb_=NULL;
                Var::append(_verb_,"%s",p_verb);
            }
            else if(v_flag==VERB && !_is_)
            {
                delete [] _verb_;
                _verb_=NULL;
                Var::append(_verb_,"%s %s",p_verb,negate);
            }
        }
        p_flag=true;
        Var::append(L," %s %s",del,l_name);
        
    }
    else if(case_no==4)
    {
        if(R!=NULL)
            Var::append(op," %s %s %s",R,del,op_name);
        else
            Var::append(op," %s %s",del,op_name);
        delete [] R;
        R=NULL;
        if(r_name!=NULL)
            Var::append(R,"%s",r_name);
        else
            R=NULL;
    }
    else if(case_no==5)
    {
        Var::append(op," %s %s",del,op_name);
    }
    else if(case_no==6)
    {
        if(r_name!=NULL)
            Var::append(R," %s %s",del,r_name);
    }
}
/*---------------------------------------------------------------------------------------------------------------------------------------
 _       _____   _____   _   _____   _____        _       _____   _____   _____   _   __   _       ___   _____   __    __  _____   _____  
| |     /  _  \ /  ___| | | /  ___| /  ___/      | |     /  _  \ /  ___| |  _  \ | | |  \ | |     /   | |  _  \  \ \  / / /  _  \ |  _  \ 
| |     | | | | | |     | | | |     | |___       | |     | | | | | |     | |_| | | | |   \| |    / /| | | |_| |   \ \/ /  | | | | | |_| | 
| |     | | | | | |  _  | | | |     \___  \      | |     | | | | | |  _  |  _  { | | | |\   |   / / | | |  _  /    \  /   | | | | |  ___/ 
| |___  | |_| | | |_| | | | | |___   ___| |      | |___  | |_| | | |_| | | |_| | | | | | \  |  / /  | | | | \ \    / /    | |_| | | |     
|_____| \_____/ \_____/ |_| \_____| /_____/      |_____| \_____/ \_____/ |_____/ |_| |_|  \_| /_/   |_| |_|  \_\  /_/     \_____/ |_|
 ---------------------------------------------------------------------------------------------------------------------------------------*/
Logics::LogBinaryOp::LogBinaryOp(Log* l,Log* r,const char* op_name,const bool _is_):
Log(op_name,_is_),left(l),right(r)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogBinaryOp::replace(Var* ref,Var* subs)
{
    left->replace(ref,subs);
    right->replace(ref,subs);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::LogBinaryOp::~LogBinaryOp()
{
    delete left;
    delete right;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogBinaryOp::print(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    left->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    del=op_name;
    right->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogBinaryOp::print(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,Var* v)
{
    left->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag,v);
    del=op_name;
    right->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag,v);
}
/*--------------------------------------------------------------------------------------------------------------
 _       _____   _____   _   _____   _____        _       _____   _____        _____   _____   __   _   _____  
| |     /  _  \ /  ___| | | /  ___| /  ___/      | |     /  _  \ /  ___|      /  ___| /  _  \ |  \ | | |  _  \ 
| |     | | | | | |     | | | |     | |___       | |     | | | | | |          | |     | | | | |   \| | | | | | 
| |     | | | | | |  _  | | | |     \___  \      | |     | | | | | |  _       | |     | | | | | |\   | | | | | 
| |___  | |_| | | |_| | | | | |___   ___| |      | |___  | |_| | | |_| |      | |___  | |_| | | | \  | | |_| | 
|_____| \_____/ \_____/ |_| \_____| /_____/      |_____| \_____/ \_____/      \_____| \_____/ |_|  \_| |_____/
 --------------------------------------------------------------------------------------------------------------*/
Logics::LogCond::LogCond(Log* l,Log* r,const char* op_name,const bool _is_):
Log(op_name,_is_),left(l),right(r)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogCond::replace(Var* ref,Var* subs)
{
    left->replace(ref,subs);
    right->replace(ref,subs);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::LogCond::~LogCond()
{
    delete left;
    delete right;
}
/*------------------------------------------------------------------------------------------------------------------------
 _       _____   _____   _   _____   _____        _       _____   _____   _____   _   __   _       ___   _____   __    __
| |     /  _  \ /  ___| | | /  ___| /  ___/      | |     /  _  \ /  ___| |  _  \ | | |  \ | |     /   | |  _  \  \ \  / /
| |     | | | | | |     | | | |     | |___       | |     | | | | | |     | |_| | | | |   \| |    / /| | | |_| |   \ \/ /
| |     | | | | | |  _  | | | |     \___  \      | |     | | | | | |  _  |  _  { | | | |\   |   / / | | |  _  /    \  /
| |___  | |_| | | |_| | | | | |___   ___| |      | |___  | |_| | | |_| | | |_| | | | | | \  |  / /  | | | | \ \    / /
|_____| \_____/ \_____/ |_| \_____| /_____/      |_____| \_____/ \_____/ |_____/ |_| |_|  \_| /_/   |_| |_|  \_\  /_/
 ------------------------------------------------------------------------------------------------------------------------*/
Logics::LogBinary::LogBinary(Var* l,Var* r,const char* op_name,const bool _is_):
Log(op_name,_is_),left(l),right(r)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogBinary::replace(Var* ref,Var* subs)
{

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogBinary::print(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    process_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag,left->name,right->name);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogBinary::print(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,Var* v)
{
    int i=0;
    if(left==NULL)
    {
        left=v;
        i++;
    }
    if(right==NULL)
    {
        right=v;
        i+=2;
    }
    process_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag,left->name,right->name);
    
    if(i>1)
    {
        right=NULL;
        i-=2;
    }
    
    if(i)
        left=NULL;
}
/*--------------------------------------------------------------------------------------------------------------------
 _       _____   _____   _   _____   _____        _       _____   _____   _   _   __   _       ___   _____   __    __ 
| |     /  _  \ /  ___| | | /  ___| /  ___/      | |     /  _  \ /  ___| | | | | |  \ | |     /   | |  _  \  \ \  / / 
| |     | | | | | |     | | | |     | |___       | |     | | | | | |     | | | | |   \| |    / /| | | |_| |   \ \/ /  
| |     | | | | | |  _  | | | |     \___  \      | |     | | | | | |  _  | | | | | |\   |   / / | | |  _  /    \  /   
| |___  | |_| | | |_| | | | | |___   ___| |      | |___  | |_| | | |_| | | |_| | | | \  |  / /  | | | | \ \    / /    
|_____| \_____/ \_____/ |_| \_____| /_____/      |_____| \_____/ \_____/ \_____/ |_|  \_| /_/   |_| |_|  \_\  /_/
 --------------------------------------------------------------------------------------------------------------------*/
Logics::LogUnary::LogUnary(Var* l,const char* op_name,const bool _is_):
Log(op_name,_is_),left(l)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogUnary::replace(Var* ref,Var* subs)
{

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogUnary::print(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    process_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag,left->name,NULL);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogUnary::print(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,Var* v)
{
    int i=0;
    if(left==NULL)
    {
        left=v;
        i++;
    }
    process_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag,left->name,NULL);
    
    if(i)
        left=NULL;
}
bool Logics::LogUnary::operator()(Var* v)
{
    int i=0;
    if(left==NULL)
    {
        left=v;
        i++;
    }
    
    bool ans =this->operator()();

    
    if(i)
        left=NULL;
    return ans;
}
/*----------------------------------
 _       _____   _____   _   _____  
| |     /  _  \ /  ___| | | /  ___/ 
| |     | | | | | |     | | | |___  
| |     | | | | | |  _  | | \___  \ 
| |___  | |_| | | |_| | | |  ___| | 
|_____| \_____/ \_____/ |_| /_____/
 ----------------------------------*/
Logics::LogIS::LogIS(Var* l,const bool _is_):LogUnary(l,"set",_is_)
{}
bool Logics::LogIS::operator() ()
{
    if(_is_)
        return (left->set);
    else
        return !(left->set);
}
Logics::Log* Logics::LogIS::clone()
{
    return new LogIS(left,_is_);
}
bool Logics::LogIS::operator()(Var* v)
{
    int i=0;
    if(left==NULL)
    {
        left=v;
        i++;
    }
    
    bool ans =this->operator()();
    
    
    if(i)
        left=NULL;
    return ans;
}
/*-------------------------------------
 _       _____   _____   _____   _____  
| |     /  _  \ /  ___| /  ___| |_   _| 
| |     | | | | | |     | |       | |   
| |     | | | | | |  _  | |  _    | |   
| |___  | |_| | | |_| | | |_| |   | |   
|_____| \_____/ \_____/ \_____/   |_|
 -------------------------------------*/
Logics::LogGT::LogGT(Var* l,Var* r,const bool _is_):LogBinary(l,r,"greater than",_is_)
{}
bool Logics::LogGT::operator() ()
{
    if(_is_)
        return (*left > *right);
    else
        return !(*left > *right);
}
Logics::Log* Logics::LogGT::clone()
{
    return new LogGT(left,right,_is_);
}
bool Logics::LogGT::operator()(Var* v)
{
    int i=0;
    if(left==NULL)
    {
        left=v;
        i++;
    }
    if(right==NULL)
    {
        right=v;
        i+=2;
    }
    bool ans =this->operator()();
    
    if(i>1)
    {
        right=NULL;
        i-=2;
    }
    
    if(i)
        left=NULL;
    return ans;
}
/*--------------------------------------
 _       _____   _____   _____   _____  
| |     /  _  \ /  ___| /  ___| | ____| 
| |     | | | | | |     | |     | |__   
| |     | | | | | |  _  | |  _  |  __|  
| |___  | |_| | | |_| | | |_| | | |___  
|_____| \_____/ \_____/ \_____/ |_____|
 --------------------------------------*/
Logics::LogGE::LogGE(Var* l,Var* r,const bool _is_):LogBinary(l,r,"greater than or equal to",_is_)
{}
bool Logics::LogGE::operator() ()
{
    if(_is_)
        return (*left >= *right);
    else
        return !(*left >= *right);
}
Logics::Log* Logics::LogGE::clone()
{
    return new LogGE(left,right,_is_);
}
bool Logics::LogGE::operator()(Var* v)
{
    int i=0;
    if(left==NULL)
    {
        left=v;
        i++;
    }
    if(right==NULL)
    {
        right=v;
        i+=2;
    }
    bool ans =this->operator()();
    
    if(i>1)
    {
        right=NULL;
        i-=2;
    }
    
    if(i)
        left=NULL;
    return ans;
}
/*----------------------------------------
 _       _____   _____   _____   _____    
| |     /  _  \ /  ___| | ____| /  _  \   
| |     | | | | | |     | |__   | | | |   
| |     | | | | | |  _  |  __|  | | | |   
| |___  | |_| | | |_| | | |___  | |_| |_  
|_____| \_____/ \_____/ |_____| \_______| 
 ----------------------------------------*/
Logics::LogEq::LogEq(Var* l,Var* r,const bool _is_):LogBinary(l,r,"equal to",_is_)
{}
bool Logics::LogEq::operator() ()
{
    if(_is_)
    {
        
        return (*left == *right);
    }
    else
        return !(*left == *right);
}
Logics::Log* Logics::LogEq::clone()
{
    return new LogEq(left,right,_is_);
}
bool Logics::LogEq::operator()(Var* v)
{
    int i=0;
    if(left==NULL)
    {
        left=v;
        i++;
    }
    if(right==NULL)
    {
        right=v;
        i+=2;
    }
    bool ans =this->operator()();
    
    if(i>1)
    {
        right=NULL;
        i-=2;
    }
    
    if(i)
        left=NULL;
    return ans;
}
/*-------------------------------------
 _       _____   _____   _       _____
| |     /  _  \ /  ___| | |     |_   _| 
| |     | | | | | |     | |       | |   
| |     | | | | | |  _  | |       | |   
| |___  | |_| | | |_| | | |___    | |   
|_____| \_____/ \_____/ |_____|   |_|
 -------------------------------------*/
Logics::LogLT::LogLT(Var* l,Var* r,const bool _is_):LogBinary(l,r,"less than",_is_)
{}
bool Logics::LogLT::operator() ()
{
    if(_is_)
        return (*left < *right);
    else
        return !(*left < *right);
}
Logics::Log* Logics::LogLT::clone()
{
    return new LogLT(left,right,_is_);
}
bool Logics::LogLT::operator()(Var* v)
{
    int i=0;
    if(left==NULL)
    {
        left=v;
        i++;
    }
    if(right==NULL)
    {
        right=v;
        i+=2;
    }
    bool ans =this->operator()();
    
    if(i>1)
    {
        right=NULL;
        i-=2;
    }
    
    if(i)
        left=NULL;
    return ans;
}
/*-------------------------------------
 _       _____   _____   _       _____  
| |     /  _  \ /  ___| | |     | ____| 
| |     | | | | | |     | |     | |__   
| |     | | | | | |  _  | |     |  __|  
| |___  | |_| | | |_| | | |___  | |___  
|_____| \_____/ \_____/ |_____| |_____| 
 -------------------------------------*/
Logics::LogLE::LogLE(Var* l,Var* r,const bool _is_):LogBinary(l,r,"less than or equal to",_is_)
{}
bool Logics::LogLE::operator() ()
{
    if(_is_)
        return (*left <= *right);
    else
        return !(*left <= *right);
}
Logics::Log* Logics::LogLE::clone()
{
    return new LogLE(left,right,_is_);
}
bool Logics::LogLE::operator()(Var* v)
{
    int i=0;
    if(left==NULL)
    {
        left=v;
        i++;
    }
    if(right==NULL)
    {
        right=v;
        i+=2;
    }
    bool ans =this->operator()();
    
    if(i>1)
    {
        right=NULL;
        i-=2;
    }
    
    if(i)
        left=NULL;
    return ans;
}
/*-------------------------------------------------
 _       _____   _____       ___   __   _   _____  
| |     /  _  \ /  ___|     /   | |  \ | | |  _  \ 
| |     | | | | | |        / /| | |   \| | | | | | 
| |     | | | | | |  _    / / | | | |\   | | | | | 
| |___  | |_| | | |_| |  / /  | | | | \  | | |_| | 
|_____| \_____/ \_____/ /_/   |_| |_|  \_| |_____/ 
 -------------------------------------------------*/
Logics::LogAnd::LogAnd(Log* l,Log* r):LogBinaryOp(l,r,"and",true)
{}
bool Logics::LogAnd::operator() ()
{return ((*left)() && (*right)());}
bool Logics::LogAnd::operator() (Var* v)
{return ((*left)(v) && (*right)(v));}
Logics::Log* Logics::LogAnd::clone()
{
    Log* l=left->clone();
    Log* r=right->clone();
    return new LogAnd(l,r);
}
/*--------------------------------------
 _       _____   _____   _____   _____   
| |     /  _  \ /  ___| /  _  \ |  _  \  
| |     | | | | | |     | | | | | |_| |  
| |     | | | | | |  _  | | | | |  _  /  
| |___  | |_| | | |_| | | |_| | | | \ \  
|_____| \_____/ \_____/ \_____/ |_|  \_\
 --------------------------------------*/
Logics::LogOr::LogOr(Log* l,Log* r):LogBinaryOp(l,r,"or",true)
{}
bool Logics::LogOr::operator() ()
{return ((*left)() || (*right)());}
bool Logics::LogOr::operator() (Var* v)
{return ((*left)(v) || (*right)(v));}
Logics::Log* Logics::LogOr::clone()
{
    Log* l=left->clone();
    Log* r=right->clone();
    return new LogOr(l,r);
}
/*---------------------------------
 _       _____   _____   _   _____  
| |     /  _  \ /  ___| | | |  ___| 
| |     | | | | | |     | | | |__   
| |     | | | | | |  _  | | |  __|  
| |___  | |_| | | |_| | | | | |     
|_____| \_____/ \_____/ |_| |_|
 ---------------------------------*/
Logics::LogIf::LogIf(Log* l,Log* r):LogCond(l,r,"if",true)
{}
bool Logics::LogIf::operator() ()
{
    if(!(*left)())
        return true;
    if((*right)())        
        return true;
    
    return false;
}
bool Logics::LogIf::operator() (Var* v)
{
    
    if(!(*left)(v))
        return true;
    
    if((*right)(v))
        return true;
    return false;
}
Logics::Log* Logics::LogIf::clone()
{
    Log* l=left->clone();
    Log* r=right->clone();
    return new LogIf(l,r);
}
void Logics::LogIf::print(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    finish_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    if(err_msg !=NULL)
        Var::append(err_msg,". ");
    
    v_flag=AUX_VERB;
    right->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    
    finish_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    
    Var::append(err_msg,", %s ","when");
    
    v_flag=VERB;
    left->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
}
void Logics::LogIf::print(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,Var* v)
{
    finish_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    if(err_msg !=NULL)
        Var::append(err_msg,". ");
    
    finish_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    if(err_msg !=NULL)
        Var::append(err_msg,". ");
    
    v_flag=AUX_VERB;
    right->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag,v);
    
    finish_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    
    Var::append(err_msg,", %s ","when");
    
    v_flag=VERB;
    left->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag,v);
}
/*-----------------------------------------
 _       _____   _____   _   _____   _____  
| |     /  _  \ /  ___| | | |  ___| |  ___| 
| |     | | | | | |     | | | |__   | |__   
| |     | | | | | |  _  | | |  __|  |  __|  
| |___  | |_| | | |_| | | | | |     | |     
|_____| \_____/ \_____/ |_| |_|     |_|
 -----------------------------------------*/
Logics::LogIff::LogIff(Log* l,Log* r):LogCond(l,r,"iff",true)
{}
bool Logics::LogIff::operator() ()
{
    if((*left)())
    {
        if((*right)())
            return true;
    }
    else
        if(!(*right)())
            return true;

    return false;
}
bool Logics::LogIff::operator() (Var* v)
{
    if((*left)(v))
    {
        if((*right)(v))
            return true;
    }
    else
        if(!(*right)(v))
            return true;
    
    return false;
}
Logics::Log* Logics::LogIff::clone()
{
    Log* l=left->clone();
    Log* r=right->clone();
    return new LogIff(l,r);
}
void Logics::LogIff::print(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    finish_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    if(err_msg !=NULL)
        Var::append(err_msg,". ");

    v_flag=AUX_VERB;
    right->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    
    finish_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    
    Var::append(err_msg,", %s ","if and only if");

    v_flag=VERB;
    left->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
}
void Logics::LogIff::print(char*& err_msg,char*& L,char*& _verb_,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,Var* v)
{
    finish_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    if(err_msg !=NULL)
        Var::append(err_msg,". ");
    
    v_flag=AUX_VERB;
    right->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag,v);
    
    finish_sntnc(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag);
    
    Var::append(err_msg,", %s ","if and only if");
    
    v_flag=VERB;
    left->print(err_msg,L,_verb_,op,R,__is__,del,v_flag,p_flag,v);
}
/*-----------------------------------------
 _       _____   _____   _   _____   _____  
| |     /  _  \ /  ___| | | /  ___| /  ___/ 
| |     | | | | | |     | | | |     | |___  
| |     | | | | | |  _  | | | |     \___  \ 
| |___  | |_| | | |_| | | | | |___   ___| | 
|_____| \_____/ \_____/ |_| \_____| /_____/
 -----------------------------------------*/
void Logics::creat_op(Var* l,const char* op_name)
{
    if(!strcmp(op_name,"set"))
        op=new LogIS(l,true);
    else if(!strcmp(op_name,"nset") || !strcmp(op_name,"!set"))
        op=new LogIS(l,false);
    
    next=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::creat_op(Var* l,const char* op_name,Var* r)
{
    if(!strcmp(op_name,"gt"))
        op=new LogGT(l,r,true);
    else if(!strcmp(op_name,"ngt") || !strcmp(op_name,"!gt"))
        op=new LogGT(l,r,false);
    else if(!strcmp(op_name,"ge"))
        op=new LogGE(l,r,true);
    else if(!strcmp(op_name,"nge") || !strcmp(op_name,"!ge"))
        op=new LogGE(l,r,false);
    else if(!strcmp(op_name,"eq"))
        op=new LogEq(l,r,true);
    else if(!strcmp(op_name,"neq") || !strcmp(op_name,"!eq"))
        op=new LogEq(l,r,false);
    else if(!strcmp(op_name,"lt"))
        op=new LogLT(l,r,true);
    else if(!strcmp(op_name,"nlt") || !strcmp(op_name,"!lt"))
        op=new LogLT(l,r,false);
    else if(!strcmp(op_name,"le"))
        op=new LogLE(l,r,true);
    else if(!strcmp(op_name,"nle") || !strcmp(op_name,"!le"))
        op=new LogLE(l,r,false);
    
    next=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Logics()
{
    op=NULL;
    next=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Logics(Var* l,const char* op_name,Var* r)
{
    creat_op(l,op_name,r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Logics(Var* l,const char* op_name)
{
    creat_op(l,op_name);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::~Logics()
{
    delete op;
    delete next;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Logics(Logics& other)
{
    if(other.op!=NULL)
        this->op=other.op->clone();
    else
        this->op=NULL;
    
    if(other.next!=NULL)
    {
        this->next=new Logics(*(other.next));
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics::Logics(Logics&& other)
{
    this->op=other.op;
    other.op=NULL;
    
    this->next=other.next;
    other.next=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator = (Logics& other)
{
    delete this->op;
    if(other.op!=NULL)
        this->op=other.op->clone();
    else
        this->op=NULL;
    
    if(other.next!=NULL)
    {
        delete this->next;
        this->next=new Logics(*(other.next));
    }
    
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator = (Logics&& other)
{
    
    delete this->op;
    this->op=other.op;
    other.op=NULL;
    
    
    delete this->next;
    this->next=other.next;
    other.next=NULL;
    
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics&& Logics::operator + (Logics&& other)
{
    if(this->op==NULL)
    {
        delete this->op;
        this->op=other.op;
        other.op=NULL;
        return std::move(*this);
    }
    Log* op_=new LogOr(this->op,other.op);
    op=op_;
    other.op=NULL;
    return std::move(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics&& Logics::operator * (Logics&& other)
{
    if(this->op==NULL)
    {
        delete this->op;
        op=other.op;
        other.op=NULL;
        return std::move(*this);
    }
    
    Log* op_=new LogAnd(this->op,other.op);
    op=op_;
    other.op=NULL;
    return std::move(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics&& Logics::operator / (Logics&& other)
{
    
    if(this->op==NULL)
    {
        delete this->op;
        this->op=other.op;
        other.op=NULL;
        return std::move(*this);
    }
    
    
    Log* op_=new LogIf(this->op,other.op);
    op=op_;
    other.op=NULL;
    return std::move(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics&& Logics::operator - (Logics&& other)
{
    if(this->op==NULL)
    {
        delete this->op;
        op=other.op;
        other.op=NULL;
        return std::move(*this);
    }
    
    Log* op_=new LogIff(this->op,other.op);
    op=op_;
    other.op=NULL;
    return std::move(*this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator += (Logics&& other)
{
    if(this->op==NULL)
    {
        delete this->op;
        op=other.op;
        other.op=NULL;
        return *this;
    }
    
    Log* op_=new LogOr(this->op,other.op);
    op=op_;
    other.op=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator *= (Logics&& other)
{
    if(this->op==NULL)
    {
        op=other.op;
        other.op=NULL;
        return *this;
    }
    
    Log* op_=new LogAnd(this->op,other.op);
    op=op_;
    other.op=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator -= (Logics&& other)
{
    if(this->op==NULL)
    {
        delete this->op;
        op=other.op;
        other.op=NULL;
        return *this;
    }
    
    Log* op_=new LogIff(this->op,other.op);
    op=op_;
    other.op=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::operator /= (Logics&& other)
{
    if(this->op==NULL)
    {
        op=other.op;
        other.op=NULL;
        return *this;
    }
    
    Log* op_=new LogIf(this->op,other.op);
    op=op_;
    other.op=NULL;
    return *this;
}
/*--------------------------------------------
 
 --------------------------------------------*/
char* Logics::print(Var* v)
{
    int v_flag=AUX_VERB;
    if(op==NULL)
        return NULL;
    
    char* err_msg=NULL;
    char* del;
    char* l=NULL;
    char* r=NULL;
    char* oper=NULL;
    char* _verb_=NULL;
    bool p_flag;
    bool __is__;
    
    op->print(err_msg,l,_verb_,oper,r,__is__,del,v_flag,p_flag,v);
    op->finish_sntnc(err_msg,l,_verb_,oper,r,__is__,del,v_flag,p_flag);
    Var::append(err_msg,".");
    
    delete [] l;
    delete [] oper;
    delete [] r;
    delete [] _verb_;
    
    if(next!=NULL)
    {
        char* next_msg=next->print(v);
        if(next_msg!=NULL)
        {
            char* err_msg_=err_msg;
            err_msg=NULL;
            Var::append(err_msg,"%s %s",err_msg_,next_msg);
            delete [] err_msg_;
        }
        delete [] next_msg;
    }
    
    return err_msg;
}
/*--------------------------------------------
 
 --------------------------------------------*/
char* Logics::operator()(Var* v)
{
    if(op==NULL)
        return NULL;
        
    if(!op->operator()(v))
    {
        return print(v);
    }
    
    if(next!=NULL)
        return next->operator()(v);
    
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
char* Logics::operator()()
{
    return operator()(NULL);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Logics::link()
{
    
    if(next==NULL)
    {
        next=new Logics();
    }
    
    return *next;
}
