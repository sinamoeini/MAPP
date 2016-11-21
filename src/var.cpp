#include "var.h"
#include "atom_types.h"

using namespace MAPP_NS;

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

