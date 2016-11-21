#include "logics.h"

using namespace MAPP_NS;
const char* Logics::Log::verb_aux ="should";
const char* Logics::Log::negate="not";
const char* Logics::Log::verb="be";
/*
const char* Logics::Log::s_verb="should be";
const char* Logics::Log::p_verb="should be";
const char* Logics::Log::s_n_verb="should not be";
const char* Logics::Log::p_n_verb="should not be";
 */
const char* Logics::Log::s_verb="is";
const char* Logics::Log::p_verb="are";
const char* Logics::Log::s_n_verb="are not";
const char* Logics::Log::p_n_verb="should not be";
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
void Logics::Log::finish_sntnc(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    if(L==NULL)
        return;
    Var::append(err_msg,"%s ",L);
    
    if(v_flag==VERB)
    {
        if(__is__ && p_flag)
            Var::append(err_msg,"%s ",p_verb);
        else if(__is__ && !p_flag)
            Var::append(err_msg,"%s ",s_verb);
    
        else if(!__is__ && p_flag)
            Var::append(err_msg,"%s ",p_n_verb);
        else
            Var::append(err_msg,"%s ",s_n_verb);
    }
    else
    {
        if(__is__)
            Var::append(err_msg,"%s %s ",verb_aux,verb);
        else
            Var::append(err_msg,"%s %s %s ",verb_aux,negate,verb);
    }
    Var::append(err_msg,"%s ",op);
    
    if(R)
        Var::append(err_msg,"%s",R);

    
    
    delete [] L;
    delete [] op;
    delete [] R;
    
    L=NULL;
    op=NULL;
    R=NULL;
    p_flag=false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::Log::start_sntnc(char*& L,char*& op,char*& R,int v_flag,const char* l_name,const char* r_name)
{
    delete [] L;
    L=NULL;
    Var::append(L,"%s",l_name);
    
    delete [] op;
    op=NULL;
    Var::append(op,"%s",op_name);
    delete [] R;
    R=NULL;
    if(r_name)
        Var::append(R,"%s",r_name);
}
/*--------------------------------------------
 3: (!L ,  O ,  R ) -> ( L0+L1 , O        , R     )
 4: ( L , !O , !R ) -> ( L     , O0+R0+O1 , R1    )
 5: ( L , !O ,  R ) -> ( L     , O0+O1    , R     )
 6: ( L ,  O , !R ) -> ( L     , O        , R0+R1 )
 
 2: (!L ,  O , !R ) ->   x
 1: (!L , !O ,  R ) ->   x
 0: (!L , !O , !R ) ->   x
 --------------------------------------------*/
void Logics::Log::process_sntnc(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,const char* l_name,const char* r_name)
{
    if(!L)
    {
        start_sntnc(L,op,R,v_flag,l_name,r_name);
        __is__=_is_;
        p_flag=false;
        return;
    }
    
    int case_no=4*eq(L,l_name)+eq(R,r_name);
    if(__is__==_is_) case_no+=2*eq(op,op_name);
    
    if(case_no<3 || case_no==7)
    {
        finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
        Var::append(err_msg,", %s ",del);
        print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
        return;
    }
 
    /*--------------------------------------------------
     3: (!L ,  O ,  R ) -> ( L0+L1 , O        , R     )
     --------------------------------------------------*/
    if(case_no==3)
    {
        p_flag=true;
        Var::append(L," %s %s",del,l_name);
    }
    /*--------------------------------------------------
     4: ( L , !O , !R ) -> ( L     , O0+R0+O1 , R1    )
     4: ( L , !O , !R ) -> ( L     , O1+R1+O0 , R0    )
     --------------------------------------------------*/
    else if(case_no==4)
    {
        if(_is_!=__is__)
        {
            if(!__is__)
            {
                char* op_=NULL;
                if(R)
                    Var::append(op_,"%s %s %s %s %s",op_name,r_name,del,negate,op);
                else
                    Var::append(op_,"%s %s %s %s",op_name,del,negate,op);
                delete [] op;
                op=op_;
                __is__=true;
                
            }
            else
            {
                
                if(R)
                    Var::append(op," %s %s %s %s",R,del,negate,op_name);
                else
                    Var::append(op," %s %s %s",del,negate,op_name);
                
                delete [] R;
                R=NULL;
                if(r_name)
                    Var::append(R,"%s",r_name);
            }
            
            
        }
        else
        {
            if(R)
                Var::append(op," %s %s %s",R,del,op_name);
            else
                Var::append(op," %s %s",del,op_name);
            delete [] R;
            R=NULL;
            if(r_name)
                Var::append(R,"%s",r_name);
        }
    }
    /*--------------------------------------------------
     5: ( L , !O ,  R ) -> ( L     , O0+O1    , R     )
     5: ( L , !O ,  R ) -> ( L     , O1+O0    , R     )
     --------------------------------------------------*/
    else if(case_no==5)
    {
        
        if(_is_!=__is__)
        {
            if(!__is__)
            {
                char* op_=NULL;
                Var::append(op_,"%s %s %s %s",op_name,del,negate,op);
                delete [] op;
                op=op_;
                __is__=true;
            }
            else
                Var::append(op," %s %s %s",del,negate,op_name);
        }
        else
            Var::append(op," %s %s",del,op_name);
    }
    /*--------------------------------------------------
     6: ( L ,  O , !R ) -> ( L     , O        , R0+R1 )
     --------------------------------------------------*/
    else if(case_no==6)
    {
        if(r_name)
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
Logics::LogBinaryOp::~LogBinaryOp()
{
    delete left;
    delete right;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogBinaryOp::print(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    del=op_name;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogBinaryOp::print(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,Var* v)
{
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
    del=op_name;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
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
void Logics::LogBinary::print(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    process_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag,left->name,right->name);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogBinary::print(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,Var* v)
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
    process_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag,left->name,right->name);
    
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
void Logics::LogUnary::print(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    process_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag,left->name,NULL);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::LogUnary::print(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,Var* v)
{
    int i=0;
    if(left==NULL)
    {
        left=v;
        i++;
    }
    process_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag,left->name,NULL);
    
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
void Logics::LogIf::print(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    if(err_msg)
        Var::append(err_msg,". ");
    
    v_flag=AUX_VERB;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    
    Var::append(err_msg,", %s ","when");
    
    v_flag=VERB;
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
}
void Logics::LogIf::print(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,Var* v)
{
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    if(err_msg)
        Var::append(err_msg,". ");
    
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    if(err_msg)
        Var::append(err_msg,". ");
    
    v_flag=AUX_VERB;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
    
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    
    Var::append(err_msg,", %s ","when");
    
    v_flag=VERB;
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
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
void Logics::LogIff::print(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag)
{
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    if(err_msg)
        Var::append(err_msg,". ");

    v_flag=AUX_VERB;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    
    Var::append(err_msg,", %s ","if and only if");

    v_flag=VERB;
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag);
}
void Logics::LogIff::print(char*& err_msg,char*& L,char*& op,char*& R,bool& __is__,char*& del,int v_flag,bool& p_flag,Var* v)
{
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    if(err_msg)
        Var::append(err_msg,". ");
    
    v_flag=AUX_VERB;
    right->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
    
    finish_sntnc(err_msg,L,op,R,__is__,del,v_flag,p_flag);
    
    Var::append(err_msg,", %s ","if and only if");
    
    v_flag=VERB;
    left->print(err_msg,L,op,R,__is__,del,v_flag,p_flag,v);
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
void Logics::assign(Var* l,const char* op_name,Var* r)
{
    creat_op(l,op_name,r);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Logics::assign(Var* l,const char* op_name)
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
    if(other.op)
        this->op=other.op->clone();
    else
        this->op=NULL;
    
    if(other.next)
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
    if(other.op)
        this->op=other.op->clone();
    else
        this->op=NULL;
    
    if(other.next)
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
    
    op->print(err_msg,l,oper,r,__is__,del,v_flag,p_flag,v);
    op->finish_sntnc(err_msg,l,oper,r,__is__,del,v_flag,p_flag);
    Var::append(err_msg,".");
    
    delete [] l;
    delete [] oper;
    delete [] r;
    delete [] _verb_;
    
    if(next)
    {
        char* next_msg=next->print(v);
        if(next_msg)
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
    
    if(next)
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
