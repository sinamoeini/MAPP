#ifndef __MAPP__var__
#define __MAPP__var__
/*--------------------------------------------
 
 --------------------------------------------*/
#include "error.h"
#include "global.h"
#include <stdio.h>
#include <string.h>
#include <type_traits>
#include <mpi.h>
enum {VERB,AUX_VERB};
namespace MAPP_NS
{
    template<class T>
    struct approp { typedef T type;};
    template<class T>
    struct approp<const T> { typedef typename approp<T>::type type; };
    template<class T>
    struct approp<T[]> { typedef typename approp<T>::type* type; };
    template<class T>
    struct approp<const T[]> { typedef typename approp<T[]>::type type; };
    template<class T, std::size_t N>
    struct approp<T[N]> { typedef typename approp<T>::type* type; };
    template<class T, std::size_t N>
    struct approp<const T[N]> { typedef typename approp<T[N]>::type type; };
    template<class T>
    struct approp<T*> { typedef typename approp<T>::type* type; };
    template<class T>
    struct approp<const T*> { typedef typename approp<T*>::type type; };
    template<class T>
    struct approp<T&> { typedef typename approp<T>::type type; };
    template<class T>
    struct approp<T&&> { typedef typename approp<T>::type type; };
    
    template<class T>
    struct approp_k {typedef T type;};
    template<class T>
    struct approp_k<T[]> {typedef typename approp_k<T>::type* type;};
    template<class T, std::size_t N>
    struct approp_k<T[N]> {typedef typename approp_k<T>::type* type;};
    template<class T>
    struct approp_k<T*> {typedef typename approp_k<T>::type* type;};
    template<class T>
    struct approp_k<T&> {typedef typename approp_k<T>::type type;};
    template<class T>
    struct approp_k<T&&> {typedef typename approp_k<T>::type type;};
    


    class FindReplace
    {
        void attach(char*&,size_t&,int);
        void attach(char*&,size_t&,char);
        char** srch;
        size_t* srch_len;
        char** rplc;
        size_t* rplc_len;
        
        int n;
        bool is_same(char*,char*);
        int look_up(char*);
    public:
        FindReplace();
        ~FindReplace();
        void reverse();
        void operator()(char*&);
        void operator()(const char*,const char*);
    };
    
    
    class Var;
    
    class VarManager
    {
    private:
        VarManager* old_vm;
    public:
        Var ** st_vars;
        Var ** dy_vars;
        int nst_vars;
        int ndy_vars;
        VarManager();
        ~VarManager();

        template<typename T0>
        Var* add_var(T0&,const char*);

        template <typename T0>
        Var* operator()(T0&&);
        template<typename T0>
        Var* operator()(bool,T0&);
    };
    
    class Logics
    {
    private:
    protected:
        class Log
        {
        private:
        protected:
        public:
            char* op_name;
            const bool _is_;
            Log(const char*,const bool);
            virtual ~Log();
            void finish_sntnc(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void start_sntnc(char*&,char*&,char*&,char*&,int,char*,char*);
            void process_sntnc(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,char*,char*);
            virtual bool operator() ()=0;
            virtual bool operator() (Var*)=0;
            virtual void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&)=0;
            virtual void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*)=0;
            int eq(const char*,const char*);
            virtual Log* clone()=0;
            virtual void replace(Var*,Var*)=0;
        };
        
        class LogUnary:public Log
        {
        public:
            Var* left;
            LogUnary(Var*,const char*,const bool);
            
            void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*);
            virtual Log* clone()=0;
            void replace(Var*,Var*);
            virtual bool operator() ()=0;
            virtual bool operator() (Var*)=0;
        };
        
        class LogBinary:public Log
        {
        public:
            Var* left;
            Var* right;
            LogBinary(Var*,Var*,const char*,const bool);
            
            void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*);
            virtual Log* clone()=0;
            void replace(Var*,Var*);
            virtual bool operator() ()=0;
            virtual bool operator() (Var*)=0;
        };
        
        class LogBinaryOp:public Log
        {
        public:
            Log* left;
            Log* right;
            LogBinaryOp(Log*,Log*,const char*,const bool);
            ~LogBinaryOp();
            
            void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*);
            virtual Log* clone()=0;
            void replace(Var*,Var*);
            virtual bool operator() ()=0;
            virtual bool operator() (Var*)=0;
        };
        
        class LogCond:public Log
        {
        public:
            Log* left;
            Log* right;
            LogCond(Log*,Log*,const char*,const bool);
            ~LogCond();
            
            virtual void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&)=0;
            virtual void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*)=0;
            virtual Log* clone()=0;
            void replace(Var*,Var*);
            virtual bool operator() ()=0;
            virtual bool operator() (Var*)=0;
        };
        
        
        class LogIS:public LogUnary
        {
        public:
            LogIS(Var*,const bool);
            bool operator() ();
            bool operator() (Var*);
            Log* clone();
        };
        
        class LogGT:public LogBinary
        {
        public:
            LogGT(Var*,Var*,const bool);
            bool operator() ();
            bool operator() (Var*);
            Log* clone();
        };
        
        class LogGE:public LogBinary
        {
        public:
            LogGE(Var*,Var*,const bool);
            bool operator() ();
            bool operator() (Var*);
            Log* clone();
        };
        
        class LogEq:public LogBinary
        {
        public:
            LogEq(Var*,Var*,const bool);
            bool operator() ();
            bool operator() (Var*);
            Log* clone();
        };
        
        class LogLT:public LogBinary
        {
        public:
            LogLT(Var*,Var*,const bool);
            bool operator() ();
            bool operator() (Var*);
            Log* clone();
        };
        
        class LogLE:public LogBinary
        {
        public:
            LogLE(Var*,Var*,const bool);
            bool operator() ();
            bool operator() (Var*);
            Log* clone();
        };
        
        class LogAnd:public LogBinaryOp
        {
        public:
            LogAnd(Log*,Log*);
            bool operator() ();
            bool operator() (Var*);            
            Log* clone();
        };
        
        class LogOr:public LogBinaryOp
        {
        public:
            LogOr(Log*,Log*);
            bool operator() ();
            bool operator() (Var*);
            Log* clone();
        };
        
        class LogIf:public LogCond
        {
        public:
            LogIf(Log*,Log*);
            bool operator() ();
            bool operator() (Var*);
            void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*);
            Log* clone();
        };
        
        class LogIff:public LogCond
        {
        public:
            LogIff(Log*,Log*);
            bool operator() ();
            bool operator() (Var*);
            void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void print(char*&,char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*);
            Log* clone();
        };
        
        void creat_op(Var*,const char*,Var*);
        void creat_op(Var*,const char*);
    public:
        Logics();
        Logics(Var*,const char*,Var*);
        Logics(Var*,const char*);
        Logics(Logics&);
        Logics(Logics&&);
        ~Logics();
        
        
        Logics& operator = (Logics&);
        Logics& operator = (Logics&&);
        Logics&& operator + (Logics&&);
        Logics&& operator * (Logics&&);
        Logics&& operator - (Logics&&);
        Logics&& operator / (Logics&&);
        Logics& operator += (Logics&&);
        Logics& operator *= (Logics&&);
        Logics& operator -= (Logics&&);
        Logics& operator /= (Logics&&);
        
        Log* op;
        Logics* next;
        
        Logics& link();
        
        char* print(Var*);
        char* operator()(Var*);
        char* operator()();

        
    };
    
    template<typename T>
    bool dynamic_static()
    {
        if(std::is_rvalue_reference<T>::value)
            return false;
        if(std::is_lvalue_reference<T>::value)
        {
            typedef typename std::remove_reference<T>::type TT;
            if(std::is_const<TT>::value)
                return false;
        }
        
        return true;
    }
    
    template<typename T>
    bool dynamic_statc()
    {
        if(std::is_rvalue_reference<T>::value)
            return false;
        if(std::is_lvalue_reference<T>::value)
        {
            typedef typename std::remove_reference<T>::type TT;
            if(std::is_const<TT>::value)
                return false;
        }
        
        return true;
    }
    
    template<typename T>
    Var* get_var(T&& v)
    {
        return MAPP_NS::g_vm->operator()(dynamic_static<T&&>(),v);
    }
    
    
    template<typename T>
    void add_var(T&& v,const char* n)
    {
        MAPP_NS::g_vm->add_var(v,n);
    }
    
    template<typename T0,typename T1>
    Logics logic(T0&& x0,const char* s,T1&& x1)
    {
        return Logics(MAPP_NS::g_vm->operator()(dynamic_static<T0&&>(),x0)
        ,s,MAPP_NS::g_vm->operator()(dynamic_static<T1&&>(),x1));
    }
    
    template<typename T0>
    Logics logic(T0&& x0,const char* s)
    {
        return Logics(MAPP_NS::g_vm->operator()(dynamic_static<T0&&>(),x0),s);
    }
    
    template<typename T0>
    Logics vlogic(const char* s,T0&& x1)
    {
        return Logics(NULL,s,MAPP_NS::g_vm->operator()(dynamic_static<T0&&>(),x1));
    }
    
    inline Logics vlogic(const char* s)
    {
        return Logics(NULL,s);
    }
    
    class Var
    {
    protected:
        char* loc;
    public:
        const void* addr;
        bool set;
        char* name;        
        
        Var(const char*,void*);
        Var();
        virtual ~Var();
        
        virtual char* scan(char*);
        virtual bool operator < (Var&);
        virtual bool operator > (Var&);
        virtual bool operator <= (Var&);
        virtual bool operator >= (Var&);
        virtual bool operator == (Var&);
        virtual Var* clone();

        void det_loc();
        void att_loc();
        virtual const char* get_type_name();
        
        template <typename... T0>
        inline static void append(char*& buff,const char* format,T0... vals)
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
        
        static inline void append(char*& buff,const char* format)
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
    
    template<typename T0>
    class var: public Var
    {
    private:
        T0 st_val;
        char* sub_scan(char*);
        T0& get_val(Var&);
    protected:
    public:
        T0& val;
        var(T0&,const char*);
        var(const T0);
        
        char* scan(char*);
        Var* clone();
        bool operator < (Var&);
        bool operator > (Var&);
        bool operator <= (Var&);
        bool operator >= (Var&);
        bool operator == (Var&);
        const char* get_type_name(){return "";}
    };
    
    template<typename T0>
    class var<T0*>: public Var
    {
    private:
        T0* st_val;
        int size;
        char* sub_scan(char*);
        T0*& get_val(Var&,int&);
    protected:
    public:
        T0*& val;
        var(T0*&,const char*);
        var(const T0*);
        template<int N>
        var(T0(&)[N]);
        ~var();
        
        char* scan(char*);
        Var* clone();
        bool operator < (Var&);
        bool operator > (Var&);
        bool operator <= (Var&);
        bool operator >= (Var&);
        bool operator == (Var&);
        const char* get_type_name(){return "";}
    };
    
    template<> var<int>::var(const int);
    template<> char* var<int>::sub_scan(char*);
    template<> const char* var<int>::get_type_name();
    
    template<> var<double>::var(const double);
    template<> char* var<double>::sub_scan(char*);
    template<> const char* var<double>::get_type_name();
    
    template<> var<char*>::var(const char*);
    template<> char* var<char*>::sub_scan(char*);
    template<> bool var<char*>::operator < (Var&);
    template<> bool var<char*>::operator > (Var&);
    template<> bool var<char*>::operator <= (Var&);
    template<> bool var<char*>::operator >= (Var&);
    template<> bool var<char*>::operator == (Var&);
    template<> const char* var<char*>::get_type_name();
    
    
}
using namespace MAPP_NS;
/*-----------------------------------------------
  _     _       ___   ____       __  _____  __
 | |   / /     /   | |  _  \    / / |_   _| \ \
 | |  / /     / /| | | |_| |   / /    | |    \ \
 | | / /     / / | | |  _  /  { {     | |     } }
 | |/ /     / /  | | | | \ \   \ \    | |    / /
 |___/     /_/   |_| |_|  \_\   \_\   |_|   /_/
 -----------------------------------------------*/
template<typename T0>
var<T0>::var(T0& val_,const char* name):
Var(name,&val_),
val(val_)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
var<T0>::var(const T0 val_):
Var(),
val(st_val)
{}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
T0& var<T0>::get_val(Var& other)
{
    var<T0>* p=dynamic_cast<var<T0>*>(&other);
    if(p==NULL)
    {
        printf("failed conversion for variable (%s)\n",other.name);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    return p->val;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool var<T0>::operator < (Var& right)
{
    return (val < get_val(right));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool var<T0>::operator > (Var& right)
{
    return (val > get_val(right));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool var<T0>::operator <= (Var& right)
{
    return (val <= get_val(right));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool var<T0>::operator >= (Var& right)
{
    return (val >= get_val(right));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool var<T0>::operator == (Var& right)
{
    return (val == get_val(right));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
char* var<T0>::scan(char* inp)
{
    char* err_msg=sub_scan(inp);
    if(err_msg!=NULL)
        return err_msg;
    
    set=true;
    return  NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
char* var<T0>::sub_scan(char*)
{
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
Var* var<T0>::clone()
{
    return new var<T0>(val,name);
}
/*-------------------------------------------------------
  _     _   ___   ____       __  _____  __   _   __ __
 | |   / / /   | |  _  \    / / |_   _| \ \ | | / / \ \
 | |  / / / /| | | |_| |   / /    | |    \ \| |/ /   \ \
 | | / / / / | | |  _  /  { {     | |     }     {     } }
 | |/ / / /  | | | | \ \   \ \    | |    / /| |\ \   / /
 |___/ /_/   |_| |_|  \_\   \_\   |_|   /_/ |_| \_\ /_/
 -------------------------------------------------------*/
template<typename T0>
var<T0*>::var(const T0*):
Var(),
val(st_val)
{
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0> template<int N>
var<T0*>::var(T0 (&v)[N]):
Var(),
val(st_val)
{
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<> template<int N>
var<char*>::var(char (&v)[N]):
Var(),
val(st_val)
{
    size=N;
    st_val=new char[size];
    name=new char[size];
    memcpy(name,v,size);
    memcpy(st_val,v,size);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
var<T0*>::var(T0*& v,const char* name_):
Var(name_,&v),
val(v)
{
    val=NULL;
    st_val=NULL;
    size=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
var<T0*>::~var()
{
    delete [] st_val;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
T0*& var<T0*>::get_val(Var& other,int& other_size)
{
    var<T0*>* p=dynamic_cast<var<T0*>*>(&other);

    if(p==NULL)
    {
        printf("failed conversion for variable (%s)\n",other.name);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    other_size=p->size;
    return p->val;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool var<T0*>::operator < (Var& right)
{
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool var<T0*>::operator > (Var& right)
{
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool var<T0*>::operator <= (Var& right)
{
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool var<T0*>::operator >= (Var& right)
{
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool var<T0*>::operator == (Var& right)
{
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
char* var<T0*>::scan(char* inp)
{
    char* err_msg=sub_scan(inp);
    if(err_msg!=NULL)
        return err_msg;
    
    set=true;
    return  NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
Var* var<T0*>::clone()
{
    if(&val==&st_val)
    {
        return new var<T0*>(st_val);
    }
    else
    {
        return new var<T0*>(val,name);
    }
    
}
/*--------------------------------------------------------------------------------------------
 _     _   ___   _____        ___  ___       ___   __   _       ___   _____   _____   _____
| |   / / /   | |  _  \      /   |/   |     /   | |  \ | |     /   | /  ___| | ____| |  _  \
| |  / / / /| | | |_| |     / /|   /| |    / /| | |   \| |    / /| | | |     | |__   | |_| |
| | / / / / | | |  _  /    / / |__/ | |   / / | | | |\   |   / / | | | |  _  |  __|  |  _  /
| |/ / / /  | | | | \ \   / /       | |  / /  | | | | \  |  / /  | | | |_| | | |___  | | \ \
|___/ /_/   |_| |_|  \_\ /_/        |_| /_/   |_| |_|  \_| /_/   |_| \_____/ |_____| |_|  \_\
 --------------------------------------------------------------------------------------------*/
template<typename T0>
Var* VarManager::add_var(T0& v,const char* var_name)
{
    var<T0>* var_=new var<T0>(v,var_name);
    Var** vars_=new Var*[ndy_vars+1];
    memcpy(vars_,dy_vars,ndy_vars*sizeof(Var*));
    delete [] dy_vars;
    dy_vars=vars_;
    dy_vars[ndy_vars]=var_;
    ndy_vars++;
    return var_;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <typename T0>
Var* VarManager::operator()(T0&& v)
{
    return operator()(dynamic_statc<T0&&>(),v);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
Var* VarManager::operator()(bool ds,T0& v)
{
    if(ds)
    {
        const void* addr=&v;
        for(int i=0;i<ndy_vars;i++)
            if(dy_vars[i]->addr==addr)
                return dy_vars[i];
    }
    else
    {
        typedef typename approp<T0>::type T1;
        
        typename approp_k<T0>::type vv=v;
        
        var<T1>* var_=new var<T1>(vv);
        
        for(int i=0;i<nst_vars;i++)
        {
            var<T1>* _var_=dynamic_cast<var<T1>*>(st_vars[i]);
            if(_var_!=NULL && *_var_==*var_)
            {
                delete var_;
                return st_vars[i];
            }
        }
        
        Var** vars_=new Var*[nst_vars+1];
        memcpy(vars_,st_vars,nst_vars*sizeof(Var*));
        delete [] st_vars;
        st_vars=vars_;
        st_vars[nst_vars]=var_;
        nst_vars++;
        return st_vars[nst_vars-1];
    }
    
    
    
    printf("error in %s at line %d\n",__PRETTY_FUNCTION__,__LINE__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
    return NULL;
}




#endif
