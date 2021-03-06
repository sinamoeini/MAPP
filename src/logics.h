#ifndef __MAPP__logics__
#define __MAPP__logics__
/*--------------------------------------------
 
 --------------------------------------------*/
#include "error.h"
#include "global.h"
#include "var.h"
#include <stdio.h>
#include <string.h>
enum {VERB,AUX_VERB};
namespace MAPP_NS
{
 
    class Logics
    {
    private:
    protected:
        class Log
        {
        private:
        protected:
        public:
            static const char* verb_aux;
            static const char* negate;
            static const char* verb;
            static const char* s_verb;
            static const char* p_verb;
            static const char* s_n_verb;
            static const char* p_n_verb;
            char* op_name;
            const bool _is_;
            Log(const char*,const bool);
            virtual ~Log();
            void finish_sntnc(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void start_sntnc(char*&,char*&,char*&,int,const char*,const char*);
            void process_sntnc(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,const char*,const char*);
            virtual bool operator() ()=0;
            virtual bool operator() (Var*)=0;
            virtual void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&)=0;
            virtual void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*)=0;
            int eq(const char*,const char*);
            virtual Log* clone()=0;
        };
        
        class LogUnary:public Log
        {
        public:
            Var* left;
            LogUnary(Var*,const char*,const bool);
            
            void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*);
            virtual Log* clone()=0;
            virtual bool operator() ()=0;
            virtual bool operator() (Var*)=0;
        };
        
        class LogBinary:public Log
        {
        public:
            Var* left;
            Var* right;
            LogBinary(Var*,Var*,const char*,const bool);
            
            void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*);
            virtual Log* clone()=0;
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
            
            void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*);
            virtual Log* clone()=0;
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
            
            virtual void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&)=0;
            virtual void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*)=0;
            virtual Log* clone()=0;
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
            void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*);
            Log* clone();
        };
        
        class LogIff:public LogCond
        {
        public:
            LogIff(Log*,Log*);
            bool operator() ();
            bool operator() (Var*);
            void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&);
            void print(char*&,char*&,char*&,char*&,bool&,char*&,int,bool&,Var*);
            Log* clone();
        };
        
        void creat_op(Var*,const char*,Var*);
        void creat_op(Var*,const char*);
        void assign(Var*,const char*,Var*);
        void assign(Var*,const char*);
    public:
        Logics();
        
        template<typename T0,typename T1>
        Logics(T0&& x0,const char* s,T1&& x1)
        {
            assign(
                   MAPP_NS::g_vm->opppp(dynamic_static<T0&&>(),x0),
                   s,
                   MAPP_NS::g_vm->opppp(dynamic_static<T1&&>(),x1));
        }
        
        template<typename T0>
        Logics(T0&& x0,const char* s)
        {
            assign(MAPP_NS::g_vm->opppp(dynamic_static<T0&&>(),x0),s);
        }
        
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

    class VLogics: public Logics
    {
    private:
    protected:
    public:
        template<typename T0>
        VLogics(const char* s,T0&& x1):
        Logics()
        {
            assign(NULL,s,MAPP_NS::g_vm->opppp(dynamic_static<T0&&>(),x1));
        }
        
        VLogics(const char* s):
        Logics()
        {
            assign(NULL,s);
        }
    };
   
    
   
}
#endif
