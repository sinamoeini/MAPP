#ifndef __MAPP__cmd__
#define __MAPP__cmd__
/*--------------------------------------------
 
 --------------------------------------------*/
#include "error.h"
#include "var.h"
#include <stdio.h>
#include <string.h>
namespace MAPP_NS
{
    class PrintStyle
    {
    private:
        void add(char* (&)[2],const char*,const char*);
    public:
        PrintStyle();
        ~PrintStyle();
        void var(const char*,const char*);
        void type(const char*,const char*);
        void cmd(const char*,const char*);
        void cmd_prop(const char*,const char*);
        void all_vars_prop(const char*,const char*);
        void all_cmd_prop(const char*,const char*);
        void vars_prop(const char*,const char*);
        void pattern_prop(const char*,const char*);
        
        char* var_brc[2];
        char* type_brc[2];
        char* cmd_brc[2];
        char* cmd_prop_brc[2];
        char* all_vars_prop_brc[2];
        char* all_cmd_prop_brc[2];
        char* vars_prop_brc[2];
        char* pattern_prop_brc[2];
    };
    
    class FileReader:protected InitPtrs
    {
    private:
        class QuantityReader
        {
        private:
            char* curr;
            bool is_del(char*);
            bool is_space(char*);
            bool is_null(char*);
            bool is_char(char*);
            bool syntax(char*);
            
            const char* dels;
            const size_t ndels;
            const char* spaces;
            const size_t nspaces;
            const char* nulls;
            const size_t nnulls;
            
            char* first;
            char* lp;
            char* rp;
            char* com;
            char* eq;
        public:
            QuantityReader();
            char* operator()(char*,char**&,int&);
        };
        
        class Ten
        {
        private:
        protected:
        public:
            Logics var_logics;
            Logics cmd_logic;
            Error*& error;
            VarManager& varmngr;
            FileReader* fr;
            char* name;
            char* format;
            int nkywrds;
            FileKeyWord* head;
            bool symmetric;
            Ten(FileReader*,const char*);
            virtual ~Ten();
            virtual bool scan(char*&,char**,int)=0;
            virtual void set_max_cmp(int&)=0;
            Logics&  add_vlog();
            Logics&  add_clog();
        };
        template<typename T0>
        class Ten_0D:public Ten
        {
        private:
            T0& v;
            FileKeyWord kywrds;
        protected:
        public:
            Ten_0D(FileReader*,const char*,T0&);
            ~Ten_0D();
            bool scan(char*&,char**,int);
            void set_max_cmp(int&);
        };
        template<typename T0>
        class Ten_1D:public Ten
        {
        private:
            int& dim0;
            T0*& v;
            FileKeyWord* kywrds;
            Logics vli0;
            Var* vi0;
            int i0;
        protected:
        public:
            Ten_1D(FileReader*,const char*,T0*&,int&);
            ~Ten_1D();
            bool scan(char*&,char**,int);
            void set_max_cmp(int&);
        };
        template<typename T0>
        class Ten_2D:public Ten
        {
        private:
            int& dim0;
            int& dim1;
            T0**& v;
            FileKeyWord** kywrds;
            Logics vli0;
            Var* vi0;
            int i0;
            Logics vli1;
            Var* vi1;
            int i1;
        protected:
        public:
            Ten_2D(FileReader*,const char*,T0**&,int&,int&);
            ~Ten_2D();
            bool scan(char*&,char**,int);
            void set_max_cmp(int&);
        };
        
        
        int& ntypes;
        Ten** tens;
        int ntens;
        FILE* fp;
        
        char* file_name;
        int* type_ref;
        int ntype_ref;
        
        void open_file();
        void close_file();
        void read_header();
        void read_body();
        int parse_line(char*,char**&,int&);
        int hash_remover(char*&);
        int read_line(char*&);
        int read_line(char*&,int&,int&);
        
        Var* curr_var;
        QuantityReader qr;
    protected:
    public:
        VarManager varmngr;
        FileReader(MAPP*);
        ~FileReader();

        void read_file(const char*);
        
        template<typename T0>
        void add_scl(const char*,T0&);
        template<typename T0>
        void add_1D(const char*,T0*&);
        template<typename T0>
        void add_2D(const char*,T0**&);
        
        
        template<typename T0>
        Logics log_u(T0&,const char*);
        template<typename T0>
        Logics log_u(T0*&,const char*);
        

        void symmetric();
        Logics&  add_vlog(){return tens[ntens-1]->add_vlog();}
        Logics&  add_clog(){return tens[ntens-1]->add_clog();}

    };
    
    
    class Pattern
    {
    private:
    protected:
    public:
        class SubPattern;
        template<typename T0>
        class SubPattern_DY;
        class SubPattern_0D;
        class SubPattern_1D;
        class SubPattern_2D;
        class SubPattern_Voigt;
        class SubPattern_2D_LT;
        
        Error*& error;
        SubPattern** sub_ptrns;
        int nsub_ptrns;
        
        void add_sp(SubPattern*);
        
        Var* curr_var;
        SubPattern* curr_sp;
        SubPattern* dy_sp;
        int nlog_tree;
        

        VarManager varmngr;
        Pattern(Error*&);
        ~Pattern();
        
       
        template<typename T0>
        void cmd_dy(const char*,T0*&,int&);
        
        void cmd(const char*);
        void cmd(bool&,const char*);
        
        void cmd_voigt(const char*,int);
        void cmd_voigt(bool*&,const char*,int);
        
        void cmd_1d(const char*,const int);
        void cmd_1d(bool*&,const char*,const int);
        
        void cmd_2d_lt(const char*,const int);
        void cmd_2d_lt(bool**&,const char*,const int);
        
        void cmd_2d(const char*,const int,const int);
        void cmd_2d(bool**&,const char*,const int,const int);
        
        bool scan(char**&,int&);
        
        template<typename T0,typename... T1>
        void add_var(T0&,T1&...);
        
        template<typename T0,typename... T1>
        void add_var(T0&,const char*,T1&...);

        template<typename T0>
        void add_var(T0&);
        
        template<typename T0>
        void add_var(T0&,const char*);
        

        void add_vdesc(int,const char*);
        void add_cdesc(const char*);
        
        Logics& add_vlog(int);
        Logics& add_clog();
        void print_info();       
        
    };
    
    class Pattern::SubPattern
    {
    private:
    protected:
        Error*& error;
        VarManager& var_mngr;

        char* vlog_chek(KeyWord&,char**&,int&);
        
        KeyWord* head_kywrd;
        int nkywrds;
        void var_adj();
        int rank;
    public:
        
        Logics cmd_logic;
        Logics* var_logics;
        char* cmd_desc;
        char** var_descs;
        char** var_types;
        Var** vars;
        int nvars;
        
        char* name;
        SubPattern(VarManager&,const char*);
        virtual ~SubPattern();
        virtual bool scan(char*&,char**& args,int& nargs)=0;
        
        template<typename T0>
        void add_var(T0&,const char*);
        template<typename T0>
        void add_var(T0*&,const char*);
        template<typename T0>
        void add_var(T0**&,const char*);
        template<typename T0,int N>
        void add_var(T0 (&x)[N],const char*);
        
        virtual char* finalize();
        
        void add_vdesc(int,const char*);
        void add_cdesc(const char*);
        
        Logics& add_vlog(int);
        Logics& add_clog();
        
        virtual void print_pattern(PrintStyle&,char*&);
        virtual void print_vars_prop(PrintStyle&,char*&,int);
        virtual void print_vars_prop(PrintStyle&,char*&);
        virtual void print_cmd_prop(PrintStyle&,char*&);
        virtual void conv_name(PrintStyle&)=0;
        void deconv_name();
    };
    

    
    class Pattern::SubPattern_0D:public Pattern::SubPattern
    {
    private:
        KeyWord kywrd;
        void init();
        void init_keywords();
        void init_keywords(bool&);
    protected:
    public:
        SubPattern_0D(VarManager&,const char*);
        SubPattern_0D(bool&,VarManager&,const char*);
        ~SubPattern_0D();
        bool scan(char*&,char**&,int&);
        
        template<typename T0>
        void add_var(T0&,const char*);
        
        void conv_name(PrintStyle&);
    };
    
    class Pattern::SubPattern_1D:public Pattern::SubPattern
    {
    private:
        int icomp0;
        Var& comp0;
        KeyWord* kywrds;
        void init();
        void init_keywords();
        void init_keywords(bool*&);
    protected:
    public:
        char* format;
        int dim;
        SubPattern_1D(VarManager&,const char*,int);
        SubPattern_1D(bool*&,VarManager&,const char*,int);
        ~SubPattern_1D();
        bool scan(char*&,char**&,int&);
        
        template<typename T0>
        void add_var(T0*&,const char*);
        
        void conv_name(PrintStyle&);
    };
    
    class Pattern::SubPattern_Voigt:public Pattern::SubPattern
    {
    private:
        int icomp0,icomp1;
        Var& comp0;
        Var& comp1;
        KeyWord* kywrds;
        void cmpt2voigt(int&,const int,const int);
        void voigt2cmpt(const int,int&,int&);
        void init();
        void init_keywords();
        void init_keywords(bool*&);
    protected:
    public:
        int dim;
        char* format;
        SubPattern_Voigt(VarManager&,const char*,int);
        SubPattern_Voigt(bool*&,VarManager&,const char*,int);
        ~SubPattern_Voigt();
        bool scan(char*&,char**&,int&);
        
        template<typename T0>
        void add_var(T0*&,const char*);
        
        void conv_name(PrintStyle&);
    };
    
    class Pattern::SubPattern_2D:public Pattern::SubPattern
    {
    private:
        int icomp0,icomp1;
        Var& comp0;
        Var& comp1;
        KeyWord** kywrds;
        void init();
        void init_keywords();
        void init_keywords(bool**&);
    protected:
    public:
        char* format;
        int dim0,dim1;
        SubPattern_2D(VarManager&,const char*,int,int);
        SubPattern_2D(bool**&,VarManager&,const char*,int,int);
        ~SubPattern_2D();
        bool scan(char*&,char**&,int&);
        
        template<typename T0>
        void add_var(T0**&,const char*);
        
        void conv_name(PrintStyle&);
    };
    
    class Pattern::SubPattern_2D_LT:public Pattern::SubPattern
    {
    private:
        int icomp0,icomp1;
        Var& comp0;
        Var& comp1;
        KeyWord** kywrds;
        void init();
        void init_keywords();
        void init_keywords(bool**&);
    protected:
    public:
        char* format;
        int dim;
        SubPattern_2D_LT(VarManager&,const char*,int);
        SubPattern_2D_LT(bool**&,VarManager&,const char*,int);
        ~SubPattern_2D_LT();
        bool scan(char*&,char**&,int&);
        
        template<typename T0>
        void add_var(T0**&,const char*);
        
        void conv_name(PrintStyle&);
    };
    
    
    template<typename T0>
    class Pattern::SubPattern_DY:public Pattern::SubPattern
    {
    private:
        T0 s_val;
        T0*& dy_vars;
        Var** v_dy_vars;
        int& ndy_vars;
    protected:
    public:
        SubPattern_DY(VarManager&,const char*,T0*&,int&);
        ~SubPattern_DY();
        bool scan(char*&,char**&,int&);
        
        void print_pattern(PrintStyle&,char*&);
        
        void conv_name(PrintStyle&){}
        char* finalize();
    };
}
/*--------------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        _____  __    __ 
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      |  _  \ \ \  / / 
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |      | | | |  \ \/ /  
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |      | | | |   \  /   
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |      | |_| |   / /    
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|      |_____/  /_/
 --------------------------------------------------------------------------------------------------------*/
template<typename T0>
Pattern::SubPattern_DY<T0>::SubPattern_DY(VarManager& var_mngr_,const char* name_,T0*& v,int& n):
SubPattern(var_mngr_,name_),
dy_vars(v),
ndy_vars(n)
{
    v_dy_vars=NULL;
    dy_vars=NULL;
    ndy_vars=0;
    char* _name=NULL;
    append(_name,"no. of %s",name_);
    var_mngr.add_var(n,_name);
    delete [] _name;
    
    var_adj();
    vars[0]=var_mngr.add_var(s_val,name);
    append(var_types[0],vars[0]->get_type_name());
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
Pattern::SubPattern_DY<T0>::~SubPattern_DY()
{
    delete [] v_dy_vars;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool Pattern::SubPattern_DY<T0>::scan(char*& err_msg,char**& args,int& nargs)
{
    
    T0* dy_vars_=new T0[ndy_vars+1];
    memcpy(dy_vars_,dy_vars,ndy_vars*sizeof(T0));
    delete [] dy_vars;
    dy_vars=dy_vars_;
    dy_vars[ndy_vars]=NULL;
    
    char* v_name=NULL;
    append(v_name,"%s_%d",name,ndy_vars);
    Var* var_=var_mngr.add_var(dy_vars[ndy_vars],v_name);
    delete [] v_name;
    
    Var** v_dy_vars_=new Var*[ndy_vars+1];
    memcpy(v_dy_vars_,v_dy_vars,ndy_vars*sizeof(Var*));
    delete [] v_dy_vars;
    v_dy_vars=v_dy_vars_;
    v_dy_vars[ndy_vars]=var_;
    ndy_vars++;
    
    err_msg=var_->scan(*args);
    if(err_msg!=NULL)
        return true;
    args++;
    nargs--;
    
    err_msg=var_logics[0](v_dy_vars[ndy_vars-1]);
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
char* Pattern::SubPattern_DY<T0>::finalize()
{return cmd_logic(NULL);}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void Pattern::SubPattern_DY<T0>::print_pattern(PrintStyle& ps,char*& buff)
{
    append(buff,ps.pattern_prop_brc[0]);
    append(buff,vars[0]->name);
    append(buff," ");
    append(buff,"...");
    append(buff,ps.pattern_prop_brc[1]);
}
/*---------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|
 ---------------------------------------------------------------------------------*/
template<typename T0>
void Pattern::SubPattern::add_var(T0& v,const char* n)
{
    
    if(dynamic_cast<SubPattern_0D*>(this)!=NULL)
        return dynamic_cast<SubPattern_0D*>(this)->add_var(v,n);
    
    printf("unsuccessful conversion at: %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0,int N>
void Pattern::SubPattern::add_var(T0 (&v)[N],const char* n)
{
    T0* vv=static_cast<T0*>(v);
    
    if(dynamic_cast<SubPattern_1D*>(this)!=NULL)
        return dynamic_cast<SubPattern_1D*>(this)->add_var(vv,n);
    if(dynamic_cast<SubPattern_Voigt*>(this)!=NULL)
        return dynamic_cast<SubPattern_Voigt*>(this)->add_var(vv,n);
    
    printf("unsuccessful conversion at: %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void Pattern::SubPattern::add_var(T0*& v,const char* n)
{
    if(dynamic_cast<SubPattern_0D*>(this)!=NULL)
        return dynamic_cast<SubPattern_0D*>(this)->add_var(v,n);
    if(dynamic_cast<SubPattern_1D*>(this)!=NULL)
        return dynamic_cast<SubPattern_1D*>(this)->add_var(v,n);
    if(dynamic_cast<SubPattern_Voigt*>(this)!=NULL)
        return dynamic_cast<SubPattern_Voigt*>(this)->add_var(v,n);
    
    printf("unsuccessful conversion at: %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void Pattern::SubPattern::add_var(T0**& v,const char* n)
{
    if(dynamic_cast<SubPattern_0D*>(this)!=NULL)
        return dynamic_cast<SubPattern_0D*>(this)->add_var(v,n);
    if(dynamic_cast<SubPattern_1D*>(this)!=NULL)
        return dynamic_cast<SubPattern_1D*>(this)->add_var(v,n);
    if(dynamic_cast<SubPattern_Voigt*>(this)!=NULL)
        return dynamic_cast<SubPattern_Voigt*>(this)->add_var(v,n);
    if(dynamic_cast<SubPattern_2D*>(this)!=NULL)
        return dynamic_cast<SubPattern_2D*>(this)->add_var(v,n);
    if(dynamic_cast<SubPattern_2D_LT*>(this)!=NULL)
        return dynamic_cast<SubPattern_2D_LT*>(this)->add_var(v,n);
    
    printf("unsuccessful conversion at: %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*------------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        _____   _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      /  _  \ |  _  \
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |      | | | | | | | |
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |      | |/| | | | | |
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |      | |_| | | |_| |
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|      \_____/ |_____/
 ------------------------------------------------------------------------------------------------------*/
template<typename T0>
void Pattern::SubPattern_0D::add_var(T0& v,const char* var_name_)
{
    kywrd.add_var(v,var_name_);
    var_adj();
    vars[nvars-1]=head_kywrd->vars[nvars-rank];
    append(var_types[nvars-1],vars[nvars-1]->get_type_name());
}
/*----------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        ___   _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      |_  | |  _  \
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |        | | | | | |
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |        | | | | | | 
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |        | | | |_| | 
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|        |_| |_____/
 ----------------------------------------------------------------------------------------------------*/
template<typename T0>
void Pattern::SubPattern_1D::add_var(T0*& v,const char* var_name_)
{
    for(int i=0;i<dim;i++)
        kywrds[i].add_var(v[i],var_name_);
    var_adj();
    vars[nvars-1]=head_kywrd->vars[nvars-rank];
    append(var_types[nvars-1],vars[nvars-1]->get_type_name());
}
/*----------------------------------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        _     _   _____   _   _____   _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      | |   / / /  _  \ | | /  ___| |_   _| 
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |      | |  / /  | | | | | | | |       | |   
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |      | | / /   | | | | | | | |  _    | |   
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |      | |/ /    | |_| | | | | |_| |   | |   
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|      |___/     \_____/ |_| \_____/   |_|
 ----------------------------------------------------------------------------------------------------------------------------*/
template<typename T0>
void Pattern::SubPattern_Voigt::add_var(T0*& v,const char* var_name_)
{
    for(int i=0;i<dim*(dim+1)/2;i++)
        kywrds[i].add_var(v[i],var_name_);
    var_adj();
    
    vars[nvars-1]=head_kywrd->vars[nvars-rank];
    append(var_types[nvars-1],vars[nvars-1]->get_type_name());
}
/*------------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        _____   _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      /___  \ |  _  \ 
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |       ___| | | | | | 
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |      /  ___/ | | | | 
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |      | |___  | |_| | 
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|      |_____| |_____/
 ------------------------------------------------------------------------------------------------------*/
template<typename T0>
void Pattern::SubPattern_2D::add_var(T0**& v,const char* var_name_)
{
    for(int i=0;i<dim0;i++)
        for(int j=0;j<dim1;j++)
            kywrds[i][j].add_var(v[i][j],var_name_);
    var_adj();
    vars[nvars-1]=head_kywrd->vars[nvars-rank];
    append(var_types[nvars-1],vars[nvars-1]->get_type_name());
}
/*---------------------------------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        _____   _____        _       _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      /___  \ |  _  \      | |     |_   _|
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |       ___| | | | | |      | |       | |
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |      /  ___/ | | | |      | |       | |
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |      | |___  | |_| |      | |___    | |
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|      |_____| |_____/      |_____|   |_|
 ---------------------------------------------------------------------------------------------------------------------------*/
template<typename T0>
void Pattern::SubPattern_2D_LT::add_var(T0**& v,const char* var_name_)
{
    for(int i=0;i<dim;i++)
        for(int j=0;j<i+1;j++)
            kywrds[i][j].add_var(v[i][j],var_name_);
    var_adj();
    vars[nvars-1]=head_kywrd->vars[nvars-rank];
    append(var_types[nvars-1],vars[nvars-1]->get_type_name());
}
/*---------------------------------------------------------
  _____       ___   _____   _____   _____   _____    __   _
 |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |
 | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |
 |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |
 | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |
 |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|
 ---------------------------------------------------------*/
template<typename T0,typename... T1>
void Pattern::add_var(T0& v,T1&... vrest)
{
    add_var(v);
    add_var(vrest...);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0,typename... T1>
void Pattern::add_var(T0& v,const char* vname,T1&... vrest)
{
    add_var(v,vname);
    add_var(vrest...);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void Pattern::add_var(T0& v)
{
    sub_ptrns[nsub_ptrns-1]->add_var(v,NULL);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void Pattern::add_var(T0& v,const char* var_name_)
{
    sub_ptrns[nsub_ptrns-1]->add_var(v,var_name_);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void Pattern::cmd_dy(const char* n,T0*& vs,int& nvs)
{
    if(dy_sp!=NULL)
    {
        error->abort("only one dynamic command");
    }
    
    dy_sp=new SubPattern_DY<T0>(varmngr,n,vs,nvs);
    curr_sp=dy_sp;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
template<typename T0>
FileReader::Ten_0D<T0>::Ten_0D(FileReader* fr,const char* name_,T0& v_):
v(v_),
Ten(fr,name_)
{
    size_t len=snprintf(NULL,0,"%s = %%s",name)+1;
    format=new char[len];
    sprintf(format,"%s = %%s",name);
    
    kywrds.init(&cmd_logic,varmngr,v,name);
    nkywrds=1;
    head=&kywrds;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
FileReader::Ten_0D<T0>::~Ten_0D()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
FileReader::Ten_1D<T0>::Ten_1D(FileReader* fr,const char* name_,T0*& v_,int& dim0_):
dim0(dim0_),
v(v_),
Ten(fr,name_)
{
    size_t len=snprintf(NULL,0,"%s ( %%d ) = %%s",name)+1;
    format=new char[len];
    sprintf(format,"%s ( %%d ) = %%s",name);

    kywrds=new FileKeyWord[dim0];
    for(int i=0;i<dim0;i++)
        kywrds[i].init(&cmd_logic,varmngr,v[i],"%s(%d)",name,i);
    
    nkywrds=dim0;
    head=&kywrds[0];
    
    char* cmp_name=NULL;
    append(cmp_name,"1st component of %s",name);
    vi0=varmngr.add_var(i0,cmp_name);
    delete [] cmp_name;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
FileReader::Ten_1D<T0>::~Ten_1D()
{
    delete [] kywrds;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
FileReader::Ten_2D<T0>::Ten_2D(FileReader* fr,const char* name_,T0**& v_,int& dim0_,int& dim1_):
dim0(dim0_),dim1(dim1_),
v(v_),
Ten(fr,name_)
{
    size_t len=snprintf(NULL,0,"%s ( %%d , %%d ) = %%s",name)+1;
    format=new char[len];
    sprintf(format,"%s ( %%d , %%d ) = %%s",name);
    
    kywrds=new FileKeyWord*[dim0];
    *kywrds=new FileKeyWord[dim0*dim1];
    for(int i=1;i<dim0;i++)
        kywrds[i]=kywrds[i-1]+dim0;
    
    for(int i=0;i<dim0;i++)
        for(int j=0;j<dim1;j++)
            kywrds[i][j].init(&cmd_logic,varmngr,v[i][j],"%s(%d,%d)",name,i,j);

    nkywrds=dim0*dim1;
    head=&kywrds[0][0];
    
    char* cmp_name=NULL;
    append(cmp_name,"1st component of %s",name);
    vi0=varmngr.add_var(i0,cmp_name);
    delete [] cmp_name;

    
    cmp_name=NULL;
    append(cmp_name,"2nd component of %s",name);
    vi1=varmngr.add_var(i1,cmp_name);
    delete [] cmp_name;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
FileReader::Ten_2D<T0>::~Ten_2D()
{
    delete [] *kywrds;
    delete [] kywrds;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void FileReader::Ten_0D<T0>::set_max_cmp(int& max_cmp)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void FileReader::Ten_1D<T0>::set_max_cmp(int& max_cmp)
{
    Var* pmin_cmp=varmngr(0);
    vli0=Logics(NULL,"ge",pmin_cmp)*Logics(NULL,"lt",varmngr(max_cmp));
    vli0=vlogic("ge",0)*Logics(NULL,"lt",varmngr(max_cmp));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void FileReader::Ten_2D<T0>::set_max_cmp(int& max_cmp)
{
    Var* pmin_cmp=varmngr(0);
    vli0=Logics(NULL,"ge",pmin_cmp)*Logics(NULL,"lt",varmngr(max_cmp));
    vli1=Logics(NULL,"ge",pmin_cmp)*Logics(NULL,"lt",varmngr(max_cmp));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool FileReader::Ten_0D<T0>::scan(char*& err_msg,char** args,int nargs)
{
    if(strcmp(name,args[0]))
        return false;
            
    if(nargs!=2)
    {
        append(err_msg,"quantity %s should be scalar",name);
        return true;
    }
        
    err_msg=kywrds.scan(args,2);
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool FileReader::Ten_1D<T0>::scan(char*& err_msg,char** args,int nargs)
{
    if(strcmp(name,args[0]))
        return false;
    
    if(nargs!=3)
    {
        append(err_msg,"quantity %s should be %d dimensional",name,1);
        return true;
    }
    
    err_msg=vi0->scan(args[2]);
    if(err_msg!=NULL)
        return true;
    
    err_msg=vli0(vi0);
    if(err_msg!=NULL)
        return true;

    int icmp=fr->type_ref[i0];
    if(icmp==-1) return true;
    
    err_msg=kywrds[icmp].scan(args,2);
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
bool FileReader::Ten_2D<T0>::scan(char*& err_msg,char** args,int nargs)
{
    if(strcmp(name,args[0]))
        return false;
    
    if(nargs!=4)
    {
        append(err_msg,"quantity %s should be %d dimensional",name,2);
        return true;
    }
    

    err_msg=vi0->scan(args[2]);
    if(err_msg!=NULL)
        return true;
    
    err_msg=vli0(vi0);
    if(err_msg!=NULL)
        return true;
    
    int icmp=fr->type_ref[i0];
    if(icmp==-1) return true;
    
    err_msg=vi1->scan(args[3]);
    if(err_msg!=NULL)
        return true;
    
    err_msg=vli1(vi1);
    if(err_msg!=NULL)
        return true;
    
    int jcmp=fr->type_ref[i1];
    if(jcmp==-1) return true;
    
    err_msg=kywrds[icmp][jcmp].scan(args,2);
    if(err_msg!=NULL)
        return true;
    
    if(symmetric)
        err_msg=kywrds[jcmp][icmp].scan(args,2);
    
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void FileReader::add_scl(const char* name,T0& val)
{
    Ten** tens_=new Ten*[ntens+1];
    memcpy(tens_,tens,ntens*sizeof(Ten*));
    delete [] tens;
    tens=tens_;
    
    tens[ntens]=new Ten_0D<T0>(this,name,val);
    ntens++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void FileReader::add_1D(const char* name,T0*& val)
{
    Ten** tens_=new Ten*[ntens+1];
    memcpy(tens_,tens,ntens*sizeof(Ten*));
    delete [] tens;
    tens=tens_;
    
    tens[ntens]=new Ten_1D<T0>(this,name,val,ntypes);
    ntens++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void FileReader::add_2D(const char* name,T0**& val)
{
    Ten** tens_=new Ten*[ntens+1];
    memcpy(tens_,tens,ntens*sizeof(Ten*));
    delete [] tens;
    tens=tens_;
    
    tens[ntens]=new Ten_2D<T0>(this,name,val,ntypes,ntypes);
    ntens++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
Logics FileReader::log_u(T0& x,const char* op)
{
    return Logics(varmngr(x),op);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
Logics FileReader::log_u(T0*& x,const char* op)
{
    return Logics(varmngr(x),op);
}

#endif

