#ifndef __MAPP__cmd__
#define __MAPP__cmd__
/*--------------------------------------------
 
 --------------------------------------------*/
#include "error.h"
#include "global.h"
#include "keyword.h"

#include <stdio.h>
#include <string.h>
namespace MAPP_NS
{
    
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
    
    class FileReader
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
        
        class AtomTypes*& atom_types;
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
        MPI_Comm& world;
    protected:
    public:
        VarManager varmngr;
        FileReader();
        ~FileReader();

        void read_file(const char*);
        
        template<typename T0>
        void add_scl(const char*,T0&);
        template<typename T0>
        void add_1D(const char*,T0*&);
        template<typename T0>
        void add_2D(const char*,T0**&);
        

        void symmetric();
        Logics&  add_vlog(){return tens[ntens-1]->add_vlog();}
        Logics&  add_clog(){return tens[ntens-1]->add_clog();}

    };
    

    
    class Pattern
    {
    private:
    protected:
    public:
        static void change_form(char*&,const char*,const char*);
        static void revert_form(char*&);
        class SubPattern;
        template<typename>
        class SubPattern_DY;
        class SubPattern_0D;
        template<size_t>
        class SubPattern_1D;
        class SubPattern_2D;
        template<size_t>
        class SubPattern_Voigt;
        template<size_t>
        class SubPattern_2D_LT;
        
        SubPattern** sub_ptrns;
        int nsub_ptrns;
        
        void add_sp(SubPattern*);
        
        Var* curr_var;
        SubPattern* curr_sp;
        SubPattern* dy_sp;
        int nlog_tree;
        

        VarManager varmngr;
        Pattern();
        ~Pattern();
        
       
        template<typename T0>
        void cmd_dy(const char*,T0*&,int&);
        
        
        void cmd(const char*);
        void cmd(bool&,const char*);
        
        template<size_t>
        void cmd_voigt(const char*);
        template<size_t>
        void cmd_voigt(bool*&,const char*);
        
        template<size_t>
        void cmd_1d(const char*);
        template<size_t>
        void cmd_1d(bool*&,const char*);
        
        template<size_t>
        void cmd_2d_lt(const char*);
        template<size_t>
        void cmd_2d_lt(bool**&,const char*);
        
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
        template<typename T0,size_t N,size_t M>
        void add_var(T0 (&)[N][M],const char*);
        template<typename T0,size_t N>
        void add_var(T0 (&)[N][N],const char*);
        template<typename T0,size_t N>
        void add_var(T0 (&)[N],const char*);
        
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
        

        static constexpr size_t det(size_t n0,size_t n1)
        {
            return n0*(n0+1) < 2*n1 ? det(n0+1,n1):n0;
        }
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
    
    template<size_t dim>
    class Pattern::SubPattern_1D:public Pattern::SubPattern
    {
    private:
        int icomp;
        var<int> comp;
        KeyWord kywrds[dim];
        void init();
        void init_keywords();
        void init_keywords(bool*&);
    protected:
    public:
        char* format;
        SubPattern_1D(VarManager&,const char*,int);
        SubPattern_1D(bool*&,VarManager&,const char*,int);
        ~SubPattern_1D();
        bool scan(char*&,char**&,int&);
        
        template<typename T0>
        void add_var(T0*&,const char*);
        
        void conv_name(PrintStyle&);
    };
    
    template<size_t dim>
    class Pattern::SubPattern_Voigt:public Pattern::SubPattern
    {
    private:
        int icomp0,icomp1;
        Var& comp0;
        Var& comp1;
        KeyWord kywrds[dim*(dim+1)/2];
        void cmpt2voigt(int&,const int,const int);
        void voigt2cmpt(const int,int&,int&);
        void init();
        void init_keywords();
        void init_keywords(bool*&);
    protected:
    public:
        char* format;
        SubPattern_Voigt(VarManager&,const char*);
        SubPattern_Voigt(bool*&,VarManager&,const char*);
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
    
    template<size_t dim>
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
        SubPattern_2D_LT(VarManager&,const char*);
        SubPattern_2D_LT(bool**&,VarManager&,const char*);
        ~SubPattern_2D_LT();
        bool scan(char*&,char**&,int&);
        
        template<typename T0>
        void add_var(T0(&)[dim][dim],const char*);
        
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
    Var::append(_name,"no. of %s",name_);
    var_mngr.adddd_var(n,_name);
    delete [] _name;
    
    var_adj();
    vars[0]=var_mngr.adddd_var(s_val,name);
    Var::append(var_types[0],vars[0]->get_type_name());
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
    Var::append(v_name,"%s_%d",name,ndy_vars);
    Var* var_=var_mngr.adddd_var(dy_vars[ndy_vars],v_name);
    delete [] v_name;
    
    Var** v_dy_vars_=new Var*[ndy_vars+1];
    memcpy(v_dy_vars_,v_dy_vars,ndy_vars*sizeof(Var*));
    delete [] v_dy_vars;
    v_dy_vars=v_dy_vars_;
    v_dy_vars[ndy_vars]=var_;
    ndy_vars++;
    
    err_msg=var_->scan(*args);
    if(err_msg)
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
    Var::append(buff,ps.pattern_prop_brc[0]);
    Var::append(buff,vars[0]->name);
    Var::append(buff," ");
    Var::append(buff,"...");
    Var::append(buff,ps.pattern_prop_brc[1]);
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
    
    if(dynamic_cast<SubPattern_0D*>(this))
        return dynamic_cast<SubPattern_0D*>(this)->add_var(v,n);
    
    printf("unsuccessful conversion at: %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0,size_t N>
void Pattern::SubPattern::add_var(T0 (&v)[N],const char* n)
{
    T0* vv=static_cast<T0*>(v);
    
    if(dynamic_cast<SubPattern_1D<N>*>(this))
        return dynamic_cast<SubPattern_1D<N>*>(this)->add_var(vv,n);
    
    if(dynamic_cast<SubPattern_Voigt<det(0,N)>*>(this))
        return dynamic_cast<SubPattern_Voigt<det(0,N)>*>(this)->add_var(vv,n);
    
    printf("unsuccessful conversion at: %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void Pattern::SubPattern::add_var(T0*& v,const char* n)
{
    if(dynamic_cast<SubPattern_0D*>(this))
        return dynamic_cast<SubPattern_0D*>(this)->add_var(v,n);
    
    printf("unsuccessful conversion at: %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0,size_t N ,size_t M>
void Pattern::SubPattern::add_var(T0 (&v)[N][M],const char* n)
{
    if(dynamic_cast<SubPattern_2D*>(this))
        return dynamic_cast<SubPattern_2D*>(this)->add_var(v,n);
    
    printf("unsuccessful conversion at: %s",__PRETTY_FUNCTION__);
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0,size_t N>
void Pattern::SubPattern::add_var(T0 (&v)[N][N],const char* n)
{
    if(dynamic_cast<SubPattern_2D*>(this))
        return dynamic_cast<SubPattern_2D*>(this)->add_var(v,n);
    if(dynamic_cast<SubPattern_2D_LT<N>*>(this))
        return dynamic_cast<SubPattern_2D_LT<N>*>(this)->add_var(v,n);
    
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
    Var::append(var_types[nvars-1],vars[nvars-1]->get_type_name());
}
/*----------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        ___   _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      |_  | |  _  \
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |        | | | | | |
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |        | | | | | | 
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |        | | | |_| | 
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|        |_| |_____/
 ----------------------------------------------------------------------------------------------------*/
template<size_t dim>
Pattern::SubPattern_1D<dim>::SubPattern_1D(VarManager& var_mngr_,const char* name_,int dim_):
SubPattern(var_mngr_,name_),
comp(icomp,"i (1st compoent of ")
{
    Var::append(comp.name,"%s)",name);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
Pattern::SubPattern_1D<dim>::SubPattern_1D(bool*& set,VarManager& var_mngr_,const char* name_,int dim_):
SubPattern_1D(var_mngr_,name_,dim_)
{
    //Var::append(comp.name,"%s)",name);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
Pattern::SubPattern_1D<dim>::~SubPattern_1D()
{
    delete [] format;
    delete [] kywrds;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_1D<dim>::init()
{
    format=NULL;
    Var::append(format,"%s[%%d]",name);
    kywrds=new KeyWord[dim];
    head_kywrd=kywrds;
    nkywrds=dim;
    rank=1;
    
    var_adj();
    vars[0]=&comp;
    const int dim_=dim;
    var_logics[0]=VLogics("ge",0)*VLogics("lt",dim_);
    
    delete [] vars[0]->name;
    vars[0]->name=NULL;
    Var::append(vars[0]->name,"i (1st compoent of %s)",name);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_1D<dim>::init_keywords()
{
    for(int i=0;i<dim;i++)
        kywrds[i].init(&cmd_logic,var_mngr,format,i);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_1D<dim>::init_keywords(bool*& set)
{
    for(int i=0;i<dim;i++)
        kywrds[i].init(&cmd_logic,var_mngr,set[i],format,i);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
bool Pattern::SubPattern_1D<dim>::scan(char*& err_msg,char**& args,int& nargs)
{
    if(sscanf(*args,format,&icomp)!=1)
        return false;
    
    err_msg=vlog_chek(kywrds[icomp],args,nargs);
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_1D<dim>::conv_name(PrintStyle& ps)
{
    char* name_=NULL;
    Var::append(name_,"%s%s%s%s",ps.cmd_brc[0],name,"[",ps.cmd_brc[1]);
    Var::append(name_,"%s",comp.name);
    Var::append(name_,"%s%s%s",ps.cmd_brc[0],"]",ps.cmd_brc[1]);
    Var::append(name_,"%c%s",'\0',name);
    delete [] name;
    name=name_;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>template<typename T0>
void Pattern::SubPattern_1D<dim>::add_var(T0*& v,const char* var_name_)
{
    for(int i=0;i<dim;i++)
        kywrds[i].add_var(v[i],var_name_);
    var_adj();
    vars[nvars-1]=head_kywrd->vars[nvars-rank];
    Var::append(var_types[nvars-1],vars[nvars-1]->get_type_name());
}
/*----------------------------------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        _     _   _____   _   _____   _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      | |   / / /  _  \ | | /  ___| |_   _| 
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |      | |  / /  | | | | | | | |       | |   
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |      | | / /   | | | | | | | |  _    | |   
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |      | |/ /    | |_| | | | | |_| |   | |   
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|      |___/     \_____/ |_| \_____/   |_|
 ----------------------------------------------------------------------------------------------------------------------------*/
template<size_t dim>
Pattern::SubPattern_Voigt<dim>::SubPattern_Voigt(VarManager& var_mngr_,const char* name_):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.adddd_var(icomp0,"")),
comp1(*var_mngr.adddd_var(icomp1,""))
{
    init();
    init_keywords();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
Pattern::SubPattern_Voigt<dim>::SubPattern_Voigt(bool*& set,VarManager& var_mngr_,const char* name_):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.adddd_var(icomp0,"")),
comp1(*var_mngr.adddd_var(icomp1,""))
{
    init();
    init_keywords(set);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
Pattern::SubPattern_Voigt<dim>::~SubPattern_Voigt()
{
    delete [] format;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_Voigt<dim>::init()
{
    format=NULL;
    Var::append(format,"%s[%%d][%%d]",name);
    
    head_kywrd=kywrds;
    nkywrds=dim*(dim+1)/2;
    rank=2;
    
    var_adj();
    var_adj();
    vars[0]=&comp0;
    vars[1]=&comp1;
    const int dim_=dim;
    var_logics[0]=VLogics("ge",0)*VLogics("lt",dim_);
    var_logics[1]=VLogics("ge",0)*VLogics("lt",dim_);
    
    delete [] vars[0]->name;
    vars[0]->name=NULL;
    Var::append(vars[0]->name,"i (1st compoent of %s)",name);
    Var::append(var_types[0],vars[0]->get_type_name());
    
    delete [] vars[1]->name;
    vars[1]->name=NULL;
    Var::append(vars[1]->name,"j (2nd compoent of %s)",name);
    Var::append(var_types[1],vars[1]->get_type_name());
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_Voigt<dim>::init_keywords()
{
    for(int i=0;i<dim*(dim+1)/2;i++)
    {
        int icmp,jcmp;
        voigt2cmpt(i,icmp,jcmp);
        kywrds[i].init(&cmd_logic,var_mngr,format,icmp,jcmp);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_Voigt<dim>::init_keywords(bool*& set)
{
    for(int i=0;i<dim*(dim+1)/2;i++)
    {
        int icmp,jcmp;
        voigt2cmpt(i,icmp,jcmp);
        kywrds[i].init(&cmd_logic,var_mngr,set[i],format,icmp,jcmp);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_Voigt<dim>::cmpt2voigt(int& icmp,const int i_,const int j_)
{
    icmp=0;
    int i,j;
    if(i_>j_)
    {
        i=i_;
        j=j_;
    }
    else
    {
        i=j_;
        j=i_;
    }
    
    for(int i0=0,j0=0;i0<dim;i0++,j0++,icmp++)
        if(i0==i && j0==j)
            return;
    
    for(int i0=dim-1;i0>-1;i0--)
    {
        if((dim-i0) %2==1)
        {
            for(int j0=i0-1;j0>-1;j0--,icmp++)
                if(i0==i && j0==j)
                    return;
        }
        else
        {
            for(int j0=0;j0<i0;j0++,icmp++)
                if(i0==i && j0==j)
                    return;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_Voigt<dim>::voigt2cmpt(const int icmp_,int& i0,int& j0)
{
    int icmp=0;
    i0=j0=0;
    for(;i0<dim;i0++,j0++,icmp++)
        if(icmp==icmp_)
            return;
    
    for(i0=dim-1;i0>-1;i0--)
    {
        if((dim-i0) %2==1)
        {
            for(j0=i0-1;j0>-1;j0--,icmp++)
                if(icmp==icmp_)
                    return;
        }
        else
        {
            for(j0=0;j0<i0;j0++,icmp++)
                if(icmp==icmp_)
                    return;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
bool Pattern::SubPattern_Voigt<dim>::scan(char*& err_msg,char**& args,int& nargs)
{
    if(sscanf(*args,format,&icomp0,&icomp1)!=2)
        return false;
    
    int ii;
    cmpt2voigt(ii,icomp0,icomp1);
    err_msg=vlog_chek(kywrds[ii],args,nargs);
    
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_Voigt<dim>::conv_name(PrintStyle& ps)
{
    char* name_=NULL;
    Var::append(name_,"%s%s%s%s",ps.cmd_brc[0],name,"[",ps.cmd_brc[1]);
    Var::append(name_,"%s",comp0.name);
    Var::append(name_,"%s%s%s",ps.cmd_brc[0],"][",ps.cmd_brc[1]);
    Var::append(name_,"%s",comp1.name);
    Var::append(name_,"%s%s%s",ps.cmd_brc[0],"]",ps.cmd_brc[1]);
    Var::append(name_,"%c%s",'\0',name);
    delete [] name;
    name=name_;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim> template<typename T0>
void Pattern::SubPattern_Voigt<dim>::add_var(T0*& v,const char* var_name_)
{
    for(int i=0;i<dim*(dim+1)/2;i++)
        kywrds[i].add_var(v[i],var_name_);
    var_adj();
    
    vars[nvars-1]=head_kywrd->vars[nvars-rank];
    Var::append(var_types[nvars-1],vars[nvars-1]->get_type_name());
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
    Var::append(var_types[nvars-1],vars[nvars-1]->get_type_name());
}
/*---------------------------------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        _____   _____        _       _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      /___  \ |  _  \      | |     |_   _|
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |       ___| | | | | |      | |       | |
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |      /  ___/ | | | |      | |       | |
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |      | |___  | |_| |      | |___    | |
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|      |_____| |_____/      |_____|   |_|
 ---------------------------------------------------------------------------------------------------------------------------*/
template<size_t dim>
Pattern::SubPattern_2D_LT<dim>::SubPattern_2D_LT(VarManager& var_mngr_,const char* name_):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.adddd_var(icomp0,"")),
comp1(*var_mngr.adddd_var(icomp1,""))
{
    init();
    init_keywords();
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
Pattern::SubPattern_2D_LT<dim>::SubPattern_2D_LT(bool**& set,VarManager& var_mngr_,const char* name_):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.adddd_var(icomp0,"")),
comp1(*var_mngr.adddd_var(icomp1,""))
{
    init();
    init_keywords(set);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
Pattern::SubPattern_2D_LT<dim>::~SubPattern_2D_LT()
{
    delete [] format;
    delete [] *kywrds;
    delete [] kywrds;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_2D_LT<dim>::init()
{
    format=NULL;
    Var::append(format,"%s[%%d][%%d]",name);
    kywrds=new KeyWord*[dim];
    *kywrds=new KeyWord[dim*(dim+1)/2];
    for(int i=1;i<dim;i++)
        kywrds[i]=kywrds[i-1]+i;
    
    rank=2;
    head_kywrd=*kywrds;
    nkywrds=dim*(dim+1)/2;
    
    var_adj();
    var_adj();
    vars[0]=&comp0;
    vars[1]=&comp1;
    const int dim_=dim;
    var_logics[0]=VLogics("ge",0)*VLogics("lt",dim_);
    var_logics[1]=VLogics("ge",0)*VLogics("le",icomp0);
    
    delete [] vars[0]->name;
    vars[0]->name=NULL;
    Var::append(vars[0]->name,"i (1st compoent of %s)",name);
    Var::append(var_types[0],vars[0]->get_type_name());
    
    delete [] vars[1]->name;
    vars[1]->name=NULL;
    Var::append(vars[1]->name,"j (2nd compoent of %s)",name);
    Var::append(var_types[1],vars[1]->get_type_name());
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_2D_LT<dim>::init_keywords()
{
    for(int i=0;i<dim;i++)
        for(int j=0;j<i+1;j++)
            kywrds[i][j].init(&cmd_logic,var_mngr,format,i,j);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_2D_LT<dim>::init_keywords(bool**& set)
{
    for(int i=0;i<dim;i++)
        for(int j=0;j<i+1;j++)
            kywrds[i][j].init(&cmd_logic,var_mngr,set[i][j],format,i,j);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
bool Pattern::SubPattern_2D_LT<dim>::scan(char*& err_msg,char**& args,int& nargs)
{
    if(sscanf(*args,format,&icomp0,&icomp1)!=2)
        return false;
    
    err_msg=vlog_chek(kywrds[icomp0][icomp1],args,nargs);
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::SubPattern_2D_LT<dim>::conv_name(PrintStyle& ps)
{
    char* name_=NULL;
    Var::append(name_,"%s%s%s%s",ps.cmd_brc[0],name,"[",ps.cmd_brc[1]);
    Var::append(name_,"%s",comp0.name);
    Var::append(name_,"%s%s%s",ps.cmd_brc[0],"][",ps.cmd_brc[1]);
    Var::append(name_,"%s",comp1.name);
    Var::append(name_,"%s%s%s",ps.cmd_brc[0],"]",ps.cmd_brc[1]);
    Var::append(name_,"%c%s",'\0',name);
    delete [] name;
    name=name_;
}
/*--------------------------------------------
     
 --------------------------------------------*/
template<size_t dim>template<typename T0>
void Pattern::SubPattern_2D_LT<dim>::add_var(T0(&v)[dim][dim],const char* var_name_)
{
    for(int i=0;i<dim;i++)
        for(int j=0;j<i+1;j++)
            kywrds[i][j].add_var(v[i][j],var_name_);
    var_adj();
    vars[nvars-1]=head_kywrd->vars[nvars-rank];
    Var::append(var_types[nvars-1],vars[nvars-1]->get_type_name());
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
    if(dy_sp)
    {
        Error::abort("only one dynamic command");
    }
    
    dy_sp=new SubPattern_DY<T0>(varmngr,n,vs,nvs);
    curr_sp=dy_sp;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::cmd_voigt(const char* n)
{
    add_sp(new SubPattern_Voigt<dim>(varmngr,n));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::cmd_voigt(bool*& set,const char* n)
{
    add_sp(new SubPattern_Voigt<dim>(set,varmngr,n));
}
/*--------------------------------------------
     
 --------------------------------------------*/
template<size_t dim>
void Pattern::cmd_1d(const char* n)
{
    add_sp(new SubPattern_1D<dim>(varmngr,n));
}
/*--------------------------------------------
     
 --------------------------------------------*/
template<size_t dim>
void Pattern::cmd_1d(bool*& set,const char* n)
{
    add_sp(new SubPattern_1D<dim>(set,varmngr,n));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::cmd_2d_lt(const char* n)
{
    add_sp(new SubPattern_2D_LT<dim>(varmngr,n));
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<size_t dim>
void Pattern::cmd_2d_lt(bool**& set,const char* n)
{
    add_sp(new SubPattern_2D_LT<dim>(set,varmngr,n));
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
    Var::append(cmp_name,"1st component of %s",name);
    vi0=varmngr.adddd_var(i0,cmp_name);
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
    Var::append(cmp_name,"1st component of %s",name);
    vi0=varmngr.adddd_var(i0,cmp_name);
    delete [] cmp_name;

    
    cmp_name=NULL;
    Var::append(cmp_name,"2nd component of %s",name);
    vi1=varmngr.adddd_var(i1,cmp_name);
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
    vli0=VLogics("ge",0)*VLogics("lt",max_cmp);
}
/*--------------------------------------------
 
 --------------------------------------------*/
template<typename T0>
void FileReader::Ten_2D<T0>::set_max_cmp(int& max_cmp)
{
    vli0=VLogics("ge",0)*VLogics("lt",max_cmp);
    vli1=VLogics("ge",0)*VLogics("lt",max_cmp);
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
        Var::append(err_msg,"quantity %s should be scalar",name);
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
        Var::append(err_msg,"quantity %s should be %d dimensional",name,1);
        return true;
    }
    
    err_msg=vi0->scan(args[2]);
    if(err_msg)
        return true;
    
    err_msg=vli0(vi0);
    if(err_msg)
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
        Var::append(err_msg,"quantity %s should be %d dimensional",name,2);
        return true;
    }
    

    err_msg=vi0->scan(args[2]);
    if(err_msg)
        return true;
    
    err_msg=vli0(vi0);
    if(err_msg)
        return true;
    
    int icmp=fr->type_ref[i0];
    if(icmp==-1) return true;
    
    err_msg=vi1->scan(args[3]);
    if(err_msg)
        return true;
    
    err_msg=vli1(vi1);
    if(err_msg)
        return true;
    
    int jcmp=fr->type_ref[i1];
    if(jcmp==-1) return true;
    
    err_msg=kywrds[icmp][jcmp].scan(args,2);
    if(err_msg)
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


#endif

