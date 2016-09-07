#include "cmd.h"
#include "atom_types.h"
#include "atoms.h"
#include "script_reader.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
PrintStyle::PrintStyle()
{
    var_brc[0]=NULL;
    var_brc[1]=NULL;
    type_brc[0]=NULL;
    type_brc[1]=NULL;
    cmd_brc[0]=NULL;
    cmd_brc[1]=NULL;
    cmd_prop_brc[0]=NULL;
    cmd_prop_brc[1]=NULL;
    all_vars_prop_brc[0]=NULL;
    all_vars_prop_brc[1]=NULL;
    all_cmd_prop_brc[0]=NULL;
    all_cmd_prop_brc[1]=NULL;
    vars_prop_brc[0]=NULL;
    vars_prop_brc[1]=NULL;
    pattern_prop_brc[0]=NULL;
    pattern_prop_brc[1]=NULL;
}
/*--------------------------------------------

 --------------------------------------------*/
PrintStyle::~PrintStyle()
{
    
    delete [] var_brc[0];
    delete [] var_brc[1];
    delete [] type_brc[0];
    delete [] type_brc[1];
    delete [] cmd_brc[0];
    delete [] cmd_brc[1];
    delete [] cmd_prop_brc[0];
    delete [] cmd_prop_brc[1];
    delete [] all_vars_prop_brc[0];
    delete [] all_vars_prop_brc[1];
    delete [] all_cmd_prop_brc[0];
    delete [] all_cmd_prop_brc[1];
    delete [] vars_prop_brc[0];
    delete [] vars_prop_brc[1];
    delete [] pattern_prop_brc[0];
    delete [] pattern_prop_brc[1];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PrintStyle::add(char* (&x)[2],const char* str_0,const char* str_1)
{
    size_t l;
    
    delete [] x[0];
    l=strlen(str_0)+1;
    x[0]=new char[l];
    memcpy(x[0],str_0,l);
    
    delete [] x[1];
    l=strlen(str_1)+1;
    x[1]=new char[l];
    memcpy(x[1],str_1,l);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void PrintStyle::var(const char* str_0,const char* str_1)
{add(var_brc,str_0,str_1);}
/*--------------------------------------------
 
 --------------------------------------------*/
void PrintStyle::type(const char* str_0,const char* str_1)
{add(type_brc,str_0,str_1);}
/*--------------------------------------------
 
 --------------------------------------------*/
void PrintStyle::cmd(const char* str_0,const char* str_1)
{add(cmd_brc,str_0,str_1);}
/*--------------------------------------------
 
 --------------------------------------------*/
void PrintStyle::cmd_prop(const char* str_0,const char* str_1)
{add(cmd_prop_brc,str_0,str_1);}
/*--------------------------------------------
 
 --------------------------------------------*/
void PrintStyle::all_vars_prop(const char* str_0,const char* str_1)
{add(all_vars_prop_brc,str_0,str_1);}
/*--------------------------------------------
 
 --------------------------------------------*/
void PrintStyle::all_cmd_prop(const char* str_0,const char* str_1)
{add(all_cmd_prop_brc,str_0,str_1);}
/*--------------------------------------------
 
 --------------------------------------------*/
void PrintStyle::vars_prop(const char* str_0,const char* str_1)
{add(vars_prop_brc,str_0,str_1);}
/*--------------------------------------------
 
 --------------------------------------------*/
void PrintStyle::pattern_prop(const char* str_0,const char* str_1)
{add(pattern_prop_brc,str_0,str_1);}

/*---------------------------------------------------------
 _____       ___   _____   _____   _____   _____    __   _  
|  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | | 
| |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| | 
|  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   | 
| |      / /  | |   | |     | |   | |___  | | \ \  | | \  | 
|_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|
 
 ---------------------------------------------------------*/
Pattern::Pattern()
{
    curr_var=NULL;
    nsub_ptrns=0;
    sub_ptrns=NULL;
    curr_sp=NULL;
    dy_sp=NULL;
    nlog_tree=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::~Pattern()
{
    for(int i=0;i<nsub_ptrns;i++)
        delete sub_ptrns[i];
    delete [] sub_ptrns;
    delete dy_sp;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::change_form(char*& str,const char* l,const char* r)
{
    if(str==NULL)
        return;
    char* _str=NULL;
    Var::append(_str,"%s%s%s%c%s",l,str,r,'\0',str);
    delete [] str;
    str=_str;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::revert_form(char*& str)
{
    if(str==NULL)
        return;
    char* pch=strchr(str,'\0');
    size_t len=strchr(pch+1,'\0')-pch;
    char* _str=new char[len];
    memcpy(_str,pch+1,len);
    delete [] str;
    str=_str;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::add_sp(SubPattern* sp)
{
    SubPattern** sub_ptrns_=new SubPattern*[nsub_ptrns+1];
    memcpy(sub_ptrns_,sub_ptrns,nsub_ptrns*sizeof(SubPattern*));
    delete [] sub_ptrns;
    sub_ptrns=sub_ptrns_;
    sub_ptrns[nsub_ptrns]=sp;
    nsub_ptrns++;
    curr_sp=sp;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::cmd(const char* n)
{
    add_sp(new SubPattern_0D(varmngr,n));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::cmd(bool& set,const char* n)
{
    add_sp(new SubPattern_0D(set,varmngr,n));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::cmd_voigt(const char* n,int d)
{
    add_sp(new SubPattern_Voigt(varmngr,n,d));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::cmd_voigt(bool*& set,const char* n,int d)
{
    add_sp(new SubPattern_Voigt(set,varmngr,n,d));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::cmd_1d(const char* n,const int d)
{
    add_sp(new SubPattern_1D(varmngr,n,d));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::cmd_1d(bool*& set,const char* n,const int d)
{
    add_sp(new SubPattern_1D(set,varmngr,n,d));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::cmd_2d_lt(const char* n,const int d)
{
    add_sp(new SubPattern_2D_LT(varmngr,n,d));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::cmd_2d_lt(bool**& set,const char* n,const int d)
{
    add_sp(new SubPattern_2D_LT(set,varmngr,n,d));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::cmd_2d(const char* n,const int d0,const int d1)
{
    add_sp(new SubPattern_2D(varmngr,n,d0,d1));
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::cmd_2d(bool**& set,const char* n,const int d0,const int d1)
{
    add_sp(new SubPattern_2D(set,varmngr,n,d0,d1));
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Pattern::scan(char**& args,int& nargs)
{
    char* err_msg=NULL;
    bool ret=sub_ptrns[0]->scan(err_msg,args,nargs);
    if(err_msg!=NULL)
        error->abort("%s",err_msg);
    if(ret==false)
        error->abort("");
    while(nargs)
    {
        ret=false;
        for(int i=1;i<nsub_ptrns && !ret;i++)
        {
            ret=sub_ptrns[i]->scan(err_msg,args,nargs);
            if(err_msg!=NULL)
                error->abort("in %s, %s",sub_ptrns[0]->name,err_msg);
        }
        
        if(!ret)
        {
            if(dy_sp==NULL)
                error->abort("in %s, unknown keyword/value %s",sub_ptrns[0]->name,*args);
            else
            {
                dy_sp->scan(err_msg,args,nargs);
                if(err_msg!=NULL)
                    error->abort("%s",err_msg);
            }
        }
    }
    
    
    for(int i=0;i<nsub_ptrns;i++)
    {
        char* msg=sub_ptrns[i]->finalize();
        if(msg!=NULL)
            error->abort("in %s, %s",sub_ptrns[0]->name,msg);
    }
    
    if(dy_sp!=NULL)
    {
        char* msg=dy_sp->finalize();
        if(msg!=NULL)
            error->abort("in %s, %s",sub_ptrns[0]->name,msg);
    }
        
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::add_vdesc(int i,const char* d)
{
    curr_sp->add_vdesc(i,d);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::add_cdesc(const char* d)
{
    curr_sp->add_cdesc(d);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Pattern::add_vlog(int i)
{
    
    return curr_sp->add_vlog(i);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Pattern::add_clog()
{
    return curr_sp->add_clog();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::print_info()
{
    if(nsub_ptrns==0)
        return;
    
    PrintStyle ps;
    ps.var("{\\color{blue}\\textit{","}}");
    ps.type("{\\color{red}$<$\\texttt{","}$>$}");
    ps.cmd("{\\color{green}\\texttt{","}}");
    ps.cmd_prop("\n","\n");
    ps.all_vars_prop("\\begin{itemize}\n","\\end{itemize}\n");
    ps.all_cmd_prop("\\begin{itemize}\n","\\end{itemize}\n");
    ps.vars_prop("\\item ","\n");
    ps.pattern_prop("\\item ","\n");
    
    /*---------------------------------------------------------------------------*/
    for(int i=0;i<varmngr.ndy_vars;i++)
    {
        varmngr.dy_vars[i]->det_loc();
        change_form(varmngr.dy_vars[i]->name,ps.var_brc[0],ps.var_brc[1]);
    }
    
    for(int i=0;i<nsub_ptrns;i++)
        sub_ptrns[i]->conv_name(ps);
    
    for(int i=0;i<nsub_ptrns;i++)
        for(int j=0;j<sub_ptrns[i]->nvars;j++)
            change_form(sub_ptrns[i]->var_types[j],ps.type_brc[0],ps.type_brc[1]);
    
    if(dy_sp!=NULL)
        change_form(dy_sp->var_types[0],ps.type_brc[0],ps.type_brc[1]);
    /*---------------------------------------------------------------------------*/
    
    char* buff=NULL;
    Var::append(buff,ps.all_cmd_prop_brc[0]);
    
    for(int i=0;i<nsub_ptrns;i++)
    {
        SubPattern* sp=sub_ptrns[i];
        sp->print_pattern(ps,buff);
        sp->print_cmd_prop(ps,buff);
        sp->print_vars_prop(ps,buff);
    }
    
    if(dy_sp!=NULL)
    {
        dy_sp->print_pattern(ps,buff);
        dy_sp->print_cmd_prop(ps,buff);
        dy_sp->print_vars_prop(ps,buff);
    }
    Var::append(buff,ps.all_cmd_prop_brc[1]);
    
    FindReplace func0;
    func0("_","\\_");
    func0(buff);
    
    printf("%s",buff);
    delete [] buff;
    
    
    /*---------------------------------------------------------------------------*/
    for(int i=0;i<varmngr.ndy_vars;i++)
    {
        revert_form(varmngr.dy_vars[i]->name);
        varmngr.dy_vars[i]->att_loc();
    }
    
    for(int i=0;i<nsub_ptrns;i++)
        sub_ptrns[i]->deconv_name();
    
    for(int i=0;i<nsub_ptrns;i++)
        for(int j=0;j<sub_ptrns[i]->nvars;j++)
            revert_form(sub_ptrns[i]->var_types[j]);
    
    if(dy_sp!=NULL)
        revert_form(dy_sp->var_types[0]);
    /*---------------------------------------------------------------------------*/
    
}
/*---------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _  
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | | 
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| | 
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   | 
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  | 
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|
 ---------------------------------------------------------------------------------*/
Pattern::SubPattern::SubPattern(VarManager& var_mngr_,const char* name_):
var_mngr(var_mngr_)
{
    name=new char[strlen(name_)+1];
    memcpy(name,name_,strlen(name_)+1);
    nvars=0;
    cmd_desc=NULL;
    var_descs=NULL;
    var_types=NULL;
    var_logics=NULL;
    head_kywrd=NULL;
    vars=NULL;
    nkywrds=0;
    nvars=0;
    rank=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::SubPattern::~SubPattern()
{    
    for(int i=0;i<nvars;i++)
    {
        delete [] var_descs[i];
        delete [] var_types[i];
    }
    delete [] cmd_desc;
    delete [] var_descs;
    delete [] var_types;
    delete [] var_logics;
    delete [] name;
    delete [] vars;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern::var_adj()
{
    Var** vars_=new Var*[nvars+1];
    memcpy(vars_,vars,nvars*sizeof(Var*));
    delete [] vars;
    vars=vars_;
    
    Logics* var_logics_=new Logics[nvars+1];
    for(int i=0;i<nvars;i++) var_logics_[i]=var_logics[i];
    delete [] var_logics;
    var_logics=var_logics_;
    
    char** var_types_=new char*[nvars+1];
    memcpy(var_types_,var_types,nvars*sizeof(char*));
    delete [] var_types;
    var_types=var_types_;
    
    char** var_descs_=new char*[nvars+1];
    memcpy(var_descs_,var_descs,nvars*sizeof(char*));
    delete [] var_descs;
    var_descs=var_descs_;
    
    var_descs[nvars]=NULL;
    var_types[nvars]=NULL;
    nvars++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern::add_vdesc(int i,const char* desc)
{
    i+=rank;
    if(i>=nvars)
        error->abort("%s",__PRETTY_FUNCTION__);

    delete [] var_descs[i];
    var_descs[i]=NULL;
    Var::append(var_descs[i],desc);
    
    size_t len=strlen(desc);
    if(var_descs[i][len-1]!='.')
        Var::append(var_descs[i],".");
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern::add_cdesc(const char* desc)
{
    delete [] cmd_desc;
    cmd_desc=NULL;
    Var::append(cmd_desc,desc);
    
    size_t len=strlen(desc);
    if(cmd_desc[len-1]!='.')
        Var::append(cmd_desc,".");
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Pattern::SubPattern::add_vlog(int i)
{
    i+=rank;
    if(i>=nvars)
        error->abort("%s",__PRETTY_FUNCTION__);
    return var_logics[i];
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& Pattern::SubPattern::add_clog()
{
    return cmd_logic;
}
/*--------------------------------------------
 
 --------------------------------------------*/
char* Pattern::SubPattern::vlog_chek(KeyWord& k,char**& args,int& nargs)
{
    char* err_msg=NULL;
    
    for(int i=0;i<rank;i++)
    {
        err_msg=var_logics[i](vars[i]);
        if(err_msg!=NULL)
            return err_msg;
    }
    
    err_msg=k.scan(args,nargs);
    if(err_msg!=NULL)
        return err_msg;
    
    for(int i=rank;i<nvars;i++)
    {
        err_msg=var_logics[i](&k[i-rank]);
        if(err_msg!=NULL)
            return err_msg;
    }
    
    
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
char* Pattern::SubPattern::finalize()
{
    for(int i=0;i<nkywrds;i++)
    {
        char* msg=head_kywrd[i].finalize();
        if(msg!=NULL)
            return msg;
    }
    
    return NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern::print_pattern(PrintStyle& ps,char*& buff)
{
    Var::append(buff,ps.pattern_prop_brc[0]);
    Var::append(buff,name);
    for(int i=rank;i<nvars;i++)
    {
        Var::append(buff," ");
        Var::append(buff,vars[i]->name);
    }
    Var::append(buff,ps.pattern_prop_brc[1]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern::print_cmd_prop(PrintStyle& ps,char*& buff)
{

    Var* v=new Var(name,NULL);
    char* msg=cmd_logic.print(v);
    delete v;
    
    
    Var::append(buff,ps.cmd_prop_brc[0]);
    
    if(cmd_desc==NULL && msg==NULL)
        return;
    
    if(cmd_desc!=NULL)
    {
        Var::append(buff,cmd_desc);
    }
    
    if(msg!=NULL)
    {
        if(cmd_desc!=NULL)
            Var::append(buff," ");
        Var::append(buff,msg);
    }
    delete [] msg;
    
    Var::append(buff,ps.cmd_prop_brc[1]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern::print_vars_prop(PrintStyle& ps,char*& buff)
{
    if(nvars==0)
        return;
    Var::append(buff,ps.all_vars_prop_brc[0]);
    
    for(int i=0;i<nvars;i++)
        print_vars_prop(ps,buff,i);
        
    Var::append(buff,ps.all_vars_prop_brc[1]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern::print_vars_prop(PrintStyle& ps,char*& buff,int i)
{
    Var::append(buff,ps.vars_prop_brc[0]);
    
    Var::append(buff,vars[i]->name);
    Var::append(buff,var_types[i]);
    
    
    char* msg=var_logics[i].print(vars[i]);
    
    if(var_descs[i]!=NULL)
    {
        Var::append(buff," ");
        Var::append(buff,var_descs[i]);
    }
    
    if(msg!=NULL)
    {
        Var::append(buff," ");
        Var::append(buff,msg);
    }
    delete [] msg;
    
    Var::append(buff,ps.vars_prop_brc[1]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern::deconv_name()
{
    if(name==NULL)
        return;
    
    char* pch=strchr(name,'\0');
    size_t len=strchr(pch+1,'\0')-pch;
    char* name_=new char[len];
    memcpy(name_,pch+1,len);
    delete [] name;
    name=name_;
}
/*------------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        _____   _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      /  _  \ |  _  \
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |      | | | | | | | |
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |      | |/| | | | | |
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |      | |_| | | |_| |
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|      \_____/ |_____/
 ------------------------------------------------------------------------------------------------------*/
Pattern::SubPattern_0D::SubPattern_0D(VarManager& var_mngr_,const char* name_):
SubPattern(var_mngr_,name_)
{
    init();
    init_keywords();
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::SubPattern_0D::SubPattern_0D(bool& set,VarManager& var_mngr_,const char* name_):
SubPattern(var_mngr_,name_)
{
    init();
    init_keywords(set);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::SubPattern_0D::~SubPattern_0D()
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_0D::init()
{
    head_kywrd=&kywrd;
    nkywrds=1;
    rank=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_0D::init_keywords()
{
    kywrd.init(&cmd_logic,var_mngr,name);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_0D::init_keywords(bool& set)
{
    kywrd.init(&cmd_logic,var_mngr,set,name);
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Pattern::SubPattern_0D::scan(char*& err_msg,char**& args,int& nargs)
{
    if(strcmp(*args,name))
        return false;

    err_msg=vlog_chek(kywrd,args,nargs);

    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_0D::conv_name(PrintStyle& ps)
{
    char* name_=NULL;
    Var::append(name_,"%s%s%s",ps.cmd_brc[0],name,ps.cmd_brc[1]);
    Var::append(name_,"%c%s",'\0',name);
    delete [] name;
    name=name_;
}
/*----------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        ___   _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      |_  | |  _  \
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |        | | | | | |
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |        | | | | | | 
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |        | | | |_| | 
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|        |_| |_____/
 ----------------------------------------------------------------------------------------------------*/
Pattern::SubPattern_1D::SubPattern_1D(VarManager& var_mngr_,const char* name_,int dim_):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.add_var(icomp0,"")),
dim(dim_)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::SubPattern_1D::SubPattern_1D(bool*& set,VarManager& var_mngr_,const char* name_,int dim_):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.add_var(icomp0,"")),
dim(dim_)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::SubPattern_1D::~SubPattern_1D()
{
    delete [] format;
    delete [] kywrds;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_1D::init()
{
    format=NULL;
    Var::append(format,"%s[%%d]",name);
    kywrds=new KeyWord[dim];
    head_kywrd=kywrds;
    nkywrds=dim;
    rank=1;
    
    var_adj();
    vars[0]=&comp0;
    const int dim_=dim;
    var_logics[0]=Logics(NULL,"ge",var_mngr(0))*Logics(NULL,"lt",var_mngr(dim_));
    
    delete [] vars[0]->name;
    vars[0]->name=NULL;
    Var::append(vars[0]->name,"i (1st compoent of %s)",name);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_1D::init_keywords()
{
    for(int i=0;i<dim;i++)
        kywrds[i].init(&cmd_logic,var_mngr,format,i);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_1D::init_keywords(bool*& set)
{
    for(int i=0;i<dim;i++)
        kywrds[i].init(&cmd_logic,var_mngr,set[i],format,i);
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Pattern::SubPattern_1D::scan(char*& err_msg,char**& args,int& nargs)
{
    if(sscanf(*args,format,&icomp0)!=1)
        return false;
    
    err_msg=vlog_chek(kywrds[icomp0],args,nargs);
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_1D::conv_name(PrintStyle& ps)
{
    char* name_=NULL;
    Var::append(name_,"%s%s%s%s",ps.cmd_brc[0],name,"[",ps.cmd_brc[1]);
    Var::append(name_,"%s",comp0.name);
    Var::append(name_,"%s%s%s",ps.cmd_brc[0],"]",ps.cmd_brc[1]);
    Var::append(name_,"%c%s",'\0',name);
    delete [] name;
    name=name_;
}
/*----------------------------------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        _     _   _____   _   _____   _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      | |   / / /  _  \ | | /  ___| |_   _| 
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |      | |  / /  | | | | | | | |       | |   
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |      | | / /   | | | | | | | |  _    | |   
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |      | |/ /    | |_| | | | | |_| |   | |   
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|      |___/     \_____/ |_| \_____/   |_|
 ----------------------------------------------------------------------------------------------------------------------------*/
Pattern::SubPattern_Voigt::SubPattern_Voigt(VarManager& var_mngr_,const char* name_,int dim_):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.add_var(icomp0,"")),
comp1(*var_mngr.add_var(icomp1,"")),
dim(dim_)
{
    init();
    init_keywords();
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::SubPattern_Voigt::SubPattern_Voigt(bool*& set,VarManager& var_mngr_,const char* name_,int dim_):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.add_var(icomp0,"")),
comp1(*var_mngr.add_var(icomp1,"")),
dim(dim_)
{
    init();
    init_keywords(set);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::SubPattern_Voigt::~SubPattern_Voigt()
{
    delete [] format;
    delete [] kywrds;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_Voigt::init()
{
    format=NULL;
    Var::append(format,"%s[%%d][%%d]",name);
    kywrds=new KeyWord[dim*(dim+1)/2];
    
    head_kywrd=kywrds;
    nkywrds=dim*(dim+1)/2;
    rank=2;
    
    var_adj();
    var_adj();
    vars[0]=&comp0;
    vars[1]=&comp1;
    const int dim_=dim;
    var_logics[0]=Logics(NULL,"ge",var_mngr(0))*Logics(NULL,"lt",var_mngr(dim_));
    var_logics[1]=Logics(NULL,"ge",var_mngr(0))*Logics(NULL,"lt",var_mngr(dim_));
    
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
void Pattern::SubPattern_Voigt::init_keywords()
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
void Pattern::SubPattern_Voigt::init_keywords(bool*& set)
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
void Pattern::SubPattern_Voigt::cmpt2voigt(int& icmp,const int i_,const int j_)
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
void Pattern::SubPattern_Voigt::voigt2cmpt(const int icmp_,int& i0,int& j0)
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
bool Pattern::SubPattern_Voigt::scan(char*& err_msg,char**& args,int& nargs)
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
void Pattern::SubPattern_Voigt::conv_name(PrintStyle& ps)
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
/*-----------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____      ___   _____   _____   _____   _____    __   _        _____   _____
/  ___/ | | | | |  _  \ |  _  \    /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      /___  \ |  _  \
| |___  | | | | | |_| | | |_| |   / /| |   | |     | |   | |__   | |_| |  |   \| |       ___| | | | | |
\___  \ | | | | |  _  { |  ___/  / / | |   | |     | |   |  __|  |  _  /  | |\   |      /  ___/ | | | |
 ___| | | |_| | | |_| | | |     / /  | |   | |     | |   | |___  | | \ \  | | \  |      | |___  | |_| |
/_____/ \_____/ |_____/ |_|    /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|      |_____| |_____/
 -----------------------------------------------------------------------------------------------------*/
Pattern::SubPattern_2D::SubPattern_2D(VarManager& var_mngr_,const char* name_,int dim_0,int dim_1):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.add_var(icomp0,"")),
comp1(*var_mngr.add_var(icomp1,"")),
dim0(dim_0),
dim1(dim_1)
{
    init();
    init_keywords();
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::SubPattern_2D::SubPattern_2D(bool**&set,VarManager& var_mngr_,const char* name_,int dim_0,int dim_1):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.add_var(icomp0,"")),
comp1(*var_mngr.add_var(icomp1,"")),
dim0(dim_0),
dim1(dim_1)
{
    init();
    init_keywords(set);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::SubPattern_2D::~SubPattern_2D()
{
    delete [] format;
    delete [] *kywrds;
    delete [] kywrds;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_2D::init()
{
    
    format=NULL;
    Var::append(format,"%s[%%d][%%d]",name);
    kywrds=new KeyWord*[dim0];
    *kywrds=new KeyWord[dim0*dim1];
    for(int i=1;i<dim0;i++)
        kywrds[i]=kywrds[i-1]+dim1;
    
    head_kywrd=*kywrds;
    nkywrds=dim0*dim1;
    rank=2;
    
    var_adj();
    var_adj();
    vars[0]=&comp0;
    vars[1]=&comp1;
    const int dim0_=dim0;
    const int dim1_=dim1;
    var_logics[0]=Logics(NULL,"ge",var_mngr(0))*Logics(NULL,"lt",var_mngr(dim0_));
    var_logics[1]=Logics(NULL,"ge",var_mngr(0))*Logics(NULL,"lt",var_mngr(dim1_));
    
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
void Pattern::SubPattern_2D::init_keywords()
{
    for(int i=0;i<dim0;i++)
        for(int j=0;j<dim1;j++)
            kywrds[i][j].init(&cmd_logic,var_mngr,format,i,j);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_2D::init_keywords(bool**& set)
{
    for(int i=0;i<dim0;i++)
        for(int j=0;j<dim1;j++)
            kywrds[i][j].init(&cmd_logic,var_mngr,set[i][j],format,i,j);
}

/*--------------------------------------------
 
 --------------------------------------------*/
bool Pattern::SubPattern_2D::scan(char*& err_msg,char**& args,int& nargs)
{

    if(sscanf(*args,format,&icomp0,&icomp1)!=2)
        return false;
    err_msg=vlog_chek(kywrds[icomp0][icomp1],args,nargs);
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_2D::conv_name(PrintStyle& ps)
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
/*---------------------------------------------------------------------------------------------------------------------------
 _____   _   _   _____   _____       ___   _____   _____   _____   _____    __   _        _____   _____        _       _____
/  ___/ | | | | |  _  \ |  _  \     /   | |_   _| |_   _| | ____| |  _  \  |  \ | |      /___  \ |  _  \      | |     |_   _|
| |___  | | | | | |_| | | |_| |    / /| |   | |     | |   | |__   | |_| |  |   \| |       ___| | | | | |      | |       | |
\___  \ | | | | |  _  { |  ___/   / / | |   | |     | |   |  __|  |  _  /  | |\   |      /  ___/ | | | |      | |       | |
 ___| | | |_| | | |_| | | |      / /  | |   | |     | |   | |___  | | \ \  | | \  |      | |___  | |_| |      | |___    | |
/_____/ \_____/ |_____/ |_|     /_/   |_|   |_|     |_|   |_____| |_|  \_\ |_|  \_|      |_____| |_____/      |_____|   |_|
 ---------------------------------------------------------------------------------------------------------------------------*/
Pattern::SubPattern_2D_LT::SubPattern_2D_LT(VarManager& var_mngr_,const char* name_,int dim_):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.add_var(icomp0,"")),
comp1(*var_mngr.add_var(icomp1,"")),
dim(dim_)
{
    init();
    init_keywords();
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::SubPattern_2D_LT::SubPattern_2D_LT(bool**& set,VarManager& var_mngr_,const char* name_,int dim_):
SubPattern(var_mngr_,name_),
comp0(*var_mngr.add_var(icomp0,"")),
comp1(*var_mngr.add_var(icomp1,"")),
dim(dim_)
{
    init();
    init_keywords(set);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Pattern::SubPattern_2D_LT::~SubPattern_2D_LT()
{
    delete [] format;
    delete [] *kywrds;
    delete [] kywrds;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_2D_LT::init()
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
    var_logics[0]=Logics(NULL,"ge",var_mngr(0))*Logics(NULL,"lt",var_mngr(dim_));
    var_logics[1]=Logics(NULL,"ge",var_mngr(0))*Logics(NULL,"le",vars[0]);
    
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
void Pattern::SubPattern_2D_LT::init_keywords()
{
    for(int i=0;i<dim;i++)
        for(int j=0;j<i+1;j++)
            kywrds[i][j].init(&cmd_logic,var_mngr,format,i,j);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_2D_LT::init_keywords(bool**& set)
{
    for(int i=0;i<dim;i++)
        for(int j=0;j<i+1;j++)
            kywrds[i][j].init(&cmd_logic,var_mngr,set[i][j],format,i,j);
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Pattern::SubPattern_2D_LT::scan(char*& err_msg,char**& args,int& nargs)
{
    if(sscanf(*args,format,&icomp0,&icomp1)!=2)
        return false;
    
    err_msg=vlog_chek(kywrds[icomp0][icomp1],args,nargs);
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Pattern::SubPattern_2D_LT::conv_name(PrintStyle& ps)
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
/*------------------------------------------------------------------------------
 _____   _   _       _____   _____    _____       ___   _____   _____   _____   
|  ___| | | | |     | ____| |  _  \  | ____|     /   | |  _  \ | ____| |  _  \  
| |__   | | | |     | |__   | |_| |  | |__      / /| | | | | | | |__   | |_| |  
|  __|  | | | |     |  __|  |  _  /  |  __|    / / | | | | | | |  __|  |  _  /  
| |     | | | |___  | |___  | | \ \  | |___   / /  | | | |_| | | |___  | | \ \  
|_|     |_| |_____| |_____| |_|  \_\ |_____| /_/   |_| |_____/ |_____| |_|  \_\
 ------------------------------------------------------------------------------*/
FileReader::FileReader():
varmngr(),
ntypes(atom_types->no_types)
{
    file_name=NULL;
    type_ref=NULL;
    ntype_ref=0;
    tens=NULL;
    ntens=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
FileReader::~FileReader()
{
    delete [] file_name;
    delete [] type_ref;
    for(int i=0;i<ntens;i++)
        delete tens[i];
    delete [] tens;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void FileReader::open_file()
{
    fp=NULL;
    int chk=1;
    if(atoms->my_p==0)
    {
        fp=fopen(file_name,"r");
        if(fp==NULL)
            chk=0;
    }
    MPI_Bcast(&chk,1,MPI_INT,0,world);
    if(!chk)
        error->abort("file %s not found",file_name);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void FileReader::close_file()
{
    if(atoms->my_p==0)
        fclose(fp);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void FileReader::read_header()
{
    
    char* line =new char[MAXCHAR];
    char** args=NULL;
    int nargs=0;
    int args_cpcty=0;
    
    while(nargs==0 && read_line(line)!=-1)
        nargs=ScriptReader::parse_line(line,args,args_cpcty);
    
    if(nargs==0)
        error->abort("%s file ended immaturely",file_name);
    
    if(nargs<ntypes)
        error->abort("the number of atoms in %s file"
        " is less than the number of atom types present in the system",file_name);
    
    ntype_ref=nargs;
    
    delete [] type_ref;
    type_ref=new int[ntype_ref];
    
    for(int ityp=0;ityp<ntype_ref;ityp++)
        type_ref[ityp]=atom_types->find_type_exist(args[ityp]);
    
    delete [] args;
    delete [] line;
    
    bool* tmp_type;
    tmp_type=new bool[atom_types->no_types];
    for(int i=0;i<atom_types->no_types;i++)
        tmp_type[i]=false;
    
    for(int i=0;i<ntype_ref;i++)
        if(type_ref[i]!=-1)
            tmp_type[type_ref[i]]=true;
    
    for(int i=0;i<atom_types->no_types;i++)
        if(!tmp_type[i])
            error->abort("element %s was not found",atom_types->atom_names[i]);
    
    for(int i=0;i<ntype_ref;i++)
        for(int j=i+1;j<ntype_ref;j++)
            if(type_ref[i]!=-1 && type_ref[i]==type_ref[j])
                error->abort("duplicate element %s in header",atom_types->atom_names[type_ref[i]]);
    
    varmngr.add_var(ntype_ref,"no. of elements in header");
    for(int i=0;i<ntens;i++)
        tens[i]->set_max_cmp(ntype_ref);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void FileReader::read_body()
{
    char** args=NULL;
    int nargs=0;
    
    char* line =new char[MAXCHAR];
    char* buff =new char[MAXCHAR];
    char* err_msg=NULL;
    while(read_line(line)!=-1)
    {
        memcpy(buff,line,MAXCHAR*sizeof(char));
        err_msg=qr(line,args,nargs);
        if(err_msg!=NULL)
            error->abort("%s",err_msg);
        if(nargs==0)
            continue;
        
        
        bool found=false;
        for(int i=0;i<ntens && !found;i++)
        {
            found=tens[i]->scan(err_msg,args,nargs);
            if(err_msg!=NULL)
                error->abort("%s, in file %s, line %s",err_msg,file_name,buff);
        }
            
        if(!found)
            error->abort("candidate keyword in file %s for %s was not found",file_name,line);
    }
    
    delete [] line;
    delete [] buff;
    delete [] args;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void FileReader::read_file(const char* file_name_)
{
    size_t len=strlen(file_name_)+1;
    delete [] file_name;
    file_name=new char[len];
    memcpy(file_name,file_name_,len);
    
    open_file();
    read_header();
    read_body();
    close_file();

    for(int i=0;i<ntens;i++)
        for(int j=0;j<tens[i]->nkywrds;j++)
        {
            char* msg=tens[i]->head[j].finalize();
            if(msg!=NULL) error->abort("%s",msg);
            delete [] msg;
        }
}
/*--------------------------------------------
 hash sign remover:
 it removes the statement after hash sign,
 replaces the withe spaces with space;
 for example:
 line=
 "  this  is a    test, #here is comment";
 will be converted to:
 newline=
 "this is a test,";
 nargs=4;
 --------------------------------------------*/
int FileReader::hash_remover(char*& line)
{
    char* last;
    last=line;
    while(*last!='#' && *last!='\0' && *last!='\n')
        last++;
    *last='\0';
    
    char* ipos=line;
    char* jpos=line;
    char* kpos=line;
    int nargs=0;
    while(ipos!=last)
    {
        while (isspace(*ipos) && ipos!=last)
            ipos++;
        if(ipos==last)
        {
            if(kpos!=line)
                kpos--;
            *kpos='\0';
            continue;
        }
        
        jpos=ipos;
        nargs++;
        
        while (!isspace(*ipos) && ipos!=last)
            ipos++;
        memmove(kpos,jpos,ipos-jpos);
        if(ipos==last)
            continue;
        
        kpos+=ipos-jpos;
        *kpos=' ';
        kpos++;
    }
    return nargs;
}
/*--------------------------------------------
 hash sign remover:
 it removes the statement after hash sign,
 replaces the withe spaces with space;
 for example:
 line=
 "  this  is a    test, #here is comment";
 will be converted to:
 newline=
 "this is a test,";
 nargs=4;
 --------------------------------------------*/
int FileReader::parse_line(char* ipos,char**& args,int& args_cpcty)
{
    char* last;
    last=ipos;
    while(*last!='#' && *last!='\0' && *last!='\n')
        last++;
    *last='\0';
    
    int nargs=0;
    while(ipos!=last)
    {
        while (isspace(*ipos) && ipos!=last)
            ipos++;
        if(ipos==last)
            continue;
        
        if(nargs+1>args_cpcty)
        {
            char** args_=new char*[nargs+1];
            memcpy(args_,args,nargs*sizeof(char*));
            delete [] args;
            args=args_;
            args_cpcty++;
        }
        args[nargs]=ipos;
        nargs++;
        
        while (!isspace(*ipos) && ipos!=last)
            ipos++;
        if(ipos==last)
            continue;
        *ipos='\0';
        ipos++;
    }
    
    return nargs;
}
/*--------------------------------------------
 read line and broadcast
 --------------------------------------------*/
int FileReader::read_line(char*& line)
{
    int lenght=0;
    if(atoms->my_p==0)
    {
        if(fgets(line,MAXCHAR,fp)==NULL)
            lenght=-1;
    }
    
    MPI_Bcast(&lenght,1,MPI_INT,0,world);
    if(lenght==-1)
        return -1;
    
    if(atoms->my_p==0)
    {
        lenght=static_cast<int>(strlen(line))+1;
    }
    MPI_Bcast(&lenght,1,MPI_INT,0,world);
    MPI_Bcast(line,lenght,MPI_CHAR,0,world);

    char* p=strchr(line,'\n');
    if(p!=NULL) *p='\0';
    
    return lenght;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void FileReader::symmetric()
{
    tens[ntens-1]->symmetric=true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
FileReader::Ten::Ten(FileReader* fr_, const char* name_):
fr(fr_),
varmngr(fr_->varmngr)
{
    size_t len=strlen(name_)+1;
    name=new char[len];
    memcpy(name,name_,len);
    format=NULL;
    symmetric=false;
    head=NULL;
    nkywrds=0;
    cmd_logic=Logics(NULL,"set");
}
/*--------------------------------------------
 
 --------------------------------------------*/
FileReader::Ten::~Ten()
{
    delete [] name;
    delete [] format;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& FileReader::Ten::add_vlog()
{
    return var_logics;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Logics& FileReader::Ten::add_clog()
{
    return cmd_logic;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool FileReader::QuantityReader::is_del(char* c)
{
    for(size_t i=0;i<ndels;i++)
        if(*c==dels[i])
            return true;
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool FileReader::QuantityReader::is_space(char* c)
{
    for(size_t i=0;i<nspaces;i++)
        if(*c==spaces[i])
            return true;
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool FileReader::QuantityReader::is_null(char* c)
{
    if(c==NULL ||*c=='\0')
        return true;
    
    for(size_t i=0;i<nnulls;i++)
        if(*c==nulls[i])
            return true;
    
    return false;
}

/*--------------------------------------------
 
 --------------------------------------------*/
bool FileReader::QuantityReader::is_char(char* c)
{
    if(is_null(c))
        return false;
    if(is_del(c))
        return false;
    if(is_space(c))
        return false;
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool FileReader::QuantityReader::syntax(char* c)
{
    if(first==NULL)
    {
        first=c;
        return true;
    }
    else if(*c=='(')
    {
        if(c==first)
            return false;
        if(lp!=NULL)
            return false;
        if(eq!=NULL)
            return false;
        lp=c;
        return true;
    }
    else if(*c==')')
    {
        
        if(rp!=NULL)
            return false;
        if(lp==NULL)
            return false;
        if(eq!=NULL)
            return false;
        if(com==NULL && c==lp+1)
            return false;
        if(com!=NULL && c==com+1)
            return false;
        rp=c;
        return true;
    }
    else if(*c==',')
    {
        
        if(lp==NULL)
            return false;
        if(rp!=NULL)
            return false;
        if(com!=NULL && c==com+1)
            return false;
        if(com==NULL && c==lp+1)
            return false;
        com=c;
        return true;
    }
    else if(*c=='=')
    {
        if(c==first)
            return false;
        if(lp!=NULL && rp==NULL)
            return false;
        eq=c;
        return true;
    }
    else if(*c=='\0')
    {
        if(eq==NULL)
            return false;
        return true;
    }
    else if (is_space(c))
    {
        return false;
    }
    
    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int FileReader::read_line(char*& line,int& line_cpcty,int& chunk)
{
    int ipos=0;
    bool file_cmplt=false;
    
    if(atoms->my_p==0)
    {
        char* mark;
        bool cmd_cmplt=false;
        bool line_cmplt=false;
        
        while(!cmd_cmplt && !file_cmplt)
        {
            line_cmplt=false;
            while(!line_cmplt && !file_cmplt)
            {
                if(ipos+chunk>line_cpcty)
                {
                    char* line_=new char[ipos+chunk];
                    memcpy(line_,line,ipos);
                    delete [] line;
                    line=line_;
                    line_cpcty=ipos+chunk;
                }
                
                if(fgets(line+ipos,chunk,fp)==NULL)
                {
                    file_cmplt=true;
                    continue;
                }
                
                if(strchr(line+ipos,'\n')!=NULL)
                    line_cmplt=true;
                else
                {
                    ipos=static_cast<int>(strlen(line));
                    chunk*=2;
                }
            }
            
            if(file_cmplt)
                continue;
            mark=strchr(line+ipos,'\\');
            if(mark==NULL)
                cmd_cmplt=true;
            else
                ipos=static_cast<int>(mark-line);
        }
        
        if(!file_cmplt)
        {
            ipos=static_cast<int>(strlen(line));
            if(ipos && line[ipos-1]=='\n')
            {
                line[ipos-1]='\0';
                ipos--;
            }
        }
    }
    
    MPI_Bcast(&file_cmplt,1,MPI_BYTE,0,world);
    
    MPI_Bcast(&ipos,1,MPI_INT,0,world);
    
    if(file_cmplt && ipos)
        error->abort("file ended unexpectedly");
    
    if(file_cmplt)
        return -1;
    
    if(ipos+1>line_cpcty)
    {
        
        
        if(line_cpcty)
            delete [] line;
        
        line =new char[ipos+1];
        
        line_cpcty=ipos+1;
        
    }
    
    MPI_Bcast(line,ipos+1,MPI_CHAR,0,world);
    
    return ipos+1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
FileReader::QuantityReader::QuantityReader():
spaces(" \t\n"),nspaces(strlen(spaces)),
dels("(,)="),ndels(strlen(dels)),
nulls("#"),nnulls(strlen(nulls))
{}
/*--------------------------------------------
 
 --------------------------------------------*/
char* FileReader::QuantityReader::operator()(char* line,char**& args,int& nargs)
{
    first=lp=rp=eq=com=NULL;
    char* p1=line;
    char* p0=NULL;
    bool stx_chk=true;
    size_t len=0;
    int n_dels=0;
    while(!is_null(p1) && stx_chk)
    {
        if(!is_space(p1) || (is_char(p0) && is_char(p1+1)))
        {
            stx_chk=syntax(p1);
            if(is_del(p1)) n_dels++;
            len++;
        }
        p0=p1;
        p1++;
        
        if(is_null(p1))
            stx_chk=syntax(p1);
    }
    
    
    if(!stx_chk)
    {
        char* err_msg=NULL;
        Var::append(err_msg,"syntax error");
        return err_msg;
    }
    
    if(n_dels==1)
        n_dels=2;
    
    delete [] args;
    nargs=n_dels;
    if(nargs)
        args=new char*[nargs];
    else
    {
        args=NULL;
        return NULL;
    }
    
    
    char** iargs=args;
    p1=line;
    p0=NULL;
    while(!is_null(p1))
    {
        if(is_char(p0) && !is_char(p1))
            *p1='\0';
        if(!is_char(p0) && is_char(p1))
        {
            *iargs=p1;
            iargs++;
        }
        p0=p1;
        p1++;
    }

    char* p=args[nargs-1];
    for(int i=nargs-1;i>1;i--)
        args[i]=args[i-1];
    args[1]=p;
    
    return NULL;
}


