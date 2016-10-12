#include "script_reader.h"
#include "macros.h"
#include "global.h"
#include "error.h"
#include "atoms.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
ScriptReader::ScriptReader(FILE*& fp_):
fp(fp_)
{
    sz_inc=100;
    line_len=MAXCHAR;
    line=new char[line_len];
}
/*--------------------------------------------
 
 --------------------------------------------*/
ScriptReader::~ScriptReader()
{
    delete [] line;
}
/*--------------------------------------------
 3 status:
 line finished
 command finished;
 file finished;
 --------------------------------------------*/
bool ScriptReader::get_line_proc0(char*& buff,int& buff_sz,int& buff_cpcty)
{
    int beg=buff_sz;
    int cpt_sz;
    char* p;
    bool line_fin;
    line_fin=false;
    
    int max_cpt_sz=buff_cpcty-buff_sz;
    p=buff+buff_sz;
    
    while(!line_fin)
    {
        p=fgets(p,max_cpt_sz,fp);
        
        if(p==NULL)
            return false;
        
        cpt_sz=static_cast<int>(strlen(p));
        buff_sz+=cpt_sz;

        // check if the line is fully captured
        if(p[cpt_sz-1]=='\n')
        {
            p[cpt_sz-1]='\0';
            buff_sz--;
            line_fin=true;
            continue;
        }
        
        
        
        char* buff_=new char[buff_cpcty+sz_inc];
        memcpy(buff_,buff,buff_cpcty*sizeof(char));
        delete [] buff;
        buff=buff_;
        
        p=buff+buff_cpcty-1;
        buff_cpcty+=sz_inc;
        max_cpt_sz=sz_inc+1;
    }
    
    char* srch=strchr(buff+beg,'#');
    if(srch!=NULL)
    {
        *srch='\0';
        buff_sz=static_cast<int>(srch-buff);
    }
    
    srch=strchr(buff+beg,'\\');
    if(srch!=NULL)
    {
        *srch='\0';
        buff_sz=static_cast<int>(srch-buff);
        return get_line_proc0(buff,buff_sz,buff_cpcty);
    }

    return true;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int ScriptReader::operator()(char**& args,int& args_cpcty)
{
    int buff_sz=0;
    if(atoms->my_p==0)
    {
        bool chk=get_line_proc0(line,buff_sz,line_len);
        if(!chk)
        {
            if(buff_sz!=0)
                buff_sz=-2;
            else
                buff_sz=-1;
        }
    }
    
    MPI_Bcast(&buff_sz,1,MPI_INT,0,world);
    if (buff_sz<0)
        return buff_sz;
    
    if(buff_sz+1>line_len)
    {
        delete [] line;
        line=new char[buff_sz+1];
        line_len=buff_sz+1;
    }
    MPI_Bcast(line,buff_sz+1,MPI_CHAR,0,world);
    
    int nargs=0;
    char* p1=line;
    char* p0=NULL;
    while(*p1!='\0')
    {
        if(is_char(p0) && !is_char(p1))
            *p1='\0';
        if(!is_char(p0) && is_char(p1))
        {
            if(nargs+1>args_cpcty)
            {
                char** args_=new char*[nargs+1];
                memcpy(args_,args,nargs*sizeof(char*));
                delete [] args;
                args=args_;
                args_cpcty=nargs+1;
            }
            args[nargs++]=p1;
        }
        p0=p1;
        p1++;
    }
    

    return nargs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool ScriptReader::is_char(char* c)
{
    if(c==NULL)
        return false;
    if(*c=='\0')
        return false;
    if(*c==' ')
        return false;
    if(*c=='\t')
        return false;
    
    return true;
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
int ScriptReader::hash_remover(char* line,char*& newline)
{
    newline = NULL;
    newline=new char[1024];
    
    int nargs=0;
    int cursor=0;
    int icursor=0;
    while (cursor < strlen(line)&& line[cursor]!='#')
    {
        if(isspace(line[cursor]))
            cursor++;
        else
        {   if(nargs!=0)
            newline[icursor++]=' ';
            while(!isspace(line[cursor])
                  && cursor < strlen(line))
                newline[icursor++]=line[cursor++];
            nargs++;
        }
    }
    newline[icursor]='\0';
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
int ScriptReader::hash_remover(char*& line)
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
 parse a command line:
 chops up a 1d array line into 2d array of
 char, also returns the number of arguments
 if the line starts with # it returns 0
 for example:
 line=
 "  this  is a    test, #here is comment";
 will be converted to:
 args[0]="this";
 args[1]="is";
 args[2]="a";
 args[3]="test,";
 nargs=4;
 --------------------------------------------*/
int ScriptReader::parse_line(const char* line,char**& args)
{
    args=NULL;
    int nargs = 0;
    int cursor = 0;
    int length= static_cast<int>(strlen(line));
    int ll[length];
    
    while (cursor < length
           && line[cursor]!='#'){
        if(isspace(line[cursor]))
            cursor++;
        else
        {
            if(nargs==0 && line[cursor]=='#')
                return 0;
            int i = 0;
            while(!isspace(line[cursor])
                  && cursor < length)
            {
                cursor++;
                i++;
            }
            ll[nargs]=i+1;
            nargs++;
        }
    }
    
    if(nargs==0) return 0;
    
    args=new char*[nargs];
    
    for(int i=0;i<nargs;i++)
        args[i]=new char[ll[i]];

    
    nargs = 0;
    cursor = 0;
    
    while (cursor<length&&line[cursor]!='#')
    {
        if( isspace(line[cursor]))
            cursor++;
        else
        {
            int i=0;
            while(!isspace(line[cursor])
                  && cursor < strlen(line))
                args[nargs][i++]=line[cursor++];
            args[nargs][i] = '\0';
            nargs++;
        }
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
int ScriptReader::parse_line(char* ipos,char**& args,int& args_cpcty)
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
int ScriptReader::concatenate(int nargs,char** args
,char*& line)
{
    int lngth=0;
    int pos=0;
    for(int i=0;i<nargs;i++)
        lngth+=static_cast<int>(strlen(args[i]));
    
    lngth+=nargs;
    line=new char[lngth];
    for(int i=0;i<nargs;i++)
    {
        lngth=static_cast<int>(strlen(args[i]));
        for(int j=0;j<lngth;j++)
        {
            line[pos]=args[i][j];
            pos++;
        }
        
        if(i!=nargs-1)
        {
            line[pos]=' ';
            pos++;
        }
        else
        {
            line[pos]='\0';
            pos++;
        }
    }
    return lngth;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ScriptReader::open_file(FILE*& fp,const char* file_name,const char* id)
{
    fp=NULL;
    int chk=1;
    if(atoms->my_p==0)
    {
        fp=fopen(file_name,id);
        if(fp==NULL)
            chk=0;
    }
    MPI_Bcast(&chk,1,MPI_INT,0,world);
    if(!chk)
        error->abort("file %s not found",file_name);
}
/*--------------------------------------------
 read line and broadcast
 --------------------------------------------*/
int ScriptReader::read_line(FILE* file,char*& line)
{
    int lenght=0;
    if(atoms->my_p==0)
    {
        if(fgets(line,MAXCHAR,file)==NULL)
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
    
    return lenght;
}
/*--------------------------------------------
 read line and broadcast
 --------------------------------------------*/
int ScriptReader::read_line(FILE* fp,char*& line,int& line_cpcty,int& chunk)
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
                    if(line_cpcty)
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

