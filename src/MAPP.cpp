/*--------------------------------------------
 Created by Sina on 05/11/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/

#include "MAPP.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "atom_types.h"
#include "neighbor.h"
#include "neighbor_md.h"
#include "thermo_dynamics.h"
#include "group.h"

#include "read_styles.h"
#include "ff_styles.h"
#include "min_styles.h"
#include "command_styles.h"
#include "md_styles.h"
#include "ls_styles.h"
#include "write_styles.h"

#include <iostream>
#include <fstream>



#include "gmres.h"


#define MAPP_VERSION "Beta"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor of the main executer
 --------------------------------------------*/
MAPP::
MAPP(int nargs,char** args,MPI_Comm communicator):
world(communicator),
atoms(new Atoms(this,3)),
x(atoms->x),
id(atoms->id)
{
    init_dubeg(true);
    
    memory=new Memory(this);
    error=new Error(this);
    timer=new Timer(this);
    atom_types = new AtomTypes(this);
    neighbor=new Neighbor_md(this);
    groups=new GroupCollection(this);
    
    forcefield=NULL;
    write=NULL;
    md=NULL;
    //min=NULL;
    dmd=NULL;

    
    type=NULL;
    ctype=NULL;
    x_d=NULL;
    c=NULL;
    c_d=NULL;
    x_dof=NULL;
    c_dof=NULL;
    
    ls=NULL;

    
    input_file=stdin;
    output=stdout;

    if(atoms->my_p==0)
        fprintf(output,"MAPP Version: %s\n",
        (char*) MAPP_VERSION);
    
    
    mode=MD_mode;
    step_no=0;
    precision=6;
    step_tally=1000;
    no_commands=0;
    
    int iarg=1;
    while(iarg<nargs)
    {
        
        if(strcmp(args[iarg],"-i")==0)
        {
            iarg++;
            if(iarg==nargs)
                error->abort("no input file");
            open_file(input_file,args[iarg],"r");
            iarg++;
        }
        else if(strcmp(args[iarg],"-o")==0)
        {
            iarg++;
            if(iarg==nargs)
                error->abort("no output file");
            open_file(output,args[iarg],"w");
            iarg++;
        }
        else
            error->abort("unknown postfix: %s",args[iarg]);
    }
    
    if(input_file==NULL)
        error->abort("input file not found");
    
    read_file();

    
    test0();
    
    if(input_file!=stdin)
        fclose(input_file);
}
/*--------------------------------------------
 destructor of the main executer
 --------------------------------------------*/
MAPP::~MAPP()
{
    if(atoms->my_p==0)
      fprintf(output,"Finito\n");
    
    if(output!=stdout)
        fclose(output);
    
    
    delete forcefield;
    delete write;
    delete md;
    
    delete groups;
    delete neighbor;
    delete atom_types;
    delete timer;
    delete memory;
    delete error;
    delete atoms;
    fin_dubeg();
}
/*--------------------------------------------
 analysing the commands
 --------------------------------------------*/
void MAPP::read_file()
{
    char* cmd;
    int cmd_cpcty=0;
    int chunk=MAXCHAR;
    char** args=NULL;
    int nargs;
    int args_cpcty=0;

    while(read_line(input_file,cmd,cmd_cpcty,chunk)!=-1)
    {
        nargs=parse_line(cmd,args,args_cpcty);
        if(nargs==0) continue;
        command_style(nargs,args);
        no_commands++;
    }
    
    if(args_cpcty) delete [] args;
    if(cmd_cpcty) delete [] cmd;
}
/*--------------------------------------------
 differnt command styles
 --------------------------------------------*/
void MAPP::command_style(int nargs,char** args)
{
    
    #define Command_Style
    #define CommandStyle(class_name,style_name) \
    else if(strcmp(args[0],#style_name)==0){    \
    class class_name* command =                 \
    new class_name(this,nargs,args);            \
    delete command;}
    
    if(strcmp(args[0],"rm")==0)
    {
        char* cmd_line;
        int lngth=concatenate(nargs,args,cmd_line);
        int shell_cmd_chk=0;
        if(atoms->my_p==0)
            shell_cmd_chk=system(cmd_line);
        MPI_Bcast(&shell_cmd_chk,1,MPI_INT,0,world);
        if(lngth) delete [] cmd_line;
    }
    #include "command_styles.h"
    else
        error->abort("unknown command:"
                     " %s",args[0]);
    #undef CommandStyle
    #undef Command_Style
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
int MAPP::parse_line(const char* line,char**& args)
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
    
    CREATE1D(args,nargs);
    
    for(int i=0;i<nargs;i++)
        CREATE1D(args[i],ll[i]);
    
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
int MAPP::hash_remover(char* line,char*& newline)
{
    newline = NULL;
    CREATE1D(newline,MAXCHAR);
    
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
int MAPP::hash_remover(char*& line)
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
int MAPP::parse_line(char* ipos,char**& args,int& args_cpcty)
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
            GROW(args,nargs,nargs+1);
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
int MAPP::concatenate(int nargs,char** args
,char*& line)
{
    int lngth=0;
    int pos=0;
    for(int i=0;i<nargs;i++)
        lngth+=static_cast<int>(strlen(args[i]));
    
    lngth+=nargs;
    CREATE1D(line,lngth);
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
 read line and broadcast
 --------------------------------------------*/
int MAPP::read_line(FILE* file,char*& line)
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
int MAPP::read_line(FILE* fp,char*& line,int& line_cpcty,int& chunk)
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
/*--------------------------------------------
 
 --------------------------------------------*/
void MAPP::open_file(FILE*& fp,const char* file_name,const char* id)
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

 --------------------------------------------*/
void MAPP::init_dubeg(bool i)
{
    if(!i)
    {
        my_debug=NULL;
        return;
    }

    char* filename=new char[MAXCHAR];
    //sprintf(filename,"debug-%d-%d",getpid(),atoms->my_p);
    sprintf(filename,"debug-%d",atoms->my_p);
    my_debug=fopen(filename,"w");
    delete [] filename;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MAPP::fin_dubeg()
{
    if(my_debug!=NULL)
        fclose(my_debug);
}
/*--------------------------------------------
 if mass unit is amu
 energy unit is eV
 t unit would be 10.1805 fs
 --------------------------------------------*/
void MAPP::test0()
{
}
/*--------------------------------------------
 test
 --------------------------------------------*/
void MAPP::test1()
{
    
    /*
    
    const char* fileName ="input_dmd";
    long byte_sz;
    char* buff;
    if(atoms->my_p==0)
    {
        fstream file(fileName,ios::in);
        file.seekg(0,ios::end);
        byte_sz=file.tellg();
    }

    MPI_Bcast(&byte_sz,1,MPI_LONG,0,world);
    byte_sz++;
    CREATE1D(buff,byte_sz);
    
    if(atoms->my_p==0)
    {
        FILE* fp=fopen(fileName,"r");
        fread (buff,1,byte_sz-1,fp);
        buff[byte_sz-1]='\0';
        fclose(fp);
    }
    
    int sz_int=static_cast<int>(byte_sz);
    long sz_long=static_cast<long>(sz_int);
    long byte_sz_=byte_sz;
    char* buff_=buff;
    while(byte_sz_)
    {
        MPI_Bcast(buff_,sz_int,MPI_CHAR,0, world);
        byte_sz_-=sz_long;
        buff_+=sz_int;
        if(byte_sz_<sz_long)
            sz_int=static_cast<int>(byte_sz_);
    }
    
    
    
    char* pch;
    char* new_buff;
    char* new_buff_;
    char* st;
    CREATE1D(new_buff,byte_sz);
    
    buff_=buff;
    new_buff_=new_buff;
    st=buff;
    pch=strchr(buff_,'#');
    while (pch!=NULL)
    {
        memcpy(new_buff_,st,pch-st);
        new_buff_+=pch-st;
        pch=strchr(pch+1,'\n');
        st=pch;
        
        pch=strchr(pch+1,'#');
    }
    memcpy(new_buff_,st,buff+byte_sz-st);
    swap(new_buff,buff);

    buff_=buff;
    new_buff_=new_buff;
    st=buff;
    pch=strchr(buff_,'\\');
    while (pch!=NULL)
    {
        memcpy(new_buff_,st,pch-st);
        new_buff_+=pch-st;
        pch=strchr(pch+1,'\n');
        st=pch+1;
        
        pch=strchr(pch+1,'\\');
    }
    memcpy(new_buff_,st,buff+byte_sz-st);
    

    int buff_length=static_cast<int>(strlen(new_buff)+1);
    if(byte_sz)
        delete [] buff;
    CREATE1D(buff,buff_length);
    memcpy(buff,new_buff,buff_length*sizeof(char));
    if(byte_sz)
        delete [] new_buff;
    
    
    int ncmds=0;
    char** cmds=NULL;
    pch=strtok(buff,"\n");
    while (pch!=NULL)
    {
        GROW(cmds,ncmds,ncmds+1);
        cmds[ncmds]=pch;
        ncmds++;
        pch=strtok(NULL,"\n");
    }
    
    for(int i=0;i<ncmds;i++)
    {
        cout << cmds[i] <<endl;
    }

    
    int* nargs;
    char*** args;
    CREATE1D(nargs,ncmds);
    CREATE1D(args,ncmds);
    for(int icmd=0;icmd<ncmds;icmd++)
    {
        nargs[icmd]=0;
        pch=strtok(cmds[icmd]," ");
        while (pch!=NULL)
        {
            GROW(args[icmd],nargs[icmd],nargs[icmd]+1);
            args[icmd][nargs[icmd]]=pch;
            nargs[icmd]++;
            pch=strtok(NULL," ");
        }
    }


    
    if(ncmds)
        delete [] cmds;
    

     */

}



