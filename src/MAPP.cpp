#include "MAPP.h"
#include "script_reader.h"
#include "atoms.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "atom_types.h"
#include "neighbor.h"
#include "neighbor_md.h"
#include "group.h"


#include "command_styles.h"
#include "md_styles.h"
#include "write_styles.h"

#include <iostream>
//#include <fstream>
#define MAPP_VERSION "Beta"
namespace MAPP_NS
{
    int step_no=0;
    int step_tally=1000;
    int mode=MD_mode;
    int precision=6;
    MPI_Comm world(MPI_COMM_WORLD);
    Neighbor* neighbor(NULL);
    ForceField* forcefield(NULL);
    Error* error(NULL);
    Memory* memory(NULL);
    Timer* timer(NULL);
    AtomTypes* atom_types(NULL);
    FILE* output(NULL);
    
    Write* write(NULL);
    MD* md(NULL);
    DMD* dmd(NULL);
    GroupCollection* groups(NULL);
    LineSearch<Min>* ls(NULL);
    Atoms<MAPP_NS::dimension>* atoms(NULL);
    MAPP* mapp(NULL);
    VarManager* g_vm(NULL);

}
using namespace MAPP_NS;
/*--------------------------------------------
 constructor of the main executer
 --------------------------------------------*/
MAPP::
MAPP(int nargs,char** args)
{
    
    delete mapp;
    mapp=this;
    atoms=new Atoms<MAPP_NS::dimension>();
    memory=new Memory();
    error=new Error();
    timer=new Timer();
    atom_types = new AtomTypes();
    neighbor=new Neighbor_md();
    groups=new GroupCollection();
    

    
    type=NULL;
    ctype=NULL;
    x_d=NULL;
    c=NULL;
    c_d=NULL;
    x_dof=NULL;
    c_dof=NULL;
    
    
    init_debug();
    
    input_file=stdin;
    output=stdout;

    if(atoms->my_p==0)
        fprintf(output,"MAPP Version: %s\n",
        (char*) MAPP_VERSION);
    
    
    
    
    no_commands=0;
    
    int iarg=1;
    while(iarg<nargs)
    {
        
        if(strcmp(args[iarg],"-i")==0)
        {
            iarg++;
            if(iarg==nargs)
                error->abort("no input file");
            ScriptReader::open_file(input_file,args[iarg],"r");
            iarg++;
        }
        else if(strcmp(args[iarg],"-o")==0)
        {
            iarg++;
            if(iarg==nargs)
                error->abort("no output file");
            ScriptReader::open_file(output,args[iarg],"w");
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
    fin_debug();
}
/*--------------------------------------------
 analysing the commands
 --------------------------------------------*/
void MAPP::read_file()
{
    char** args=NULL;
    int args_cpcty=0;
    int nargs=0;
    ScriptReader sr(input_file);
    while(nargs>=0)
    {
        nargs=sr(args,args_cpcty);
        if(nargs<=0)
            continue;
        command_style(nargs,args);
        no_commands++;
    }
    delete [] args;
    if(nargs==-2)
        error->abort("file ended unexpectedly");
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
    new class_name(nargs,args);                 \
    delete command;}
    
    if(strcmp(args[0],"rm")==0)
    {
        char* cmd_line;
        int lngth=ScriptReader::concatenate(nargs,args,cmd_line);
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

 --------------------------------------------*/
#ifdef MAPP_DEBUG
#include <unistd.h>
#endif
void MAPP::init_debug()
{
#ifdef MAPP_DEBUG
    if(!i)
    {
        my_debug=NULL;
        return;
    }

    char* filename=new char[MAXCHAR];
    sprintf(filename,"debug-%d",atoms->my_p);
    my_debug=fopen(filename,"w");
    delete [] filename;
    fprintf(my_debug,"my pid is %d\n",getpid());
#endif
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MAPP::fin_debug()
{
#ifdef MAPP_DEBUG
    fclose(my_debug);
#endif
}
/*--------------------------------------------
 if mass unit is amu
 energy unit is eV
 t unit would be 10.1805 fs
 boltzmann 8.617332478e-5
 planck 0.40623423822503807
 Commmand cmd
 
 command: md
 
 list of styles: nh
 
 Cmd St_0 St_1 ... {cmd^i st^i_0 st^i_1 ...} {cmd^j st^j_0 st^j_1 ...}
 
 for example for md
 cmd_name: md
 style_0: nh
 style_1: ntaut|npt|navt
 0
 1
 0,0
 --------------------------------------------*/
void MAPP::test0()
{
}
