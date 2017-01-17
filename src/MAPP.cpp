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
#include "comm.h"
#include "ff.h"

#include "command_styles.h"
#include "md_styles.h"
#include "write_styles.h"

#include <iostream>
#define MAPP_VERSION "Beta"
namespace MAPP_NS
{
    int step_no=0;
    int step_tally=1000;
    int mode=MD_mode;
    int precision=6;
    MPI_Comm __world__(MPI_COMM_WORLD);
    ForceField* forcefield(NULL);
    Error* error(NULL);
    int Error::my_p;
    FILE* output(NULL);
    
    Write* write(NULL);
    MD* md(NULL);
    DMD* dmd(NULL);
    Dynamic* dynamic(NULL);
    LineSearch<Min>* ls(NULL);
    Atoms* atoms(NULL);
    MAPP* mapp(NULL);
    Communication* comm(NULL);
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
    comm=new Communication(__world__);
    atoms=new Atoms(comm);
    Error::my_p=comm->my_p;
    timer=new Timer();
    atom_types = new AtomTypes();
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
                Error::abort("no input file");
            ScriptReader::open_file(input_file,args[iarg],"r");
            iarg++;
        }
        else if(strcmp(args[iarg],"-o")==0)
        {
            iarg++;
            if(iarg==nargs)
                Error::abort("no output file");
            ScriptReader::open_file(output,args[iarg],"w");
            iarg++;
        }
        else
            Error::abort("unknown postfix: %s",args[iarg]);
    }
    
    if(!input_file)
        Error::abort("input file not found");
    
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
    delete atom_types;
    delete timer;
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
        Error::abort("file ended unexpectedly");
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
        MPI_Bcast(&shell_cmd_chk,1,MPI_INT,0,comm->world);
        if(lngth) delete [] cmd_line;
    }
    #include "command_styles.h"
    else
        Error::abort("unknown command:"
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
//#include "logics.h"
//#define mul(A,B) A.print(NULL),B.print(NULL),(A*B).print(NULL)
#include <iomanip>
using namespace std;
void MAPP::test0()
{
    /*
    VarManager vm;
    int l0,l1,r0,r1;
    var<int>L0(l0,"L0");
    var<int>L1(l1,"L1");
    var<int>R0(r0,"R0");
    var<int>R1(r1,"R1");
    
    const char* op1="eq";
    const char* op0="!gt";
    const char* form="case %s\n\t%s\n+\t%s\n=\t%s\n"
    "--------------------------------------------------------------------------\n";
    const char* c1="1: (!L , !O ,  R ) ->   x";
    const char* c2="2: (!L ,  O , !R ) ->   x";
    const char* c3="3: (!L ,  O ,  R ) -> ( L0+L1 , O        , R     )";
    const char* c4="4: ( L , !O , !R ) -> ( L     , O0+R0+O1 , R1    )";
    const char* c5="5: ( L , !O ,  R ) -> ( L     , O0+O1    , R     )";
    const char* c6="6: ( L ,  O , !R ) -> ( L     , O        , R0+R1 )";

    
    printf(form,c4,mul(Logics(l0,op0,r0),Logics(l0,op1,r1)));
    printf(form,c2,mul(Logics(l0,op0,r0),Logics(l1,op0,r1)));
    printf(form,c1,mul(Logics(l0,op0,r0),Logics(l1,op1,r0)));
    printf(form,c3,mul(Logics(l0,op0,r0),Logics(l1,op0,r0)));
    printf(form,c5,mul(Logics(l0,op0,r0),Logics(l0,op1,r0)));
    printf(form,c6,mul(Logics(l0,op0,r0),Logics(l0,op0,r1)));
    
    const char* op2="!eq";
    printf(form,c2,mul(Logics(l0,op0,r0),Logics(l1,op2,r1)));
    printf(form,c3,mul(Logics(l0,op0,r0),Logics(l1,op2,r0)));
    printf(form,c6,mul(Logics(l0,op0,r0),Logics(l0,op2,r1)));
    */
}
