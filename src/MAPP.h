/*--------------------------------------------
 Created by Sina on 05/11/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/

#ifndef __MAPP__MAPP__
#define __MAPP__MAPP__
#include <mpi.h>
#include "atoms.h"
#include "ls.h"
namespace MAPP_NS
{
    enum {MD_mode,DMD_mode};
    
    class MAPP
    {
    private:
    protected:
    public:
        MAPP(int,char**,MPI_Comm);
        ~MAPP();
        MPI_Comm world;
        int mode;
        int step_no;
        int step_tally;
        int no_commands;
        int precision;

        /* 
         * classes that need to be initialized in constructor
         */
        class Atoms* atoms;
        class Memory* memory;
        class Error* error;
        class Timer* timer;
        class AtomTypes* atom_types;
        class Neighbor* neighbor;
        
        /* 
         * classes that might be initialized mid simulation
         */
        class ForceField* forcefield;
        class Write* write;
        class MD* md;
        class Min* min;
        class Clock* clock;


        /* 
         * atomic vectors
         */
        Vec<type0>*& x;
        Vec<int>*& id;
        Vec<md_type>* type;
        Vec<dmd_type>* ctype;
        Vec<type0>* x_d;
        Vec<type0>* c;
        Vec<type0>* c_d;
        Vec<type0>* f;
        Vec<byte>* dof;
        Vec<byte>* cdof;
        
        

        /* 
         * linesearch crap
         */
        LineSearch<Min>* ls;
        
        
        FILE* output;
        FILE* input_file;
        
        void read_file();
        void min_style(int,char**);
        void ff_style(int,char**);
        void md_style(int,char**);
        void ls_style(int,char**);
        void clock_style(int,char**);
        void read_style(int,char**);
        void write_style(int,char**);
        void command_style(int,char**);
        
        //extra functions that can be used by other files
        int read_line(FILE*,char*&);
        int parse_line(const char*,char**&);
        int parse_line(char*,char**&,int&);
        int hash_remover(char*,char*&);
        int hash_remover(char*&);
        int concatenate(int,char**,char*&);
        void test();
        
        int read_line(FILE*,char*&,int&,int&);
        void init_dubeg(int);
        void fin_dubeg();
        FILE* my_debug;
    };
}

#endif
