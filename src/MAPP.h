/*--------------------------------------------
 Created by Sina on 05/11/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/

#ifndef __MAPP__MAPP__
#define __MAPP__MAPP__
#include "type_def.h"
#include <stdio.h>
namespace MAPP_NS
{
    template<typename> class Vec;
    class MAPP
    {
    private:
    protected:
    public:
        MAPP(int,char**);
        ~MAPP();
        int no_commands;
        Vec<atom_type>* type;
        Vec<atom_type>* ctype;
        Vec<type0>* x_d;
        Vec<type0>* c;
        Vec<type0>* c_d;
        Vec<bool>* x_dof;
        Vec<bool>* c_dof;        
        
        FILE* input_file;
        FILE* my_debug;
        
        class Timer* timer;
        class GroupCollection* groups;
        class AtomTypes* atom_types;
        
        
        void read_file();
        void command_style(int,char**);

        void init_debug();
        void fin_debug();
        void test0();
    };
}


#endif
