/*--------------------------------------------
 Created by Sina on 05/12/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef MAPP_init_h
#define MAPP_init_h
#include "type_def.h"
#include <mpi.h>
namespace MAPP_NS
{
    enum {MD_mode,DMD_mode};

    template<const int> class Atoms;
    template<typename> class Vec;
    template<class Func> class LineSearch;
    class Min;

    constexpr int dimension=3;
    extern int mode;
    extern int precision;
    extern class Neighbor* neighbor;
    extern class ForceField* forcefield;
    extern class Error* error;
    extern class Timer* timer;
    extern FILE* output;
    extern MPI_Comm world;
    extern class Write* write;
    extern class MD* md;
    extern class DMD* dmd;
    extern class GroupCollection* groups;
    extern class Memory* memory;
    extern LineSearch<Min>* ls;
    extern int step_no;
    extern int step_tally;
    extern class AtomTypes* atom_types;
    extern class Atoms<MAPP_NS::dimension>* atoms;
    extern class MAPP* mapp;
    extern class VarManager* g_vm;
}
#endif








