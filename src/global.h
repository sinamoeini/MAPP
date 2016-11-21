#ifndef MAPP_global_h
#define MAPP_global_h
#include "type_def.h"
#include <mpi.h>
namespace MAPP_NS
{
    enum {MD_mode,DMD_mode};

    class Atoms;
    template<typename> class Vec;
    template<class Func> class LineSearch;
    class Min;

    constexpr int __dim__=3;
    
    
    extern class MAPP* mapp;
        extern int mode;
        extern int precision;
        extern int step_no;
        extern int step_tally;
        extern MPI_Comm world;
    
        extern class Atoms* atoms;
        extern class ForceField* forcefield;
            extern class Neighbor* neighbor;
    
    
        extern class Timer* timer;
        extern class Error* error;
    
        extern FILE* output;
        extern class Write* write;

        extern class GroupCollection* groups;
    
        extern class AtomTypes* atom_types;
    
    
        extern LineSearch<Min>* ls;
        extern class MD* md;
        extern class DMD* dmd;
    
    
    
    
    
    extern class VarManager* g_vm;
}
#endif








