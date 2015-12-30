/*--------------------------------------------
 Created by Sina on 05/12/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef MAPP_init_h
#define MAPP_init_h
#include <mpi.h>
#include "MAPP.h"
namespace MAPP_NS
{
    class InitPtrs
    {
    private:
    public:
        InitPtrs(MAPP* _mapp):
        mapp(_mapp),
        memory(_mapp->memory),
        error(_mapp->error),
        timer(_mapp->timer),
        forcefield(_mapp->forcefield),
        world(_mapp->world),
        atoms(_mapp->atoms),
        groups(_mapp->groups),
        neighbor(_mapp->neighbor),
        ls(_mapp->ls),
        atom_types(_mapp->atom_types),
        output(_mapp->output),
        write(_mapp->write),
        step_no(_mapp->step_no),
        step_tally(_mapp->step_tally)
        {}
        virtual ~InitPtrs(){}
    protected:
        MAPP* mapp;
        Memory*& memory;
        Error*& error;
        Timer*& timer;
        ForceField*& forcefield;
        MPI_Comm &world;
        Neighbor*& neighbor;
        GroupCollection*& groups;
        Atoms*& atoms;
        LineSearch<Min>*& ls;
        AtomTypes*& atom_types;
        FILE*& output;
        Write*& write;
        int& step_no;
        int& step_tally;
    };
    
}
#endif








