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
        InitPtrs(MAPP* ptr):
        mapp(ptr),
        memory(ptr->memory),
        error(ptr->error),
        timer(ptr->timer),
        forcefield(ptr->forcefield),
        world(ptr->world),
        atoms(ptr->atoms),
        neighbor(ptr->neighbor),
        ls(mapp->ls),
        atom_types(ptr->atom_types),
        output(ptr->output),
        write(ptr->write),
        step_no(ptr->step_no),
        step_tally(ptr->step_tally)
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








