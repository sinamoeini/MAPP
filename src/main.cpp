/*--------------------------------------------
 Created by Sina on 05/12/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "MAPP.h"
#include <stdlib.h>
#include <mpi.h>
using namespace MAPP_NS;
int main(int nargs, char** args)
{
    MPI_Init(&nargs,&args);
    MAPP* mapp= new MAPP(nargs,args);
    delete mapp;
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}


