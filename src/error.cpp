/*--------------------------------------------
 Created by Sina on 05/11/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <cstdarg>
#include <cstdlib>
#include <mpi.h>
#include "error.h"
#include "atoms.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor of the error handler:
 --------------------------------------------*/
Error::Error(MAPP* mapp):InitPtrs(mapp)
{
    MPI_Comm_rank(world,&my_no);
}
/*--------------------------------------------
 destructor of the error handler:
 --------------------------------------------*/
Error::~Error()
{
}
/*--------------------------------------------
 output the error line and abort the code
 --------------------------------------------*/
void Error::abort(const char* msg,...)
{
    char err_msg[MAXCHAR];
    
    va_list args;
    va_start (args, msg);
    vsprintf (err_msg,msg, args);
    
    if (my_no==0)
        printf("ABORTED! %s \n"
               ,err_msg);
    va_end (args);
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*--------------------------------------------
 output the error line and abort the code
 --------------------------------------------*/
/*
void Error::abort(int line,const char* file
,const char* msg,...)
{
    char err_msg[MAXCHAR];
    
    va_list args;
    va_start (args, msg);
    vsprintf (err_msg,msg, args);
    
    if (my_no==0)
        printf("ABORTED! %s \n"
               ,err_msg);
    
    
    va_end (args);
    
    char* pch=NULL;
    pch=strchr(file,'/');
    int start=0;
    while (pch!=NULL)
    {
        start=static_cast<int>(pch-file);
        pch=strchr(pch+1,'/');
    }
    start++;
    int fin=static_cast<int>(strlen(file));
    
    if (my_no==0)
    {
        printf("For more details see ");
        for(int i=start;i<fin;i++)
            printf("%c",file[i]);
        printf(":%d\n",line);
    }
    
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
*/
/*--------------------------------------------
 output the warning msg 
 --------------------------------------------*/
void Error::warning(const char *msg,...)
{
    char war_msg[MAXCHAR];
    
    va_list args;
    va_start (args, msg);
    vsprintf (war_msg,msg, args);
    
    if (my_no==0)
        printf("WARNING: %s \n"
               ,war_msg);
    va_end (args);
}
