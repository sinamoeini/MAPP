/*--------------------------------------------
 Created by Sina on 05/11/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <stdlib.h>
#include <stdarg.h>
#include "error.h"
#include "atoms.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor of the error handler:
 --------------------------------------------*/
Error::Error()
{
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
    
    if(my_p==0)
        fprintf(output,"ABORTED! %s \n",err_msg);
    va_end (args);

    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*--------------------------------------------
 output the error line and abort the code
 --------------------------------------------*/
void Error::abort_sing(const char* msg,...)
{
    char err_msg[MAXCHAR];
    
    va_list args;
    va_start (args, msg);
    vsprintf (err_msg,msg, args);

    fprintf(output,"ABORTED (%d)! %s \n",atoms->my_p,err_msg);
    va_end (args);
    MPI_Abort(__world__,911);
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*--------------------------------------------
 output the error line and abort the code
 --------------------------------------------*/
void Error::print(const char* msg,...)
{
    char err_msg[MAXCHAR];
    
    va_list args;
    va_start (args, msg);
    vsprintf (err_msg,msg, args);
    
    if(my_p==0)
        fprintf(output,"%s",err_msg);
    va_end (args);
}
/*--------------------------------------------
 output the error line and abort the code
 --------------------------------------------*/
void Error::abort()
{
    MPI_Finalize();
    exit(EXIT_FAILURE);
}
/*--------------------------------------------
 output the warning msg 
 --------------------------------------------*/
void Error::warning(const char *msg,...)
{
    char war_msg[MAXCHAR];
    
    va_list args;
    va_start (args, msg);
    vsprintf (war_msg,msg, args);
    
    if(my_p==0)
        fprintf(output,"WARNING: %s \n"
               ,war_msg);
    va_end (args);
}
