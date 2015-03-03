#ifndef __MAPP__memory__
#define __MAPP__memory__
#include <iostream>
#include <exception>
#include <stdio.h>
#include "init.h"
#include "error.h"

#define GROW(A,oldsize,newsize)     \
memory->grow(A,oldsize,newsize,#A,  \
__LINE__,__FILE__,__FUNCTION__)
#define CREATE1D(A,d0)              \
memory->create(A,d0,#A,__LINE__,    \
__FILE__,__FUNCTION__)
#define CREATE2D(A,d0,d1)           \
memory->create(A,d0,d1,#A,__LINE__, \
__FILE__,__FUNCTION__)

using namespace std;
namespace MAPP_NS {
    class Memory : protected InitPtrs{
    private:
    protected:
    public:
        Memory(MAPP *);
        ~Memory();
        
        template <typename TYPE>
        TYPE* create(TYPE*&,int,const char*
        ,int,const char*,const char*);

        template <typename TYPE>
        TYPE** create(TYPE**&,int,int,const char*
        ,int,const char*,const char*);
        
        template <typename TYPE>
        TYPE* grow(TYPE*&,int,int,const char*
        ,int,const char*,const char*);

    };
}

using namespace MAPP_NS;
/*--------------------------------------------
 create 1d vector
 --------------------------------------------*/
template<typename TYPE>
TYPE* Memory::create(TYPE*& array,int d0
,const char* name,int line_no,const char* file
,const char* function)
{
    try
    {
        array = new TYPE [d0];
    }
    catch(bad_alloc&)
    {
        error->abort("memory allocation failure "
        "in file  %s, function %s, line: %d for "
        "variable: %s",file,function,line_no,name);
    }
    return array;
}
/*--------------------------------------------
 create 2d vector
 --------------------------------------------*/
template <typename TYPE>
TYPE** Memory::create(TYPE**& array,int d0
,int d1,const char* name,int line_no,
const char* file,const char* function)
{
    create(array,d0,name,line_no,file,function);
    for(int i=0;i<d0;i++)
        create(array[i],d1,name,line_no,file,function);
    return array;
}
/*--------------------------------------------
 grow 1d vector
 --------------------------------------------*/
template <typename TYPE>
TYPE* Memory::grow(TYPE*& array,int oldsize,
int newsize,const char* name,int line_no,
const char* file,const char* function)
{
    if (oldsize==0)
    {
        return create(array,newsize,name,line_no,file,function);
    }
    else if (oldsize==newsize)
    {
        return array;
    }
    else
    {
        
        TYPE* newarray=array;
        try
        {
            int size1=newsize;
            int size2=oldsize;
            int size=MIN(size1,size2);
            newarray = new TYPE[newsize];
            memcpy(newarray,array,size*sizeof(TYPE));
            delete [] array;
            array=newarray;
        }
        catch (bad_alloc&)
        {
            error->abort("reallocation "
            "failure in file  %s, function %s, "
            "line: %d for variable: %s",file,function,line_no,name);
        }
        return array;
    }
}

#endif
