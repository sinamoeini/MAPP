#ifndef __MAPP__memory__
#define __MAPP__memory__
#include <exception>
#include <new>
#include "macros.h"
namespace MAPP_NS
{
    class Memory {
    private:
    protected:
    public:
        
        template <typename TYPE>
        static TYPE* create(TYPE*&,long,const char*
        ,int,const char*,const char*);

        template <typename TYPE>
        static TYPE** create(TYPE**&,long,long,const char*
        ,int,const char*,const char*);
        
        template <typename TYPE>
        static TYPE** create_2d(TYPE**&,long,long,const char*
        ,int,const char*,const char*);
        
        template <typename TYPE>
        static void del_2d(TYPE**&);
        
        template <typename TYPE>
        static TYPE* grow(TYPE*&,long,long,const char*
        ,int,const char*,const char*);

    };
}
using namespace MAPP_NS;
#include "error.h"
#include "global.h"
/*--------------------------------------------
 create 1d vector
 --------------------------------------------*/
template<typename TYPE>
TYPE* Memory::create(TYPE*& array,long d0
,const char* name,int line_no,const char* file
,const char* function)
{
    try
    {
        array = new TYPE [d0];
    }
    catch(std::bad_alloc&)
    {
        Error::abort("memory allocation failure "
        "in file  %s, function %s, line: %d for "
        "variable: %s",file,function,line_no,name);
    }
    return array;
}
/*--------------------------------------------
 create 2d vector
 --------------------------------------------*/
template <typename TYPE>
TYPE** Memory::create(TYPE**& array,long d0
,long d1,const char* name,int line_no,
const char* file,const char* function)
{
    create(array,d0,name,line_no,file,function);
    for(int i=0;i<d0;i++)
        create(array[i],d1,name,line_no,file,function);
    return array;
}
/*--------------------------------------------
 create 2d vector
 --------------------------------------------*/
template <typename TYPE>
TYPE** Memory::create_2d(TYPE**& array,long d0
,long d1,const char* name,int line_no,
const char* file,const char* function)
{
    create(array,d0,name,line_no,file,function);
    create(*array,d0*d1,name,line_no,file,function);
    
    for(int i=1;i<d0;i++)
        array[i]=array[i-1]+d1;
    return array;
}
/*--------------------------------------------
 create 2d vector
 --------------------------------------------*/
template <typename TYPE>
void Memory::del_2d(TYPE**& array)
{
    if(array==NULL)
    {
        delete [] array;
        return;
    }
    delete [] *array;
    delete [] array;    
}
/*--------------------------------------------
 grow 1d vector
 --------------------------------------------*/
template <typename TYPE>
TYPE* Memory::grow(TYPE*& array,long oldsize,
long newsize,const char* name,int line_no,
const char* file,const char* function)
{
    if(oldsize==0)
    {
        return create(array,newsize,name,line_no,file,function);
    }
    else if(oldsize==newsize)
    {
        return array;
    }
    else
    {
        
        TYPE* newarray=array;
        try
        {
            long size1=newsize;
            long size2=oldsize;
            long size=MIN(size1,size2);
            newarray = new TYPE[newsize];
            memcpy(newarray,array,size*sizeof(TYPE));
            delete [] array;
            array=newarray;
        }
        catch (std::bad_alloc&)
        {
            Error::abort("reallocation "
            "failure in file  %s, function %s, "
            "line: %d for variable: %s",file,function,line_no,name);
        }
        return array;
    }
}

#endif
