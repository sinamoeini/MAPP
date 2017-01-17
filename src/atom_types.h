/*--------------------------------------------
 Created by Sina on 06/27/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__atom_types__
#define __MAPP__atom_types__
#include "type_def.h"
namespace MAPP_NS
{
    class AtomTypes 
    {
    private:
        void assign_color_rad(char*,type0*);
    protected:
    public:
        AtomTypes();
        ~AtomTypes();
        
        int no_types;
        char** atom_names;
        type0** clr_rad;
        type0* mass;
        
        int add_type(type0,char*);
        int find_type(const char*);
        int find_type_exist(char*);
        
        atom_type get_atom_type(int,atom_type*,type0*);
    };
}
//using namespace MAPP_NS;
//namespace MAPP_NS
//{
//    template<class T>
//    class PerTypeAttr
//    {
//    private:
//    protected:
//    public:
//        static unsigned int get_rank();
//        static T* allocate(atom_type,size_t);
//        static void deallocate(T&);
//        T& var;
//        AtomTypes*& atom_types;
//        
//        
//    };
//}
///*--------------------------------------------
// 
// --------------------------------------------*/
//template<class T>
//unsigned int PerTypeAttr<T>::get_rank()
//{
//    return 0;
//}
///*--------------------------------------------
// 
// --------------------------------------------*/
//template<class T>
//T* PerTypeAttr<T>::allocate(atom_type sz0,size_t sz1)
//{
//    return new T[sz1];
//}
///*--------------------------------------------
// 
// --------------------------------------------*/
//template<class T>
//void PerTypeAttr<T>::deallocate(T&)
//{
//}
///*--------------------------------------------
// 
// --------------------------------------------*/
//namespace MAPP_NS
//{
//    template<class T>
//    class PerTypeAttr<T*>
//    {
//    private:
//    protected:
//    public:
//        static unsigned int get_rank();
//        static T** allocate(atom_type,size_t);
//        static void deallocate(T*&);
//        
//        T*& var;
//        AtomTypes*& atom_types;
//        PerTypeAttr<T*>(AtomTypes*&,T*&);
//        
//    };
//}
///*--------------------------------------------
// 
// --------------------------------------------*/
//template<class T>
//unsigned int PerTypeAttr<T*>::get_rank()
//{
//    return 1+PerTypeAttr<T>::get_rank();
//}
///*--------------------------------------------
// 
// --------------------------------------------*/
//template<class T>
//T** PerTypeAttr<T*>::allocate(atom_type sz0,size_t sz1)
//{
//    
//    T** ans=new T*[sz1];
//    *ans=PerTypeAttr<T>::allocate(sz0,sz0*sz1);
//    for(atom_type i=1;i<sz1;i++)
//        ans[i]=ans[i-1]+sz0;
//    return ans;
//}
///*--------------------------------------------
// 
// --------------------------------------------*/
//template<class T>
//void PerTypeAttr<T*>::deallocate(T*& v)
//{
//    if(v)
//        PerTypeAttr<T>::deallcate(*v);
//    delete [] v;
//    v=NULL;
//}
///*--------------------------------------------
// 
// --------------------------------------------*/
//template<class T>
//PerTypeAttr<T*>::PerTypeAttr(AtomTypes*& atom_types_,T*& var_):
//atom_types(atom_types_),
//var(var_)
//{
//    
//}
//


#endif 
