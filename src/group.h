#ifndef __MAPP__group__
#define __MAPP__group__
#include "init.h"
namespace MAPP_NS
{
    class Group:protected InitPtrs
    {
    private:
    protected:
        virtual void update()=0;
        int* grp_idx;
        int grp_sz;
        int grp_sz_tot;
    public:
        Group(MAPP*,const char*);
        ~Group();
        void get_idx(int&,int*&);
        char* name;
    };
    
    class GroupCollection:protected InitPtrs
    {
    private:
    protected:
    public:
        GroupCollection(MAPP*);
        ~GroupCollection();
        Group** grps;
        int ngrps;
        Group* find_grp(const char*);
        void del_grp(Group*);
        Group* add_grp_file(const char*,int,char**);
        Group* add_grp_type(const char*,int,char**);
    };
    
    class Group_general:public Group
    {
    private:
        int* grp_id;
        unsigned long xchng_id;
    protected:
        void update();
    public:
        Group_general(MAPP*,const char*,int,char**);
        ~Group_general();

    };
    
    class Group_type:public Group
    {
    private:
        bool* type;
        int ntypes;
    protected:
        void update();
    public:
        Group_type(MAPP*,const char*,int,char**);
        ~Group_type();
    };
}
#endif
