#ifndef __MAPP__group__
#define __MAPP__group__
#include <mpi.h>
namespace MAPP_NS
{
    class Group
    {
    private:
    protected:
        virtual void update()=0;
        int* grp_idx;
        int grp_sz;
        int grp_sz_tot;
        class AtomTypes*& atom_types;
        MPI_Comm& world;
    public:
        Group(const char*);
        virtual ~Group();
        void get_idx(int&,int*&);
        char* name;
    };
    
    class GroupCollection
    {
    private:
    protected:
    public:
        GroupCollection();
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
        Group_general(const char*,int,char**);
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
        Group_type(const char*,int,char**);
        ~Group_type();
    };
}
#endif
