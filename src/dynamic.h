#ifndef __MAPP__dynamic__
#define __MAPP__dynamic__

#include <mpi.h>
#include "global.h"
#include "exchange.h"

namespace MAPP_NS
{
    class vec;
    template<typename> class Vec;
    class Dynamic
    {
    private:
        //move these to simulation
        
        void store_x0();
        bool decide();
        void store_arch_vecs();
        void restore_arch_vecs_();
        void restore_arch_vecs();
        
        
        bool box_chng;
        
        vec** arch_vecs;
        int narch_vecs;
        int nxchng_vecs;
        int nupdt_vecs;
        Vec<type0>* x0;
        Vec<int>* id_arch;
    protected:
    public:
        Dynamic(class Atoms*&,class Communication*&,vec* const *,int,vec* const *,int,vec* const *,int);
        
        Dynamic(class Atoms*&,class Communication*&,std::initializer_list<vec*>,std::initializer_list<vec*>,std::initializer_list<vec*>);
        Dynamic(class Atoms*&,class Communication*&,std::initializer_list<vec*>,std::initializer_list<vec*>);
        ~Dynamic();
        void add_xchng(vec*);
        void add_updt(vec*);
        
        Atoms*& atoms;
        Communication*& comm;
        Exchange* xchng;
        Update* updt;
        
        void update(vec**,int);
        void update(vec*);
        void init_xchng();
        void fin_xchng();
        void init(bool);
        void fin();
    };
}



#endif
