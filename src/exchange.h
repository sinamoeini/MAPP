//
//  exchange.hpp
//  MAPP
//
//  Created by Sina on 11/21/16.
//  Copyright © 2016 MIT. All rights reserved.
//

#ifndef __MAPP__exchange__
#define __MAPP__exchange__
#include "global.h"
/*------------------------------------------------------------------
 _____  __    __  _____   _   _       ___   __   _   _____   _____  
| ____| \ \  / / /  ___| | | | |     /   | |  \ | | /  ___| | ____| 
| |__    \ \/ /  | |     | |_| |    / /| | |   \| | | |     | |__   
|  __|    }  {   | |     |  _  |   / / | | | |\   | | |  _  |  __|  
| |___   / /\ \  | |___  | | | |  / /  | | | | \  | | |_| | | |___  
|_____| /_/  \_\ \_____| |_| |_| /_/   |_| |_|  \_| \_____/ |_____|
 ------------------------------------------------------------------*/
namespace MAPP_NS
{
    template<typename> class Vec;
    class vec;
    class Exchange
    {
    private:
        
        /*things that reference cannot be removed*/
        int& natms;
        Vec<type0>*& x;
        unsigned long& xchng_id;
        
        /*things that reference cannot be removed*/
        const int my_p;
        int (&neigh_p)[__dim__][2];
        type0 (&s_lo)[__dim__];
        type0 (&s_hi)[__dim__];
        
        
        int buff_grw;
        byte* snd_buff[2];
        int snd_buff_sz[2];
        int snd_buff_cpcty[2];
        
        byte* rcv_buff;
        int rcv_buff_sz;
        int rcv_buff_cpcty;
        
        void load(int&,int);
        void load(byte*&,int&);
        int xchng_buff(int,int);
        
        vec**& vecs;
        int& nxchng_vecs;
        int tot_xchng_sz;
        
#ifndef BLOCK_COMM
        MPI_Request request[2];
#endif
        MPI_Comm& world;
    protected:
    public:        
        Exchange(class Atoms*,class Communication*,int&);
        ~Exchange();
        void full_xchng();
    };
}
using namespace MAPP_NS;
/*------------------------------------------------
 _   _   _____   _____       ___   _____   _____  
| | | | |  _  \ |  _  \     /   | |_   _| | ____| 
| | | | | |_| | | | | |    / /| |   | |   | |__   
| | | | |  ___/ | | | |   / / | |   | |   |  __|  
| |_| | | |     | |_| |  / /  | |   | |   | |___  
\_____/ |_|     |_____/ /_/   |_|   |_|   |_____|
 ------------------------------------------------*/
namespace MAPP_NS
{
    
    class Update
    {
    private:
        /*things that reference cannot be removed*/
        int& natms;
        int& natms_ph;
        
        /*things that reference can be removed*/
        type0 (&H)[__dim__][__dim__];
        type0 (&B)[__dim__][__dim__];
        type0 (&max_cut_s)[__dim__];
        const int& my_p;
        int (&neigh_p)[__dim__][2];
        type0 (&s_lo)[__dim__];
        type0 (&s_hi)[__dim__];
        type0& max_cut;
        Vec<type0>*& x;
        
        vec**& vecs;
        int& nvecs;
        int& nxchng_vecs;
        int& nupdt_vecs;
        int tot_updt_vecs_sz;
        
        
        
        
        
        int tot_ncomms;
        int ncomms[2*__dim__];
        bool pbc_correction[2*__dim__];
        type0 s_bnd[2*__dim__];
        
        int** snd_atms_lst;
        int* snd_atms_lst_sz;
        int* snd_atms_lst_cpcty;
        int max_snd_atms_lst_sz;
        int snd_atms_lst_grw;

        int* rcv_atms_lst_sz;
        int max_rcv_atms_lst_sz;

        byte* snd_buff;        
        int snd_buff_sz;
        int snd_buff_cpcty;
        int snd_buff_grw;

        byte* rcv_buff;
        int rcv_buff_sz;
        int rcv_buff_cpcty;
        int rcv_buff_grw;
        
        
        
        
        
        void add_to_snd_lst(int&,int&);
        void reserve_snd_buff(int);
        void reserve_rcv_buff(int);
        
        class LoadUnLoadUpdate;
        class LoadUnLoadUpdateComm;
        class LoadUnLoadUpdateSelfComm;
        LoadUnLoadUpdate* comm_manager[__dim__];
        
    protected:
    public:
        Update(Atoms*,Communication*,int&,int&);
        ~Update();
        void reset();
        void update(vec*,bool);
        void update(vec**,int,bool);
        void list();
        void rm_rdndncy();
    };

}

/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class Update::LoadUnLoadUpdate
    {
    private:
    protected:
    public:
        LoadUnLoadUpdate(){};
        virtual ~LoadUnLoadUpdate(){};
        virtual void load(int&)=0;
        virtual void unload(int&,int&,int&)=0;
        virtual void load_unload(int&,int&,int&)=0;
        virtual void update_mult(int&,int&,int&,vec**&,int&,int&)=0;
        virtual void update_sing(int&,int&,int&,vec*&)=0;
        virtual void xchng_buff(int&,int&,byte*&,int&,int&,byte*&)=0;
    };
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    
    class Update::LoadUnLoadUpdateComm:
    public Update::LoadUnLoadUpdate
    {
    private:
        
        int**& snd_atms_lst;
        int*& snd_atms_lst_sz;
        
        int*& rcv_atms_lst_sz;
        
        byte*& snd_buff;
        int& snd_buff_sz;
        int& snd_buff_cpcty;
        int& snd_buff_grw;
        
        byte*& rcv_buff;
        int& rcv_buff_sz;
        int& rcv_buff_cpcty;
        int& rcv_buff_grw;
        
        vec**& vecs;
        int& nupdt_vecs;
        int& tot_updt_vecs_sz;
        MPI_Comm& world;
#ifndef BLOCK_COMM
        MPI_Request request[2];
#endif
    protected:
    public:
        LoadUnLoadUpdateComm(Update*,MPI_Comm&);
        void load(int&);
        void unload(int&,int&,int&);
        void load_unload(int&,int&,int&);
        void update_mult(int&,int&,int&,vec**&,int&,int&);
        void update_sing(int&,int&,int&,vec*&);
        void xchng_buff(int&,int&,byte*&,int&,int&,byte*&);
    };
}
/*----------------------------------------------------------------------------------------
 
 ----------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class Update::LoadUnLoadUpdateSelfComm:
    public Update::LoadUnLoadUpdate
    {
    private:
        int**& snd_atms_lst;
        int*& snd_atms_lst_sz;
        
        int*& rcv_atms_lst_sz;
        
        vec**& vecs;
        int& nupdt_vecs;
    protected:
    public:
        LoadUnLoadUpdateSelfComm(Update*);
        void load(int&);
        void unload(int&,int&,int&);
        void load_unload(int&,int&,int&);
        void update_mult(int&,int&,int&,vec**&,int&,int&);
        void update_sing(int&,int&,int&,vec*&);
        void xchng_buff(int&,int&,byte*&,int&,int&,byte*&);
    };
}

#endif /* exchange_hpp */
