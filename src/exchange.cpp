#include "exchange.h"
#include "atoms.h"
#include "comm.h"
#include "ff.h"
#include "neighbor.h"
using namespace MAPP_NS;
/*------------------------------------------------------------------
 _____  __    __  _____   _   _       ___   __   _   _____   _____  
| ____| \ \  / / /  ___| | | | |     /   | |  \ | | /  ___| | ____| 
| |__    \ \/ /  | |     | |_| |    / /| | |   \| | | |     | |__   
|  __|    }  {   | |     |  _  |   / / | | | |\   | | |  _  |  __|  
| |___   / /\ \  | |___  | | | |  / /  | | | | \  | | |_| | | |___  
|_____| /_/  \_\ \_____| |_| |_| /_/   |_| |_|  \_| \_____/ |_____|
 ------------------------------------------------------------------*/
/*--------------------------------------------
 
 --------------------------------------------*/
Exchange::Exchange(Atoms* atoms,Communication* comm,int& nxchng_vecs_):
natms(atoms->natms),
x(atoms->x),

world(comm->world),
my_p(comm->my_p),
neigh_p(comm->neigh_p),
s_lo(comm->s_lo),
s_hi(comm->s_hi),
xchng_id(comm->xchng_id),

vecs(atoms->vecs),
nxchng_vecs(nxchng_vecs_)
{
    
    
    buff_grw=1024;
    
    snd_buff[0]=snd_buff[1]=NULL;
    snd_buff_sz[0]=snd_buff_sz[1]=0;
    snd_buff_cpcty[0]=snd_buff_cpcty[1]=0;
    
    rcv_buff=NULL;
    rcv_buff_sz=0;
    rcv_buff_cpcty=0;
    
    tot_xchng_sz=0;
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        tot_xchng_sz+=vecs[ivec]->byte_sz;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/

Exchange::~Exchange()
{

    delete [] snd_buff[0];
    delete [] snd_buff[1];
    delete [] rcv_buff;
}
/*--------------------------------------------
 
 --------------------------------------------*/

inline void Exchange::load(int& iatm,int idir)
{
    if(snd_buff_cpcty[idir]<snd_buff_sz[idir]+tot_xchng_sz)
    {
        byte* tmp_buff=new byte[snd_buff_sz[idir]+tot_xchng_sz+buff_grw];
        memcpy(tmp_buff,snd_buff[idir],snd_buff_sz[idir]);
        delete [] snd_buff[idir];
        snd_buff[idir]=tmp_buff;
        snd_buff_cpcty[idir]=snd_buff_sz[idir]+tot_xchng_sz+buff_grw;
    }
    byte* tmp_buff=&snd_buff[idir][snd_buff_sz[idir]];
    
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        vecs[ivec]->pop_out(tmp_buff,iatm);
    
    snd_buff_sz[idir]+=tot_xchng_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/

inline void Exchange::load(byte*& buff,int& idir)
{
    if(snd_buff_cpcty[idir]<snd_buff_sz[idir]+tot_xchng_sz)
    {
        byte* tmp_buff=new byte[snd_buff_sz[idir]+tot_xchng_sz+buff_grw];
        memcpy(tmp_buff,snd_buff[idir],snd_buff_sz[idir]);
        delete [] snd_buff[idir];
        snd_buff[idir]=tmp_buff;
        snd_buff_cpcty[idir]=snd_buff_sz[idir]+tot_xchng_sz+buff_grw;
    }
    memcpy(&snd_buff[idir][snd_buff_sz[idir]],buff,tot_xchng_sz);
    snd_buff_sz[idir]+=tot_xchng_sz;
    buff+=tot_xchng_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/

int Exchange::xchng_buff(int idim,int idir)
{
    rcv_buff_sz=0;
    int max_snd_sz;
    MPI_Allreduce(&snd_buff_sz[idir],&max_snd_sz,1,MPI_INT,MPI_MAX,world);
    if(max_snd_sz==0)
        return 0;
    
    MPI_Sendrecv(&snd_buff_sz[idir],1,MPI_INT,neigh_p[idim][idir],0,
                 &rcv_buff_sz,1,MPI_INT,neigh_p[idim][1-idir],0,
                 world,MPI_STATUS_IGNORE);
    
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+buff_grw];
        rcv_buff_cpcty=rcv_buff_sz+buff_grw;
    }
    
#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff[idir],snd_buff_sz[idir],MPI_BYTE,neigh_p[idim][idir],0,
                 rcv_buff,rcv_buff_sz,MPI_BYTE,neigh_p[idim][1-idir],0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_buff_sz)
        MPI_Irecv(rcv_buff,rcv_buff_sz,
        MPI_BYTE,neigh_p[idim][1-idir],
        0,
        world,&request[0]);
    
    if(snd_buff_sz[idir])
        MPI_Isend(snd_buff[idir],snd_buff_sz[idir],
        MPI_BYTE,neigh_p[idim][idir],
        0,
        world,&request[1]);

    if(rcv_buff_sz && snd_buff_sz[idir])
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_buff_sz && snd_buff_sz[idir]==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_buff_sz==0 && snd_buff_sz[idir])
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    
    snd_buff_sz[idir]=0;
    
    /* 
     * here make sure, that we have enough space
     * for the potential atoms to be inserted in
     * xchng vectors
     */
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        vecs[ivec]->reserve(rcv_buff_sz/tot_xchng_sz);
    
    return 1;
}
/*--------------------------------------------
 
 --------------------------------------------*/

void Exchange::full_xchng()
{
    
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        vecs[ivec]->resize(natms);
    int disp;
    type0 s,ds_lo,ds_hi;
    int iatm;
    int xchng_lcl=0;
    
    for(int idim=0;idim<__dim__;idim++)
    {
        if(my_p==neigh_p[idim][0] && my_p==neigh_p[idim][1])
            continue;
        
        disp=idim*sizeof(type0);
        
        snd_buff_sz[0]=snd_buff_sz[1]=0;
        iatm=0;
        while (iatm<x->vec_sz)
        {
            s=(*x)(iatm,idim);
            
            if(s_lo[idim]<=s && s<s_hi[idim])
                iatm++;
            else
            {
                if(s>=s_hi[idim])
                {
                    ds_hi=s-s_hi[idim];
                    ds_lo=1.0+s_lo[idim]-s;
                }
                else
                {
                    ds_hi=1.0+s-s_hi[idim];
                    ds_lo=s_lo[idim]-s;
                }
                
                if(ds_hi<ds_lo)
                    load(iatm,snd_to_frnt);
                else
                    load(iatm,snd_to_bhnd);                
              
            }
        }
        
        if(snd_buff_sz[0]+snd_buff_sz[1])
            xchng_lcl=1;
        
        byte* buff_tmp;
        for(int idir=0;idir<2;idir++)
        {
            while(xchng_buff(idim,idir))
            {
                buff_tmp=rcv_buff;
                while(rcv_buff_sz)
                {
                    s=*(type0*)(buff_tmp+disp);
                    if(s_lo[idim]<=s && s<s_hi[idim])
                        for(int ivec=0;ivec<nxchng_vecs;ivec++)
                            vecs[ivec]->pop_in(buff_tmp);
                    else
                        load(buff_tmp,idir);
                    
                    rcv_buff_sz-=tot_xchng_sz;
                }
            }
            
        }
    }
    int xchng;
    MPI_Allreduce(&xchng_lcl,&xchng,1,MPI_INT,MPI_MAX,world);
    if(xchng) xchng_id++;
    natms=x->vec_sz;
}
/*------------------------------------------------
 _   _   _____   _____       ___   _____   _____  
| | | | |  _  \ |  _  \     /   | |_   _| | ____| 
| | | | | |_| | | | | |    / /| |   | |   | |__   
| | | | |  ___/ | | | |   / / | |   | |   |  __|  
| |_| | | |     | |_| |  / /  | |   | |   | |___  
\_____/ |_|     |_____/ /_/   |_|   |_|   |_____|
 ------------------------------------------------*/
/*--------------------------------------------
 
 --------------------------------------------*/
Update::Update(Atoms* atoms,Communication* comm,int& nupdt_vecs_,int& nxchng_vecs_):
natms(atoms->natms),
natms_ph(atoms->natms_ph),
H(atoms->H),
B(atoms->B),
max_cut_s(atoms->max_cut_s),

my_p(comm->my_p),
neigh_p(comm->neigh_p),
s_lo(comm->s_lo),
s_hi(comm->s_hi),

max_cut(atoms->max_cut),
x(atoms->x),

vecs(atoms->vecs),
nvecs(atoms->nvecs),
nupdt_vecs(nupdt_vecs_),
nxchng_vecs(nxchng_vecs_)
{

    snd_buff=NULL;
    snd_buff_cpcty=0;
    snd_buff_grw=1024;
    
    rcv_buff=NULL;
    rcv_buff_cpcty=0;
    rcv_buff_grw=1024;
    snd_atms_lst_grw=8;
    
    tot_updt_vecs_sz=0;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        tot_updt_vecs_sz+=vecs[ivec]->byte_sz;
    
    int icurs=0;
    for(int idim=0;idim<__dim__;idim++)
    {
        if(my_p==neigh_p[idim][0] && my_p==neigh_p[idim][1])
            comm_manager[idim]
            =new LoadUnLoadUpdateSelfComm(this);
        else
            comm_manager[idim]=new LoadUnLoadUpdateComm(this,comm->world);
        
        // snd_to_bhnd && rcv_fm_frnt
        if(comm->my_loc[idim]==comm->tot_p_grid[idim]-1)
            pbc_correction[icurs]=true;
        else
            pbc_correction[icurs]=false;
        icurs++;
        
        // snd_to_frnt && rcv_fm_bhnd
        if(comm->my_loc[idim]==0)
            pbc_correction[icurs]=true;
        else
            pbc_correction[icurs]=false;
        icurs++;
    }

    tot_ncomms=0;
    snd_atms_lst=NULL;
    snd_atms_lst_cpcty=snd_atms_lst_sz=rcv_atms_lst_sz=NULL;
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
Update::~Update()
{
    delete [] rcv_buff;
    delete [] snd_buff;
    
    for(int i=0;i<tot_ncomms;i++)
        delete [] snd_atms_lst[i];
    
    delete [] snd_atms_lst;
    delete [] snd_atms_lst_cpcty;
    delete [] snd_atms_lst_sz;
    delete [] rcv_atms_lst_sz;
    
    for(int idim=0;idim<__dim__;idim++)
        delete comm_manager[idim];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Update::reset()
{
    for(int idim=0;idim<__dim__;idim++)
    {
        type0 tmp=0.0;
        for(int jdim=idim;jdim<__dim__;jdim++)
            tmp+=B[jdim][idim]*B[jdim][idim];
        
        max_cut_s[idim]=sqrt(tmp);
    }
    for(int idim=0;idim<__dim__;idim++)
        max_cut_s[idim]*=max_cut;
    
    int icurs=0;
    int tot_ncomms_=0;
    for(int idim=0;idim<__dim__;idim++)
    {
        // snd_to_bhnd && rcv_fm_frnt
        s_bnd[icurs]=s_lo[idim]+max_cut_s[idim];
        tot_ncomms_+=static_cast<int>(max_cut_s[idim]/(s_hi[idim]-s_lo[idim]))+1;
        ncomms[icurs++]=tot_ncomms_;
        
        // snd_to_frnt && rcv_fm_bhnd
        s_bnd[icurs]=s_hi[idim]-max_cut_s[idim];
        tot_ncomms_+=static_cast<int>(max_cut_s[idim]/(s_hi[idim]-s_lo[idim]))+1;
        ncomms[icurs++]=tot_ncomms_;
    }
    
    if(tot_ncomms_==tot_ncomms)
        return;
    
    for(int i=0;i<tot_ncomms;i++)
        delete [] snd_atms_lst[i];
    delete [] snd_atms_lst;
    delete [] snd_atms_lst_cpcty;
    delete [] snd_atms_lst_sz;
    delete [] rcv_atms_lst_sz;
    
    tot_ncomms=tot_ncomms_;
    snd_atms_lst=new int*[tot_ncomms];
    snd_atms_lst_cpcty=new int[tot_ncomms];
    snd_atms_lst_sz=new int[tot_ncomms];
    rcv_atms_lst_sz=new int[tot_ncomms];
    for(int i=0;i<tot_ncomms;i++)
    {
        snd_atms_lst_cpcty[i]=0;
        snd_atms_lst[i]=NULL;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::add_to_snd_lst(int& icomm,int& iatm)
{
    if(snd_atms_lst_sz[icomm]+1>snd_atms_lst_cpcty[icomm])
    {
        int* tmp_lst=new int[snd_atms_lst_sz[icomm]+1+snd_atms_lst_grw];
        memcpy(tmp_lst,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]*sizeof(int));
        delete [] snd_atms_lst[icomm];
        snd_atms_lst[icomm]=tmp_lst;
        snd_atms_lst_cpcty[icomm]=snd_atms_lst_sz[icomm]+1+snd_atms_lst_grw;
    }
    snd_atms_lst[icomm][snd_atms_lst_sz[icomm]]=iatm;
    snd_atms_lst_sz[icomm]++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::reserve_rcv_buff(int xtra)
{
    if(rcv_buff_cpcty<xtra+rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[xtra+rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=xtra+rcv_buff_sz+rcv_buff_grw;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::reserve_snd_buff(int xtra)
{
    if(snd_buff_cpcty<xtra+snd_buff_sz)
    {
        delete [] snd_buff;
        snd_buff=new byte[xtra+snd_buff_sz+snd_buff_grw];
        snd_buff_cpcty=xtra+snd_buff_sz+snd_buff_grw;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Update::update(vec** updt_vecs,int nupdt_vecs,bool x_xst)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    int tot_byte_sz=0;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        updt_vecs[ivec]->vec_sz=natms;
        tot_byte_sz+=updt_vecs[ivec]->byte_sz;
    }
    snd_buff_sz=rcv_buff_sz=0;
    reserve_snd_buff(tot_byte_sz*max_snd_atms_lst_sz);
    reserve_rcv_buff(tot_byte_sz*max_rcv_atms_lst_sz);
    

    int icurs=0;
    int icomm=0;
    for(int idim=0;idim<__dim__;idim++)
    {
        for(int idir=0;idir<2;idir++)
        {
            while(icomm<ncomms[icurs])
            {
                comm_manager[idim]->update_mult(icomm
                ,neigh_p[idim][idir]
                ,neigh_p[idim][1-idir]
                ,updt_vecs,nupdt_vecs,tot_byte_sz);
                
                if(x_xst && pbc_correction[icurs])
                {
                    if(idir)
                    {
                        for(int iatm=x->vec_sz-1;iatm>x->vec_sz-rcv_atms_lst_sz[icomm]-1;iatm--)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                x_vec[iatm*x_dim+jdim]-=H[idim][jdim];
                    }
                    else
                    {
                        for(int iatm=x->vec_sz-1;iatm>x->vec_sz-rcv_atms_lst_sz[icomm]-1;iatm--)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                x_vec[iatm*x_dim+jdim]+=H[idim][jdim];
                    }
                }
                icomm++;
            }
            icurs++;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Update::update(vec* updt_vec,bool x_xst)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    snd_buff_sz=0;
    reserve_snd_buff(updt_vec->byte_sz*max_snd_atms_lst_sz);
    updt_vec->vec_sz=natms;
    
    int icurs=0;
    int icomm=0;
    for(int idim=0;idim<__dim__;idim++)
    {
        for(int idir=0;idir<2;idir++)
        {
            while(icomm<ncomms[icurs])
            {
                comm_manager[idim]->update_sing(icomm
                ,neigh_p[idim][idir]
                ,neigh_p[idim][1-idir]
                ,updt_vec);
                
                if(x_xst && pbc_correction[icurs])
                {
                    if(idir)
                    {
                        for(int iatm=x->vec_sz-1;iatm>x->vec_sz-rcv_atms_lst_sz[icomm]-1;iatm--)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                x_vec[iatm*x_dim+jdim]-=H[idim][jdim];
                    }
                    else
                    {
                        for(int iatm=x->vec_sz-1;iatm>x->vec_sz-rcv_atms_lst_sz[icomm]-1;iatm--)
                            for(int jdim=0;jdim<idim+1;jdim++)
                                x_vec[iatm*x_dim+jdim]+=H[idim][jdim];
                    }
                }
                icomm++;
            }
            icurs++;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Update::list()
{
    natms_ph=0;
    type0*& x_vec=x->begin();
    int x_dim=x->dim;
    int icurs=0;
    int icomm=0;
    int lo_atm,hi_atm;
    int last_atm;
    bool dir;
    type0 inc;
    max_snd_atms_lst_sz=max_rcv_atms_lst_sz=0;
    
    for(int idim=0;idim<__dim__;idim++)
    {
        
        last_atm=x->vec_sz;
        
        inc=1.0;
        dir=true;
        for(int idir=0;idir<2;idir++)
        {
            lo_atm=0;
            hi_atm=last_atm;
            while(icomm<ncomms[icurs])
            {
                snd_buff_sz=0;
                rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm]=0;
                
                for(int iatm=lo_atm;iatm<hi_atm;iatm++)
                    if((x_vec[iatm*x_dim+idim]<s_bnd[icurs])==dir)
                        add_to_snd_lst(icomm,iatm);
                
                comm_manager[idim]->load_unload(icomm
                ,neigh_p[idim][idir]
                ,neigh_p[idim][1-idir]);
                
                lo_atm=x->vec_sz-rcv_atms_lst_sz[icomm];
                hi_atm=x->vec_sz;
                if(pbc_correction[icurs])
                    for(int iatm=lo_atm;iatm<hi_atm;iatm++)
                        x_vec[iatm*x_dim+idim]+=inc;
                
                max_snd_atms_lst_sz=MAX(max_snd_atms_lst_sz,snd_atms_lst_sz[icomm]);
                max_rcv_atms_lst_sz=MAX(max_rcv_atms_lst_sz,rcv_atms_lst_sz[icomm]);
                icomm++;
            }
            icurs++;
            inc-=2.0;
            dir=!dir;
        }
    }
    natms_ph=x->vec_sz-natms;
    for(int ivec=nxchng_vecs;ivec<nvecs;ivec++)
    {
        vecs[ivec]->vec_sz=0;
        vecs[ivec]->resize(x->vec_sz);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Update::rm_rdndncy()
{
    snd_buff_sz=rcv_buff_sz=0;
    reserve_rcv_buff(max_snd_atms_lst_sz);
    reserve_snd_buff(natms_ph);
    
    byte* mark=snd_buff;
    forcefield->neighbor->mark_redndnt_ph(mark);
    
    int rcv_atms_lst_sz_;
    int snd_atms_lst_sz_=0;
    int snd_atms_lst_cpcty_=max_snd_atms_lst_sz;
    int* snd_atms_lst_=NULL;
    if(snd_atms_lst_cpcty_) snd_atms_lst_=new int[snd_atms_lst_cpcty_];
    
    int nlocomm;
    byte* mark_=mark+natms_ph;
    int icurs=2*__dim__-1;
    int jcomm=tot_ncomms-1;
    max_snd_atms_lst_sz=max_rcv_atms_lst_sz=0;
    for(int idim=__dim__-1;idim>-1;idim--)
    {
        for(int idir=1;idir>-1;idir--)
        {
            if(icurs!=0)
                nlocomm=ncomms[icurs-1];
            else
                nlocomm=0;
            
            while(jcomm>nlocomm-1)
            {
                mark_-=rcv_atms_lst_sz[jcomm];
                comm_manager[idim]->xchng_buff(
                neigh_p[idim][1-idir],rcv_atms_lst_sz[jcomm],mark_,
                neigh_p[idim][idir],snd_atms_lst_sz[jcomm],rcv_buff);
                
                snd_atms_lst_sz_=0;
                for(int i=0; i<snd_atms_lst_sz[jcomm];i++)
                {
                    if(rcv_buff[i]=='1')
                    {
                        if(snd_atms_lst[jcomm][i]>=natms)
                            mark[snd_atms_lst[jcomm][i]-natms]='1';
                        snd_atms_lst_[snd_atms_lst_sz_++]=snd_atms_lst[jcomm][i];
                    }
                }
                memcpy(snd_atms_lst[jcomm],snd_atms_lst_,snd_atms_lst_sz_*sizeof(int));
                snd_atms_lst_sz[jcomm]=snd_atms_lst_sz_;
                
                rcv_atms_lst_sz_=0;
                for(int i=0; i<rcv_atms_lst_sz[jcomm];i++)
                    if(mark_[i]=='1')
                        rcv_atms_lst_sz_++;
                rcv_atms_lst_sz[jcomm]=rcv_atms_lst_sz_;
                
                max_snd_atms_lst_sz=MAX(max_snd_atms_lst_sz,snd_atms_lst_sz[jcomm]);
                max_rcv_atms_lst_sz=MAX(max_rcv_atms_lst_sz,rcv_atms_lst_sz[jcomm]);
                
                jcomm--;
            }
            icurs--;
        }
    }
    
    delete [] snd_atms_lst_;
    
    int old_2_new_cpcty=natms+natms_ph;
    int* old_2_new=NULL;
    if(old_2_new_cpcty) old_2_new=new int[old_2_new_cpcty];
    
    int list_sz=0;
    int list_cpcty=natms_ph;
    int* list=NULL;
    if(list_cpcty) list=new int[list_cpcty];
    
    for(int iatm=0;iatm<natms;iatm++)
        old_2_new[iatm]=iatm;
    
    icurs=natms;
    for(int iatm=natms;iatm<natms+natms_ph;iatm++)
        if(mark[iatm-natms]=='1')
        {
            old_2_new[iatm]=icurs++;
            list[list_sz++]=iatm;
        }
    
    int new_natms_ph=list_sz;

    for(int icomm=0;icomm<tot_ncomms;icomm++)
        for(int i=0; i<snd_atms_lst_sz[icomm];i++)
            snd_atms_lst[icomm][i]=old_2_new[snd_atms_lst[icomm][i]];

    forcefield->neighbor->rename_atoms(old_2_new);
    delete [] old_2_new;
    
    int* list_=list;
    
    int vec_sz=natms;
    while (*list_==natms+icurs)
    {
        list_++;
        vec_sz++;
        list_sz--;
    }
    
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        vecs[ivec]->vec_sz=vec_sz;
        vecs[ivec]->cpy_pst(list_,list_sz);
    }
    
    delete [] list;
    
    natms_ph=new_natms_ph;
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
Update::LoadUnLoadUpdateComm::
LoadUnLoadUpdateComm(Update* updt,MPI_Comm& world_):
LoadUnLoadUpdate(),
world(world_),

snd_atms_lst(updt->snd_atms_lst),
snd_atms_lst_sz(updt->snd_atms_lst_sz),

rcv_atms_lst_sz(updt->rcv_atms_lst_sz),

snd_buff(updt->snd_buff),
snd_buff_sz(updt->snd_buff_sz),
snd_buff_cpcty(updt->snd_buff_cpcty),
snd_buff_grw(updt->snd_buff_grw),

rcv_buff(updt->rcv_buff),
rcv_buff_sz(updt->rcv_buff_sz),
rcv_buff_cpcty(updt->rcv_buff_cpcty),
rcv_buff_grw(updt->rcv_buff_grw),

vecs(updt->vecs),
nupdt_vecs(updt->nupdt_vecs),
tot_updt_vecs_sz(updt->tot_updt_vecs_sz)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateComm::load
(int& iatm)
{
    if(snd_buff_cpcty<snd_buff_sz+tot_updt_vecs_sz)
    {
        byte* tmp_snd_buff=new byte[snd_buff_sz+tot_updt_vecs_sz+snd_buff_grw];
        memcpy(tmp_snd_buff,snd_buff,snd_buff_sz);
        delete [] snd_buff;
        snd_buff=tmp_snd_buff;
        snd_buff_cpcty=snd_buff_sz+tot_updt_vecs_sz+snd_buff_grw;
    }
    
    byte* tmp_snd_buff=snd_buff+snd_buff_sz;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        vecs[ivec]->cpy(tmp_snd_buff,iatm);
    snd_buff_sz+=tot_updt_vecs_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateComm::unload
(int& icomm,int& snd_p,int& rcv_p)
{
    
    MPI_Sendrecv(&snd_atms_lst_sz[icomm],1,MPI_INT,snd_p,0,
                 &rcv_atms_lst_sz[icomm],1,MPI_INT,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
    
    rcv_buff_sz=rcv_atms_lst_sz[icomm]*tot_updt_vecs_sz;
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=rcv_atms_lst_sz[icomm]*tot_updt_vecs_sz+rcv_buff_grw;
    }

#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff,snd_buff_sz,MPI_BYTE,snd_p,0,
                 rcv_buff,rcv_buff_sz,MPI_BYTE,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_atms_lst_sz[icomm])
        MPI_Irecv(rcv_buff,rcv_buff_sz,
        MPI_BYTE,rcv_p,
        0,
        world,&request[0]);
    
    if(snd_atms_lst_sz[icomm])
        MPI_Isend(snd_buff,snd_buff_sz,
        MPI_BYTE,snd_p,
        0,
        world,&request[1]);
    
    if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm])
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    
    byte* tmp_rcv_buff0=rcv_buff;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        vecs[ivec]->pst(tmp_rcv_buff0,tot_updt_vecs_sz,rcv_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateComm::load_unload
(int& icomm,int& snd_p,int& rcv_p)
{
    snd_buff_sz=snd_atms_lst_sz[icomm]*tot_updt_vecs_sz;
    if(snd_buff_cpcty<snd_buff_sz)
    {
        delete [] snd_buff;
        snd_buff=new byte[snd_buff_sz+snd_buff_grw];
        snd_buff_cpcty=snd_buff_sz+snd_buff_grw;
    }
    
    byte* tmp_snd_buff=snd_buff;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        vecs[ivec]->cpy(tmp_snd_buff,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
    
    MPI_Sendrecv(&snd_atms_lst_sz[icomm],1,MPI_INT,snd_p,0,
                 &rcv_atms_lst_sz[icomm],1,MPI_INT,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
    
    rcv_buff_sz=rcv_atms_lst_sz[icomm]*tot_updt_vecs_sz;
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=rcv_atms_lst_sz[icomm]*tot_updt_vecs_sz+rcv_buff_grw;
    }

#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff,snd_buff_sz,MPI_BYTE,snd_p,0,
                 rcv_buff,rcv_buff_sz,MPI_BYTE,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_atms_lst_sz[icomm])
        MPI_Irecv(rcv_buff,rcv_buff_sz,
        MPI_BYTE,rcv_p,
        0,
        world,&request[0]);
    
    if(snd_atms_lst_sz[icomm])
        MPI_Isend(snd_buff,snd_buff_sz,
        MPI_BYTE,snd_p,
        0,
        world,&request[1]);
    
    if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm])
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    
    byte* tmp_rcv_buff=rcv_buff;
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        vecs[ivec]->reserve(rcv_atms_lst_sz[icomm]);
        vecs[ivec]->pst(tmp_rcv_buff,rcv_atms_lst_sz[icomm]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateComm::update_mult
(int& icomm,int& snd_p,int& rcv_p,vec**& vecs
,int& nvecs,int& vecs_byte_sz)
{
    byte* tmp_snd_buff=snd_buff;
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->cpy(tmp_snd_buff,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff,snd_atms_lst_sz[icomm]*vecs_byte_sz,MPI_BYTE,snd_p,0,
                 rcv_buff,rcv_atms_lst_sz[icomm]*vecs_byte_sz,MPI_BYTE,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_atms_lst_sz[icomm])
        MPI_Irecv(rcv_buff,rcv_atms_lst_sz[icomm]*vecs_byte_sz,
        MPI_BYTE,rcv_p,
        0,
        world,&request[0]);
    
    if(snd_atms_lst_sz[icomm])
        MPI_Isend(snd_buff,snd_atms_lst_sz[icomm]*vecs_byte_sz,
        MPI_BYTE,snd_p,
        0,
        world,&request[1]);
    
    if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm])
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    byte* tmp_rcv_buff=rcv_buff;
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->pst(tmp_rcv_buff,rcv_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateComm::update_sing
(int& icomm,int& snd_p,int& rcv_p,vec*& v)
{
    byte* tmp_snd_buff=snd_buff;
    v->cpy(tmp_snd_buff,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff,snd_atms_lst_sz[icomm]*v->byte_sz,MPI_BYTE,snd_p,0,
                 v->end(),rcv_atms_lst_sz[icomm]*v->byte_sz,MPI_BYTE,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_atms_lst_sz[icomm])
        MPI_Irecv(v->end(),rcv_atms_lst_sz[icomm]*v->byte_sz,
        MPI_BYTE,rcv_p,
        0,
        world,&request[0]);
    
    if(snd_atms_lst_sz[icomm])
        MPI_Isend(snd_buff,snd_atms_lst_sz[icomm]*v->byte_sz,
        MPI_BYTE,snd_p,
        0,
        world,&request[1]);
    
    if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm])
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    v->vec_sz+=rcv_atms_lst_sz[icomm];
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateComm::xchng_buff
(int& snd_p,int& snd_buff_sz_,byte*& snd_buff_
,int& rcv_p,int& rcv_buff_sz_,byte*& rcv_buff_)
{
#ifdef BLOCK_COMM
    MPI_Sendrecv(snd_buff_,snd_buff_sz_,MPI_BYTE,snd_p,0,
                 rcv_buff_,rcv_buff_sz_,MPI_BYTE,rcv_p,0,
                 world,MPI_STATUS_IGNORE);
#else
    if(rcv_buff_sz_)
        MPI_Irecv(rcv_buff_,rcv_buff_sz_,
        MPI_BYTE,rcv_p,
        0,
        world,&request[0]);
    
    if(snd_buff_sz_)
        MPI_Isend(snd_buff_,snd_buff_sz_,
        MPI_BYTE,snd_p,
        0,
        world,&request[1]);
    
    if(rcv_buff_sz_ && snd_buff_sz_)
        MPI_Waitall(2,request,MPI_STATUS_IGNORE);
    else if(rcv_buff_sz_ && snd_buff_sz_==0)
        MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
    else if(rcv_buff_sz_==0 && snd_buff_sz_)
        MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
}
/*------------------------------------------------------------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------------------------------------------------------------*/
Update::LoadUnLoadUpdateSelfComm::
LoadUnLoadUpdateSelfComm(Update* updt):
LoadUnLoadUpdate(),

snd_atms_lst(updt->snd_atms_lst),
snd_atms_lst_sz(updt->snd_atms_lst_sz),

rcv_atms_lst_sz(updt->rcv_atms_lst_sz),

vecs(updt->vecs),
nupdt_vecs(updt->nupdt_vecs)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateSelfComm::
load(int& iatm)
{
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
        vecs[ivec]->cpy_pst(iatm);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateSelfComm::unload
(int& icomm,int&,int&)
{
    rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm];
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateSelfComm::load_unload
(int& icomm,int&,int&)
{
    rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm];
    
    for(int ivec=0;ivec<nupdt_vecs;ivec++)
    {
        vecs[ivec]->reserve(rcv_atms_lst_sz[icomm]);
        vecs[ivec]->cpy_pst(snd_atms_lst[icomm],rcv_atms_lst_sz[icomm]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateSelfComm::update_mult
(int& icomm,int&,int&,vec**& vecs
,int& nvecs,int&)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->cpy_pst(snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateSelfComm::update_sing
(int& icomm,int&,int&,vec*& v)
{
    v->cpy_pst(snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Update::LoadUnLoadUpdateSelfComm::xchng_buff
(int&,int& snd_buff_sz_,byte*& snd_buff_
,int&,int& rcv_buff_sz_,byte*& rcv_buff_)
{
    memcpy(rcv_buff_,snd_buff_,rcv_buff_sz_);
}


