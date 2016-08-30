#include "atoms.h"
#include "xmath.h"
#include "timer.h"
#include "ff.h"
#include "neighbor.h"
#include "memory.h"
using namespace MAPP_NS;
/*--------------------------------------
 _____   _          __      ___   _____
/  ___/ | |        / /     /   | |  _  \
| |___  | |  __   / /     / /| | | |_| |
\___  \ | | /  | / /     / / | | |  ___/
 ___| | | |/   |/ /     / /  | | | |
/_____/ |___/|___/     /_/   |_| |_|
 
 --------------------------------------*/
Atoms::Swap::Swap(Atoms* atoms_,vec** vecs,int nvecs):
atoms(atoms_),
world(atoms->world),
swap_vecs(vecs),
nswap_vecs(nvecs),
dimension(atoms->dimension),
neigh_p(atoms->comm->neigh_p),
my_p(atoms->comm->my_p),
H(atoms->H),
B(atoms->B),
max_cut_s(atoms->max_cut_s),
s_lo(atoms->s_lo),
s_hi(atoms->s_hi),
x(atoms->x)
{

    snd_buff=NULL;
    snd_buff_cpcty=0;
    snd_buff_grw=1024;
    
    rcv_buff=NULL;
    rcv_buff_cpcty=0;
    rcv_buff_grw=1024;
    snd_atms_lst_grw=8;
    
    tot_swap_vecs_sz=0;
    for(int ivec=0;ivec<nswap_vecs;ivec++)
        tot_swap_vecs_sz+=swap_vecs[ivec]->byte_sz;
    
    
    
    
    ncomms=new int[2*dimension];
    pbc_correction=new int[2*dimension];
    s_bnd=new type0[2*dimension];
    comm_manager=new LoadUnLoadUpdate*[dimension];
    int icurs=0;
    for(int idim=0;idim<dimension;idim++)
    {
        if(my_p==neigh_p[idim][0] && my_p==neigh_p[idim][1])
            comm_manager[idim]=new LoadUnLoadUpdateSelfComm(this);
        else
            comm_manager[idim]=new LoadUnLoadUpdateComm(this);
        
        // snd_to_bhnd && rcv_fm_frnt
        if(atoms->comm->my_loc[idim]==atoms->comm->tot_p_grid[idim]-1)
            pbc_correction[icurs]=1;
        else
            pbc_correction[icurs]=0;
        icurs++;
        
        // snd_to_frnt && rcv_fm_bhnd
        if(atoms->comm->my_loc[idim]==0)
            pbc_correction[icurs]=1;
        else
            pbc_correction[icurs]=0;
        icurs++;
    }

    tot_ncomms=0;
    snd_atms_lst=NULL;
    snd_atms_lst_cpcty=snd_atms_lst_sz=rcv_atms_lst_sz=NULL;
    reset();
}
/*--------------------------------------------
 
 --------------------------------------------*/
Atoms::Swap::~Swap()
{
    delete [] rcv_buff;
    delete [] snd_buff;
    
    for(int i=0;i<tot_ncomms;i++)
        delete [] snd_atms_lst[i];
    
    delete [] snd_atms_lst;
    delete [] snd_atms_lst_cpcty;
    delete [] snd_atms_lst_sz;
    delete [] rcv_atms_lst_sz;
    
    delete [] s_bnd;
    delete [] pbc_correction;
    delete [] ncomms;
    
    for(int idim=0;idim<dimension;idim++)
        delete comm_manager[idim];
    delete [] comm_manager;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::Swap::reset()
{
    for(int idim=0;idim<dimension;idim++)
    {
        type0 tmp=0.0;
        for(int jdim=idim;jdim<dimension;jdim++)
            tmp+=B[jdim][idim]*B[jdim][idim];
        
        max_cut_s[idim]=sqrt(tmp);
    }
    for(int idim=0;idim<dimension;idim++)
        max_cut_s[idim]*=atoms->max_cut;
    
    int icurs=0;
    int tot_ncomms_=0;
    for(int idim=0;idim<dimension;idim++)
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
inline void Atoms::Swap::add_to_snd_lst(int& icomm,int& iatm)
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
inline void Atoms::Swap::reserve_rcv_buff(int xtra)
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
inline void Atoms::Swap::reserve_snd_buff(int xtra)
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
void Atoms::Swap::update(vec** vecs,int nvecs,bool x_xst)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    int tot_byte_sz=0;
    for(int ivec=0;ivec<nvecs;ivec++)
    {
        vecs[ivec]->vec_sz=atoms->natms;
        tot_byte_sz+=vecs[ivec]->byte_sz;
    }
    snd_buff_sz=rcv_buff_sz=0;
    reserve_snd_buff(tot_byte_sz*max_snd_atms_lst_sz);
    reserve_rcv_buff(tot_byte_sz*max_rcv_atms_lst_sz);
    

    int icurs=0;
    int icomm=0;
    for(int idim=0;idim<dimension;idim++)
    {
        for(int idir=0;idir<2;idir++)
        {
            while(icomm<ncomms[icurs])
            {
                comm_manager[idim]->update_mult(icomm
                ,neigh_p[idim][idir]
                ,neigh_p[idim][1-idir]
                ,vecs,nvecs,tot_byte_sz);
                
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
void Atoms::Swap::update(vec* v,bool x_xst)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    snd_buff_sz=0;
    reserve_snd_buff(v->byte_sz*max_snd_atms_lst_sz);
    v->vec_sz=atoms->natms;
    
    int icurs=0;
    int icomm=0;
    for(int idim=0;idim<dimension;idim++)
    {
        for(int idir=0;idir<2;idir++)
        {
            while(icomm<ncomms[icurs])
            {
                comm_manager[idim]->update_sing(icomm
                ,neigh_p[idim][idir]
                ,neigh_p[idim][1-idir]
                ,v);
                
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
void Atoms::Swap::list()
{
    type0*& x_vec=x->begin();
    int x_dim=x->dim;
    int icurs=0;
    int icomm=0;
    int lo_atm,hi_atm;
    int last_atm;
    bool dir;
    type0 inc;
    max_snd_atms_lst_sz=max_rcv_atms_lst_sz=0;
    
    for(int idim=0;idim<dimension;idim++)
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
    atoms->natms_ph=x->vec_sz-atoms->natms;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::Swap::eliminate_redundancy()
{
    int natms=atoms->natms,natms_ph=atoms->natms_ph;
    snd_buff_sz=rcv_buff_sz=0;
    reserve_rcv_buff(max_snd_atms_lst_sz);
    reserve_snd_buff(natms_ph);
    
    byte* mark=snd_buff;
    memset(mark,'0',natms_ph);

    int** neighbor_list=atoms->neighbor->neighbor_list;
    int* neighbor_list_sz=atoms->neighbor->neighbor_list_size;
    for(int iatm=0;iatm<natms;iatm++)
        for(int j=0;j<neighbor_list_sz[iatm];j++)
            if(neighbor_list[iatm][j]>=natms)
                mark[neighbor_list[iatm][j]-natms]='1';
    
    int rcv_atms_lst_sz_;
    int snd_atms_lst_sz_=0;
    int snd_atms_lst_cpcty_=max_snd_atms_lst_sz;
    int* snd_atms_lst_=new int[snd_atms_lst_cpcty_];
    
    int nlocomm;
    byte* mark_=mark+natms_ph;
    int icurs=2*dimension-1;
    int jcomm=tot_ncomms-1;
    max_snd_atms_lst_sz=max_rcv_atms_lst_sz=0;
    for(int idim=dimension-1;idim>-1;idim--)
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
    
    if(snd_atms_lst_cpcty_)
        delete [] snd_atms_lst_;
    
    int old_2_new_cpcty=natms+natms_ph;
    int* old_2_new=new int[old_2_new_cpcty];
    
    int list_sz=0;
    int list_cpcty=natms_ph;
    int* list=new int[list_cpcty];
    
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

    for(int iatm=0;iatm<natms;iatm++)
        for(int j=0;j<neighbor_list_sz[iatm];j++)
            neighbor_list[iatm][j]=old_2_new[neighbor_list[iatm][j]];

    if(old_2_new_cpcty)
        delete [] old_2_new;
    
    int* list_=list;
    
    int vec_sz=natms;
    while (*list_==natms+icurs)
    {
        list_++;
        vec_sz++;
        list_sz--;
    }
    
    for(int ivec=0;ivec<nswap_vecs;ivec++)
    {
        swap_vecs[ivec]->vec_sz=vec_sz;
        swap_vecs[ivec]->cpy_pst(list_,list_sz);
    }
    
    if(list_cpcty)
        delete [] list;
    
    atoms->natms_ph=new_natms_ph;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Atoms::Swap::LoadUnLoadUpdateComm::
LoadUnLoadUpdateComm(Swap* swap):
LoadUnLoadUpdate(),
world(swap->world),

snd_atms_lst(swap->snd_atms_lst),
snd_atms_lst_sz(swap->snd_atms_lst_sz),

rcv_atms_lst_sz(swap->rcv_atms_lst_sz),

snd_buff(swap->snd_buff),
snd_buff_sz(swap->snd_buff_sz),
snd_buff_cpcty(swap->snd_buff_cpcty),
snd_buff_grw(swap->snd_buff_grw),

rcv_buff(swap->rcv_buff),
rcv_buff_sz(swap->rcv_buff_sz),
rcv_buff_cpcty(swap->rcv_buff_cpcty),
rcv_buff_grw(swap->rcv_buff_grw),

swap_vecs(swap->swap_vecs),
nswap_vecs(swap->nswap_vecs),
tot_swap_vecs_sz(swap->tot_swap_vecs_sz)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateComm::load
(int& iatm)
{
    if(snd_buff_cpcty<snd_buff_sz+tot_swap_vecs_sz)
    {
        byte* tmp_snd_buff=new byte[snd_buff_sz+tot_swap_vecs_sz+snd_buff_grw];
        memcpy(tmp_snd_buff,snd_buff,snd_buff_sz);
        delete [] snd_buff;
        snd_buff=tmp_snd_buff;
        snd_buff_cpcty=snd_buff_sz+tot_swap_vecs_sz+snd_buff_grw;
    }
    
    byte* tmp_snd_buff=snd_buff+snd_buff_sz;
    for(int ivec=0;ivec<nswap_vecs;ivec++)
        swap_vecs[ivec]->cpy(tmp_snd_buff,iatm);
    snd_buff_sz+=tot_swap_vecs_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateComm::unload
(int& icomm,int& snd_p,int& rcv_p)
{
    MPI_Status status[2];
    MPI_Request request[2];
    
    MPI_Sendrecv(&snd_atms_lst_sz[icomm],1,MPI_INT,snd_p,0
    ,&rcv_atms_lst_sz[icomm],1,MPI_INT,rcv_p,0,world,&status[0]);
    
    rcv_buff_sz=rcv_atms_lst_sz[icomm]*tot_swap_vecs_sz;
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=rcv_atms_lst_sz[icomm]*tot_swap_vecs_sz+rcv_buff_grw;
    }
    
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
        MPI_Waitall(2,request,status);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],&status[0]);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],&status[1]);
    
    byte* tmp_rcv_buff0=rcv_buff;
    for(int ivec=0;ivec<nswap_vecs;ivec++)
        swap_vecs[ivec]->pst(tmp_rcv_buff0,tot_swap_vecs_sz,rcv_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateComm::load_unload
(int& icomm,int& snd_p,int& rcv_p)
{
    snd_buff_sz=snd_atms_lst_sz[icomm]*tot_swap_vecs_sz;
    if(snd_buff_cpcty<snd_buff_sz)
    {
        delete [] snd_buff;
        snd_buff=new byte[snd_buff_sz+snd_buff_grw];
        snd_buff_cpcty=snd_buff_sz+snd_buff_grw;
    }
    
    byte* tmp_snd_buff=snd_buff;
    for(int ivec=0;ivec<nswap_vecs;ivec++)
        swap_vecs[ivec]->cpy(tmp_snd_buff,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
    
    
    MPI_Status status[2];
    MPI_Request request[2];
    
    MPI_Sendrecv(&snd_atms_lst_sz[icomm],1,MPI_INT,snd_p,0
    ,&rcv_atms_lst_sz[icomm],1,MPI_INT,rcv_p,0,world,&status[0]);
    
    rcv_buff_sz=rcv_atms_lst_sz[icomm]*tot_swap_vecs_sz;
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+rcv_buff_grw];
        rcv_buff_cpcty=rcv_atms_lst_sz[icomm]*tot_swap_vecs_sz+rcv_buff_grw;
    }
    
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
        MPI_Waitall(2,request,status);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],&status[0]);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],&status[1]);
    
    byte* tmp_rcv_buff=rcv_buff;
    for(int ivec=0;ivec<nswap_vecs;ivec++)
    {
        swap_vecs[ivec]->reserve(rcv_atms_lst_sz[icomm]);
        swap_vecs[ivec]->pst(tmp_rcv_buff,rcv_atms_lst_sz[icomm]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateComm::update_mult
(int& icomm,int& snd_p,int& rcv_p,vec**& vecs
,int& nvecs,int& vecs_byte_sz)
{
    byte* tmp_snd_buff=snd_buff;
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->cpy(tmp_snd_buff,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
    
    MPI_Status status[2];
    MPI_Request request[2];
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
        MPI_Waitall(2,request,status);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],&status[0]);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],&status[1]);
    
    byte* tmp_rcv_buff=rcv_buff;
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->pst(tmp_rcv_buff,rcv_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateComm::update_sing
(int& icomm,int& snd_p,int& rcv_p,vec*& v)
{
    byte* tmp_snd_buff=snd_buff;

    v->cpy(tmp_snd_buff,snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
    
    MPI_Status status[2];
    MPI_Request request[2];
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
        MPI_Waitall(2,request,status);
    else if(rcv_atms_lst_sz[icomm] && snd_atms_lst_sz[icomm]==0)
        MPI_Waitall(1,&request[0],&status[0]);
    else if(rcv_atms_lst_sz[icomm]==0 && snd_atms_lst_sz[icomm])
        MPI_Waitall(1,&request[1],&status[1]);
    
    v->vec_sz+=rcv_atms_lst_sz[icomm];
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateComm::xchng_buff
(int& snd_p,int& snd_buff_sz_,byte*& snd_buff_
,int& rcv_p,int& rcv_buff_sz_,byte*& rcv_buff_)
{
    MPI_Status status[2];
    MPI_Request request[2];
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
        MPI_Waitall(2,request,status);
    else if(rcv_buff_sz_ && snd_buff_sz_==0)
        MPI_Waitall(1,&request[0],&status[0]);
    else if(rcv_buff_sz_==0 && snd_buff_sz_)
        MPI_Waitall(1,&request[1],&status[1]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
Atoms::Swap::LoadUnLoadUpdateSelfComm::
LoadUnLoadUpdateSelfComm(Swap* swap):
LoadUnLoadUpdate(),

snd_atms_lst(swap->snd_atms_lst),
snd_atms_lst_sz(swap->snd_atms_lst_sz),

rcv_atms_lst_sz(swap->rcv_atms_lst_sz),

swap_vecs(swap->swap_vecs),
nswap_vecs(swap->nswap_vecs)
{
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateSelfComm::
load(int& iatm)
{
    for(int ivec=0;ivec<nswap_vecs;ivec++)
        swap_vecs[ivec]->cpy_pst(iatm);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateSelfComm::unload
(int& icomm,int& ,int& )
{
    rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm];
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateSelfComm::load_unload
(int& icomm,int&,int&)
{
    rcv_atms_lst_sz[icomm]=snd_atms_lst_sz[icomm];
    
    for(int ivec=0;ivec<nswap_vecs;ivec++)
    {
        swap_vecs[ivec]->reserve(rcv_atms_lst_sz[icomm]);
        swap_vecs[ivec]->cpy_pst(snd_atms_lst[icomm],rcv_atms_lst_sz[icomm]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateSelfComm::update_mult
(int& icomm,int&,int&,vec**& vecs
,int& nvecs,int&)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs[ivec]->cpy_pst(snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateSelfComm::update_sing
(int& icomm,int&,int&,vec*& v)
{
    v->cpy_pst(snd_atms_lst[icomm],snd_atms_lst_sz[icomm]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Swap::LoadUnLoadUpdateSelfComm::xchng_buff
(int&,int& snd_buff_sz_,byte*& snd_buff_
,int&,int& rcv_buff_sz_,byte*& rcv_buff_)
{
    memcpy(rcv_buff_,snd_buff_,rcv_buff_sz_);
}
/*---------------------------------------
__    __  _____   _   _   __   _   _____  
\ \  / / /  ___| | | | | |  \ | | /  ___| 
 \ \/ /  | |     | |_| | |   \| | | |     
  }  {   | |     |  _  | | |\   | | |  _
 / /\ \  | |___  | | | | | | \  | | |_| | 
/_/  \_\ \_____| |_| |_| |_|  \_| \_____/

 ---------------------------------------*/
Atoms::Xchng::Xchng(Atoms* atoms
,vec** x_vecs,int nx_vecs):
world(atoms->world),
xchng_vecs(x_vecs),
nxchng_vecs(nx_vecs),
neigh_p(atoms->comm->neigh_p),
my_p(atoms->comm->my_p),
s_lo(atoms->comm->s_lo),
s_hi(atoms->comm->s_hi),
dimension(atoms->dimension),
x(atoms->x),
xchng_id(atoms->xchng_id)
{
    buff_grw=1024;
    
    snd_buff=new byte*[2];
    snd_buff_sz=new int[2];
    snd_buff_cpcty=new int[2];
    snd_buff_sz[0]=snd_buff_sz[1]=0;
    snd_buff_cpcty[0]=snd_buff_cpcty[1]=0;
    
    rcv_buff_sz=0;
    rcv_buff_cpcty=0;
    
    tot_xchng_sz=0;
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        tot_xchng_sz+=xchng_vecs[ivec]->byte_sz;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Atoms::Xchng::~Xchng()
{
    for(int i=0;i<2;i++)
        if(snd_buff_cpcty[i])
            delete [] snd_buff[i];
    
    delete [] snd_buff_cpcty;
    delete [] snd_buff_sz;
    
    if(rcv_buff_cpcty)
        delete [] rcv_buff;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Xchng::load(int& iatm,int idir)
{
    if(snd_buff_cpcty[idir]<snd_buff_sz[idir]+tot_xchng_sz)
    {
        byte* tmp_buff=new byte[snd_buff_sz[idir]+tot_xchng_sz+buff_grw];
        memcpy(tmp_buff,snd_buff[idir],snd_buff_sz[idir]);
        if(snd_buff_cpcty[idir])
            delete [] snd_buff[idir];
        snd_buff[idir]=tmp_buff;
        snd_buff_cpcty[idir]=snd_buff_sz[idir]+tot_xchng_sz+buff_grw;
    }
    byte* tmp_buff=&snd_buff[idir][snd_buff_sz[idir]];
    
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        xchng_vecs[ivec]->pop_out(tmp_buff,iatm);
    
    snd_buff_sz[idir]+=tot_xchng_sz;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void Atoms::Xchng::load(byte*& buff,int& idir)
{
    if(snd_buff_cpcty[idir]<snd_buff_sz[idir]+tot_xchng_sz)
    {
        byte* tmp_buff=new byte[snd_buff_sz[idir]+tot_xchng_sz+buff_grw];
        memcpy(tmp_buff,snd_buff[idir],snd_buff_sz[idir]);
        if(snd_buff_cpcty[idir])
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
int Atoms::Xchng::xchng_buff(int idim,int idir)
{
    rcv_buff_sz=0;
    MPI_Request request[2];
    MPI_Status status[2];
    int max_snd_sz;
    MPI_Allreduce(&snd_buff_sz[idir],&max_snd_sz,1,MPI_INT,MPI_MAX,world);
    if(max_snd_sz==0)
        return 0;
    
    MPI_Sendrecv(&snd_buff_sz[idir],1,MPI_INT,neigh_p[idim][idir],0,&rcv_buff_sz,1,MPI_INT,neigh_p[idim][1-idir],0,world,&status[0]);
    if(rcv_buff_cpcty<rcv_buff_sz)
    {
        if(rcv_buff_cpcty)
            delete [] rcv_buff;
        rcv_buff=new byte[rcv_buff_sz+buff_grw];
        rcv_buff_cpcty=rcv_buff_sz+buff_grw;
    }
    
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
        MPI_Waitall(2,request,status);
    else if(rcv_buff_sz && snd_buff_sz[idir]==0)
        MPI_Waitall(1,&request[0],&status[0]);
    else if(rcv_buff_sz==0 && snd_buff_sz[idir])
        MPI_Waitall(1,&request[1],&status[1]);
    
    snd_buff_sz[idir]=0;
    
    /* 
     * here make sure, that we have enough space
     * for the potential atoms to be inserted in
     * xchng vectors
     */
    for(int ivec=0;ivec<nxchng_vecs;ivec++)
        xchng_vecs[ivec]->reserve(rcv_buff_sz/tot_xchng_sz);
    
    return 1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::Xchng::full_xchng()
{
    int disp;
    type0 s,ds_lo,ds_hi;
    int iatm;
    int xchng_lcl=0;
    
    for(int idim=0;idim<dimension;idim++)
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
                            xchng_vecs[ivec]->pop_in(buff_tmp);
                    else
                        load(buff_tmp,idir);
                    
                    rcv_buff_sz-=tot_xchng_sz;
                }
            }
            
        }
    }
    int xchng;
    MPI_Allreduce(&xchng_lcl,&xchng,1,MPI_INT,MPI_MAX,world);
    if(xchng)
        xchng_id++;
}
/*----------------------------------------------------------------------------------------------------------------
 _____   _____       ___  ___       ___  ___   _   _   __   _   _   _____       ___   _____   _   _____   __   _  
/  ___| /  _  \     /   |/   |     /   |/   | | | | | |  \ | | | | /  ___|     /   | |_   _| | | /  _  \ |  \ | | 
| |     | | | |    / /|   /| |    / /|   /| | | | | | |   \| | | | | |        / /| |   | |   | | | | | | |   \| | 
| |     | | | |   / / |__/ | |   / / |__/ | | | | | | | |\   | | | | |       / / | |   | |   | | | | | | | |\   | 
| |___  | |_| |  / /       | |  / /       | | | |_| | | | \  | | | | |___   / /  | |   | |   | | | |_| | | | \  | 
\_____| \_____/ /_/        |_| /_/        |_| \_____/ |_|  \_| |_| \_____| /_/   |_|   |_|   |_| \_____/ |_|  \_|
 
 ----------------------------------------------------------------------------------------------------------------*/
Atoms::Communincation::Communincation(Atoms* atoms):
world(atoms->world),
dimension(atoms->dimension)
{
    /* beginning of communication related parameters */
    
    MPI_Comm_rank(world,&my_p);
    MPI_Comm_size(world,&tot_p);
    
    int name_length;
    char* node_name=new char[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(node_name,&name_length);
    if(MPI_MAX_PROCESSOR_NAME)
        node_name[name_length]='\0';
    name_length++;
    

    int* name_lenghts=new int[tot_p];
    for(int i=0;i<tot_p;i++)
        name_lenghts[i]=0;
    
    name_lenghts[my_p]=name_length;
    
    
    int* all_name_lenghts=new int[tot_p];
    MPI_Allreduce(name_lenghts,all_name_lenghts,tot_p, MPI_INT,MPI_SUM,world);
    if(tot_p)
        delete [] name_lenghts;
    
    char** all_names=new char*[tot_p];
    for(int i=0;i<tot_p;i++)
        all_names[i]=new char[all_name_lenghts[i]];
    
    
    
    for(int i=0;i<name_length;i++)
        all_names[my_p][i]=node_name[i];
    if(MPI_MAX_PROCESSOR_NAME)
        delete [] node_name;
    
    
    for(int i=0;i<tot_p;i++)
        MPI_Bcast(all_names[i],all_name_lenghts[i],MPI_CHAR,i,world);
    
    int* node_no=new int[tot_p];
    
    for(int i=0;i<tot_p;i++)
        node_no[i]=-1;  
    
    tot_n=0;
    for(int i=0;i<tot_p;i++)
    {
        if(node_no[i]==-1)
        {
            node_no[i]=tot_n;
            tot_n++;
        }
        for(int j=i+1;j<tot_p;j++)
            if(strcmp(all_names[i],all_names[j])==0)
                node_no[j]=node_no[i];
    }
    
    
    for(int i=0;i<tot_p;i++)
        if(all_name_lenghts[i])
            delete [] all_names[i];
    
    if(tot_p)
    {
        delete [] all_names;
        delete [] all_name_lenghts;
    }
    
    p_per_n=new int[tot_n];
    
    for(int i=0;i<tot_n;i++)
        p_per_n[i]=0;
    
    for(int i=0;i<tot_p;i++)
        p_per_n[node_no[i]]++;
    
    my_n=node_no[my_p];

    
    n_p_grid=new int*[tot_n];
    for(int i=0;i<tot_n;i++)
        n_p_grid[i]=new int[p_per_n[i]];
    
    
    int pos;
    for(int inode=0;inode<tot_n;inode++)
    {
        pos=0;
        for(int iproc=0;iproc<tot_p;iproc++)
        {
            if(node_no[iproc]==inode)
            {
                n_p_grid[inode][pos]=iproc;
                pos++;
            }
        }
    }
    if(tot_p)
        delete [] node_no;
    
    /* end of communication related parameters */
    
    tot_p_grid=new int[dimension];
    my_loc=new int[dimension];
    neigh_p=new int*[dimension];
    for(int idim=0;idim<dimension;idim++)
        neigh_p[idim]=new int[2];
    s_lo=new type0[dimension];
    s_hi=new type0[dimension];
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Atoms::Communincation::~Communincation()
{
    if(dimension)
    {
        for(int i=0;i<dimension;i++)
        {
            delete [] neigh_p[i];
        }
        delete [] neigh_p;
        delete [] my_loc;
        delete [] tot_p_grid;
        delete [] s_lo;
        delete [] s_hi;
    }
    
    for(int i=0;i<tot_n;i++)
        if(p_per_n[i])
            delete [] n_p_grid[i];
    if(tot_n)
    {
        delete [] n_p_grid;
        delete [] p_per_n;
    }
}
/*--------------------------------------------
 automatically grid the domain
 --------------------------------------------*/
void Atoms::Communincation::auto_grid(type0** H)
{
    int eq_p_per_n=1;
    for(int i=0;i<tot_n;i++)
        if(p_per_n[i]!=p_per_n[0])
            eq_p_per_n=0;
    
    type0 vol;
    type0* area=new type0[dimension];
    int prin_dimension=0;
    
    if(eq_p_per_n)
    {
        for(int i=0;i<dimension;i++)
        {
            area[i]=1.0;
            for(int j=0;j<dimension;j++)
                if(i!=j)
                    area[i]*=H[j][j];
        }
        
        type0 min_area=-1.0;
        for(int i=0;i<dimension;i++)
            if(min_area<0.0 || area[i]<min_area)
            {
                min_area=area[i];
                prin_dimension=i;
            }
        
    }
    
    XMath* xmath= new XMath();
    int** fac_list;
    int fac_list_size=xmath->fac(tot_p,dimension,fac_list);
    delete  xmath;
    
    type0 ratio=-1.0;
    type0 tmp_ratio;
    
    if(eq_p_per_n && tot_n>1)
    {
        for(int ifac=0;ifac<fac_list_size;ifac++)
        {
            if(fac_list[ifac][prin_dimension]==tot_n)
            {
                vol=1.0;
                for(int i=0;i<dimension;i++)
                {
                    vol*=H[i][i]/static_cast<type0>(fac_list[ifac][i]);
                    area[i]=1.0;
                    for(int j=0;j<dimension;j++)
                        if(i!=j)
                            area[i]*=H[j][j]/static_cast<type0>(fac_list[ifac][j]);
                }
                
                tmp_ratio=0.0;
                for(int i=0;i<dimension;i++)
                    tmp_ratio+=2.0*area[i];
                
                tmp_ratio=tmp_ratio/vol;
                if(tmp_ratio<ratio||ratio<0.0)
                {
                    ratio=tmp_ratio;
                    for(int i=0;i<dimension;i++)
                        tot_p_grid[i]=fac_list[ifac][i];
                }
            }
            
        }
        
        
        for(int i=0;i<fac_list_size;i++)
            delete [] fac_list[i];
        delete [] fac_list;
        
        
        if(dimension>1)
        {
            int* tmp_tot_p_grid=new int[dimension-1];
            int* dim_indx=new int[dimension-1];
            int* my_tmp_loc=new int[dimension-1];
            int* my_tmp_denom=new int[dimension-1];
            
            
            int pos=0;
            
            for(int i=0;i<dimension;i++)
            {
                if(i!=prin_dimension)
                {
                    tmp_tot_p_grid[pos]=tot_p_grid[i];
                    dim_indx[pos]=i;
                    pos++;
                }
            }
            
            int my_p_in_my_node=-1;
            for(int i=0;i<p_per_n[my_n];i++)
                if(n_p_grid[my_n][i]==my_p)
                    my_p_in_my_node=i;
            
            for(int i=0;i<dimension-1;i++)
            {
                
                int no=1;
                for(int j=0;j<i;j++)
                    no*=tmp_tot_p_grid[j];
                my_tmp_denom[i]=no;
            }
            
            
            int ttmp=my_p_in_my_node;
            for(int i=dimension-2;i>-1;i--)
            {
                my_tmp_loc[i]=ttmp/my_tmp_denom[i];
                ttmp-=my_tmp_loc[i]*my_tmp_denom[i];
            }
            
            
            for(int i=0;i<dimension-1;i++)
            {
                my_loc[dim_indx[i]]=my_tmp_loc[i];
            }
            my_loc[prin_dimension]=my_n;
            
            
            
            int nop;
            int nom;
            int tmp0,tmp1;
            for(int i=0;i<dimension-1;i++)
            {
                
                tmp0=my_tmp_loc[i]-1;
                if(tmp0==-1)
                    tmp0=tmp_tot_p_grid[i]-1;
                
                tmp1=my_tmp_loc[i]+1;
                if(tmp1==tmp_tot_p_grid[i])
                    tmp1=0;
                
                nop=0;
                nom=0;
                
                for(int j=0;j<dimension-1;j++)
                {
                    if(i==j)
                    {
                        nom+=tmp0*my_tmp_denom[j];
                        nop+=tmp1*my_tmp_denom[j];
                    }
                    else
                    {
                        nom+=my_tmp_loc[j]*my_tmp_denom[j];
                        nop+=my_tmp_loc[j]*my_tmp_denom[j];
                    }
                }
                
                neigh_p[dim_indx[i]][1]=n_p_grid[my_n][nop];
                neigh_p[dim_indx[i]][0]=n_p_grid[my_n][nom];
            }
            
            tmp0=my_n-1;
            if(tmp0==-1)
                tmp0=tot_n-1;
            neigh_p[prin_dimension][0]=n_p_grid[tmp0][my_p_in_my_node];
            
            tmp1=my_n+1;
            if(tmp1==tot_n)
                tmp1=0;
            neigh_p[prin_dimension][1]=n_p_grid[tmp1][my_p_in_my_node];
            
            
            delete [] dim_indx;
            delete [] tmp_tot_p_grid;
            delete [] my_tmp_loc;
            delete [] my_tmp_denom;
            
            
        }
        else
        {
            tot_p_grid[0]=tot_p;
            
            int my_p_in_my_node=-1;
            for(int i=0;i<p_per_n[my_n];i++)
                if(n_p_grid[my_n][i]==my_p)
                    my_p_in_my_node=i;
            
            
            my_loc[0]=p_per_n[my_n]*my_n;
            my_loc[0]+=my_p_in_my_node;
            
            int tmp0p,tmp1p,tmp0n,tmp1n;
            tmp0p=my_p_in_my_node;
            tmp1p=my_p_in_my_node;
            tmp0n=my_n;
            tmp1n=my_n;
            
            tmp0p--;
            if(tmp0p==-1)
            {
                tmp0p=p_per_n[my_n]-1;
                tmp0n--;
                if(tmp0n==-1)
                    tmp0n=tot_n-1;
            }
            
            tmp1p++;
            if(tmp1p==p_per_n[my_n])
            {
                tmp1p=0;
                tmp1n++;
                if(tmp1n==tot_n)
                    tmp1n=0;
            }
            
            neigh_p[0][0]=n_p_grid[tmp0n][tmp0p];
            neigh_p[0][1]=n_p_grid[tmp1n][tmp1p];
            
        }
        
    }
    else
    {
        for(int ifac=0;ifac<fac_list_size;ifac++)
        {
            vol=1.0;
            for(int i=0;i<dimension;i++)
            {
                vol*=H[i][i]/static_cast<type0>(fac_list[ifac][i]);
                area[i]=1.0;
                for(int j=0;j<dimension;j++)
                    if(i!=j)
                        area[i]*=H[j][j]/static_cast<type0>(fac_list[ifac][j]);
            }
            
            tmp_ratio=0.0;
            for(int i=0;i<dimension;i++)
                tmp_ratio+=2.0*area[i];
            
            tmp_ratio=tmp_ratio/vol;
            if(tmp_ratio<ratio||ratio<0.0)
            {
                ratio=tmp_ratio;
                for(int i=0;i<dimension;i++)
                    tot_p_grid[i]=fac_list[ifac][i];
            }
            
        }
        
        for(int i=0;i<fac_list_size;i++)
            delete [] fac_list[i];
        delete [] fac_list;
        
        int* list=new int[dimension];
        
        MPI_Comm cartesian;
        for(int i=0;i<dimension;i++)
            list[i]=1;
        MPI_Cart_create(world,dimension,tot_p_grid,list,1,&cartesian);
        MPI_Cart_get(cartesian,dimension,tot_p_grid,list,my_loc);
        for(int i=0;i<dimension;i++)
            MPI_Cart_shift(cartesian,i,1,&neigh_p[i][0],&neigh_p[i][1]);
        
        MPI_Comm_free(&cartesian);
        
        delete [] list;
    }
    
    delete [] area;
    
    
    for(int i=0;i<dimension;i++)
    {
        s_lo[i]=
        static_cast<type0>(my_loc[i])
        /static_cast<type0>(tot_p_grid[i]);
        
        s_hi[i]=
        (static_cast<type0>(my_loc[i])+1.0)
        /static_cast<type0>(tot_p_grid[i]);
    }
}
/*--------------------------------------------
 manually grid the domain
 --------------------------------------------*/
void Atoms::Communincation::man_grid(int* n)
{
    for(int i=0;i<dimension;i++)
        tot_p_grid[i]=n[i];
    
    int* list=new int[dimension];
    for(int i=0;i<dimension;i++)
        list[i]=1;
    
    MPI_Comm cartesian;
    MPI_Cart_create(world,dimension,tot_p_grid,list,1,&cartesian);
    MPI_Cart_get(cartesian,dimension,tot_p_grid,list,my_loc);
    for(int i=0;i<dimension;i++)
        MPI_Cart_shift(cartesian,i,1,&neigh_p[i][0],&neigh_p[i][1]);
    
    MPI_Comm_free(&cartesian);
    
    if(dimension)
        delete [] list;
    
    for(int i=0;i<dimension;i++)
    {
        s_lo[i]=
        static_cast<type0>(my_loc[i])
        /static_cast<type0>(tot_p_grid[i]);
        
        s_hi[i]=
        (static_cast<type0>(my_loc[i])+1.0)
        /static_cast<type0>(tot_p_grid[i]);
    }
}
/*---------------------------------------------------------------------------
      ___   _____   _____       ___  ___   _____
     /   | |_   _| /  _  \     /   |/   | /  ___/
    / /| |   | |   | | | |    / /|   /| | | |___
   / / | |   | |   | | | |   / / |__/ | | \___  \
  / /  | |   | |   | |_| |  / /       | |  ___| |
 /_/   |_|   |_|   \_____/ /_/        |_| /_____/
 
 ---------------------------------------------------------------------------*/
Atoms::Atoms(MAPP* mapp,int dim):
world(mapp->world),
neighbor(mapp->neighbor),
forcefield(mapp->forcefield),
error(mapp->error),
timer(mapp->timer),
output(mapp->output),
dimension(dim),
nvecs(0)
{
    
    comm=new Communincation(this);
    tot_p=comm->tot_p;
    my_p=comm->my_p;
    s_lo=comm->s_lo;
    s_hi=comm->s_hi;
    
    grid_established=false;
    box_chng=false;
    skin=0.5;
    tot_natms=0;
    natms=0;
    natms_ph=0;
    xchng_id=0;
    
    x=new Vec<type0>(this,dimension);
    id= new Vec<int>(this,1);
    

    max_cut_s=new type0[dimension];
    H=new type0*[dimension];
    B=new type0*[dimension];
    *H=new type0[dimension*dimension];
    *B=new type0[dimension*dimension];
    
    for(int idim=1;idim<dimension;idim++)
    {
        B[idim]=B[idim-1]+dimension;
        H[idim]=H[idim-1]+dimension;
    }
    
    for(int idim=0;idim<dimension;idim++)
    {
        max_cut_s[idim]=0.0;
        for(int jdim=0;jdim<dimension;jdim++)
            H[idim][jdim]=B[idim][jdim]=0.0;
    }
        
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Atoms::~Atoms()
{
    delete comm;
    if(dimension)
    {
        delete [] *B;
        delete [] *H;
        delete [] B;
        delete [] H;
        delete [] max_cut_s;
    }

    while(nvecs)
        delete vecs[0];
}
/*--------------------------------------------
 add a new vec
 --------------------------------------------*/
void Atoms::add_vec(vec* v)
{
    vec** vecs_=new vec*[nvecs+1];
    memcpy(vecs_,vecs,nvecs*sizeof(vec*));
    if(nvecs)
        delete [] vecs;
    vecs=vecs_;
    vecs[nvecs]=v;
    nvecs++;
}
/*--------------------------------------------
 add a new vec with name
 --------------------------------------------*/
void Atoms::add_vec(vec* v,const char* name)
{
    // first see if the name is duplicate
    for(int ivec=0;ivec<nvecs;ivec++)
        if(strcmp(name,vecs[ivec]->name)==0)
            error->abort("vector %s already exist",name);
    
    
    vec** vecs_=new vec*[nvecs+1];
    memcpy(vecs_,vecs,nvecs*sizeof(vec*));
    if(nvecs)
        delete [] vecs;
    vecs=vecs_;
    vecs[nvecs]=v;
    nvecs++;
}
/*--------------------------------------------
 add a new vec with name
 --------------------------------------------*/
vec* Atoms::find_vec(const char* name)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        if(vecs[ivec]->name!=NULL && strcmp(name,vecs[ivec]->name)==0)
            return vecs[ivec];

    error->abort("vector %s does not exist",name);
    
    return NULL;
}
/*--------------------------------------------
 remove a vector
 --------------------------------------------*/
void Atoms::del_vec(vec* v)
{
    if(nvecs-1)
    {
        vec** vecs_=new vec*[nvecs-1];
        int ivec=0;
        for(int jvec=0;jvec<nvecs;jvec++)
            if(vecs[jvec]!=v)
                vecs_[ivec++]=vecs[jvec];
            
        if(nvecs)
            delete [] vecs;
        vecs=vecs_;
        nvecs--;
    }
    else
    {
        delete [] vecs;
        nvecs=0;
    }
}
/*--------------------------------------------
 automatically grid the domain
 --------------------------------------------*/
void Atoms::auto_grid()
{
    type0 det=1.0;
    for(int idim=0;idim<dimension;idim++)
        det*=H[idim][idim]*H[idim][idim];
    if(det==0.0)
        error->abort("cannot perform automatic "
        "domain griding, the domain has not been established");
        
    comm->auto_grid(H);
    grid_established=true;
    
    if(my_p==0)
    {
        fprintf(output,"automatic grid performed: ");
        
        for(int i=0;i<dimension-1;i++)
            fprintf(output,"%dx",comm->tot_p_grid[i]);
        
        fprintf(output,"%d",comm->tot_p_grid[dimension-1]);
        
        fprintf(output,"\n");
    }
}
/*--------------------------------------------
 manually grid the domain
 --------------------------------------------*/
void Atoms::man_grid(int* n)
{
    comm->man_grid(n);
    xchng=new Xchng(this,vecs,nvecs);
    x2s(natms);
    xchng->full_xchng();
    natms=x->vec_sz;
    s2x(natms);
    delete xchng;
    grid_established=true;
    
    if(my_p==0)
    {
        fprintf(output,"manual grid performed: ");
        
        for(int i=0;i<dimension-1;i++)
            fprintf(output,"%dx",comm->tot_p_grid[i]);
        
        fprintf(output,"%d",comm->tot_p_grid[dimension-1]);
        
        fprintf(output,"\n");
    }
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s(int no)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    
    for(int icmp=0;icmp<no*x_dim;icmp+=x_dim)
    {
        for(int idim=0;idim<dimension;idim++)
        {
            x_vec[icmp+idim]=x_vec[icmp+idim]*B[idim][idim];
            for(int jdim=idim+1;jdim<dimension;jdim++)
                x_vec[icmp+idim]+=x_vec[icmp+jdim]*B[jdim][idim];
        }
        for(int idim=0;idim<dimension;idim++)
        {
            while(x_vec[icmp+idim]<0.0)
                x_vec[icmp+idim]++;
            while(x_vec[icmp+idim]>=1.0)
                x_vec[icmp+idim]--;
        }
    }
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x(int no)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    
    for(int icmp=0;icmp<no*x_dim;icmp+=x_dim)
        for(int idim=0;idim<dimension;idim++)
        {
            x_vec[icmp+idim]=x_vec[icmp+idim]*H[idim][idim];
            for(int jdim=idim+1;jdim<dimension;jdim++)
                x_vec[icmp+idim]+=x_vec[icmp+jdim]*H[jdim][idim];
        }

}
/*--------------------------------------------
 insert a number of atoms
 --------------------------------------------*/
void Atoms::insert(byte* buff,vec** vecs_,int nvecs_,int natms_)
{
    int stride=0;
    for(int ivec=0;ivec<nvecs_;ivec++) stride+=vecs_[ivec]->byte_sz;
    

    for(int ivec=0;ivec<nvecs_;ivec++)
        vecs_[ivec]->pst(buff,stride,natms_);
    
    natms+=natms_;
}
/*--------------------------------------------
 make room for some local atoms and phantom 
 atoms; this is used for grand canocical monte 
 carlo, when a successfull insertion trial has
 occured. using this function we make room for 
 the new entries. The new entries are insrerted 
 manually by GCMC
 
 *** we might need a better name for this 
 function
 --------------------------------------------*/
void Atoms::add()
{
    for(int i=0;i<nvecs;i++)
        vecs[i]->add();
    natms++;
}
/*--------------------------------------------
 delete some local atoms and phantom atoms; 
 this is used for grand canocical monte carlo, 
 when a successfull deletion trial has occured.
 it takes a the list of local atoms and phantoms
 
 !! it is assumed that both lists are ascending
 
 *** we might need a better name for this
 function
 --------------------------------------------*/
void Atoms::del(int& del_idx)
{
    for(int i=0;i<nvecs;i++)
        vecs[i]->del(del_idx);
    natms--;
}
/*--------------------------------------------
 restart
 --------------------------------------------*/
void Atoms::restart()
{
    tot_natms=0;
    natms=0;
    natms_ph=0;
    while(nvecs)
        delete vecs[0];
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool Atoms::xchng_chk(unsigned long& xchng_id_)
{
    if(xchng_id_==xchng_id)
        return true;
    else
    {
        xchng_id_=xchng_id;
        return false;
    }
}
/*--------------------------------------------
 init a simulation
 --------------------------------------------*/
void Atoms::init(VecLst* vec_list_,bool box_chng_)
{
    type0 det=1.0;
    for(int idim=0;idim<dimension;idim++)
        det*=H[idim][idim]*H[idim][idim];
    if(det==0.0)
        error->abort("cannot start the simulation "
        "domain griding, the domain is not defined");
    if(!grid_established)
        error->abort("cannot start the simulation "
        "domain grid has not been established");
    
    /*
     ?. check wether the forcefield is NULL
     */
    if(forcefield==NULL)
        error->abort("cannot start the simulation, force field has not been established");
    /*
     ?. check wether the neighbor is NULL
     */
    if(neighbor==NULL)
        error->abort("cannot start the simulation, neighbior list has not been established");
    
    
    timer->init();
    timer->start(COMM_TIME_mode);
    box_chng=box_chng_;
    /*
     ?. reorder vecs so that vec_list->xchng_vecs
     are the first vec_list->nxchng_vecs vectors
     */
    
    int* vec_mapping=new int[nvecs];
    for(int ivec=0;ivec<nvecs;ivec++) vec_mapping[ivec]=-1;
    int vecs_mapped=0;
    for(int ivec=0;ivec<nvecs;ivec++)
    {
        int jvec=0;
        while (jvec<vec_list_->nxchng_vecs && vec_list_->xchng_vecs[jvec]!=vecs[ivec])
            jvec++;
        if(jvec!=vec_list_->nxchng_vecs)
        {
            vec_mapping[ivec]=jvec;
            vecs_mapped++;
        }
    }
    
    if(vecs_mapped!=vec_list_->nxchng_vecs)
        error->abort("not all xchng vectors were found");

    for(int ivec=0;ivec<nvecs;ivec++)
        if(vec_mapping[ivec]==-1)
        {
            vec_mapping[ivec]=vecs_mapped;
            vecs_mapped++;
        }
    
    vec** vecs_= new vec*[nvecs];
    for(int ivec=0;ivec<nvecs;ivec++)
        vecs_[vec_mapping[ivec]]=vecs[ivec];
    if(nvecs)
    {
        delete [] vec_mapping;
        delete [] vecs;
    }
    vecs=vecs_;
    
    /*id_arch=new Vec<int>(this,1);
    int* id_0=id->begin();
    int* id_1=id_arch->begin();
    memcpy(id_1,id_0,natms*sizeof(int));*/
    /*
     ?. set vec_list to vec_list_
     */
    vec_list=vec_list_;
    
    if(vec_list->narch_vecs)
    {
        id_arch=new Vec<int>(this,1);
        int* id_0=id->begin();
        int* id_1=id_arch->begin();
        memcpy(id_1,id_0,natms*sizeof(int));
        del_vec(id_arch);
        for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
            del_vec(vec_list->arch_vecs[ivec]);
        
    }
    

    /*
     ?. add the x0 vector
     */
    x0=new Vec<type0>(this,dimension);
    
    /*
     ?. check wether the grid is created
     */
    // not done yet
    

    

    
    /*
     ?. make sure that all vec_list->xchng_vecs sizes are
     natms if not increase or decrease them.
     */
    for(int ivec=0;ivec<vec_list->nxchng_vecs;ivec++)
        vecs[ivec]->resize(natms);

    /*
     ?. create a new xchng
     */
    xchng=new Xchng(this,vec_list->xchng_vecs,vec_list->nxchng_vecs);
    /*
     ?. do x2s()
     */
    

    x2s(natms);

    /*
     ?. do the xchng
     */
    xchng->full_xchng();

    /* 
     ?. update natms
     */
    natms=x->vec_sz;
    /*
     ?. initiate forcefield,
     and obtain max_cut
     */
    forcefield->init();
    max_cut=forcefield->max_cut();
    /*
     ?. create the new swaps (it does calculate max_cut_s)
     */
    swap=new Swap(this,vec_list->updt_vecs,vec_list->nupdt_vecs);

    /*
     ?. set natms_ph to 0
     */
    natms_ph=0;
    /*
     ?. setup and find the phantom atoms and updating simultaniously
     */
    swap->list();
    /*
     ?. build the neighbor list and store x0
     */
    neighbor->init();
    swap->eliminate_redundancy();

    for(int ivec=vec_list->nxchng_vecs;ivec<nvecs;ivec++)
    {
        vecs[ivec]->vec_sz=0;
        vecs[ivec]->resize(natms+natms_ph);
    }
    
    type0* x_vec=x->begin();
    type0* x0_vec=x0->begin();
    int x_dim=x->dim;
    if(box_chng)
    {
        for(int iatm=0;iatm<natms+natms_ph;iatm++)
            for(int idim=0;idim<dimension;idim++)
                x0_vec[iatm*dimension+idim]=x_vec[iatm*x_dim+idim];
    }
    else
    {
        for(int iatm=0;iatm<natms;iatm++)
            for(int idim=0;idim<dimension;idim++)
                x0_vec[iatm*dimension+idim]=x_vec[iatm*x_dim+idim];
    }
    timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void Atoms::fin()
{
    delete xchng;
    delete swap;
    delete x0;

    if(vec_list->narch_vecs)
    {
        re_arrange(vec_list->arch_vecs,vec_list->narch_vecs);
        for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
            add_vec(vec_list->arch_vecs[ivec]);
        add_vec(id_arch);
        delete id_arch;
    }
    
    for(int ivec=0;ivec<nvecs;ivec++)
    {
        vecs[ivec]->vec_sz=natms;
        vecs[ivec]->shrink_to_fit();
    }
    
    neighbor->fin();
    forcefield->fin();
    timer->fin();
}
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void Atoms::update(vec* updt_vec)
{
    update(&updt_vec,1);
}
/*--------------------------------------------
 update a number of vectors
 --------------------------------------------*/
void Atoms::update(vec** updt_vecs,int nupdt_vecs)
{
    timer->start(COMM_TIME_mode);
    bool x_xst=false;
    for(int ivec=0;x_xst==false && ivec<nupdt_vecs;ivec++)
        if(updt_vecs[ivec]==x)
            x_xst=true;
    if(x_xst==false)
    {
        if(nupdt_vecs==1)
            swap->update(updt_vecs[0],false);
        else
            swap->update(updt_vecs,nupdt_vecs,false);
        timer->stop(COMM_TIME_mode);
        return;
    }
    
    
    type0 skin_sq=0.25*skin*skin;
    if(box_chng)
    {
        if(nupdt_vecs==1)
            swap->update(x,true);
        else
            swap->update(updt_vecs,nupdt_vecs,true);

        type0 dx_sq;
        int succ,succ_lcl=1;
        type0* x_vec=x->begin();
        type0* x0_vec=x0->begin();
        int x_dim=x->dim;
        for(int iatm=0;succ_lcl && iatm<natms+natms_ph;iatm++)
        {
            dx_sq=0.0;
            for(int idim=0;idim<dimension;idim++)
            {
                dx_sq+=(x0_vec[iatm*dimension+idim]-x_vec[iatm*x_dim+idim])
                *(x0_vec[iatm*dimension+idim]-x_vec[iatm*x_dim+idim]);
            }
            if(dx_sq>skin_sq)
                succ_lcl=0;
        }
        MPI_Allreduce(&succ_lcl,&succ,1,MPI_INT,MPI_MIN,world);
        if(succ)
        {
            timer->stop(COMM_TIME_mode);
            return;
        }
        
        x2s(natms);
        for(int ivec=0;ivec<vec_list->nxchng_vecs;ivec++)
            vecs[ivec]->resize(natms);
        xchng->full_xchng();
        natms=x->vec_sz;
        natms_ph=0;
        
        swap->reset();
        swap->list();
        
        neighbor->create_list(box_chng);
        
        natms_ph=x->vec_sz-natms;
        for(int ivec=vec_list->nxchng_vecs;ivec<nvecs;ivec++)
        {
            vecs[ivec]->vec_sz=0;
            vecs[ivec]->resize(natms+natms_ph);
        }
        
        x_vec=x->begin();
        x0_vec=x0->begin();
        for(int iatm=0;iatm<natms+natms_ph;iatm++)
            for(int idim=0;idim<dimension;idim++)
                x0_vec[iatm*dimension+idim]=x_vec[iatm*x_dim+idim];
        
    }
    else
    {   type0 dx_sq;
        int succ,succ_lcl=1;
        type0* x_vec=x->begin();
        type0* x0_vec=x0->begin();
        int x_dim=x->dim;
        for(int iatm=0;succ_lcl && iatm<natms;iatm++)
        {
            dx_sq=0.0;
            for(int idim=0;idim<dimension;idim++)
            {
                dx_sq+=(x0_vec[iatm*dimension+idim]-x_vec[iatm*x_dim+idim])
                *(x0_vec[iatm*dimension+idim]-x_vec[iatm*x_dim+idim]);
            }
            if(dx_sq>skin_sq)
                succ_lcl=0;
        }
        MPI_Allreduce(&succ_lcl,&succ,1,MPI_INT,MPI_MIN,world);
        
        if(succ)
        {
            if(nupdt_vecs==1)
                swap->update(x,true);
            else
                swap->update(updt_vecs,nupdt_vecs,true);
            
            timer->stop(COMM_TIME_mode);
            return;
        }
        
        x2s(natms);
        for(int ivec=0;ivec<vec_list->nxchng_vecs;ivec++)
            vecs[ivec]->resize(natms);
        xchng->full_xchng();
        natms=x->vec_sz;
        natms_ph=0;
        
        swap->list();
        
        neighbor->create_list(box_chng);
        
        for(int ivec=vec_list->nxchng_vecs;ivec<nvecs;ivec++)
        {
            vecs[ivec]->vec_sz=0;
            vecs[ivec]->resize(natms+natms_ph);
        }
        
        x_vec=x->begin();
        x0_vec=x0->begin();
        for(int iatm=0;iatm<natms;iatm++)
            for(int idim=0;idim<dimension;idim++)
                x0_vec[iatm*dimension+idim]=x_vec[iatm*x_dim+idim];
    }

    timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::init_xchng()
{
    timer->start(COMM_TIME_mode);
    x2s(natms);
    for(int ivec=0;ivec<vec_list->nxchng_vecs;ivec++)
        vecs[ivec]->resize(natms);
    xchng->full_xchng();
    natms=x->vec_sz;
    natms_ph=0;
    swap->list();
    for(int ivec=vec_list->nxchng_vecs;ivec<nvecs;ivec++)
    {
        vecs[ivec]->vec_sz=0;
        vecs[ivec]->resize(natms+natms_ph);
    }
    timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::fin_xchng()
{
    timer->start(COMM_TIME_mode);
    swap->list();
    neighbor->create_list(box_chng);
    
    for(int ivec=vec_list->nxchng_vecs;ivec<nvecs;ivec++)
    {
        vecs[ivec]->vec_sz=0;
        vecs[ivec]->resize(natms+natms_ph);
    }
    
    type0* x_vec=x->begin();
    type0* x0_vec=x0->begin();
    int x_dim=x->dim;

    x_vec=x->begin();
    x0_vec=x0->begin();
    for(int iatm=0;iatm<natms;iatm++)
        for(int idim=0;idim<dimension;idim++)
            x0_vec[iatm*dimension+idim]=x_vec[iatm*x_dim+idim];
    timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::set_skin(type0 skin_)
{
    skin=skin_;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 Atoms::get_skin()
{
    return skin;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::reset()
{
    xchng=new Xchng(this,vecs,nvecs);
    x2s(natms);
    xchng->full_xchng();
    natms=x->vec_sz;
    s2x(natms);
    delete xchng;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::re_arrange(vec** arch_vecs,int narch_vecs)
{
    auto fit_2_shrk=
    [](int*& arr,int sz,int cpcty)->void
    {
        if(cpcty==sz)
            return;
        if(sz==0)
        {
            delete [] arr;
            arr=NULL;
        }
        
        int* _arr=new int[sz];
        memcpy(_arr,arr,sz*sizeof(int));
        delete [] arr;
        arr=_arr;
    };
    
    auto sort_by_p=
    [] (int*& comm_size,MPI_Comm& world,int tot_p,int my_p,
        int* lst,int lst_sz,
        int* slst,int* slst_ref,int slst_sz)->void
    {
        comm_size=new int[tot_p];
        for(int ip=0;ip<tot_p;ip++)
            comm_size[ip]=0;
        
        int* lst_sz_per_p=new int[tot_p];
        int* _lst_=NULL;
        MPI_Allgather(&lst_sz,1,MPI_INT,lst_sz_per_p,1,MPI_INT,world);
        int _lst_sz_;
        int max_lst_sz=0;
        for(int ip=0;ip<tot_p;ip++)
            max_lst_sz=MAX(max_lst_sz,lst_sz_per_p[ip]);
        if(max_lst_sz)
            _lst_=new int[max_lst_sz];
        
        
        int* fnd_slst=NULL;
        int* ufnd_slst=NULL;
        int* fnd_slst_ref=NULL;
        int* ufnd_slst_ref=NULL;
        
        
        if(slst_sz)
        {
            fnd_slst=new int[slst_sz];
            ufnd_slst=new int[slst_sz];
            fnd_slst_ref=new int[slst_sz];
            ufnd_slst_ref=new int[slst_sz];
        }
        
        for(int ip=0;ip<tot_p;ip++)
        {
            _lst_sz_=lst_sz_per_p[ip];
            if(my_p==ip)
                memcpy(_lst_,lst,sizeof(int)*_lst_sz_);
            MPI_Bcast(_lst_,_lst_sz_,MPI_INT,ip,world);
            if(my_p==ip)
                continue;
            if(_lst_sz_ && slst_sz)
            {
                bool incmplt=true;
                int slst_pos=0,_lst_pos_=0;
                int nfnd=0;
                int nufnd=0;

                
                while (incmplt)
                {
                    while(incmplt && _lst_[_lst_pos_]!=slst[slst_pos])
                    {
                        if(_lst_[_lst_pos_]<slst[slst_pos])
                        {
                            while(incmplt && _lst_[_lst_pos_]<slst[slst_pos])
                            {
                                _lst_pos_++;
                                if(_lst_pos_==_lst_sz_)
                                    incmplt=false;
                            }
                        }
                        else
                        {
                            while(incmplt && slst[slst_pos]<_lst_[_lst_pos_])
                            {
                                ufnd_slst[nufnd]=slst[slst_pos];
                                ufnd_slst_ref[nufnd]=slst_ref[slst_pos];
                                nufnd++;
                                slst_pos++;
                                if(slst_pos==slst_sz)
                                    incmplt=false;
                            }
                        }
                    }
                    
                    while(incmplt && _lst_[_lst_pos_]==slst[slst_pos])
                    {
                        fnd_slst[nfnd]=slst[slst_pos];
                        fnd_slst_ref[nfnd]=slst_ref[slst_pos];
                        
                        nfnd++;
                        _lst_pos_++;
                        slst_pos++;
                        
                        if(_lst_pos_==_lst_sz_ || slst_pos==slst_sz)
                            incmplt=false;
                    }
                }
                
                memcpy(slst,fnd_slst,nfnd*sizeof(int));
                memcpy(slst+nfnd,ufnd_slst,nufnd*sizeof(int));
                
                memcpy(slst_ref,fnd_slst_ref,nfnd*sizeof(int));
                memcpy(slst_ref+nfnd,ufnd_slst_ref,nufnd*sizeof(int));
                
                slst+=nfnd;
                slst_ref+=nfnd;
                slst_sz-=nfnd;
                comm_size[ip]=nfnd;
            }
        }
        
        delete [] fnd_slst_ref;
        delete [] ufnd_slst_ref;
        delete [] fnd_slst;
        delete [] ufnd_slst;
        delete [] _lst_;
        delete [] lst_sz_per_p;
        
    };
    
    class act
    {
    public:
        virtual void eq(int)=0;
        virtual void neq(int)=0;
    };
    
    class act_arch: public act
    {
        /*
         for iact to manage: it manages the arch_id
         the arch ids are catagorized into 2 groups:
         0. ids that are found here and we need the map for their new place: orig_inds,natms_to_kep
         1. ids that are not found here and therefore must be sent to
         another proc: ids_to_snd,inds_to_snd,natms_to_snd
         */
    public:
        int*& idx_lst_arch;
        int*& id_lst_snd;
        int*& idx_lst_snd;
        
        int& nkep;
        int& nsnd;
        
        int* id_lst_arch;
        
        act_arch(
        int*& _id_lst_snd,
        int*& _idx_lst_snd,
        int*& _idx_lst_arch,
        int& _nkep,
        int& _nsnd,
        int* _id_lst_arch):
        id_lst_snd(_id_lst_snd),
        idx_lst_snd(_idx_lst_snd),
        idx_lst_arch(_idx_lst_arch),
        nkep(_nkep),
        nsnd(_nsnd),
        id_lst_arch(_id_lst_arch)
        {}
        void eq(int idx)
        {
            idx_lst_arch[nkep++]=idx;
        }
        void neq(int idx)
        {
            idx_lst_snd[nsnd]=idx;
            id_lst_snd[nsnd++]=id_lst_arch[idx];
        }
    };
    
    class act_curr: public act
    {
        /*
         for jact to manage: it manages the orig_id
         the current ids are catagorized into 2 groups:
         0. ids that are found here and we need the map for their new place: arch_inds,natms_to_kep
         1. ids that are not found here and therefore must be recieved from
         another proc: ids_to_rcv,inds_to_rcv,natms_to_rcv
         */
        
    public:
        int*& idx_lst_curr;
        int*& id_lst_rcv;
        int*& idx_lst_rcv;
        
        int& nkep;
        int& nrcv;
        
        int* id_lst_curr;
        
        act_curr(
        int*& _id_lst_rcv,
        int*& _idx_lst_rcv,
        int*& _idx_lst_curr,
        int& _nkep,
        int& _nrcv,
        int* _id_lst_curr):
        id_lst_rcv(_id_lst_rcv),
        idx_lst_rcv(_idx_lst_rcv),
        idx_lst_curr(_idx_lst_curr),
        nkep(_nkep),
        nrcv(_nrcv),
        id_lst_curr(_id_lst_curr)
        {}
        void eq(int idx)
        {
            idx_lst_curr[nkep++]=idx;
        }
        void neq(int idx)
        {
            idx_lst_rcv[nrcv]=idx;
            id_lst_rcv[nrcv++]=id_lst_curr[idx];
        }
    };
    
    // input
    int* id_lst_arch=id_arch->begin();
    int* id_lst_curr=id->begin();
    int natms_arch=id_arch->vec_sz;
    int natms_curr=id->vec_sz;
    int byte_sz=0;
    for(int ivec=0;ivec<narch_vecs;ivec++)
        byte_sz+=arch_vecs[ivec]->byte_sz;
    
    //output
    int* idx_lst_arch=NULL;
    int* idx_lst_snd=NULL;
    int* id_lst_snd=NULL;
    
    int* idx_lst_curr=NULL;
    int* idx_lst_rcv=NULL;
    int* id_lst_rcv=NULL;
    
    if(natms_arch)
    {
        idx_lst_arch=new int[natms_arch];
        idx_lst_snd=new int[natms_arch];
        id_lst_snd=new int[natms_arch];
    }
    
    if(natms_curr)
    {
        idx_lst_curr=new int[natms_curr];
        idx_lst_rcv=new int[natms_curr];
        id_lst_rcv=new int[natms_curr];
    }
    
    
    int nkep,nkep_arch=0,nkep_curr=0,nrcv=0,nsnd=0;
    
    act_arch* _act_arch_=new act_arch(id_lst_snd,idx_lst_snd,idx_lst_arch,nkep_arch,nsnd,id_lst_arch);
    act_curr* _act_curr_=new act_curr(id_lst_rcv,idx_lst_rcv,idx_lst_curr,nkep_curr,nrcv,id_lst_curr);
    act* _act_arch=_act_arch_;
    act* _act_curr=_act_curr_;
    XMath* xmath=new XMath();
    xmath->srch_lst_lst(id_lst_arch,natms_arch,_act_arch,id_lst_curr,natms_curr,_act_curr);
    delete _act_curr_;
    delete _act_arch_;
    delete xmath;
    
    nkep=nkep_curr;
    
    fit_2_shrk(idx_lst_arch,nkep_arch,natms_arch);
    fit_2_shrk(idx_lst_snd,nsnd,natms_arch);
    fit_2_shrk(id_lst_snd,nsnd,natms_arch);
    
    fit_2_shrk(idx_lst_curr,nkep_curr,natms_curr);
    fit_2_shrk(idx_lst_rcv,nrcv,natms_curr);
    fit_2_shrk(id_lst_rcv,nrcv,natms_curr);
    
    
    int* nsnd_comm;
    int* nrcv_comm;
    sort_by_p(nsnd_comm,world,tot_p,my_p,id_lst_rcv,nrcv,id_lst_snd,idx_lst_snd,nsnd);
    sort_by_p(nrcv_comm,world,tot_p,my_p,id_lst_snd,nsnd,id_lst_rcv,idx_lst_rcv,nrcv);
    
    delete [] id_lst_snd;
    delete [] id_lst_rcv;
    
    byte** snd_buff=new byte*[tot_p];
    byte** rcv_buff=new byte*[tot_p];
    *snd_buff=NULL;
    *rcv_buff=NULL;
    
    if (nsnd)
    {
        *snd_buff=new byte[byte_sz*nsnd];
        for(int ip=1;ip<tot_p;ip++)
            snd_buff[ip]=snd_buff[ip-1]+byte_sz*nsnd_comm[ip-1];
    }
    if(nrcv)
    {
        *rcv_buff=new byte[byte_sz*nrcv];
        for(int ip=1;ip<tot_p;ip++)
            rcv_buff[ip]=rcv_buff[ip-1]+byte_sz*nrcv_comm[ip-1];
    }
    
    byte* _snd_buff=*snd_buff;
    for(int ivec=0;ivec<narch_vecs;ivec++)
        for(int i=0;i<nsnd;i++)
            arch_vecs[ivec]->cpy(_snd_buff,idx_lst_snd[i]);
    
    delete [] idx_lst_snd;
    
    for(int ivec=0;ivec<narch_vecs;ivec++)
        arch_vecs[ivec]->rearrange(idx_lst_arch,idx_lst_curr,nkep,natms_curr);
    
    delete [] idx_lst_arch;
    delete [] idx_lst_curr;
    
    
    MPI_Status status[2];
    MPI_Request request[2];
    
    for(int idisp=1;idisp<tot_p;idisp++)
    {
        int rcv_p=my_p-idisp;
        if(rcv_p<0) rcv_p+=tot_p;
        int snd_p=my_p+idisp;
        if(snd_p>=tot_p) snd_p-=tot_p;
        
        if(nrcv_comm[rcv_p])
        {
            MPI_Irecv(rcv_buff[rcv_p],nrcv_comm[rcv_p]*byte_sz,
                      MPI_BYTE,rcv_p,
                      0,
                      world,&request[0]);
        }
        
        if(nsnd_comm[snd_p])
        {
            MPI_Isend(snd_buff[snd_p],nsnd_comm[snd_p]*byte_sz,
                      MPI_BYTE,snd_p,
                      0,
                      world,&request[1]);
        }
        
        if(nrcv_comm[rcv_p] && nsnd_comm[snd_p])
            MPI_Waitall(2,request,status);
        else if(nrcv_comm[rcv_p] && nsnd_comm[snd_p]==0)
            MPI_Waitall(1,&request[0],&status[0]);
        else if(nrcv_comm[rcv_p]==0 && nsnd_comm[snd_p])
            MPI_Waitall(1,&request[1],&status[1]);
        
    }
    
    delete [] *snd_buff;
    delete [] snd_buff;
    delete [] nsnd_comm;
    delete [] nrcv_comm;
    
    byte* _rcv_buff=*rcv_buff;
    for(int ivec=0;ivec<narch_vecs;ivec++)
        for(int i=0;i<nrcv;i++)
            arch_vecs[ivec]->cpy(_rcv_buff,idx_lst_rcv[i]);
    
    delete [] idx_lst_rcv;
    delete [] *rcv_buff;
    delete [] rcv_buff;
}
/*----------------------------------------------
  _     _   _____   _____   _       _____   _____
 | |   / / | ____| /  ___| | |     /  ___/ |_   _|
 | |  / /  | |__   | |     | |     | |___    | |
 | | / /   |  __|  | |     | |     \___  \   | |
 | |/ /    | |___  | |___  | |___   ___| |   | |
 |___/     |_____| \_____| |_____| /_____/   |_|
 
 ----------------------------------------------*/
VecLst::VecLst(Atoms* atoms):
nxchng_vecs(2),
nupdt_vecs(1)
{
    xchng_vecs=new vec*[nxchng_vecs];
    xchng_vecs[0]=atoms->x;
    xchng_vecs[1]=atoms->id;
    tot_xchng_byte_sz=atoms->x->byte_sz
    +atoms->id->byte_sz;
    
    updt_vecs=new vec*[nupdt_vecs];
    updt_vecs[0]=atoms->x;
    tot_updt_byte_sz=atoms->x->byte_sz;
    
    narch_vecs=0;
    tot_arch_byte_sz=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
VecLst::~VecLst()
{
    if(nxchng_vecs)
        delete [] xchng_vecs;
    if(nupdt_vecs)
        delete [] updt_vecs;
    if(narch_vecs)
        delete [] arch_vecs;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::add_updt(vec* v)
{
    int ivec=0;
    while(ivec<nupdt_vecs && updt_vecs[ivec]!=v)
        ivec++;
    
    if(ivec!=nupdt_vecs)
        return;
    
    add_xchng(v);
    
    vec** updt_vecs_=new vec*[nupdt_vecs+1];
    memcpy(updt_vecs_,updt_vecs,nupdt_vecs*sizeof(vec*));
    if(nupdt_vecs)
        delete [] updt_vecs;
    updt_vecs=updt_vecs_;
    
    updt_vecs[nupdt_vecs]=v;
    tot_updt_byte_sz+=v->byte_sz;
    nupdt_vecs++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::del_updt(vec* v)
{
    if(v==updt_vecs[0])
        return;
    
    int ivec=0;
    while(ivec<nupdt_vecs && updt_vecs[ivec]!=v)
        ivec++;
    if(ivec==nupdt_vecs)
        return;
    
    vec** updt_vecs_=new vec*[nupdt_vecs-1];
    for(int jvec=0;jvec<ivec;jvec++)
        updt_vecs_[jvec]=updt_vecs[jvec];
    for(int jvec=ivec+1;jvec<nupdt_vecs;jvec++)
        updt_vecs_[jvec-1]=updt_vecs[jvec];
    if(nupdt_vecs)
        delete [] updt_vecs;
    updt_vecs=updt_vecs_;
    
    tot_updt_byte_sz-=v->byte_sz;
    nupdt_vecs--;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::add_xchng(vec* v)
{
    int ivec=0;
    while(ivec<nxchng_vecs && xchng_vecs[ivec]!=v)
        ivec++;
    
    if(ivec!=nxchng_vecs)
        return;
    
    vec** xchng_vecs_=new vec*[nxchng_vecs+1];
    memcpy(xchng_vecs_,xchng_vecs,nxchng_vecs*sizeof(vec*));
    if(nxchng_vecs)
        delete [] xchng_vecs;
    xchng_vecs=xchng_vecs_;
    
    xchng_vecs[nxchng_vecs]=v;
    tot_xchng_byte_sz+=v->byte_sz;
    nxchng_vecs++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::del_xchng(vec* v)
{
    del_updt(v);
    if(v==xchng_vecs[0] || v==xchng_vecs[1])
        return;
    
    int ivec=0;
    while(ivec<nxchng_vecs && xchng_vecs[ivec]!=v)
        ivec++;
    if(ivec==nxchng_vecs)
        return;
    
    vec** xchng_vecs_=new vec*[nxchng_vecs-1];
    for(int jvec=0;jvec<ivec;jvec++)
        xchng_vecs_[jvec]=xchng_vecs[jvec];
    for(int jvec=ivec+1;jvec<nxchng_vecs;jvec++)
        xchng_vecs_[jvec-1]=xchng_vecs[jvec];
    if(nxchng_vecs)
        delete [] xchng_vecs;
    xchng_vecs=xchng_vecs_;
    
    tot_xchng_byte_sz-=v->byte_sz;
    nxchng_vecs--;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::add_arch(vec* v)
{
    int ivec=0;
    try
    {
        
        while(ivec<nxchng_vecs && xchng_vecs[ivec]!=v)
            ivec++;
        if(ivec!=nxchng_vecs)
            throw 1;
    }
    catch(int i)
    {
        printf("archive vector cannot be exchange vector\n");
    }

    
    ivec=0;
    while(ivec<narch_vecs && arch_vecs[ivec]!=v)
        ivec++;
    if(ivec!=narch_vecs)
        return;
    
    vec** arch_vecs_=new vec*[narch_vecs+1];
    memcpy(arch_vecs_,arch_vecs,narch_vecs*sizeof(vec*));
    if(narch_vecs)
        delete [] arch_vecs;
    arch_vecs=arch_vecs_;
    
    arch_vecs[narch_vecs]=v;
    tot_arch_byte_sz+=v->byte_sz;
    narch_vecs++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::del_arch(vec* v)
{
    int ivec=0;
    while(ivec<narch_vecs && arch_vecs[ivec]!=v)
        ivec++;
    if(ivec==narch_vecs)
        return;
    
    vec** arch_vecs_=new vec*[narch_vecs-1];
    for(int jvec=0;jvec<ivec;jvec++)
        arch_vecs_[jvec]=arch_vecs[jvec];
    for(int jvec=ivec+1;jvec<narch_vecs;jvec++)
        arch_vecs_[jvec-1]=arch_vecs[jvec];
    if(narch_vecs)
        delete [] arch_vecs;
    arch_vecs=arch_vecs_;
    
    tot_arch_byte_sz-=v->byte_sz;
    narch_vecs--;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::add_updt(vec** vecs,int nvecs)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        add_updt(vecs[ivec]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::del_updt(vec** vecs,int nvecs)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        del_updt(vecs[ivec]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::add_xchng(vec** vecs,int nvecs)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        add_xchng(vecs[ivec]);
}/*--------------------------------------------
  
  --------------------------------------------*/
void VecLst::del_xchng(vec** vecs,int nvecs)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        del_xchng(vecs[ivec]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::add_arch(vec** vecs,int nvecs)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        add_arch(vecs[ivec]);
}/*--------------------------------------------
  
  --------------------------------------------*/
void VecLst::del_arch(vec** vecs,int nvecs)
{
    for(int ivec=0;ivec<nvecs;ivec++)
        del_arch(vecs[ivec]);
}

