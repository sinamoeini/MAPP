#include "dynamic.h"
#include "comm.h"
#include "atoms.h"
#include "xmath.h"
#include "error.h"
#include "timer.h"
#include "MAPP.h"
#include "ff.h"
#include "neighbor.h"
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
Dynamic::Dynamic(Atoms*& atoms_,
Communication*& comm_,
vec* const * updt_vecs_,int nupdt_vecs_,
vec* const * xchng_vecs_,int nxchng_vecs_,
vec* const * arch_vecs_,int narch_vecs_):
comm(comm_),
atoms(atoms_)
{
    auto is_in=[](const vec* v,vec* const * vs,int nvs)->bool
    {
        for(int i=0;i<nvs;i++)
            if(vs[i]==v) return true;
        return false;
    };
    
    int new_nvecs=atoms->nvecs-narch_vecs_;
    vec** new_vecs=new vec*[new_nvecs];
    new_vecs[0]=atoms->x;
    memcpy(new_vecs+1,updt_vecs_,nupdt_vecs_*sizeof(vec*));
    new_vecs[nupdt_vecs_+1]=atoms->id;
    memcpy(new_vecs+nupdt_vecs_+2,xchng_vecs_,nxchng_vecs_*sizeof(vec*));
    int ivec=new_nvecs-1;
    vec** vecs=atoms->vecs;
    
    
    for(int i=0;i<atoms->nvecs;i++)
    {
        if(is_in(vecs[i],new_vecs,new_nvecs) || is_in(vecs[i],arch_vecs_,narch_vecs_))
            continue;
        new_vecs[ivec--]=vecs[i];
    }
    
    nxchng_vecs=nxchng_vecs_+nupdt_vecs_+2;
    nupdt_vecs=nupdt_vecs_+1;
    delete [] atoms->vecs;
    atoms->vecs=new_vecs;
    atoms->nvecs=new_nvecs;
    

    
    narch_vecs=narch_vecs_;
    arch_vecs=NULL;
    if(!narch_vecs) return;
    
    arch_vecs=new vec*[narch_vecs];
    memcpy(arch_vecs,arch_vecs_,narch_vecs*sizeof(vec*));
}
/*--------------------------------------------
 
 --------------------------------------------*/
Dynamic::Dynamic(Atoms*& atoms_,Communication*& comm_,
std::initializer_list<vec*> v0,std::initializer_list<vec*> v1
,std::initializer_list<vec*> v2):
Dynamic(atoms_,comm_,v0.begin(),static_cast<int>(v0.size()),v1.begin(),static_cast<int>(v1.size()),v2.begin(),static_cast<int>(v2.size()))
{
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
Dynamic::Dynamic(Atoms*& atoms_,Communication*& comm_,
std::initializer_list<vec*> v0,std::initializer_list<vec*> v1):
Dynamic(atoms_,comm_,v0.begin(),static_cast<int>(v0.size()),v1.begin(),static_cast<int>(v1.size()),NULL,0)
{
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::add_xchng(vec* v)
{
    vec** vecs=atoms->vecs;
    int nvces=atoms->nvecs;
    int ivec=nxchng_vecs;
    for(;vecs[ivec]!=v && ivec<nvces;ivec++){}
    vecs[ivec]=vecs[nxchng_vecs];
    vecs[nxchng_vecs++]=v;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::add_updt(vec* v)
{
    vec** vecs=atoms->vecs;
    int nvces=atoms->nvecs;
    int ivec=nupdt_vecs;
    for(;vecs[ivec]!=v && ivec<nvces;ivec++){}
    vecs[ivec]=vecs[nupdt_vecs];
    vecs[nupdt_vecs++]=v;
    if(ivec<nxchng_vecs)
        nxchng_vecs++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
Dynamic::~Dynamic()
{
    delete [] arch_vecs;
}
/*--------------------------------------------
 init a simulation
 --------------------------------------------*/
void Dynamic::init(bool box_chng_)
{
    type0 (&H)[__dim__][__dim__]=atoms->H;
    type0 det=1.0;
    for(int idim=0;idim<__dim__;idim++)
        det*=H[idim][idim]*H[idim][idim];
    if(det==0.0)
        Error::abort("cannot start the simulation "
        "domain griding, the domain is not defined");
    if(!comm->grid_established)
        Error::abort("cannot start the simulation "
        "domain grid has not been established");
    

    if(!forcefield)
        Error::abort("cannot start the simulation, force field has not been established");

    
    mapp->timer->init();
    mapp->timer->start(COMM_TIME_mode);

    
    box_chng=box_chng_;
    
    store_arch_vecs();
    
    x0=new Vec<type0>(atoms,__dim__);
    xchng=new Exchange(atoms,comm,nxchng_vecs);
    atoms->x2s_lcl();
    xchng->full_xchng();

    
    forcefield->init();
    atoms->max_cut=forcefield->max_cut();
    updt=new Update(atoms,comm,nupdt_vecs,nxchng_vecs);
    
    updt->list();
    forcefield->neighbor->init();
    updt->rm_rdndncy();

    store_x0();
    mapp->timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void Dynamic::fin()
{

    delete x0;

    restore_arch_vecs();
    vec** vecs=atoms->vecs;
    int nvecs=atoms->nvecs;
    int natms=atoms->natms;
    for(int ivec=0;ivec<nvecs;ivec++)
    {
        vecs[ivec]->vec_sz=natms;
        vecs[ivec]->shrink_to_fit();
    }

    forcefield->neighbor->fin();
    forcefield->fin();
    mapp->timer->fin();
    delete xchng;
    delete updt;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::store_x0()
{
    type0* x_vec=atoms->x->begin();
    type0* x0_vec=x0->begin();
    int x_dim=atoms->x->dim;
    int last_atm=atoms->natms;
    if(box_chng) last_atm+=atoms->natms_ph;
    if(x_dim==__dim__)
    {
        memcpy(x0_vec,x_vec,last_atm*__dim__*sizeof(type0));
        return;
    }
    for(int iatm=0;iatm<last_atm;iatm++,x0_vec+=__dim__,x_vec+=x_dim)
        memcpy(x0_vec,x_vec,__dim__*sizeof(type0));
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline bool Dynamic::decide()
{
    type0 skin_sq=0.25*comm->skin*comm->skin;
    type0 dx_sq;
    int succ,succ_lcl=1;
    type0* x_vec=atoms->x->begin();
    type0* x0_vec=x0->begin();
    int x_dim=atoms->x->dim;
    int last_atm=atoms->natms;
    if(box_chng) last_atm+=atoms->natms_ph;
    
    for(int iatm=0;succ_lcl && iatm<last_atm;iatm++,x0_vec+=__dim__,x_vec+=x_dim)
    {
        dx_sq=0.0;
        for(int idim=0;idim<__dim__;idim++)
            dx_sq+=(x0_vec[idim]-x_vec[idim])*(x0_vec[idim]-x_vec[idim]);
        if(dx_sq>skin_sq)
            succ_lcl=0;
    }

    MPI_Allreduce(&succ_lcl,&succ,1,MPI_INT,MPI_MIN,comm->world);
    if(succ) return true;
    return false;
}
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void Dynamic::update(vec* updt_vec)
{
    update(&updt_vec,1);
}
/*--------------------------------------------
 update a number of vectors
 --------------------------------------------*/
void Dynamic::update(vec** updt_vecs,int nupdt_vecs)
{
    mapp->timer->start(COMM_TIME_mode);
    bool x_xst=false;
    for(int ivec=0;x_xst==false && ivec<nupdt_vecs;ivec++)
        if(updt_vecs[ivec]==atoms->x)
            x_xst=true;
    if(x_xst==false)
    {
        if(nupdt_vecs==1)
            updt->update(updt_vecs[0],false);
        else
            updt->update(updt_vecs,nupdt_vecs,false);
        mapp->timer->stop(COMM_TIME_mode);
        return;
    }
    
    
    if(box_chng)
    {
        if(nupdt_vecs==1)
            updt->update(atoms->x,true);
        else
            updt->update(updt_vecs,nupdt_vecs,true);

        if(decide())
        {
            mapp->timer->stop(COMM_TIME_mode);
            return;
        }
        
        atoms->x2s_lcl();
        xchng->full_xchng();
        
        updt->reset();
        updt->list();
        forcefield->neighbor->create_list(box_chng);
        store_x0();
    }
    else
    {
        if(decide())
        {
            
            if(nupdt_vecs==1)
                updt->update(atoms->x,true);
            else
                updt->update(updt_vecs,nupdt_vecs,true);
            mapp->timer->stop(COMM_TIME_mode);
            return;
        }
        
        atoms->x2s_lcl();
        xchng->full_xchng();
        
        updt->list();
        forcefield->neighbor->create_list(box_chng);
        
        store_x0();
    }

    mapp->timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::init_xchng()
{
    mapp->timer->start(COMM_TIME_mode);
    atoms->x2s_lcl();
    xchng->full_xchng();
    updt->list();

    mapp->timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::fin_xchng()
{
    mapp->timer->start(COMM_TIME_mode);

    updt->list();
    forcefield->neighbor->create_list(box_chng);
    store_x0();
    mapp->timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::store_arch_vecs()
{
    if(!narch_vecs) return;
    id_arch=new Vec<int>(atoms,1);
    int* id_0=atoms->id->begin();
    int* id_1=id_arch->begin();
    memcpy(id_1,id_0,atoms->natms*sizeof(int));
    atoms->del_vec(id_arch);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::restore_arch_vecs()
{
    if(!narch_vecs) return;
    
    const int tot_p=comm->tot_p;
    const int my_p=comm->my_p;
    MPI_Comm& world=comm->world;
    
    int byte_sz=0;
    for(int ivec=0;ivec<narch_vecs;ivec++)
        byte_sz+=arch_vecs[ivec]->byte_sz;
    
    int natms_old=id_arch->vec_sz;
    int* id_old=id_arch->begin();
    int* key_old=NULL;
    if(natms_old) key_old=new int[natms_old];
    for(int i=0;i<natms_old;i++) key_old[i]=i;
    XMath::quicksort(key_old,key_old+natms_old
    ,[&id_old](int* rank_i,int* rank_j){return (id_old[*rank_i]<id_old[*rank_j]);}
    ,[](int* rank_i,int* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    int natms_new=atoms->id->vec_sz;
    int* id_new=atoms->id->begin();
    int* key_new=NULL;
    if(natms_new) key_new=new int[natms_new];
    for(int i=0;i<natms_new;i++) key_new[i]=i;
    XMath::quicksort(key_new,key_new+natms_new
    ,[&id_new](int* rank_i,int* rank_j){return (id_new[*rank_i]<id_new[*rank_j]);}
    ,[](int* rank_i,int* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    
    int nsnd=0;
    int* idx_lst_snd=NULL;
    int* id_lst_snd=NULL;
    if(natms_old)
    {
        idx_lst_snd=new int[natms_old];
        id_lst_snd=new int[natms_old];
    }
    
    int nrcv=0;
    int* idx_lst_rcv=NULL;
    int* id_lst_rcv=NULL;
    if(natms_new)
    {
        idx_lst_rcv=new int[natms_new];
        id_lst_rcv=new int[natms_new];
    }
    
    int nkeep=0;
    int* idx_keep_old=NULL;
    int* idx_keep_new=NULL;
    if(MIN(natms_old,natms_new))
    {
        idx_keep_old=new int[MIN(natms_old,natms_new)];
        idx_keep_new=new int[MIN(natms_old,natms_new)];
    }
    
    int i_old=0,i_new=0;
    for(;i_old<natms_old && i_new<natms_new;)
    {
        if(id_old[key_old[i_old]]<id_new[key_new[i_new]])
        {
            id_lst_snd[nsnd]=id_old[key_old[i_old]];
            idx_lst_snd[nsnd++]=key_old[i_old++];
        }
        else if(id_old[key_old[i_old]]>id_new[key_new[i_new]])
        {
            id_lst_rcv[nrcv]=id_new[key_new[i_new]];
            idx_lst_rcv[nrcv++]=key_new[i_new++];
        }
        else
        {
            idx_keep_old[nkeep]=key_old[i_old++];
            idx_keep_new[nkeep++]=key_new[i_new++];
        }
    }
    
    for(;i_old<natms_old;)
    {
        id_lst_snd[nsnd]=id_old[key_old[i_old]];
        idx_lst_snd[nsnd++]=key_old[i_old++];
    }
    for(;i_new<natms_new;)
    {
        id_lst_rcv[nrcv]=id_new[key_new[i_new]];
        idx_lst_rcv[nrcv++]=key_new[i_new++];
    }
    
    delete [] key_old;
    delete [] key_new;
    
    /*-------------------------------------------------*/
    
    auto shrnk_2_fit=
    [](int*& arr,const int sz,const int cpcty)->void
    {
        if(cpcty==0 || cpcty==sz)
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
    
    /*---------------------------------------------------*/
    
    shrnk_2_fit(id_lst_snd,nsnd,natms_old);
    shrnk_2_fit(idx_lst_snd,nsnd,natms_old);
    shrnk_2_fit(id_lst_rcv,nrcv,natms_new);
    shrnk_2_fit(idx_lst_rcv,nrcv,natms_new);
    shrnk_2_fit(idx_keep_old,nkeep,MIN(natms_old,natms_new));
    shrnk_2_fit(idx_keep_new,nkeep,MIN(natms_old,natms_new));

#ifdef ARCH_RESTORE_DEBUG
    int tot_nsnd;
    int tot_nrcv;
    MPI_Allreduce(&nsnd,&tot_nsnd,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&nrcv,&tot_nrcv,1,MPI_INT,MPI_SUM,world);
    if(tot_nsnd!=tot_nrcv)
        Error::abort("sanity check failed number of atoms to be sent(%d) "
        "is unequal to number of atoms to be recieved(%d)",tot_nsnd,tot_nrcv);
#endif
    
    auto sort_idx_by_p=
    [&world,&my_p,&tot_p] (const int* lst,int lst_sz,const int* sub_lst,int sub_lst_sz,
            int*& sub_lst_idx,int* comm_size)->void
    {
        for(int ip=0;ip<tot_p;ip++)
            comm_size[ip]=0;
        
        int max_lst_sz;
        MPI_Allreduce(&lst_sz,&max_lst_sz,1,MPI_INT,MPI_MAX,world);
        int mother_lst_sz;
        int* mother_lst=NULL;
        if(max_lst_sz) mother_lst=new int[max_lst_sz];
        
        int fnd_sz,ufnd_sz;
        
        int* ufnd=NULL;
        int* ufnd_idx=NULL;
        int* fnd_idx=NULL;
        if(sub_lst_sz)
        {
            ufnd=new int[sub_lst_sz];
            ufnd_idx=new int[sub_lst_sz];
            fnd_idx=new int[sub_lst_sz];
        }
        
        fnd_sz=0;
        ufnd_sz=sub_lst_sz;
        memcpy(ufnd,sub_lst,sub_lst_sz*sizeof(int));
        memcpy(ufnd_idx,sub_lst_idx,sub_lst_sz*sizeof(int));
        for(int imother_lst,iufnd,ufnd_tmp_sz,ip=0;ip<tot_p;ip++)
        {
            if(ip==my_p)
            {
                mother_lst_sz=lst_sz;
                memcpy(mother_lst,lst,mother_lst_sz*sizeof(int));
            }
            MPI_Bcast(&mother_lst_sz,1,MPI_INT,ip,world);
            MPI_Bcast(mother_lst,mother_lst_sz,MPI_INT,ip,world);
            
            if(ip==my_p || ufnd_sz==0 || mother_lst_sz==0) continue;
            
            ufnd_tmp_sz=0;
            for(imother_lst=0,iufnd=0;imother_lst<mother_lst_sz && iufnd<ufnd_sz;)
            {
                if(mother_lst[imother_lst]<ufnd[iufnd])
                    imother_lst++;
                else if(mother_lst[imother_lst]>ufnd[iufnd])
                {
                    ufnd_idx[ufnd_tmp_sz]=ufnd_idx[iufnd];
                    ufnd[ufnd_tmp_sz++]=ufnd[iufnd++];
                }
                else
                {
                    fnd_idx[fnd_sz++]=ufnd_idx[iufnd++];
                    imother_lst++;
                    comm_size[ip]++;
                }
            }
            
            for(;iufnd<ufnd_sz;)
            {
                ufnd_idx[ufnd_tmp_sz]=ufnd_idx[iufnd];
                ufnd[ufnd_tmp_sz++]=ufnd[iufnd++];
            }
            
            ufnd_sz=ufnd_tmp_sz;
            
        }
        
        delete [] sub_lst_idx;
        sub_lst_idx=fnd_idx;
        delete [] ufnd;
        delete [] ufnd_idx;
        delete [] mother_lst;
    };
   
    int* nsnd_comm=new int[tot_p];
    int* nrcv_comm=new int[tot_p];
    sort_idx_by_p(id_lst_rcv,nrcv,id_lst_snd,nsnd,idx_lst_snd,nsnd_comm);
    sort_idx_by_p(id_lst_snd,nsnd,id_lst_rcv,nrcv,idx_lst_rcv,nrcv_comm);
    delete [] id_lst_snd;
    delete [] id_lst_rcv;
 
#ifdef ARCH_RESTORE_DEBUG
    int comm_sz_snd=0,comm_sz_rcv=0;
    for(int i=0;i<tot_p;i++)
    {
        comm_sz_snd+=nsnd_comm[i];
        comm_sz_rcv+=nrcv_comm[i];
    }
    if(nsnd!=comm_sz_snd)
        Error::abort_sing("sanity check failed just found %d out of %d atoms that should be sent");
    if(nrcv!=comm_sz_rcv)
        Error::abort_sing("sanity check failed just found %d out of %d atoms that should be recieved");
#endif
    
    byte** snd_buff=new byte*[tot_p];
    *snd_buff=NULL;
    if (nsnd)
    {
        *snd_buff=new byte[byte_sz*nsnd];
        for(int ip=1;ip<tot_p;ip++)
            snd_buff[ip]=snd_buff[ip-1]+byte_sz*nsnd_comm[ip-1];
    }
    
    byte** rcv_buff=new byte*[tot_p];
    *rcv_buff=NULL;
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
        arch_vecs[ivec]->rearrange(idx_keep_old,idx_keep_new,nkeep,natms_new);
    delete [] idx_keep_old;
    delete [] idx_keep_new;

#ifndef BLOCK_COMM
    MPI_Request request[2];
#endif
    for(int idisp=1;idisp<tot_p;idisp++)
    {
        int rcv_p=my_p-idisp;
        if(rcv_p<0) rcv_p+=tot_p;
        int snd_p=my_p+idisp;
        if(snd_p>=tot_p) snd_p-=tot_p;
#ifdef BLOCK_COMM
        MPI_Sendrecv(snd_buff[snd_p],nsnd_comm[snd_p]*byte_sz,MPI_BYTE,snd_p,0,
                     rcv_buff[rcv_p],nrcv_comm[rcv_p]*byte_sz,MPI_BYTE,rcv_p,0,
                     world,MPI_STATUS_IGNORE);
#else
        if(nrcv_comm[rcv_p])
            MPI_Irecv(rcv_buff[rcv_p],nrcv_comm[rcv_p]*byte_sz,
                      MPI_BYTE,rcv_p,
                      0,
                      world,&request[0]);
        
        if(nsnd_comm[snd_p])
            MPI_Isend(snd_buff[snd_p],nsnd_comm[snd_p]*byte_sz,
                      MPI_BYTE,snd_p,
                      0,
                      world,&request[1]);
        
        if(nrcv_comm[rcv_p] && nsnd_comm[snd_p])
            MPI_Waitall(2,request,MPI_STATUS_IGNORE);
        else if(nrcv_comm[rcv_p] && nsnd_comm[snd_p]==0)
            MPI_Waitall(1,&request[0],MPI_STATUS_IGNORE);
        else if(nrcv_comm[rcv_p]==0 && nsnd_comm[snd_p])
            MPI_Waitall(1,&request[1],MPI_STATUS_IGNORE);
#endif
    }
    delete [] nsnd_comm;
    delete [] nrcv_comm;
    
    delete [] *snd_buff;
    delete [] snd_buff;
    
    byte* _rcv_buff=*rcv_buff;
    for(int ivec=0;ivec<narch_vecs;ivec++)
        for(int i=0;i<nrcv;i++)
            arch_vecs[ivec]->pst_to(_rcv_buff,idx_lst_rcv[i]);
    
    delete [] idx_lst_rcv;
    delete [] *rcv_buff;
    delete [] rcv_buff;
    
    for(int ivec=0;ivec<narch_vecs;ivec++)
        atoms->add_vec(arch_vecs[ivec]);
    atoms->add_vec(id_arch);
    delete id_arch;
    id_arch=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Dynamic::restore_arch_vecs_()
{
    if(!narch_vecs) return;
    
    const int tot_p=comm->tot_p;
    const int my_p=comm->my_p;
    MPI_Comm& world=comm->world;
    
    auto shrnk_2_fit=
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
    
    auto sort_idx_by_p=
    [&world,&my_p,&tot_p] (const int* lst,int lst_sz,const int* sub_lst,int sub_lst_sz,
            int*& sub_lst_idx,int* comm_size)->void
    {
        for(int ip=0;ip<tot_p;ip++)
            comm_size[ip]=0;
        
        int max_lst_sz;
        MPI_Allreduce(&lst_sz,&max_lst_sz,1,MPI_INT,MPI_MAX,world);
        int* mother_lst=NULL;
        int mother_lst_sz;
        if(max_lst_sz)
            mother_lst=new int[max_lst_sz];
        
        int fnd_sz=0,ufnd_sz=sub_lst_sz,ufnd_tmp_sz;
        int* ufnd=NULL;
        int* ufnd_tmp=NULL;
        int* fnd_idx=NULL;
        int* ufnd_idx=NULL;
        int* ufnd_idx_tmp=NULL;
        if(sub_lst_sz)
        {
            ufnd=new int[sub_lst_sz];
            ufnd_tmp=new int[sub_lst_sz];
            fnd_idx=new int[sub_lst_sz];
            ufnd_idx=new int[sub_lst_sz];
            ufnd_idx_tmp=new int[sub_lst_sz];
        }
        
        memcpy(ufnd,sub_lst,sub_lst_sz*sizeof(int));
        memcpy(ufnd_idx,sub_lst_idx,sub_lst_sz*sizeof(int));
        
        
        int imother_lst,iufnd;
        for(int ip=0;ip<tot_p;ip++)
        {
            if(ip==my_p)
            {
                mother_lst_sz=lst_sz;
                memcpy(mother_lst,lst,mother_lst_sz*sizeof(int));
            }
            MPI_Bcast(&mother_lst_sz,1,MPI_INT,ip,world);
            MPI_Bcast(mother_lst,mother_lst_sz,MPI_INT,ip,world);
            
            if(ip==my_p || ufnd_sz==0 || mother_lst_sz==0) continue;
            
            ufnd_tmp_sz=0;
            for(imother_lst=0,iufnd=0;imother_lst<mother_lst_sz && iufnd<ufnd_sz;)
            {
                if(mother_lst[imother_lst]<ufnd[iufnd])
                    imother_lst++;
                else if(mother_lst[imother_lst]>ufnd[iufnd])
                {
                    ufnd_idx_tmp[ufnd_tmp_sz]=ufnd_idx[iufnd];
                    ufnd_tmp[ufnd_tmp_sz++]=ufnd[iufnd++];
                }
                else
                {
                    fnd_idx[fnd_sz++]=ufnd_idx[iufnd++];
                    imother_lst++;
                    comm_size[ip]++;
                }
            }
            
            for(;iufnd<ufnd_sz;)
            {
                ufnd_idx_tmp[ufnd_tmp_sz]=ufnd_idx[iufnd];
                ufnd_tmp[ufnd_tmp_sz++]=ufnd[iufnd++];
            }
            
            std::swap(ufnd_tmp,ufnd);
            std::swap(ufnd_idx_tmp,ufnd_idx);
            ufnd_sz=ufnd_tmp_sz;
            
        }
        
        delete [] sub_lst_idx;
        sub_lst_idx=fnd_idx;
        
        delete [] ufnd;
        delete [] ufnd_tmp;
        delete [] ufnd_idx;
        delete [] ufnd_idx_tmp;
        delete [] mother_lst;
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
    int* id_lst_curr=atoms->id->begin();
    int natms_arch=id_arch->vec_sz;
    int natms_curr=atoms->id->vec_sz;
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
    XMath::srch_lst_lst(id_lst_arch,natms_arch,_act_arch,id_lst_curr,natms_curr,_act_curr);
    delete _act_curr_;
    delete _act_arch_;
    
    nkep=nkep_curr;
    
    shrnk_2_fit(idx_lst_arch,nkep_arch,natms_arch);
    shrnk_2_fit(idx_lst_snd,nsnd,natms_arch);
    shrnk_2_fit(id_lst_snd,nsnd,natms_arch);
    
    shrnk_2_fit(idx_lst_curr,nkep_curr,natms_curr);
    shrnk_2_fit(idx_lst_rcv,nrcv,natms_curr);
    shrnk_2_fit(id_lst_rcv,nrcv,natms_curr);
    
    
    int* nsnd_comm=new int[tot_p];
    int* nrcv_comm=new int[tot_p];
    sort_idx_by_p(id_lst_rcv,nrcv,id_lst_snd,nsnd,idx_lst_snd,nsnd_comm);
    sort_idx_by_p(id_lst_snd,nsnd,id_lst_rcv,nrcv,idx_lst_rcv,nrcv_comm);
    delete [] id_lst_snd;
    delete [] id_lst_rcv;
    
    byte** snd_buff=new byte*[tot_p];
    *snd_buff=NULL;
    if (nsnd)
    {
        *snd_buff=new byte[byte_sz*nsnd];
        for(int ip=1;ip<tot_p;ip++)
            snd_buff[ip]=snd_buff[ip-1]+byte_sz*nsnd_comm[ip-1];
    }
    
    byte** rcv_buff=new byte*[tot_p];
    *rcv_buff=NULL;
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
    
    for(int idisp=1;idisp<tot_p;idisp++)
    {
        int rcv_p=my_p-idisp;
        if(rcv_p<0) rcv_p+=tot_p;
        int snd_p=my_p+idisp;
        if(snd_p>=tot_p) snd_p-=tot_p;
        MPI_Sendrecv(snd_buff[snd_p],nsnd_comm[snd_p]*byte_sz,MPI_BYTE,snd_p,0,
                     rcv_buff[rcv_p],nrcv_comm[rcv_p]*byte_sz,MPI_BYTE,rcv_p,0,
                     world,MPI_STATUS_IGNORE);
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
    
    
    for(int ivec=0;ivec<narch_vecs;ivec++)
        atoms->add_vec(arch_vecs[ivec]);
    atoms->add_vec(id_arch);
    delete id_arch;
    id_arch=NULL;
}
