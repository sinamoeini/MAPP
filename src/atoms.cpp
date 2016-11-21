#include "atoms.h"
using namespace MAPP_NS;
/*---------------------------------------------------------------------------
      ___   _____   _____       ___  ___   _____
     /   | |_   _| /  _  \     /   |/   | /  ___/
    / /| |   | |   | | | |    / /|   /| | | |___
   / / | |   | |   | | | |   / / |__/ | | \___  \
  / /  | |   | |   | |_| |  / /       | |  ___| |
 /_/   |_|   |_|   \_____/ /_/        |_| /_____/
 
 ---------------------------------------------------------------------------*/
Atoms::Atoms():
nvecs(0)
{
    vecs=NULL;
    id_arch=NULL;
    
    grid_established=false;
    box_chng=false;
    skin=0.5;
    tot_natms=0;
    natms=0;
    natms_ph=0;
    xchng_id=0;
    
    x=new Vec<type0>(this,__dim__);
    id= new Vec<int>(this,1);

    
    
    for(int idim=0;idim<__dim__;idim++)
    {
        max_cut_s[idim]=0.0;
        for(int jdim=0;jdim<__dim__;jdim++)
            H[idim][jdim]=B[idim][jdim]=0.0;
    }
    
    comm=new Communincation<__dim__>(this);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Atoms::~Atoms()
{
    delete comm;


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
        if(vecs[ivec]->name && strcmp(name,vecs[ivec]->name)==0)
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
        
        delete [] vecs;
        vecs=vecs_;
        nvecs--;
    }
    else
    {
        delete [] vecs;
        vecs=NULL;
        nvecs=0;
    }
}
/*--------------------------------------------
 automatically grid the domain
 --------------------------------------------*/
void Atoms::auto_grid()
{
    type0 det=1.0;
    for(int idim=0;idim<__dim__;idim++)
        det*=H[idim][idim]*H[idim][idim];
    if(det==0.0)
        error->abort("cannot perform automatic "
        "domain griding, the domain has not been established");
        
    comm->auto_grid(H);
    grid_established=true;
    
    if(my_p==0)
    {
        fprintf(output,"automatic grid performed: ");
        
        for(int i=0;i<__dim__-1;i++)
            fprintf(output,"%dx",comm->tot_p_grid[i]);
        
        fprintf(output,"%d",comm->tot_p_grid[__dim__-1]);
        
        fprintf(output,"\n");
    }
}
/*--------------------------------------------
 manually grid the domain
 --------------------------------------------*/
void Atoms::man_grid(int* n)
{
    comm->man_grid(n);
    xchng=new Exchange<__dim__>(this,vecs,nvecs);
    x2s_lcl();
    xchng->full_xchng();
    s2x_lcl();
    delete xchng;
    grid_established=true;
    
    if(my_p==0)
    {
        fprintf(output,"manual grid performed: ");
        
        for(int i=0;i<__dim__-1;i++)
            fprintf(output,"%dx",comm->tot_p_grid[i]);
        
        fprintf(output,"%d",comm->tot_p_grid[__dim__-1]);
        
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
#ifdef UNROLL
    for(int i=0;i<no;i++)
        XMatrixVector::x2s<__dim__>(x_vec+i*x_dim,B);
#else
    for(int i=0;i<no;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*B[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
                x_vec[idim]+=x_vec[jdim]*B[jdim][idim];
            while(x_vec[idim]<0.0)
                x_vec[idim]++;
            while(x_vec[idim]>=1.0)
                x_vec[idim]--;
        }
#endif
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x(int no)
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
#ifdef UNROLL
    for(int i=0;i<no;i++)
        XMatrixVector::s2x<__dim__>(x_vec+i*x_dim,H);
#else
    for(int i=0;i<no;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*H[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
                x_vec[idim]+=x_vec[jdim]*H[jdim][idim];
        }
#endif
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s_lcl()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
#ifdef UNROLL
    for(int i=0;i<natms;i++)
        XMatrixVector::x2s<__dim__>(x_vec+i*x_dim,B);
#else
    for(int i=0;i<natms;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*B[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
                x_vec[idim]+=x_vec[jdim]*B[jdim][idim];
            
            while(x_vec[idim]<0.0)
                x_vec[idim]++;
            while(x_vec[idim]>=1.0)
                x_vec[idim]--;
        }
#endif
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x_lcl()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
#ifdef UNROLL
    for(int i=0;i<natms;i++)
        XMatrixVector::s2x<__dim__>(x_vec+i*x_dim,H);
#else
    for(int i=0;i<natms;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*H[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
                x_vec[idim]+=x_vec[jdim]*H[jdim][idim];
        }
#endif
}
/*--------------------------------------------
 x2s
 --------------------------------------------*/
void Atoms::x2s_all()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    int nall=natms+natms_ph;
#ifdef UNROLL
    for(int i=0;i<nall;i++)
        XMatrixVector::x2s<__dim__>(x_vec+i*x_dim,B);
#else
    for(int i=0;i<nall;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*B[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
                x_vec[idim]+=x_vec[jdim]*B[jdim][idim];
            
            while(x_vec[idim]<0.0)
                x_vec[idim]++;
            while(x_vec[idim]>=1.0)
                x_vec[idim]--;
        }
#endif
}
/*--------------------------------------------
 s2x
 --------------------------------------------*/
void Atoms::s2x_all()
{
    type0* x_vec=x->begin();
    int x_dim=x->dim;
    int nall=natms+natms_ph;
#ifdef UNROLL
    for(int i=0;i<nall;i++)
        XMatrixVector::s2x<__dim__>(x_vec+i*x_dim,H);
#else
    for(int i=0;i<nall;i++,x_vec+=x_dim)
        for(int idim=0;idim<_dimension_;idim++)
        {
            x_vec[idim]=x_vec[idim]*H[idim][idim];
            for(int jdim=idim+1;jdim<_dimension_;jdim++)
                x_vec[idim]+=x_vec[jdim]*H[jdim][idim];
        }
#endif
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
 
 --------------------------------------------*/
void Atoms::reorder_vecs(VecLst* vec_list_)
{
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
}
/*--------------------------------------------
 init a simulation
 --------------------------------------------*/
void Atoms::init(VecLst* vec_list_,bool box_chng_)
{
    
    
    type0 det=1.0;
    for(int idim=0;idim<__dim__;idim++)
        det*=H[idim][idim]*H[idim][idim];
    if(det==0.0)
        error->abort("cannot start the simulation "
        "domain griding, the domain is not defined");
    if(!grid_established)
        error->abort("cannot start the simulation "
        "domain grid has not been established");
    

    if(!forcefield)
        error->abort("cannot start the simulation, force field has not been established");
    
    if(!neighbor)
        error->abort("cannot start the simulation, neighbior list has not been established");
    
    
    timer->init();
    timer->start(COMM_TIME_mode);
    

    
    reorder_vecs(vec_list_);
    box_chng=box_chng_;
    vec_list=vec_list_;
    
    store_arch_vecs();
    
    x0=new Vec<type0>(this,__dim__);
    


    xchng=new Exchange<__dim__>(this,vec_list->xchng_vecs,vec_list->nxchng_vecs);

    x2s_lcl();
    xchng->full_xchng();
    
    forcefield->init();
    max_cut=forcefield->max_cut();

    updt=new Update<__dim__>(this,vec_list->updt_vecs,vec_list->nupdt_vecs);
    updt->list();
    
    neighbor->init();
    updt->rm_rdndncy();
    store_x0();
    
    timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 update one vectors
 --------------------------------------------*/
void Atoms::fin()
{
    delete xchng;
    delete updt;
    delete x0;

    restore_arch_vecs();

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
 
 --------------------------------------------*/
void Atoms::store_x0()
{
    type0* x_vec=x->begin();
    type0* x0_vec=x0->begin();
    int x_dim=x->dim;
    int last_atm=natms;
    if(box_chng) last_atm+=natms_ph;
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
inline bool Atoms::decide()
{
    type0 skin_sq=0.25*skin*skin;
    type0 dx_sq;
    int succ,succ_lcl=1;
    type0* x_vec=x->begin();
    type0* x0_vec=x0->begin();
    int x_dim=x->dim;
    int last_atm=natms;
    if(box_chng) last_atm+=natms_ph;
    
    for(int iatm=0;succ_lcl && iatm<last_atm;iatm++,x0_vec+=__dim__,x_vec+=x_dim)
    {
        dx_sq=0.0;
        for(int idim=0;idim<__dim__;idim++)
            dx_sq+=(x0_vec[idim]-x_vec[idim])*(x0_vec[idim]-x_vec[idim]);
        if(dx_sq>skin_sq)
            succ_lcl=0;
    }
    
    MPI_Allreduce(&succ_lcl,&succ,1,MPI_INT,MPI_MIN,world);
    if(succ) return true;
    return false;
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
            updt->update(updt_vecs[0],false);
        else
            updt->update(updt_vecs,nupdt_vecs,false);
        timer->stop(COMM_TIME_mode);
        return;
    }
    
    
    if(box_chng)
    {
        if(nupdt_vecs==1)
            updt->update(x,true);
        else
            updt->update(updt_vecs,nupdt_vecs,true);

        if(decide())
        {
            timer->stop(COMM_TIME_mode);
            return;
        }
        
        
        x2s_lcl();
        xchng->full_xchng();
        
        updt->reset();
        updt->list();
        
        neighbor->create_list(box_chng);

        store_x0();
    }
    else
    {
        if(decide())
        {
            if(nupdt_vecs==1)
                updt->update(x,true);
            else
                updt->update(updt_vecs,nupdt_vecs,true);
            
            timer->stop(COMM_TIME_mode);
            return;
        }
        
        x2s_lcl();
        xchng->full_xchng();
        
        updt->list();
        
        neighbor->create_list(box_chng);
        
        
        store_x0();
    }

    timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::init_xchng()
{
    timer->start(COMM_TIME_mode);

    x2s_lcl();
    xchng->full_xchng();
    updt->list();

    timer->stop(COMM_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::fin_xchng()
{
    timer->start(COMM_TIME_mode);
    updt->list();
    neighbor->create_list(box_chng);
    store_x0();
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
    xchng=new Exchange<__dim__>(this,vecs,nvecs);
    x2s_lcl();
    xchng->full_xchng();
    s2x_lcl();
    delete xchng;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::store_arch_vecs()
{
    if(vec_list->narch_vecs==0) return;
    id_arch=new Vec<int>(this,1);
    int* id_0=id->begin();
    int* id_1=id_arch->begin();
    memcpy(id_1,id_0,natms*sizeof(int));
    del_vec(id_arch);
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        del_vec(vec_list->arch_vecs[ivec]);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::restore_arch_vecs()
{
    if(vec_list->narch_vecs==0) return;
    
    int byte_sz=0;
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        byte_sz+=vec_list->arch_vecs[ivec]->byte_sz;
    
    int natms_old=id_arch->vec_sz;
    int* id_old=id_arch->begin();
    int* key_old=NULL;
    if(natms_old) key_old=new int[natms_old];
    for(int i=0;i<natms_old;i++) key_old[i]=i;
    XMath::quicksort(key_old,key_old+natms_old
    ,[&id_old](int* rank_i,int* rank_j){return (id_old[*rank_i]<id_old[*rank_j]);}
    ,[](int* rank_i,int* rank_j){std::swap(*rank_i,*rank_j);}
    );
    
    int natms_new=id->vec_sz;
    int* id_new=id->begin();
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
        error->abort("sanity check failed number of atoms to be sent(%d) "
        "is unequal to number of atoms to be recieved(%d)",tot_nsnd,tot_nrcv);
#endif
    
    auto sort_idx_by_p=
    [this] (const int* lst,int lst_sz,const int* sub_lst,int sub_lst_sz,
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
        error->abort_sing("sanity check failed just found %d out of %d atoms that should be sent");
    if(nrcv!=comm_sz_rcv)
        error->abort_sing("sanity check failed just found %d out of %d atoms that should be recieved");
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
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        for(int i=0;i<nsnd;i++)
            vec_list->arch_vecs[ivec]->cpy(_snd_buff,idx_lst_snd[i]);
    delete [] idx_lst_snd;
    
    
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        vec_list->arch_vecs[ivec]->rearrange(idx_keep_old,idx_keep_new,nkeep,natms_new);
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
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        for(int i=0;i<nrcv;i++)
            vec_list->arch_vecs[ivec]->pst_to(_rcv_buff,idx_lst_rcv[i]);
    
    delete [] idx_lst_rcv;
    delete [] *rcv_buff;
    delete [] rcv_buff;
    
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        add_vec(vec_list->arch_vecs[ivec]);
    add_vec(id_arch);
    delete id_arch;
    id_arch=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::restore_arch_vecs_()
{
    if(vec_list->narch_vecs==0) return;
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
    [this] (const int* lst,int lst_sz,const int* sub_lst,int sub_lst_sz,
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
    int* id_lst_curr=id->begin();
    int natms_arch=id_arch->vec_sz;
    int natms_curr=id->vec_sz;
    int byte_sz=0;
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        byte_sz+=vec_list->arch_vecs[ivec]->byte_sz;
    
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
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        for(int i=0;i<nsnd;i++)
            vec_list->arch_vecs[ivec]->cpy(_snd_buff,idx_lst_snd[i]);
    
    delete [] idx_lst_snd;
    
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        vec_list->arch_vecs[ivec]->rearrange(idx_lst_arch,idx_lst_curr,nkep,natms_curr);
    
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
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        for(int i=0;i<nrcv;i++)
            vec_list->arch_vecs[ivec]->cpy(_rcv_buff,idx_lst_rcv[i]);
    
    delete [] idx_lst_rcv;
    delete [] *rcv_buff;
    delete [] rcv_buff;
    
    
    for(int ivec=0;ivec<vec_list->narch_vecs;ivec++)
        add_vec(vec_list->arch_vecs[ivec]);
    add_vec(id_arch);
    delete id_arch;
    id_arch=NULL;
}
/*----------------------------------------------
  _     _   _____   _____   _       _____   _____
 | |   / / | ____| /  ___| | |     /  ___/ |_   _|
 | |  / /  | |__   | |     | |     | |___    | |
 | | / /   |  __|  | |     | |     \___  \   | |
 | |/ /    | |___  | |___  | |___   ___| |   | |
 |___/     |_____| \_____| |_____| /_____/   |_|
 
 ----------------------------------------------*/
VecLst::VecLst(Atoms* atoms)
{
    tot_xchng_byte_sz=tot_updt_byte_sz=tot_arch_byte_sz=0;
    nxchng_vecs=nupdt_vecs=narch_vecs=0;
    xchng_vecs=updt_vecs=arch_vecs=NULL;
    add_xchng(atoms->x);
    add_xchng(atoms->id);
    add_updt(atoms->x);
}
/*--------------------------------------------
 
 --------------------------------------------*/
VecLst::~VecLst()
{
    delete [] xchng_vecs;
    delete [] updt_vecs;
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
    
    vec** updt_vecs_=NULL;
    if(nupdt_vecs-1)
        updt_vecs_=new vec*[nupdt_vecs-1];
    
    for(int jvec=0;jvec<ivec;jvec++)
        updt_vecs_[jvec]=updt_vecs[jvec];
    for(int jvec=ivec+1;jvec<nupdt_vecs;jvec++)
        updt_vecs_[jvec-1]=updt_vecs[jvec];
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
    
    vec** xchng_vecs_=NULL;;
    if(nxchng_vecs-1)
        xchng_vecs_=new vec*[nxchng_vecs-1];
    
    for(int jvec=0;jvec<ivec;jvec++)
        xchng_vecs_[jvec]=xchng_vecs[jvec];
    for(int jvec=ivec+1;jvec<nxchng_vecs;jvec++)
        xchng_vecs_[jvec-1]=xchng_vecs[jvec];
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
    
    vec** arch_vecs_=NULL;;
    if(narch_vecs-1)
        arch_vecs_=new vec*[narch_vecs-1];
    
    for(int jvec=0;jvec<ivec;jvec++)
        arch_vecs_[jvec]=arch_vecs[jvec];
    for(int jvec=ivec+1;jvec<narch_vecs;jvec++)
        arch_vecs_[jvec-1]=arch_vecs[jvec];
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

