/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "sgcmc.h"
#include "memory.h"
#include "random.h"
#include "atom_types.h"
#include "neighbor.h"
#include "ff.h"
#include "md.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
SGCMC::SGCMC(MAPP* mapp,int m_,dmd_type gas_type_,type0 mu_,type0 T_,int seed):
GCMC(mapp,gas_type_,mu_,T_,seed),
m(m_)
{
    s_x_buff=NULL;

    cell_size=new type0[dim];
    ncells_per_dim=new int[dim];
    cell_denom=new int[dim];
    icell_coord=new int[dim];
    jcell_coord=new int[dim];
    
    head_atm=NULL;
    cell_coord_buff=NULL;
    ncells=0;
    
    /*--------------------------------------------------
     find the relative neighbor list for cells here we 
     figure out the number of neighboring cells and also 
     allocate the memory for rel_neigh_lst and 
     rel_neigh_lst_coord. values for rel_neigh_lst_coord
     is assigned here,but finding the values for 
     rel_neigh_lst requires knowledge of the box and 
     domain, it will be differed to create()
     --------------------------------------------------*/
    int countr[dim];
    for(int i=0;i<dim;i++)
        countr[i]=-m;
    int max_no_neighs=1;
    for(int i=0;i<dim;i++)
        max_no_neighs*=2*m+1;
    nneighs=0;
    rel_neigh_lst_coord=new int[max_no_neighs*dim];
    int* rel_neigh_lst_coord_=rel_neigh_lst_coord;
    int sum;
    int rc_sq=m*m;
    for(int i=0;i<max_no_neighs;i++)
    {
        sum=0;
        for(int j=0;j<dim;j++)
        {
            sum+=countr[j]*countr[j];
            if(countr[j]!=0)
                sum+=1-2*fabs(countr[j]);
        }
        
        if(sum<rc_sq)
        {
            for(int j=0;j<dim;j++)
                rel_neigh_lst_coord_[j]=countr[j];
            rel_neigh_lst_coord_+=dim;
            nneighs++;
        }
        
        countr[0]++;
        for(int j=0;j<dim-1;j++)
            if(countr[j]==m+1)
            {
                countr[j]=-m;
                countr[j+1]++;
            }
    }
    rel_neigh_lst_coord_=new int[nneighs*dim];
    memcpy(rel_neigh_lst_coord_,rel_neigh_lst_coord,nneighs*dim*sizeof(int));
    delete [] rel_neigh_lst_coord;
    rel_neigh_lst_coord=rel_neigh_lst_coord_;
    
    

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
SGCMC::~SGCMC()
{
    delete [] rel_neigh_lst_coord;
    delete [] jcell_coord;
    delete [] icell_coord;
    delete [] cell_denom;
    delete [] ncells_per_dim;
    delete [] cell_size;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::init()
{
    GCMC::init();
    MPI_Scan(&ngas,&ngas_before,1,MPI_INT,MPI_SUM,world);
    ngas_before-=ngas;
    box_setup();
    vars=new type0[forcefield->gcmc_n_vars];
    lcl_vars=new type0[forcefield->gcmc_n_vars];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::fin()
{
    delete [] vars;
    delete [] lcl_vars;
    box_dismantle();
    GCMC::fin();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::box_setup()
{
    GCMC::box_setup();
    
    ncells=1;
    for(int i=0;i<dim;i++)
    {
        cell_size[i]=cut_s[i]/static_cast<type0>(m);
        ncells_per_dim[i]=static_cast<int>((s_hi[i]-s_lo[i])/cell_size[i])+1;
        
        cell_denom[i]=ncells;
        ncells*=ncells_per_dim[i];
    }
    
    head_atm=new int[ncells];
    cell_coord_buff=new int[max_ntrial_atms*dim];
    s_x_buff=new type0[dim*max_ntrial_atms];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::box_dismantle()
{
    delete [] head_atm;
    delete [] cell_coord_buff;
    delete [] s_x_buff;
    
    head_atm=NULL;
    cell_coord_buff=NULL;
    s_x_buff=NULL;
    ncells=0;
}
/*--------------------------------------------
 construct the bin list
 --------------------------------------------*/
void SGCMC::xchng(bool box_chng,int nattmpts)
{
    curr_comm=&world;
    curr_root=0;
        
    
    atoms->init_xchng();
    if(box_chng)
        box_setup();
    
    /*--------------------------------------------------
     here we allocate the memory for cell_vec & next_vec
     --------------------------------------------------*/
    if(forcefield->gcmc_tag_enabled) tag_vec_p=new Vec<int>(atoms,1);
    else tag_vec_p=NULL;
    cell_vec_p=new Vec<int>(atoms,1);
    next_vec_p=new Vec<int>(atoms,1);
    s_vec_p=new Vec<type0>(atoms,dim);
    memcpy(s_vec_p->begin(),mapp->x->begin(),sizeof(type0)*natms*dim);
    
    /*--------------------------------------------------
     here we reset head_atm
     --------------------------------------------------*/
    for(int i=0;i<ncells;i++) head_atm[i]=-1;
    
    /*--------------------------------------------------
     here we assign values for cell_vec, next_vec &
     headt_atm; it is supposed that we have fractional
     coordinates at this point
     --------------------------------------------------*/
    int* next_vec=next_vec_p->begin();
    int* cell_vec=cell_vec_p->begin();
    type0* s=mapp->x->begin()+(natms-1)*dim;
    for(int i=natms-1;i>-1;i--,s-=dim)
    {
        find_cell_no(s,cell_vec[i]);
        next_vec[i]=head_atm[cell_vec[i]];
        head_atm[cell_vec[i]]=i;
    }
    
    neighbor->create_list(box_chng);

    forcefield->init_xchng();
    for(int i=0;i<atoms->nvecs;i++)
        atoms->vecs[i]->resize(natms);
        
    ngas=0;
    md_type* type=mapp->type->begin();
    for(int i=0;i<natms;i++)
        if(type[i]==gas_type) ngas++;
    MPI_Scan(&ngas,&ngas_before,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&ngas,&tot_ngas,1,MPI_INT,MPI_SUM,world);
    ngas_before-=ngas;
    

    next_jatm_p=&SGCMC::next_jatm_reg;
    dof_diff=0;
    for(int i=0;i<nattmpts;i++)
        attmpt();

    
    forcefield->fin_xchng();
    memcpy(mapp->x->begin(),s_vec_p->begin(),sizeof(type0)*natms*dim);
    
    delete tag_vec_p;
    delete s_vec_p;
    delete next_vec_p;
    delete cell_vec_p;
    atoms->fin_xchng();


}
/*--------------------------------------------
 this must be used only for local atoms
 --------------------------------------------*/
inline void SGCMC::find_cell_no(type0*& s,int& cell_no)
{
    cell_no=0;
    for(int i=0;i<dim;i++)
        cell_no+=cell_denom[i]*MIN(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),ncells_per_dim[i]-1);
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void SGCMC::find_cell_coord(type0*& s,int*& cell_coord)
{
    for(int i=0;i<dim;i++)
    {
        if(s[i]<s_lo[i])
            cell_coord[i]=-static_cast<int>((s_lo[i]-s[i])/cell_size[i])-1;
        else if(s_hi[i]<=s[i])
            cell_coord[i]=MAX(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),ncells_per_dim[i]-1);
        else
            cell_coord[i]=MIN(static_cast<int>((s[i]-s_lo[i])/cell_size[i]),ncells_per_dim[i]-1);
    }
}
/*--------------------------------------------
 attempt an insertion
 --------------------------------------------*/
void SGCMC::attmpt()
{
    if(random->uniform()<0.5)
    {
        xchng_mode=INS_MODE;
        im_root=true;
        int iproc=-1;
        for(int i=0;i<dim;i++)
        {
            s_buff[i]=random->uniform();
            if(s_buff[i]<s_lo[i] || s_buff[i]>=s_hi[i])
                im_root=false;
        }
        if(im_root)
            iproc=atoms->my_p;
        MPI_Allreduce(&iproc,&curr_root,1,MPI_INT,MPI_MAX,world);
    }
    else
    {
        if(tot_ngas==0)
            return;
        xchng_mode=DEL_MODE;
        igas=static_cast<int>(tot_ngas*random->uniform());
        int iproc=-1;
        im_root=false;
        if(ngas_before<=igas && igas<ngas_before+ngas)
        {
            im_root=true;
            iproc=atoms->my_p;
            int n=igas-ngas_before;
            int icount=-1;
            md_type* type=mapp->type->begin();
            
            del_idx=0;
            for(;icount!=n;del_idx++)
                if(type[del_idx]==gas_type) icount++;
            del_idx--;
            gas_id=mapp->id->begin()[del_idx];
            memcpy(s_buff,s_vec_p->begin()+del_idx*dim,dim*sizeof(type0));
        }
        MPI_Allreduce(&iproc,&curr_root,1,MPI_INT,MPI_MAX,world);
        MPI_Bcast(s_buff,dim,MPI_TYPE0,curr_root,world);
        MPI_Bcast(&gas_id,1,MPI_INT,curr_root,world);
    }
    
    prep_s_x_buff();
    if(tag_vec_p) reset_tag();
    forcefield->pre_gcmc_energy(this);
    
    delta_u=forcefield->gcmc_energy(this);
    MPI_Bcast(&delta_u,1,MPI_TYPE0,curr_root,world);
    
    type0 fac;
    root_succ=false;
    if(xchng_mode==INS_MODE)
    {
        fac=z_fac*vol/((static_cast<type0>(tot_ngas)+1.0)*exp(beta*delta_u));
        if(random->uniform()<fac)
        {
            root_succ=true;
            ins_succ();
            forcefield->post_gcmc_energy(this);
        }
    }
    else
    {
        fac=static_cast<type0>(tot_ngas)*exp(beta*delta_u)/(z_fac*vol);
        if(random->uniform()<fac)
        {
            root_succ=true;
            del_succ();
            forcefield->post_gcmc_energy(this);
        }
    }
    
}
/*--------------------------------------------
 attempt an insertion
 --------------------------------------------*/
void SGCMC::ins_attmpt()
{
    xchng_mode=INS_MODE;
    for(int i=0;i<dim;i++)
        s_buff[i]=random->uniform();
    
    prep_s_x_buff();
    if(tag_vec_p) reset_tag();
    forcefield->xchng_energy_timer(this);
}
/*--------------------------------------------
 things to do after a successful insertion
 --------------------------------------------*/
void SGCMC::ins_succ()
{

    int new_id=get_new_id();
    for(int i=0;i<dim;i++) vel_buff[i]=random->gaussian()*sigma;
    if(im_root)
    {
        atoms->add();
        
        memcpy(mapp->x->begin()+(natms-1)*dim,s_x_buff,dim*sizeof(type0));
        memcpy(mapp->x_d->begin()+(natms-1)*dim,vel_buff,dim*sizeof(type0));
        mapp->type->begin()[natms-1]=gas_type;
        mapp->id->begin()[natms-1]=new_id;
        if(tag_vec_p) tag_vec_p->begin()[natms-1]=-1;
        if(mapp->x_dof)
        {
            bool* dof=mapp->x_dof->begin()+(natms-1)*dim;
            for(int i=0;i<dim;i++) dof[i]=true;
        }
        
        
        memcpy(s_vec_p->begin()+(natms-1)*dim,s_buff,dim*sizeof(type0));
        
        int cell_=0;
        for(int i=0;i<dim;i++) cell_+=cell_denom[i]*cell_coord_buff[i];
        cell_vec_p->begin()[natms-1]=cell_;
        
        
        int* nxt_p=head_atm+cell_;
        while(*nxt_p!=-1)
            nxt_p=next_vec_p->begin()+*nxt_p;
        *nxt_p=natms-1;
        next_vec_p->begin()[natms-1]=-1;
        ngas++;
    }
    else
    {
        if(atoms->my_p>curr_root)
            ngas_before++;
    }
    
    dof_diff+=dim;
    tot_ngas++;
    atoms->tot_natms++;    
}
/*--------------------------------------------
 attempt a deletion
 --------------------------------------------*/
void SGCMC::del_attmpt()
{
    if(tot_ngas==0) return;
    xchng_mode=DEL_MODE;
    
    
    igas=static_cast<int>(tot_ngas*random->uniform());
    int iproc,iproc_=-1;
    if(ngas_before<=igas && igas<ngas_before+ngas)
    {
        iproc_=atoms->my_p;
        int n=igas-ngas_before;
        int icount=-1;
        md_type* type=mapp->type->begin();
        
        del_idx=0;
        for(;icount!=n;del_idx++)
            if(type[del_idx]==gas_type) icount++;
        del_idx--;
        gas_id=mapp->id->begin()[del_idx];
        memcpy(s_buff,s_vec_p->begin()+del_idx*dim,dim*sizeof(type0));
    }
    MPI_Allreduce(&iproc_,&iproc,1,MPI_INT,MPI_MAX,world);
    MPI_Bcast(s_buff,dim,MPI_TYPE0,iproc,world);
    MPI_Bcast(&gas_id,1,MPI_INT,iproc,world);

    prep_s_x_buff();
    if(tag_vec_p) reset_tag();
    forcefield->xchng_energy_timer(this);
}
/*--------------------------------------------
 things to do after a successful deletion
 --------------------------------------------*/
void SGCMC::del_succ()
{

    add_del_id();

    if(im_root)
    {
        /*--------------------------------------------------
         0.0 find the link to del_idx
         --------------------------------------------------*/
        int* p_2_idx=head_atm+cell_vec_p->begin()[del_idx];
        while(*p_2_idx!=del_idx)
            p_2_idx=next_vec_p->begin()+*p_2_idx;
        /*--------------------------------------------------
         0.1 replace link to del_idx with link from del_idx
         --------------------------------------------------*/
        *p_2_idx=next_vec_p->begin()[del_idx];
        
        if(del_idx!=natms-1)
        {
            /*--------------------------------------------------
             1.0 find the first link to after del_idx
             --------------------------------------------------*/
            int* p_2_first_aft_del=head_atm+cell_vec_p->begin()[natms-1];
            while(*p_2_first_aft_del<del_idx)
                p_2_first_aft_del=next_vec_p->begin()+*p_2_first_aft_del;
            
            /*--------------------------------------------------
             1.1 find the link to natms-1
             --------------------------------------------------*/
            int* p_2_last=p_2_first_aft_del;
            while(*p_2_last!=natms-1)
                p_2_last=next_vec_p->begin()+*p_2_last;
            
            
            /*--------------------------------------------------
             1.2 remove the link to natms-1 and end it ther
             --------------------------------------------------*/
            *p_2_last=-1;
            
            /*--------------------------------------------------
             1.3 insert the new link at del_idx, but for now it 
             is at natms-1 after the move operation which 
             happens next it will go del_idx
             --------------------------------------------------*/
            next_vec_p->begin()[natms-1]=*p_2_first_aft_del;
            *p_2_first_aft_del=del_idx;
            
        }
        
        
        atoms->del(del_idx);
        ngas--;
    }
    else
    {
        if(igas<ngas_before)
            ngas_before--;
    }
    dof_diff-=dim;
    tot_ngas--;
    atoms->tot_natms--;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::prep_s_x_buff()
{
    
    im_root=true;
    for(int i=0;i<dim;i++)
        if(s_buff[i]<s_lo[i] || s_buff[i]>=s_hi[i])
            im_root=false;
    
    
    int n_per_dim[dim];
    type0 s;
    int no;
    ntrial_atms=1;
    for(int i=0;i<dim;i++)
    {
        no=0;
        s=s_buff[i];
        
        for(type0 s_=s;s_<s_hi_ph[i];s_++)
            if(s_lo_ph[i]<=s_ && s_<s_hi_ph[i])
                s_trials[i][no++]=s_;
        
        for(type0 s_=s-1.0;s_lo_ph[i]<=s_;s_--)
            if(s_lo_ph[i]<=s_ && s_<s_hi_ph[i])
                s_trials[i][no++]=s_;
        
        ntrial_atms*=no;
        n_per_dim[i]=no;
    }
    
    
    int count[dim];
    type0* buff=s_x_buff;
    int* cell_coord=cell_coord_buff;
    type0** H=atoms->H;
    for(int i=0;i<dim;i++) count[i]=0;
    for(int i=0;i<ntrial_atms;i++)
    {
        for(int j=0;j<dim;j++)
            buff[j]=s_trials[j][count[j]];
        
        find_cell_coord(buff,cell_coord);
        
        for(int j=0;j<dim;j++)
        {
            buff[j]*=H[j][j];
            for(int k=j+1;k<dim;k++)
                buff[j]+=buff[k]*H[k][j];
        }
        
        count[0]++;
        for(int j=0;j<dim-1;j++)
            if(count[j]==n_per_dim[j])
            {
                count[j]=0;
                count[j+1]++;
            }
        cell_coord+=dim;
        buff+=dim;
    }
}
/*--------------------------------------------
 this is used for insertion trial
 --------------------------------------------*/
void SGCMC::reset_iatm()
{
    itrial_atm=-1;
    next_iatm();
}
/*--------------------------------------------
 this is used for insertion trial
 --------------------------------------------*/
void SGCMC::next_iatm()
{
    itrial_atm++;
    
    if(itrial_atm==ntrial_atms)
    {
        /*--------------------------------------------------
         if we have reached the end of the list
         --------------------------------------------------*/
        iatm=-1;
        return;
    }
    
    /*--------------------------------------------------
     give the atom a number that cannot be the
     same as the existing ones:
     
     iatm=natms+natms_ph+itrial_atm;
     
     natms+natms_ph <= iatm
     iatm < natms+natms_ph+ntrial_atms
     --------------------------------------------------*/
    if(itrial_atm==0 && im_root && xchng_mode==DEL_MODE)
        iatm=del_idx;
    else
        iatm=natms+itrial_atm;
    
    /*--------------------------------------------------
     assign the cell number and the already calculated
     cell coordinates
     --------------------------------------------------*/
    for(int i=0;i<dim;i++)
        icell_coord[i]=cell_coord_buff[dim*itrial_atm+i];
    
    /*--------------------------------------------------
     assign the position of iatm
     --------------------------------------------------*/
    ix=s_x_buff+itrial_atm*dim;
    
}
/*--------------------------------------------
 find the next cell that probably
 contains local atoms
 --------------------------------------------*/
inline void SGCMC::next_jcell()
{
    ineigh++;
    if(ineigh==nneighs)
    {
        jcell=-1;
        return;
    }
    jcell=0;
    for(int i=0;i<dim;i++)
    {
        jcell_coord[i]=icell_coord[i]+rel_neigh_lst_coord[ineigh*dim+i];
        if(jcell_coord[i]<0 || jcell_coord[i]>ncells_per_dim[i]-1)
            return next_jcell();
        jcell+=cell_denom[i]*jcell_coord[i];
    }
    
    if(head_atm[jcell]==-1)
        return next_jcell();
    
    jatm=head_atm[jcell];
    
}
/*--------------------------------------------
 this is used for insertion trial
 --------------------------------------------*/
void SGCMC::reset_jatm()
{
    ineigh=0;
    jatm_next=-1;
    next_jatm_reg();
    
    /*
    ineigh=-1;
    
    next_jcell();
    if(jcell==-1)
    {
        if(itrial_atm==0 && im_root)
        {
            iself=0;
            next_jatm_p=&SGCMC::next_jatm_self;
            return next_jatm_self();
        }
        else
        {
            jatm=-1;
            return;
        }
        
    }
    
    if(jatm==iatm)
        return next_jatm_reg();
    
    jx=mapp->x->begin()+dim*jatm;
    jtype=mapp->type->begin()[jatm];
    rsq=0.0;
    for(int i=0;i<dim;i++)
        rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
    
    if(rsq>=cut_sq[itype][jtype])
        return next_jatm_reg();
    if(tag_vec_p) tag_vec_p->begin()[jatm]=icomm;
    */
}
/*--------------------------------------------
 find the next lcl atom
 --------------------------------------------*/
inline void SGCMC::next_jatm_reg()
{
 
    while(true)
    {
        while(jatm_next==-1 && ineigh<nneighs)
        {
            bool lcl=true;
            jcell=0;
            for(int i=0;i<dim && lcl;i++)
            {
                jcell_coord[i]=icell_coord[i]+rel_neigh_lst_coord[ineigh*dim+i];
                jcell+=cell_denom[i]*jcell_coord[i];

                if(jcell_coord[i]<0 || jcell_coord[i]>ncells_per_dim[i]-1)
                    lcl=false;
            }
            
            if(lcl)
                jatm_next=head_atm[jcell];
            ineigh++;
        }
        
        jatm=jatm_next;
        if(jatm==-1)
        {
            
            if(itrial_atm==0 && im_root)
            {
                iself=0;
                next_jatm_p=&SGCMC::next_jatm_self;
                return next_jatm_self();
            }
            else return;
        }
        
        jatm_next=next_vec_p->begin()[jatm];
        
        if(jatm==iatm) continue;
        
        jx=mapp->x->begin()+dim*jatm;
        jtype=mapp->type->begin()[jatm];
        
        rsq=0.0;
        for(int i=0;i<dim;i++) rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
        
        if(rsq>=cut_sq[itype][jtype]) continue;
        
        if(tag_vec_p) tag_vec_p->begin()[jatm]=icomm;
        return;
    }
  
 
    /*
    if(next_vec_p->begin()[jatm]==-1)
    {
        next_jcell();
        if(jcell==-1)
        {
            if(itrial_atm==0 && im_root)
            {
                iself=0;
                next_jatm_p=&SGCMC::next_jatm_self;
                return next_jatm_self();
            }
            else
            {
                jatm=-1;
                return;
            }
            
        }
        else
        {
            if(jatm==iatm)
                return next_jatm_reg();

            jx=mapp->x->begin()+dim*jatm;
            jtype=mapp->type->begin()[jatm];
            rsq=0.0;
            for(int i=0;i<dim;i++)
                rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
            if(rsq>=cut_sq[itype][jtype])
                return next_jatm_reg();
            if(tag_vec_p) tag_vec_p->begin()[jatm]=icomm;
            return;
        }
    }
    
    jatm=next_vec_p->begin()[jatm];
    if(jatm==iatm)
        return next_jatm_reg();
    
    jx=mapp->x->begin()+dim*jatm;
    jtype=mapp->type->begin()[jatm];
    rsq=0.0;
    for(int i=0;i<dim;i++)
        rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
    if(rsq>=cut_sq[itype][jtype])
        return next_jatm_reg();
    if(tag_vec_p) tag_vec_p->begin()[jatm]=icomm;
    */
}
/*--------------------------------------------
 find the next interactin image atom
 --------------------------------------------*/
inline void SGCMC::next_jatm_self()
{
    
    while(true)
    {
        iself++;
        if(iself==ntrial_atms)
        {
            next_jatm_p=&SGCMC::next_jatm_reg;
            jatm=-1;
            return;
        }
        
        jatm=natms+iself;
        jx=s_x_buff+iself*dim;
        jtype=gas_type;
        rsq=0.0;
        for(int i=0;i<dim;i++)
            rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
        
        if(rsq>=cut_sq[itype][jtype]) continue;
        
        return;
    }
    
    /*
    iself++;
    if(iself==ntrial_atms)
    {
        next_jatm_p=&SGCMC::next_jatm_reg;
        jatm=-1;
        return;
    }
    jatm=natms+iself;
    jx=s_x_buff+iself*dim;
    jtype=gas_type;
    rsq=0.0;
    for(int i=0;i<dim;i++)
        rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
    if(rsq>=cut_sq[itype][jtype])
        next_jatm_self();
     */
}
/*--------------------------------------------
 find the next interactin image atom
 --------------------------------------------*/
void SGCMC::next_jatm()
{
    (this->*next_jatm_p)();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::next_icomm()
{
    icomm=-1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::reset_icomm()
{
    icomm=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void SGCMC::refresh()
{
    for(int i=0;i<ncells;i++)
        head_atm[i]=-1;
    int* next_vec=next_vec_p->begin()+natms-1;
    int* cell_vec=cell_vec_p->begin()+natms-1;
    for(int i=natms-1;i>-1;i--,next_vec--,cell_vec--)
    {
        *next_vec=head_atm[*cell_vec];
        head_atm[*cell_vec]=i;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void SGCMC::reset_tag()
{
    int* tag=tag_vec_p->begin();
    for(int i=0;i<natms;i++)
        tag[i]=-1;
}






