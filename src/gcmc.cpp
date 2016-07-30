/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "gcmc.h"
#include "memory.h"
#include "random.h"
#include "atom_types.h"
#include "neighbor.h"
#include "ff.h"
#include "md.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
GCMC::GCMC(MAPP* mapp,int m_,dmd_type gas_type_,type0 mu_,type0 T_,int seed):InitPtrs(mapp),
dim(atoms->dimension),
x_dim(mapp->x->dim),
gas_type(gas_type_),
mu(mu_),
T(T_),
cut_sq(forcefield->cut_sq),
m(m_),
s_hi(atoms->s_hi),
s_lo(atoms->s_lo),
natms(atoms->natms),
natms_ph(atoms->natms_ph),
itype(gas_type)
{
    random=new Random(mapp,seed);
    vel_buff=new type0[dim];
    cut_s=new type0[dim];
    s_lo_ph=new type0[dim];
    s_hi_ph=new type0[dim];
    cell_size=new type0[dim];
    ncells_per_dim=new int[dim];
    cell_denom=new int[dim];
    icell_coord=new int[dim];
    jcell_coord=new int[dim];
    nimages_per_dim=new int*[dim];
    *nimages_per_dim=new int[dim*2];
    for(int i=1;i<dim;i++) nimages_per_dim[i]=nimages_per_dim[i-1]+2;
    ins_s_trials=new type0*[dim];
    
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
    rel_neigh_lst=new int[nneighs];
    
    
    head_atm=NULL;
    *ins_s_trials=NULL;
    ins_cell=NULL;
    ins_cell_coord=NULL;
    ins_buff=NULL;
    ins_cell=NULL;
    del_lst=NULL;
    ncells=0;
    
    del_ids=NULL;
    del_ids_sz=del_ids_cpcty=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
GCMC::~GCMC()
{
    delete [] rel_neigh_lst;
    delete [] rel_neigh_lst_coord;
    delete [] ins_s_trials;
    delete [] *nimages_per_dim;
    delete [] nimages_per_dim;
    delete [] jcell_coord;
    delete [] icell_coord;
    delete [] cell_denom;
    delete [] ncells_per_dim;
    delete [] cell_size;
    delete [] s_hi_ph;
    delete [] s_lo_ph;
    delete [] cut_s;
    delete [] vel_buff;
    delete random;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void GCMC::add_del_id()
{
    if(del_ids_sz+1>del_ids_cpcty)
    {
        int* del_ids_=new int[del_ids_cpcty+1];
        memcpy(del_ids_,del_ids,del_ids_sz*sizeof(int));
        del_ids_cpcty++;
        delete [] del_ids;
        del_ids=del_ids_;
    }
    
    del_ids[del_ids_sz++]=gas_id;
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline int GCMC::get_new_id()
{
    if(del_ids_sz)
    {
        del_ids_sz--;
        return del_ids[del_ids_sz];
    }
    else
    {
        max_id++;
        return max_id;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::init()
{
    cut=forcefield->cut[itype][0];
    for(int i=1;i<atom_types->no_types;i++)
        cut=MAX(cut,forcefield->cut[itype][i]);
    
    gas_mass=atom_types->mass[gas_type];
    kbT=mapp->md->boltz*T;
    beta=1.0/kbT;
    lambda=mapp->md->hplanck/sqrt(2.0*M_PI*kbT*gas_mass);
    sigma=sqrt(kbT/gas_mass);
    z_fac=1.0;
    for(int i=0;i<dim;i++) z_fac/=lambda;
    zz_fac=-kbT*log(zz_fac);
    z_fac*=exp(beta*mu);
    
    box_setup();
    
    ngas=0;
    md_type* type=mapp->type->begin();
    for(int i=0;i<natms;i++)
        if(type[i]==gas_type) ngas++;
    MPI_Scan(&ngas,&ngas_before,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&ngas,&tot_ngas,1,MPI_INT,MPI_SUM,world);
    ngas_before-=ngas;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::fin()
{
    delete [] head_atm;
    delete [] *ins_s_trials;
    delete [] ins_cell;
    delete [] ins_cell_coord;
    delete [] ins_buff;
    delete [] del_lst;
    
    head_atm=NULL;
    *ins_s_trials=NULL;
    ins_cell=NULL;
    ins_cell_coord=NULL;
    ins_buff=NULL;
    ins_cell=NULL;
    del_lst=NULL;
    ncells=0;
    
    delete [] del_ids;
    del_ids=NULL;
    del_ids_sz=del_ids_cpcty=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::box_setup()
{
    int sz=0;
    ncells=1;
    vol=1.0;
    int max_ntrial_atms=1;
    for(int i=0;i<dim;i++)
    {
        vol*=atoms->H[i][i];
        type0 tmp=0.0;
        for(int j=i;j<dim;j++)
            tmp+=atoms->B[j][i]*atoms->B[j][i];
        
        cut_s[i]=sqrt(tmp)*cut;
        
        cell_size[i]=cut_s[i]/static_cast<type0>(m);
        
        s_lo_ph[i]=s_lo[i]-cut_s[i];
        s_hi_ph[i]=s_hi[i]+cut_s[i];
        nimages_per_dim[i][0]=static_cast<int>(floor(s_hi_ph[i]));
        nimages_per_dim[i][1]=-static_cast<int>(floor(s_lo_ph[i]));
        max_ntrial_atms*=1+nimages_per_dim[i][0]+nimages_per_dim[i][1];
        sz+=1+nimages_per_dim[i][0]+nimages_per_dim[i][1];
        
        ncells_per_dim[i]=static_cast<int>
        ((s_hi[i]-s_lo[i])/cell_size[i])+1+2*m;
        
        cell_denom[i]=ncells;
        ncells*=ncells_per_dim[i];
    }
    
    /*--------------------------------------------------
     here we assign values for rel_neigh_lst using the
     values determined for rel_neigh_lst_coord and
     cell_denom
     --------------------------------------------------*/
    int* rel_neigh_lst_coord_=rel_neigh_lst_coord;
    for(int i=0;i<nneighs;i++,rel_neigh_lst_coord_+=dim)
    {
        rel_neigh_lst[i]=0;
        for(int j=0;j<dim;j++)
            rel_neigh_lst[i]+=
            rel_neigh_lst_coord_[j]*cell_denom[j];
    }
    
    
    delete [] head_atm;
    delete [] *ins_s_trials;
    delete [] ins_cell;
    delete [] ins_cell_coord;
    delete [] ins_buff;
    delete [] ins_cell;
    delete [] del_lst;
    
    
    /*--------------------------------------------------
     here we allocate the memory for head_atm
     --------------------------------------------------*/
    head_atm=new int[ncells];
    
    /*--------------------------------------------------
     here we allocate the memory for ins_s_trials
     --------------------------------------------------*/
    *ins_s_trials=new type0[sz];
    for(int i=1;i<dim;i++)
        ins_s_trials[i]=ins_s_trials[i-1]+1+nimages_per_dim[i-1][0]+nimages_per_dim[i-1][1];
    
    /*--------------------------------------------------
     here we allocate the memory for ins_buff & del_buff
     --------------------------------------------------*/
    ins_buff=new type0[dim*max_ntrial_atms];
    ins_cell=new int[max_ntrial_atms];
    ins_cell_coord=new int[max_ntrial_atms*dim];
    del_lst=new int[max_ntrial_atms];
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void GCMC::refresh()
{
    for(int i=0;i<ncells;i++)
        head_atm[i]=-1;
    int nall=natms+natms_ph;
    int* next_vec=next_vec_p->begin()+nall-1;
    int* cell_vec=cell_vec_p->begin()+nall-1;
    for(int i=nall-1;i>-1;i--,next_vec--,cell_vec--)
        if(*cell_vec!=-1)
        {
            *next_vec=head_atm[*cell_vec];
            head_atm[*cell_vec]=i;
        }

}
/*--------------------------------------------
 
 --------------------------------------------*/
inline void GCMC::find_cell_no(type0*& s,int& bin_no)
{
    bin_no=0;
    for(int i=0;i<dim;i++)
    {
        if(s[i]>=s_lo_ph[i] && s[i]<s_hi_ph[i])
            bin_no+=cell_denom[i]*(m+
            static_cast<int>(floor((s[i]-s_lo[i])/cell_size[i])));
        else
        {
            bin_no=-1;
            return;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
inline bool GCMC::lcl(type0*& s)
{
    for(int i=0;i<dim;i++)
        if(s[i]<s_lo[i] || s[i]>=s_hi[i])
            return false;
    return true;
}
/*--------------------------------------------
 this is used for insertion trial
 --------------------------------------------*/
inline void GCMC::next_iatm_ins()
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
    iatm=natms+natms_ph+itrial_atm;
    
    /*--------------------------------------------------
     assign the cell number and the already calculated 
     cell coordinates
     --------------------------------------------------*/
    icell=ins_cell[itrial_atm];
    for(int i=0;i<dim;i++)
        icell_coord[i]=ins_cell_coord[dim*itrial_atm+i];
    
    /*--------------------------------------------------
     assign the position of iatm
     --------------------------------------------------*/
    ix=ins_buff+itrial_atm*dim;
    
    /*--------------------------------------------------
     reset the neighboring cell number
     --------------------------------------------------*/
    ineigh=-1;
    
    if(itrial_atm==0 && first_atm_lcl)
    {
        /*--------------------------------------------------
         if this atom does belong to us we need to consider
         all neighbors including the phantom ones
         --------------------------------------------------*/
        next_jatm_p=&GCMC::next_jatm_all;
        /*--------------------------------------------------
         find the first non-empty cell
         --------------------------------------------------*/
        next_jcell_all();
        if(jcell==-1)
        {
            /*--------------------------------------------------
             if we did not find a non-emty neighboring cell we 
             need to go through interaction of this atom with
             its images
             --------------------------------------------------*/
            iself=0;
            next_jatm_p=&GCMC::next_jatm_self;
            next_jatm_self();
        }
        else
        {
            jx=mapp->x->begin()+x_dim*jatm;
            jtype=mapp->type->begin()[jatm];
            rsq=0.0;
            for(int i=0;i<dim;i++)
                rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
            /*--------------------------------------------------
             if we get an atom that is non-interacting we need
             to go the next local atom
             --------------------------------------------------*/
            if(rsq>=cut_sq[itype][jtype])
                next_jatm_all();
        }
    }
    else
    {
        /*--------------------------------------------------
         if this atom does not belong to us we need to
         consider only local neighbors excluding the phantom
         ones
         --------------------------------------------------*/
        next_jatm_p=&GCMC::next_jatm_lcl;
        /*--------------------------------------------------
         find the first non-empty local cell
         --------------------------------------------------*/
        next_jcell_lcl();
        if(jcell==-1)
        {
            /*--------------------------------------------------
             if we did not find a non-emty neighboring cell 
             we have reached the end of the list
             --------------------------------------------------*/
            jatm=-1;
        }
        else
        {
            if(jatm>=natms)
            {
                /*--------------------------------------------------
                 if we get a phantom atom we need to go to the next
                 local atom
                 --------------------------------------------------*/
                next_jatm_lcl();
                return;
            }
            jx=mapp->x->begin()+x_dim*jatm;
            jtype=mapp->type->begin()[jatm];
            rsq=0.0;
            for(int i=0;i<dim;i++)
                rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
            /*--------------------------------------------------
             if we get an atom that is non-interacting we need 
             to go the next local atom
             --------------------------------------------------*/
            if(rsq>=cut_sq[itype][jtype])
                next_jatm_lcl();
        }
    }
    
}
/*--------------------------------------------
 this is used for deletion trial
 --------------------------------------------*/
inline void GCMC::next_iatm_del()
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
     get the atom number from deletion list
     --------------------------------------------------*/
    iatm=del_lst[itrial_atm];
    
    /*--------------------------------------------------
     assign the cell number and calculate the cell 
     coordinates
     --------------------------------------------------*/
    icell=cell_vec_p->begin()[iatm];
    int tmp=icell;
    for(int i=dim-1;i>-1;i--)
    {
        icell_coord[i]=tmp/cell_denom[i];
        tmp-=icell_coord[i]*cell_denom[i];
    }
    
    /*--------------------------------------------------
     assign the position of iatm
     --------------------------------------------------*/
    ix=mapp->x->begin()+iatm*x_dim;
    
    /*--------------------------------------------------
     reset the neighboring cell number
     --------------------------------------------------*/
    ineigh=-1;
    
    if(iatm<natms)
    {
        /*--------------------------------------------------
         if this atom does belong to us we need to consider
         all neighbors including the phantom ones
         --------------------------------------------------*/
        next_jatm_p=&GCMC::next_jatm_all;
        /*--------------------------------------------------
         find the first non-empty cell
         --------------------------------------------------*/
        next_jcell_all();
        if(jcell==-1)
        {
            /*--------------------------------------------------
             if we did not find a non-emty neighboring cell
             we have reached the end of the list
             --------------------------------------------------*/
            jatm=-1;
        }
        else
        {
            if(jatm==iatm)
            {
                /*--------------------------------------------------
                 if we get the same atom we need to go to the next
                 atom
                 --------------------------------------------------*/
                next_jatm_all();
                return;
            }
                
            jx=mapp->x->begin()+x_dim*jatm;
            jtype=mapp->type->begin()[jatm];
            rsq=0.0;
            for(int i=0;i<dim;i++)
                rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
            /*--------------------------------------------------
             if we get an atom that is non-interacting we need
             to go the next atom
             --------------------------------------------------*/
            if(rsq>=cut_sq[itype][jtype])
                next_jatm_all();
        }
    }
    else
    {
        /*--------------------------------------------------
         if this atom does not belong to us we need to
         consider only local neighbors excluding the phantom
         ones
         --------------------------------------------------*/
        next_jatm_p=&GCMC::next_jatm_lcl;
        /*--------------------------------------------------
         find the first non-empty local cell
         --------------------------------------------------*/
        next_jcell_lcl();
        if(jcell==-1)
        {
            /*--------------------------------------------------
             if we did not find a non-emty neighboring cell
             we have reached the end of the list
             --------------------------------------------------*/
            jatm=-1;
        }
        else
        {
            /*--------------------------------------------------
             if we get the same atom or a phantom atom we need 
             to go to the next local atom
             --------------------------------------------------*/
            if(jatm>=natms || jatm==iatm)
            {
                next_jatm_lcl();
                return;
            }
            jx=mapp->x->begin()+x_dim*jatm;
            jtype=mapp->type->begin()[jatm];
            rsq=0.0;
            for(int i=0;i<dim;i++)
                rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
            /*--------------------------------------------------
             if we get an atom that is non-interacting we need
             to go the next local atom
             --------------------------------------------------*/
            if(rsq>=cut_sq[itype][jtype])
                next_jatm_lcl();
        }
    }
}
/*--------------------------------------------
 find the next cell that probably
 contains local atoms
 --------------------------------------------*/
inline void GCMC::next_jcell_lcl()
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
        if(jcell_coord[i]<m || jcell_coord[i]>=ncells_per_dim[i]-m)
        {
            next_jcell_lcl();
            return;
        }
        jcell+=cell_denom[i]*jcell_coord[i];
    }
    
    if(head_atm[jcell]==-1)
    {
        next_jcell_lcl();
        return;
    }
    
    jatm=head_atm[jcell];

}
/*--------------------------------------------
 find the next cell
 --------------------------------------------*/
inline void GCMC::next_jcell_all()
{
    ineigh++;
    if(ineigh==nneighs)
    {
        jcell=-1;
        return;
    }
    
    jcell=icell+rel_neigh_lst[ineigh];
    
    if(head_atm[jcell]==-1)
    {
        next_jcell_all();
        return;
    }
    
    jatm=head_atm[jcell];
}
/*--------------------------------------------
 find the next atom
 --------------------------------------------*/
inline void GCMC::next_jatm_all()
{
    if(next_vec_p->begin()[jatm]==-1)
    {
        next_jcell_all();
        if(jcell==-1)
        {
            if(itrial_atm==0 && first_atm_lcl && xchng_mode==INS_MODE)
            {
                iself=0;
                next_jatm_p=&GCMC::next_jatm_self;
                next_jatm_self();
            }
            else
                jatm=-1;
        }
        else
        {
            if(jatm==iatm)
            {
                next_jatm_all();
                return;
            }
            jx=mapp->x->begin()+x_dim*jatm;
            jtype=mapp->type->begin()[jatm];
            rsq=0.0;
            for(int i=0;i<dim;i++)
                rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
            if(rsq>=cut_sq[itype][jtype])
                next_jatm_all();
        }
        return;
    }
    
    
    
    jatm=next_vec_p->begin()[jatm];
    if(jatm==iatm)
    {
        next_jatm_all();
        return;
    }
    jx=mapp->x->begin()+x_dim*jatm;
    jtype=mapp->type->begin()[jatm];
    rsq=0.0;
    for(int i=0;i<dim;i++)
        rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
    if(rsq>=cut_sq[itype][jtype])
        next_jatm_all();
}
/*--------------------------------------------
 find the next lcl atom
 --------------------------------------------*/
inline void GCMC::next_jatm_lcl()
{
    
    if(next_vec_p->begin()[jatm]==-1)
    {
        next_jcell_lcl();
        if(jcell==-1)
        {
            jatm=-1;
        }
        else
        {
            if(jatm>=natms || jatm==iatm)
            {
                next_jatm_lcl();
                return;
            }
            jx=mapp->x->begin()+x_dim*jatm;
            jtype=mapp->type->begin()[jatm];
            rsq=0.0;
            for(int i=0;i<dim;i++)
                rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
            if(rsq>=cut_sq[itype][jtype])
                next_jatm_lcl();
        }
        return;
    }
    
    jatm=next_vec_p->begin()[jatm];
    if(jatm>=natms || jatm==iatm)
    {
        next_jatm_lcl();
        return;
    }
    jx=mapp->x->begin()+x_dim*jatm;
    jtype=mapp->type->begin()[jatm];
    rsq=0.0;
    for(int i=0;i<dim;i++)
        rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
    if(rsq>=cut_sq[itype][jtype])
        next_jatm_lcl();
}
/*--------------------------------------------
 find the next interactin image atom
 --------------------------------------------*/
inline void GCMC::next_jatm_self()
{
    iself++;
    if(iself==ntrial_atms)
    {
        iself=0;
        jatm=-1;
        return;
    }
    jatm=natms+natms_ph+iself;
    jx=ins_buff+iself*dim;
    jtype=gas_type;
    rsq=0.0;
    for(int i=0;i<dim;i++)
        rsq+=(ix[i]-jx[i])*(ix[i]-jx[i]);
    if(rsq>=cut_sq[itype][jtype])
        next_jatm_self();
}
/*--------------------------------------------
 random position for insertion
 --------------------------------------------*/
void GCMC::ins_attmpt()
{
    itrial_atm=-1;
    iatm=-1;
    xchng_mode=INS_MODE;
    first_atm_lcl=true;
    ntrial_atms=1;
    int n_per_dim[dim];
    type0 s,tmp;
    for(int i=0;i<dim;i++)
    {
        n_per_dim[i]=0;
        s=random->uniform();
        if(s <s_lo[i] || s>= s_hi[i])
            first_atm_lcl=false;
        
        tmp=s;
        if(s_lo_ph[i]<=tmp && tmp< s_hi_ph[i])
            ins_s_trials[i][n_per_dim[i]++]=tmp++;
        
        for(int j=0;j<nimages_per_dim[i][0];j++)
            if(s_lo_ph[i]<=tmp && tmp< s_hi_ph[i])
                ins_s_trials[i][n_per_dim[i]++]=tmp++;
        
        tmp=s-1.0;
        for(int j=0;j<nimages_per_dim[i][1];j++)
            if(s_lo_ph[i]<=tmp && tmp< s_hi_ph[i])
                ins_s_trials[i][n_per_dim[i]++]=tmp--;
        ntrial_atms*=n_per_dim[i];
        
    }
    
    int icurs[dim];
    for(int i=0;i<dim;i++)
        icurs[i]=0;
    
    type0* buff=ins_buff;
    int* cell_coord=ins_cell_coord;
    type0** H=atoms->H;
    
    
#ifdef DEBUG_GCMC
    int chk=1;
#endif
    
    for(int i=0;i<ntrial_atms;i++)
    {
        ins_cell[i]=0;

#ifdef DEBUG_GCMC
        chk=1;
#endif
        for(int j=0;j<dim;j++)
        {
            buff[j]=ins_s_trials[j][icurs[j]];
            cell_coord[j]=static_cast<int>(floor((buff[j]-s_lo[j])/cell_size[j]))+m;
            ins_cell[i]+=cell_coord[j]*cell_denom[j];
            
            
            
            
#ifdef DEBUG_GCMC
            if(i==0 && first_atm_lcl)
            {
                if(cell_coord[j]<m || cell_coord[j]>=ncells_per_dim[j]-m)
                    error->abort("this does not work in inserton");
            }
            else
            {
                if(cell_coord[j]<0 || cell_coord[j]>=ncells_per_dim[j])
                    error->abort("this does not work in inserton for ghost");
                
                
                if(cell_coord[j]<m || cell_coord[j]>=ncells_per_dim[j]-m-1)
                    chk*=0;
            }
#endif

            buff[j]=buff[j]*H[j][j];
            for(int k=j+1;k<dim;k++)
                buff[j]+=buff[k]*H[k][j];
        }
        
        
#ifdef DEBUG_GCMC
        if(i!=0 || !first_atm_lcl)
        {
            if(chk)
            {
                error->abort("this does not work in inserton for ghost inside the boundry");
            }
        }
#endif
        
        
        buff+=dim;
        cell_coord+=dim;
        
        icurs[0]++;
        for(int j=0;j<dim-1;j++)
            if(icurs[j]==n_per_dim[j])
            {
                icurs[j]=0;
                icurs[j+1]++;
            }
    }
    
    
    
    if(ntrial_atms)
        next_iatm_ins();
    forcefield->xchng_energy_timer(this);
}
/*--------------------------------------------
 find random atom for deletion
 --------------------------------------------*/
void GCMC::del_attmpt()
{
    iatm=-1;
    itrial_atm=-1;
    xchng_mode=DEL_MODE;
    first_atm_lcl=false;
    ntrial_atms=0;
    if(tot_ngas==0) return;
    int atm_id_=0;
    int idx=-1;
    igas=static_cast<int>(tot_ngas*random->uniform());
    int* id=mapp->id->begin();
    if(ngas_before<=igas && igas<ngas_before+ngas)
    {
        first_atm_lcl=true;
        int n=igas-ngas_before;
        int icount=-1;
        md_type* type=mapp->type->begin();
        
        for(idx=0;icount!=n;idx++)
            if(type[idx]==gas_type) icount++;
        idx--;
        atm_id_=id[idx];
        del_lst[ntrial_atms++]=idx;
    }
    MPI_Allreduce(&atm_id_,&gas_id,1,MPI_INT,MPI_SUM,world);

    
    for(int i=natms;i<natms+natms_ph;i++)
        if(id[i]==gas_id)
            del_lst[ntrial_atms++]=i;
    if(ntrial_atms)
        next_iatm_del();
    forcefield->xchng_energy_timer(this);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::ins_succ()
{
    for(int i=0;i<dim;i++) vel_buff[i]=random->gaussian()*sigma;
    int new_id=get_new_id();
    int iproc_=0,iproc;
    type0* buff=ins_buff;
    int* cell_buff=ins_cell;
    if(first_atm_lcl)
        iproc_=atoms->my_p;
    MPI_Allreduce(&iproc_,&iproc,1,MPI_INT,MPI_SUM,world);
    if(first_atm_lcl)
    {
        atoms->add(1,ntrial_atms-1);
        ngas++;
        
        if(mapp->x_dof)
        {
            bool* dof=mapp->x_dof->begin()+(natms-1)*x_dim;
            for(int i=0;i<dim;i++)
                dof[i]=true;
        }
        memcpy(mapp->x->begin()+x_dim*(natms-1),buff,dim*sizeof(type0));
        memcpy(mapp->x_d->begin()+dim*(natms-1),vel_buff,dim*sizeof(type0));
        mapp->type->begin()[natms-1]=gas_type;
        mapp->id->begin()[natms-1]=new_id;
        cell_vec_p->begin()[natms-1]=ins_cell[0];
        buff+=dim;
        cell_buff++;
        ntrial_atms--;
    }
    else
    {
        if(ntrial_atms)
            atoms->add(0,ntrial_atms);
        if(iproc<atoms->my_p)
            ngas_before++;
    }
    
    type0* x=mapp->x->begin()+x_dim*(natms+natms_ph-ntrial_atms);
    md_type* type=mapp->type->begin()+natms+natms_ph-ntrial_atms;
    int* id=mapp->id->begin()+natms+natms_ph-ntrial_atms;
    int* cell=cell_vec_p->begin()+natms+natms_ph-ntrial_atms;
    for(int i=0;i<ntrial_atms;i++,x+=x_dim,buff+=dim)
    {
        memcpy(x,buff,dim*sizeof(type0));
        type[i]=gas_type;
        id[i]=new_id;
        cell[i]=cell_buff[i];
    }
    
    dof_diff+=dim;
    tot_ngas++;
    atoms->tot_natms++;
    refresh();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::del_succ()
{
    add_del_id();
    if(first_atm_lcl)
    {
        atoms->del(del_lst,1,del_lst+1,ntrial_atms-1);
        ngas--;
    }
    else
    {
        if(ntrial_atms)
            atoms->del(del_lst,0,del_lst,ntrial_atms);
        if(igas<ngas_before)
            ngas_before--;
    }
    
    dof_diff-=dim;
    tot_ngas--;
    atoms->tot_natms--;
    if(ntrial_atms)
        refresh();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::next_iatm()
{
    if(xchng_mode==INS_MODE)
        next_iatm_ins();
    else
        next_iatm_del();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void GCMC::next_jatm()
{
    (this->*next_jatm_p)();
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool GCMC::decide(type0& en)
{
    type0 fac;
    if(xchng_mode==INS_MODE)
    {
        fac=z_fac*vol/((static_cast<type0>(tot_ngas)+1.0)*exp(beta*en));
        if(random->uniform()<fac)
        {
            ins_succ();
            return true;
        }
    }
    else
    {
        fac=static_cast<type0>(tot_ngas)*exp(beta*en)/(z_fac*vol);
        if(random->uniform()<fac)
        {
            del_succ();
            return true;
        }
    }
    return false;
}
/*--------------------------------------------
 
 --------------------------------------------*/
bool GCMC::calc_mu(type0& en)
{
    if(xchng_mode==INS_MODE)
    {
        //en-kbT*log(vol/static_cast<type0>(tot_ngas+1))+zz_fac;
    }
    else
    {
        //en-kbT*log(vol/static_cast<type0>(tot_ngas))+zz_fac;
    }
    return false;
}
/*--------------------------------------------
 construct the bin list
 --------------------------------------------*/
void GCMC::xchng(bool box_chng,int nattmpts)
{
    if(box_chng)
        box_setup();
    
    atoms->init_xchng();
    
    /*--------------------------------------------------
     here we allocate the memory for cell_vec & next_vec
     --------------------------------------------------*/
    cell_vec_p=new Vec<int>(atoms,1);
    cell_vec_p->resize(0);
    cell_vec_p->resize(natms+natms_ph);
    
    next_vec_p=new Vec<int>(atoms,1);
    next_vec_p->resize(0);
    next_vec_p->resize(natms+natms_ph);
    
    
    /*--------------------------------------------------
     resize id and update so we have the ids of phantom
     atoms (useful for deletion)
     --------------------------------------------------*/
    mapp->id->resize(natms+natms_ph);
    atoms->update(mapp->id);
    
    
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
    int nall=natms+natms_ph;
    type0* s=mapp->x->begin()+(nall-1)*x_dim;
    
#ifdef DEBUG_GCMC
    int y,chk=1;
#endif
    
    for(int i=nall-1;i>-1;i--,s-=x_dim)
    {
#ifndef DEBUG_GCMC
        find_cell_no(s,cell_vec[i]);
#endif
        
#ifdef DEBUG_GCMC
        chk=1;
        cell_vec[i]=0;
        for(int j=0;j<dim && cell_vec[i]!=-1;j++)
        {
            if(s[j]>=s_lo_ph[j] && s[j]<s_hi_ph[j])
            {
                y=m+static_cast<int>(floor((s[j]-s_lo[j])/cell_size[j]));
                if(i<natms)
                {
                    if(y<m || y>=ncells_per_dim[j]-m)
                        error->abort("this does not work");
                }
                else
                {
                    if(y<0 || y>=ncells_per_dim[j])
                        error->abort("this does not work for ghost");
                    if(y<m || y>=ncells_per_dim[j]-m-1)
                        chk*=0;
                }
                
                cell_vec[i]+=cell_denom[j]*y;
            }
            else
            {
                cell_vec[i]=-1;
            }
        }
        if(cell_vec[i]!=-1 && (i>=natms && chk))
        {
            error->abort("this does not work for ghost 2nd");
        }
#endif
        
        
        
        if(cell_vec[i]!=-1)
        {
            next_vec[i]=head_atm[cell_vec[i]];
            head_atm[cell_vec[i]]=i;
        }
    }
    
    neighbor->create_list(box_chng);
    forcefield->init_xchng();
    
    ngas=0;
    md_type* type=mapp->type->begin();
    for(int i=0;i<natms;i++)
        if(type[i]==gas_type) ngas++;
    MPI_Scan(&ngas,&ngas_before,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&ngas,&tot_ngas,1,MPI_INT,MPI_SUM,world);
    ngas_before-=ngas;
    
    int max_id_=0;
    int* id=mapp->id->begin();
    for(int i=0;i<natms;i++)
        max_id_=MAX(id[i],max_id_);
    MPI_Allreduce(&max_id_,&max_id,1,MPI_INT,MPI_MAX,world);
    
    dof_diff=0;
    for(int i=0;i<nattmpts;i++)
    {
        if(random->uniform()<0.5)
            ins_attmpt();
        else
            del_attmpt();
    }

    mapp->id->resize(natms);
    delete next_vec_p;
    delete cell_vec_p;
    
    forcefield->fin_xchng();
    atoms->fin_xchng();
}






