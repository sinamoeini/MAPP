#include "neighbor.h"
#include "atoms.h"
#include "ff.h"
#include "xmath.h"
#include "atom_types.h"
using namespace MAPP_NS;
/*---------------------------------------------------------------------------
      ___   _____   _____       ___  ___   _____
     /   | |_   _| /  _  \     /   |/   | /  ___/
    / /| |   | |   | | | |    / /|   /| | | |___
   / / | |   | |   | | | |   / / |__/ | | \___  \
  / /  | |   | |   | |_| |  / /       | |  ___| |
 /_/   |_|   |_|   \_____/ /_/        |_| /_____/
 
 ---------------------------------------------------------------------------*/
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Atoms::Atoms(MAPP* mapp)
:InitPtrs(mapp)
{
    
    //world=communicator;
    dimension=3;
    tot_natms=0;
    natms=0;
    natms_ph=0;
    atm_vec_ph_size=0;
    atm_vec_size=0;
    tot_cut_ph=0.0;
    skin=0.5;
    
    CREATE2D(H,dimension,dimension);
    CREATE2D(B,dimension,dimension);
    
    CREATE1D(s_lo,dimension);
    CREATE1D(s_hi,dimension);
    CREATE1D(s_ph_lo,dimension);
    CREATE1D(s_ph_hi,dimension);
    CREATE1D(cut_ph_s,dimension);
    CREATE2D(comm_need,dimension,2);
    
    CREATE1D(tot_p_grid,dimension);
    CREATE1D(my_loc,dimension);
    CREATE2D(neigh_p,dimension,2);
    
    
    //communication related parameters
    MPI_Comm_rank(world,&my_p_no);
    MPI_Comm_size(world,&tot_p);
    ph_lst= new SwapLst(mapp);
    
    snd_buff_0_capacity=0;
    snd_buff_1_capacity=0;
    rcv_buff_capacity=0;
    snd_ph_buff_capacity=0;
    rcv_ph_buff_capacity=0;
    
    no_vecs=0;
    tot_byte_size=0;
    
    
    int name_length;
    char* node_name;
    CREATE1D(node_name,MPI_MAX_PROCESSOR_NAME);
    MPI_Get_processor_name(node_name,&name_length);
    node_name[name_length] = '\0';
    name_length++;
    

    int* name_lenghts;
    CREATE1D(name_lenghts,tot_p);
    for(int i=0;i<tot_p;i++)
        name_lenghts[i]=0;
    
    
    name_lenghts[my_p_no]=name_length;
    int* all_name_lenghts;
    CREATE1D(all_name_lenghts,tot_p);
    MPI_Allreduce(name_lenghts,all_name_lenghts,tot_p, MPI_INT,MPI_SUM,world);
    delete [] name_lenghts;
    
    
    
    
    char** all_names;
    CREATE1D(all_names,tot_p);
    for(int i=0;i<tot_p;i++)
        CREATE1D(all_names[i],all_name_lenghts[i]);
    
    
    
    for(int i=0;i<name_length;i++)
        all_names[my_p_no][i]=node_name[i];
    delete [] node_name;
    
    
    for(int i=0;i<tot_p;i++)
        MPI_Bcast(all_names[i],all_name_lenghts[i],MPI_CHAR,i,world);
    
    
    int* node_no;
    CREATE1D(node_no,tot_p);
    
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
    delete [] all_names;
    
    
    CREATE1D(p_per_n,tot_n);
    
    for(int i=0;i<tot_n;i++)
        p_per_n[i]=0;
    
    for(int i=0;i<tot_p;i++)
        p_per_n[node_no[i]]++;
    
    my_n_no=node_no[my_p_no];

    
    CREATE1D(n_p_grid,tot_n);
    for(int i=0;i<tot_n;i++)
        CREATE1D(n_p_grid[i],p_per_n[i]);
    
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
    delete [] node_no;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Atoms::~Atoms()
{

    if(dimension)
    {
        for(int i=0;i<dimension;i++)
        {
            delete [] H[i];
            delete [] B[i];
            delete [] neigh_p[i];
            delete [] comm_need[i];
        }
        delete [] H;
        delete [] B;
        delete [] neigh_p;
        delete [] comm_need;
        
        delete [] s_lo;
        delete [] s_hi;
        delete [] s_ph_lo;
        delete [] s_ph_hi;
        delete [] cut_ph_s;
        
        delete [] tot_p_grid;
        delete [] my_loc;
        
    }
    
    delete ph_lst;
    
    
    if(snd_buff_0_capacity)
        delete [] snd_buff_0;
    if(snd_buff_1_capacity)
        delete [] snd_buff_1;
    if(rcv_buff_capacity)
        delete [] rcv_buff;
    if(snd_ph_buff_capacity)
        delete [] snd_ph_buff;
    if(rcv_ph_buff_capacity)
        delete [] rcv_ph_buff;
    
    
    for(int i=0;i<no_vecs;i++)
        delete vectors[i];
    if(no_vecs)
        delete [] vectors;
    
    for(int i=0;i<tot_n;i++)
        delete [] n_p_grid[i];
    delete [] n_p_grid;
    delete [] p_per_n;
    
}
/*--------------------------------------------
 clean up the and reset
 --------------------------------------------*/
inline void Atoms::setup_comm_need()
{
    
    for (int i=0;i<dimension;i++)
        comm_need[i][0]=comm_need[i][1]=static_cast<int>
        (cut_ph_s[i]/(s_hi[i]-s_lo[i]))+1;
    /*
    int chk,chk_lcl,snd_p,rcv_p,icomm;
    MPI_Status status;
    type0 snd_frac,rcv_frac,my_s;
    for(int idim=0;idim<dimension;idim++)
        for(int idir=0;idir<2;idir++)
        {
            
            snd_p=neigh_p[idim][idir];
            rcv_p=neigh_p[idim][1-idir];
            
            snd_frac=s_hi[idim]-s_lo[idim];
            my_s=0.0;
            icomm=0;
            chk=0;
            chk_lcl=0;
            while (chk==0)
            {
                MPI_Sendrecv(&snd_frac,1,
                MPI_TYPE0,snd_p,0,&rcv_frac,
                1,MPI_TYPE0,rcv_p,0,world,
                &status);
                
                my_s+=rcv_frac;
                snd_frac=rcv_frac;
                
                if(cut_ph_s[idim]<=my_s)
                    chk_lcl=1;
    
                MPI_Allreduce(&chk_lcl,&chk,1,MPI_INT,MPI_MIN,world);
                
                
                icomm++;
            }
            comm_need[idim][idir]=icomm;
        }
    */
}
/*--------------------------------------------
 clean up the and reset
 --------------------------------------------*/
void Atoms::restart()
{
    
    
    delete ph_lst;
    
    
    if(snd_buff_0_capacity)
        delete [] snd_buff_0;
    if(snd_buff_1_capacity)
        delete [] snd_buff_1;
    if(rcv_buff_capacity)
        delete [] rcv_buff;
    if(snd_ph_buff_capacity)
        delete [] snd_ph_buff;
    if(rcv_ph_buff_capacity)
        delete [] rcv_ph_buff;
    
    if(no_vecs)
        delete [] vectors;
    
    ph_lst= new SwapLst(mapp);
    
    snd_buff_0_capacity=0;
    snd_buff_1_capacity=0;
    rcv_buff_capacity=0;
    snd_ph_buff_capacity=0;
    rcv_ph_buff_capacity=0;
    
    no_vecs=0;
    tot_byte_size=0;
    
    tot_natms=0;
    natms=0;
    natms_ph=0;
    atm_vec_ph_size=0;
    atm_vec_size=0;
    tot_cut_ph=0.0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::chng_dim(int dim)
{
    if(dim<1)
        error->abort("box dimension should be more than 0");
    
    if(dimension)
    {
        for(int i=0;i<dimension;i++)
        {
            delete [] H[i];
            delete [] B[i];
            delete [] neigh_p[i];
            delete [] comm_need[i];
        }
        delete [] H;
        delete [] B;
        delete [] neigh_p;
        delete [] comm_need;
        
        delete [] s_lo;
        delete [] s_hi;
        delete [] s_ph_lo;
        delete [] s_ph_hi;
        delete [] cut_ph_s;
        
        delete [] tot_p_grid;
        delete [] my_loc;
        
    }
    
    
    dimension=dim;
    
    CREATE2D(H,dimension,dimension);
    CREATE2D(B,dimension,dimension);
    
    CREATE1D(s_lo,dimension);
    CREATE1D(s_hi,dimension);
    CREATE1D(s_ph_lo,dimension);
    CREATE1D(s_ph_hi,dimension);
    CREATE1D(cut_ph_s,dimension);
    CREATE2D(comm_need,dimension,2);
    
    CREATE1D(tot_p_grid,dimension);
    CREATE1D(my_loc,dimension);
    CREATE2D(neigh_p,dimension,2);
}
/*--------------------------------------------
 this fuction:
 
 0. sets maximum thickness of phantom area
 (tot_cut_ph) using max cutoff & skin
 1. converts tot_cut_ph to fractional 
 coordinates (cut_ph_s[i]) 
 2. finds the boundary of our local box 
 + phantom area, in terms of fractional 
 coordinates (s_ph_lo[i] & s_ph_hi[i])
 3. creates a new swap list, in case the new 
 number of swaps is not the same as before
 --------------------------------------------*/
void Atoms::set_max_cutoff(type0 cut)
{
    tot_cut_ph=cut+skin;
    type0 tmp;
    for (int i=0;i<dimension;i++)
    {
        tmp=0.0;
        for(int j=i;j<dimension;j++)
            tmp+=B[j][i]*B[j][i];
        
        cut_ph_s[i]=tot_cut_ph*sqrt(tmp);
    }
       
    for (int i=0;i<dimension;i++)
    {
        s_ph_hi[i]=s_hi[i]-cut_ph_s[i];
        s_ph_lo[i]=s_lo[i]+cut_ph_s[i];
    }
    /*
    for (int i=0;i<dimension;i++)
        comm_need[i][0]=comm_need[i][1]=static_cast<int>
        (cut_ph_s[i]/(s_hi[i]-s_lo[i]))+1;
     */
    setup_comm_need();
    
    int nswap=0;
    for(int i=0;i<dimension;i++)
        nswap+=comm_need[i][0]+comm_need[i][1];
    
    ph_lst->newlist(nswap);
}
/*--------------------------------------------
 this function:
 0. finds the swap lists for each dimension
 & direction
 1. uses xchng_ph(int,int,int*,int
 ,class VecLst*) to transfer said lists to
 neighboring processors
 
 caveat: it is assumed that x is already is
 converted to s, i.e. the coordinates are
 frational.
 --------------------------------------------*/
void Atoms::setup_ph(int box_chng,class VecLst* list)
{
    
    /*
     if the box size has changed
     we need to:
     0. find cut_ph_s
     1. assign s_ph_lo & s_ph_hi
     2. find all the number of needed communications
     3. and make a new swap list
     */
    if(box_chng)
    {
        type0 tmp;
        for (int i=0;i<dimension;i++)
        {
            tmp=0.0;
            for(int j=i;j<dimension;j++)
                tmp+=B[j][i]*B[j][i];
            
            cut_ph_s[i]=tot_cut_ph*sqrt(tmp);
        }
        
        for (int i=0;i<dimension;i++)
        {
            s_ph_hi[i]=s_hi[i]-cut_ph_s[i];
            s_ph_lo[i]=s_lo[i]+cut_ph_s[i];
        }

        setup_comm_need();
        
        int nswap=0;
        for(int i=0;i<dimension;i++)
            nswap+=comm_need[i][0]+comm_need[i][1];
        
        ph_lst->newlist(nswap);
    }
    
    
    
    
    natms_ph=0;
    ph_lst->reset();
    int lo_r_dir,hi_r_dir,lo_r_dim,hi_r_dim,tmp_r;
    int icomm,iatm,idim;
    type0* s;
    int s_dim=vectors[0]->dim;
    int iswap=0;
    
    lo_r_dim=0;
    hi_r_dim=natms+natms_ph;
    for(idim=0;idim<dimension;idim++)
    {
        lo_r_dir=lo_r_dim;
        hi_r_dir=hi_r_dim;
        
        /*
         send to back & reciev from front dir=0
         */
        for(icomm=0;icomm<comm_need[idim][0];icomm++)
        {
            tmp_r=natms+natms_ph;
            vectors[0]->ret(s);
            for(iatm=lo_r_dir;iatm<hi_r_dir;iatm++)
                if(s[iatm*s_dim+idim]<=s_ph_lo[idim])
                    ph_lst->add(iswap,iatm);
            /*
             find the snd_p[iswap and rcv_p[iswap]
             int snd_proc=neigh_p[idim][0];
             int rcv_proc=neigh_p[idim][1];
             */
            ph_lst->snd_p[iswap]=neigh_p[idim][0];
            ph_lst->rcv_p[iswap]=neigh_p[idim][1];
            /*
             pbc_correction[iswap]=0 or (1+idim) or -(1+idim)
             
             if((dir==0&&my_loc[dim]==tot_p_grid[dim]-1) ||
             (dir==1&&my_loc[dim]==0))
             {
             int last_atm_no=natms+natms_ph;
             int strt_atm_no=last_atm_no-no_new_atms;
             type0* s;
             vectors[0]->ret(s);
             int s_dim=vectors[0]->dim;
             if(dir==0)
             for (int iatm=strt_atm_no;iatm<last_atm_no;iatm++)
             s[iatm*s_dim+dim]++;
             else if(dir==1)
             for (int iatm=strt_atm_no;iatm<last_atm_no;iatm++)
             s[iatm*s_dim+dim]--;
             }
             
             */
            if(my_loc[idim]==tot_p_grid[idim]-1)
                ph_lst->pbc_correction[iswap]=1+idim;
            else
                ph_lst->pbc_correction[iswap]=0;
            
            xchng_ph(iswap,list);
            
            hi_r_dir=natms+natms_ph;
            lo_r_dir=tmp_r;
            iswap++;
        }
        
        lo_r_dir=lo_r_dim;
        hi_r_dir=hi_r_dim;
        
        /*
         send to front & recieve from back dir=1
         */
        for(icomm=0;icomm<comm_need[idim][1];icomm++)
        {
            tmp_r=natms+natms_ph;
            vectors[0]->ret(s);
            for(iatm=lo_r_dir;iatm<hi_r_dir;iatm++)
                if(s_ph_hi[idim]<=s[iatm*s_dim+idim])
                    ph_lst->add(iswap,iatm);
            
            /*
             find the snd_p[iswap and rcv_p[iswap]
             int snd_proc=neigh_p[idim][0];
             int rcv_proc=neigh_p[idim][1];
             */
            ph_lst->snd_p[iswap]=neigh_p[idim][1];
            ph_lst->rcv_p[iswap]=neigh_p[idim][0];
            /*
             pbc_correction[iswap]=0 or (1+idim) or -(1+idim)
             
             if((dir==0&&my_loc[dim]==tot_p_grid[dim]-1) ||
             (dir==1&&my_loc[dim]==0))
             {
             int last_atm_no=natms+natms_ph;
             int strt_atm_no=last_atm_no-no_new_atms;
             type0* s;
             vectors[0]->ret(s);
             int s_dim=vectors[0]->dim;
             if(dir==0)
             for (int iatm=strt_atm_no;iatm<last_atm_no;iatm++)
             s[iatm*s_dim+dim]++;
             else if(dir==1)
             for (int iatm=strt_atm_no;iatm<last_atm_no;iatm++)
             s[iatm*s_dim+dim]--;
             }
             
             */
            if(my_loc[idim]==0)
                ph_lst->pbc_correction[iswap]=-(1+idim);
            else
                ph_lst->pbc_correction[iswap]=0;
            
            xchng_ph(iswap,list);
            
            hi_r_dir=natms+natms_ph;
            lo_r_dir=tmp_r;
            iswap++;
        }
        
        lo_r_dim=0;
        hi_r_dim=natms+natms_ph;
    }
    
}
/*--------------------------------------------
 this function:
 0. sends/recieves phantom atoms to/from
 neighboring processors.
 1. if needed makes correction to positions (s)
 due to periodic boundary conditions
 --------------------------------------------*/
void Atoms::xchng_ph(int iswap,class VecLst* list)
{
    if(my_p_no==ph_lst->rcv_p[iswap])
    {
        
        ph_lst->rcv_size[iswap]=ph_lst->snd_size[iswap];
        
        /*
         now all of the values of this swap is set
         */
        
        /*
         first make sure that there is enough space for the new
         phantom atoms
         */
        if(ph_lst->rcv_size[iswap]+natms+natms_ph>atm_vec_ph_size)
        {
            grow(1,ph_lst->rcv_size[iswap],list);
        }
        
        /*
         now do the copying
         maybe later you can make the copying much more efficient
         keep the old copy & add a new one to take care of the lists
         */
        
        for(int ivec=0;ivec<list->ph_no_vecs;ivec++)
            vectors[list->ph_vec_list[ivec]]->copy2bottom(ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
        
        /*
         in order to be consistant with unpack_ph above
         */
        natms_ph+=ph_lst->rcv_size[iswap];
    }
    else
    {
        MPI_Request request[2];
        MPI_Status status[2];
        
        /*
         first find the number of atoms to rcv in this swap, 
         i.e. ph_lst->rcv_size[iswap]
         */
        MPI_Sendrecv(&ph_lst->snd_size[iswap],1,MPI_INT,ph_lst->snd_p[iswap],
        0,
        &ph_lst->rcv_size[iswap],1,MPI_INT,ph_lst->rcv_p[iswap],
        0,
        world,&status[0]);
        /*
         now all of the values of this swap is set
         */
        
        
        
        if(list->ph_no_vecs!=1)
        {
            /*
             if there is nothing to send & something to recieve, just recieve
             */
            if (ph_lst->rcv_size[iswap]&&ph_lst->snd_size[iswap]==0)
            {
                /*
                 adjust the size of the buffers
                 */
                int rcv_buff_size=list->ph_byte_size*ph_lst->rcv_size[iswap];
                buff_size_management(rcv_ph_buff,rcv_ph_buff_capacity,rcv_buff_size);
                
                /*
                 prepare for recieve
                 */
                MPI_Irecv(rcv_ph_buff,rcv_buff_size,
                MPI_BYTE,ph_lst->rcv_p[iswap],
                0,
                world,&request[0]);
                
                /*
                 do the recieving
                 */
                MPI_Wait(&request[0],&status[0]);
                
                /*
                 unpack the recieved buffer
                 */
                unpack_ph(rcv_ph_buff,list,ph_lst->rcv_size[iswap]);
            }
            /*
             if there is something to send & nothing to recieve, just send
             */
            else if (ph_lst->rcv_size[iswap]==0&&ph_lst->snd_size[iswap])
            {
                /*
                 adjust the size of the buffers
                 */
                int snd_buff_size=list->ph_byte_size*ph_lst->snd_size[iswap];
                buff_size_management(snd_ph_buff,snd_ph_buff_capacity,snd_buff_size);
                
                /*
                 now pack up the list and put them in snd_ph_buff
                 */
                pack_ph(snd_ph_buff,list,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
                
                /*
                 prepare for send
                 */
                MPI_Isend(snd_ph_buff,snd_buff_size,
                MPI_BYTE,ph_lst->snd_p[iswap],
                0,
                world,&request[1]);
                
                /*
                 do the sending
                 */
                MPI_Wait(&request[1],&status[1]);
            }
            /*
             if there is something to send & something to recieve, send and recieve
             */
            else if (ph_lst->rcv_size[iswap]&&ph_lst->snd_size[iswap])
            {
                /*
                 adjust the size of the buffers
                 */
                int snd_buff_size=list->ph_byte_size*ph_lst->snd_size[iswap];
                int rcv_buff_size=list->ph_byte_size*ph_lst->rcv_size[iswap];
                buff_size_management(snd_ph_buff,snd_ph_buff_capacity,snd_buff_size);
                buff_size_management(rcv_ph_buff,rcv_ph_buff_capacity,rcv_buff_size);
                
                
                /*
                 now pack up the list and put them in snd_ph_buff
                 */
                pack_ph(snd_ph_buff,list,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
                
                
                /*
                 prepare for recieve
                 */
                MPI_Irecv(rcv_ph_buff,rcv_buff_size,
                MPI_BYTE,ph_lst->rcv_p[iswap],
                0,
                world,&request[0]);
                
                /*
                 prepare for send
                 */
                MPI_Isend(snd_ph_buff,snd_buff_size,
                MPI_BYTE,ph_lst->snd_p[iswap],
                0,
                world,&request[1]);
                
                /*
                 do the sending & recieving
                 */
                MPI_Waitall(2,request,status);
                
                /*
                 unpack the recieved buffer
                 */
                unpack_ph(rcv_ph_buff,list,ph_lst->rcv_size[iswap]);
            }
            
            /*
             remember that unpack_ph will take care of increasing the size of 
             the atomic vectors if necessary
             */
            
            /*
             if there was something to recieve natms_ph are updated by now
             */
        }
        else
        {
            /*
             first make sure that there is enough space for the new
             phantom atoms
             */
            if(ph_lst->rcv_size[iswap]+natms+natms_ph>atm_vec_ph_size)
            {
                grow(1,ph_lst->rcv_size[iswap],list);
            }
            
            /*
             if there is nothing to send & something to recieve, just recieve
             */
            if (ph_lst->rcv_size[iswap]&&ph_lst->snd_size[iswap]==0)
            {
                
                /*
                 prepare for recieve
                 */
                MPI_Irecv(vectors[list->ph_vec_list[0]]->ret(natms+natms_ph),ph_lst->rcv_size[iswap]*list->ph_byte_size,
                MPI_BYTE,ph_lst->rcv_p[iswap],
                0,
                world,&request[0]);
                
                /*
                 do the recieving
                 */
                MPI_Wait(&request[0],&status[0]);
            }
            /*
             if there is something to send & nothing to recieve, just send
             */
            else if (ph_lst->rcv_size[iswap]==0&&ph_lst->snd_size[iswap])
            {
                /*
                 adjust the size of the buffers
                 */
                int snd_buff_size=list->ph_byte_size*ph_lst->snd_size[iswap];
                buff_size_management(snd_ph_buff,snd_ph_buff_capacity,snd_buff_size);
                
                /*
                 now pack up the list and put them in snd_ph_buff
                 */
                pack_ph(snd_ph_buff,list,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
                
                /*
                 prepare for send
                 */
                MPI_Isend(snd_ph_buff,snd_buff_size,
                MPI_BYTE,ph_lst->snd_p[iswap],
                0,
                world,&request[1]);
                
                /*
                 do the sending
                 */
                MPI_Wait(&request[1],&status[1]);
            }
            /*
             if there is something to send & something to recieve, send and recieve
             */
            else if (ph_lst->rcv_size[iswap]&&ph_lst->snd_size[iswap])
            {
                /*
                 adjust the size of the buffers
                 */
                int snd_buff_size=list->ph_byte_size*ph_lst->snd_size[iswap];
                buff_size_management(snd_ph_buff,snd_ph_buff_capacity,snd_buff_size);
                
                /*
                 now pack up the list and put them in snd_ph_buff
                 */
                pack_ph(snd_ph_buff,list,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
                
                
                /*
                 prepare for recieve
                 */
                MPI_Irecv(vectors[list->ph_vec_list[0]]->ret(natms+natms_ph),ph_lst->rcv_size[iswap]*list->ph_byte_size,
                MPI_BYTE,ph_lst->rcv_p[iswap],
                0,
                world,&request[0]);
                
                /*
                 prepare for send
                 */
                MPI_Isend(snd_ph_buff,snd_buff_size,
                MPI_BYTE,ph_lst->snd_p[iswap],
                0,
                world,&request[1]);
                
                /*
                 do the sending & recieving
                 */
                MPI_Waitall(2,request,status);
                
            }
            
            /*
             in order to be consistant with unpack_ph
             */
            natms_ph+=ph_lst->rcv_size[iswap];
        }

    }

    /*
     now it is time for correction due to periodic boundary conditions 
     remember that we are in s mode now
     deifnitely x (or s) is included in this list, so there is no need to 
     check to see wether it is there
     */
    
    if(ph_lst->pbc_correction[iswap])
    {
        type0* s;
        vectors[0]->ret(s);
        int s_dim=vectors[0]->dim;
        
        int lo_atm=natms+natms_ph-ph_lst->rcv_size[iswap];
        
        if(ph_lst->pbc_correction[iswap]>0)
        {
            int idim=ph_lst->pbc_correction[iswap]-1;
            
            for(int iatm=lo_atm;iatm<natms+natms_ph;iatm++)
                s[iatm*s_dim+idim]++;
        }
        else
        {
            int idim=-ph_lst->pbc_correction[iswap]-1;
            for(int iatm=lo_atm;iatm<natms+natms_ph;iatm++)
                s[iatm*s_dim+idim]--;
            
        }
    }


    
}
/*--------------------------------------------
 updates the phantom atoms of a group of vectors
 --------------------------------------------*/
void Atoms::xchng_ph(int iswap,int* vec_list
,int no_vecs,int vec_byte_size)
{
    
    
    
    if(my_p_no==ph_lst->rcv_p[iswap])
    {
        
        for(int ivec=0;ivec<no_vecs;ivec++)
            vectors[vec_list[ivec]]->copy2bottom(ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
         
        /*
         in order to be consistant with unpack_ph above
         */
        natms_ph+=ph_lst->rcv_size[iswap];
    }
    else
    {
        if(no_vecs!=1)
        {
            MPI_Request request[2];
            MPI_Status status[2];
            if (ph_lst->rcv_size[iswap]&&ph_lst->snd_size[iswap]==0)
            {
                int rcv_buff_size=vec_byte_size*ph_lst->rcv_size[iswap];
                buff_size_management(rcv_ph_buff,rcv_ph_buff_capacity,rcv_buff_size);
                
                MPI_Irecv(rcv_ph_buff,rcv_buff_size,
                          MPI_BYTE,ph_lst->rcv_p[iswap],
                          0,
                          world,&request[0]);
                MPI_Wait(&request[0],&status[0]);
                
                unpack_ph(rcv_ph_buff,vec_list,no_vecs,ph_lst->rcv_size[iswap]);
            }
            else if (ph_lst->rcv_size[iswap]==0&&ph_lst->snd_size[iswap])
            {
                int snd_buff_size=vec_byte_size*ph_lst->snd_size[iswap];
                buff_size_management(snd_ph_buff,snd_ph_buff_capacity,snd_buff_size);
                
                pack_ph(snd_ph_buff,vec_list,no_vecs,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
                
                MPI_Isend(snd_ph_buff,snd_buff_size,
                          MPI_BYTE,ph_lst->snd_p[iswap],
                          0,
                          world,&request[1]);
                MPI_Wait(&request[1],&status[1]);
            }
            else if (ph_lst->rcv_size[iswap]&&ph_lst->snd_size[iswap])
            {
                
                int snd_buff_size=vec_byte_size*ph_lst->snd_size[iswap];
                int rcv_buff_size=vec_byte_size*ph_lst->rcv_size[iswap];
                buff_size_management(snd_ph_buff,snd_ph_buff_capacity,snd_buff_size);
                buff_size_management(rcv_ph_buff,rcv_ph_buff_capacity,rcv_buff_size);
                
                pack_ph(snd_ph_buff,vec_list,no_vecs,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
                
                MPI_Irecv(rcv_ph_buff,rcv_buff_size,
                          MPI_BYTE,ph_lst->rcv_p[iswap],
                          0,
                          world,&request[0]);
                MPI_Isend(snd_ph_buff,snd_buff_size,
                          MPI_BYTE,ph_lst->snd_p[iswap],
                          0,
                          world,&request[1]);
                MPI_Waitall(2,request,status);
                
                unpack_ph(rcv_ph_buff,vec_list,no_vecs,ph_lst->rcv_size[iswap]);
            }
        }
        else
        {
            MPI_Request request[2];
            MPI_Status status[2];
            if (ph_lst->rcv_size[iswap]&&ph_lst->snd_size[iswap]==0)
            {
                MPI_Irecv(vectors[vec_list[0]]->ret(natms+natms_ph),ph_lst->rcv_size[iswap]*vec_byte_size,
                          MPI_BYTE,ph_lst->rcv_p[iswap],
                          0,
                          world,&request[0]);
                MPI_Wait(&request[0],&status[0]);
            }
            else if (ph_lst->rcv_size[iswap]==0&&ph_lst->snd_size[iswap])
            {
                int snd_buff_size=vec_byte_size*ph_lst->snd_size[iswap];
                buff_size_management(snd_ph_buff,snd_ph_buff_capacity,snd_buff_size);
                
                pack_ph(snd_ph_buff,vec_list,no_vecs,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
                
                MPI_Isend(snd_ph_buff,snd_buff_size,
                          MPI_BYTE,ph_lst->snd_p[iswap],
                          0,
                          world,&request[1]);
                MPI_Wait(&request[1],&status[1]);
            }
            else if (ph_lst->rcv_size[iswap]&&ph_lst->snd_size[iswap])
            {
                
                int snd_buff_size=vec_byte_size*ph_lst->snd_size[iswap];
                buff_size_management(snd_ph_buff,snd_ph_buff_capacity,snd_buff_size);
                
                pack_ph(snd_ph_buff,vec_list,no_vecs,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
                
                MPI_Irecv(vectors[vec_list[0]]->ret(natms+natms_ph),ph_lst->rcv_size[iswap]*vec_byte_size,
                          MPI_BYTE,ph_lst->rcv_p[iswap],
                          0,
                          world,&request[0]);
                MPI_Isend(snd_ph_buff,snd_buff_size,
                          MPI_BYTE,ph_lst->snd_p[iswap],
                          0,
                          world,&request[1]);
                MPI_Waitall(2,request,status);
            }
            
            /*
             in order to be consistant with unpack_ph above
             */
            natms_ph+=ph_lst->rcv_size[iswap];
            
        }
    }
    
    
    /*
     now it is time for correction due to periodic boundary conditions
     remember that we are in x mode now
     
     first we need to know if x is included in these vectors
     */
    if(vec_list[0]==0)
    {
        if(ph_lst->pbc_correction[iswap])
        {
            type0* x;
            vectors[0]->ret(x);
            int x_dim=vectors[0]->dim;
            
            int lo_atm=natms+natms_ph-ph_lst->rcv_size[iswap];
            
            if(ph_lst->pbc_correction[iswap]>0)
            {
                int idim=ph_lst->pbc_correction[iswap]-1;
                int jdim;
                for(int iatm=lo_atm;iatm<natms+natms_ph;iatm++)
                    for(jdim=0;jdim<idim+1;jdim++)
                        x[iatm*x_dim+jdim]+=H[idim][jdim];
            }
            else
            {
                int idim=-ph_lst->pbc_correction[iswap]-1;
                int jdim;
                for(int iatm=lo_atm;iatm<natms+natms_ph;iatm++)
                    for(jdim=0;jdim<idim+1;jdim++)
                        x[iatm*x_dim+jdim]-=H[idim][jdim];
                
            }
        }
    }
    
}
/*--------------------------------------------
 this function:
 
 0. exchanges own atoms with neighboring 
 processors so that each processor has their
 own proper atoms (limited to their local 
 domain)
 
 caveat: it is assumed that x is already is
 converted to s, i.e. the coordinates are
 frational.
 --------------------------------------------*/
void Atoms::xchng_cmplt(class VecLst* list)
{
    for(int idim=0;idim<dimension;idim++)
        xchng_cmplt(idim,list);

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::xchng_cmplt(int idim,class VecLst* list)
{
    if(tot_p_grid[idim]==1)
        return;
    int icomp;

    type0* s;
    vectors[0]->ret(s);
    int s_dim=vectors[0]->dim;
    
    int snd_buff_0_size=0;
    int snd_buff_1_size=0;
    
    int iatm=0;
    
    while (iatm<natms)
    {
        icomp=s_dim*iatm;
        if(s[icomp+idim]<s_lo[idim])
        {
            if(my_loc[idim]!=tot_p_grid[idim]-1)
                pack_cmplt_prtl(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
            else
            {
                if(0.5<=s[icomp+idim])
                    pack_cmplt_prtl(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
                else
                    pack_cmplt_prtl(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
            }
            
        }
        else if(s[icomp+idim]>=s_hi[idim])
        {
            if(my_loc[idim]!=0)
                pack_cmplt_prtl(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
            else
            {
                if (s[icomp+idim]<0.5)
                    pack_cmplt_prtl(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
                else
                    pack_cmplt_prtl(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
            }
            
        }
        else iatm++;
    }
    
    int rcv_buff_size;
    MPI_Request request;
    MPI_Status status;
    
    int no_atms;
    int d_size=sizeof(type0);
    
    int max_snd_buff_0_size;
    MPI_Allreduce(&snd_buff_0_size,&max_snd_buff_0_size,1,
                  MPI_INT,MPI_MAX,world);
    while(max_snd_buff_0_size)
    {
        int rcv_p=neigh_p[idim][1];
        int snd_p=neigh_p[idim][0];
        MPI_Sendrecv(&snd_buff_0_size,1,
                     MPI_INT,snd_p,0,&rcv_buff_size,
                     1,MPI_INT,rcv_p,0,world,
                     &status);
        
        buff_size_management(rcv_buff,rcv_buff_capacity,rcv_buff_size);
        
        if(rcv_buff_size)
            MPI_Irecv(rcv_buff,rcv_buff_size,
                      MPI_BYTE,rcv_p,0,world,
                      &request);
        if(snd_buff_0_size)
            MPI_Send(snd_buff_0,snd_buff_0_size,
                     MPI_BYTE,snd_p,0,world);
        
        snd_buff_0_size=0;
        if (rcv_buff_size)
        {
            MPI_Wait(&request,&status);
            
            no_atms=rcv_buff_size/list->byte_size;
            
            for(int i=0;i<no_atms;i++)
            {
                type0 s_i;
                memcpy(&s_i,&rcv_buff[i*list->byte_size+idim*d_size],d_size);
                if(s_i<s_lo[idim])
                {
                    if(snd_buff_0_capacity<snd_buff_0_size+list->byte_size)
                    {
                        GROW(snd_buff_0,snd_buff_0_capacity,snd_buff_0_size+list->byte_size);
                        snd_buff_0_capacity=snd_buff_0_size+list->byte_size;
                    }
                    memcpy(&(snd_buff_0[snd_buff_0_size]),&(rcv_buff[i*list->byte_size]),list->byte_size);
                    snd_buff_0_size+=list->byte_size;
                }
                else
                    unpack_cmplt(rcv_buff,i*list->byte_size,list);
            }
        }
        MPI_Allreduce(&snd_buff_0_size,&max_snd_buff_0_size,1,
                      MPI_INT,MPI_MAX,world);
    }
    
    int max_snd_buff_1_size;
    MPI_Allreduce(&snd_buff_1_size,&max_snd_buff_1_size,1,
                  MPI_INT,MPI_MAX,world);
    while(max_snd_buff_1_size)
    {
        int rcv_p=neigh_p[idim][0];
        int snd_p=neigh_p[idim][1];
        MPI_Sendrecv(&snd_buff_1_size,1,
                     MPI_INT,snd_p,0,&rcv_buff_size,
                     1,MPI_INT,rcv_p,0,world,
                     &status);
        
        buff_size_management(rcv_buff,rcv_buff_capacity,rcv_buff_size);

        if(rcv_buff_size)
            MPI_Irecv(rcv_buff,rcv_buff_size,
                      MPI_BYTE,rcv_p,0,world,
                      &request);
        if(snd_buff_1_size)
            MPI_Send(snd_buff_1,snd_buff_1_size,
                     MPI_BYTE,snd_p,0,world);
        
        snd_buff_1_size=0;
        if (rcv_buff_size)
        {
            MPI_Wait(&request,&status);
            
            no_atms=rcv_buff_size/list->byte_size;
            for(int i=0;i<no_atms;i++)
            {
                type0 s_i;
                memcpy(&s_i,&rcv_buff[i*list->byte_size+idim*d_size],d_size);
                if(s_hi[idim]<=s_i)
                {
                    if(snd_buff_1_capacity<snd_buff_1_size+list->byte_size)
                    {
                        GROW(snd_buff_1,snd_buff_1_capacity,snd_buff_1_size+list->byte_size);
                        snd_buff_1_capacity=snd_buff_1_size+list->byte_size;
                    }
                    memcpy(&snd_buff_1[snd_buff_1_size],&rcv_buff[i*list->byte_size],list->byte_size);
                    snd_buff_1_size+=list->byte_size;
                }
                else
                    unpack_cmplt(rcv_buff,i*list->byte_size,list);
            }
        }
        MPI_Allreduce(&snd_buff_1_size,&max_snd_buff_1_size,1,
                      MPI_INT,MPI_MAX,world);
    }
}
/*--------------------------------------------
 
 0. exchanges own atoms with neighboring
 processors so that each processor has their
 own proper atoms (limited to their local
 domain). unlike xchng_cmplt(class VecLst*)it 
 assumes that atoms have not moved more than 1 
 processor & performs only 1 communication in 
 each dimension and direction. hence the name
 partial. if an atom moves more than one 
 processors it produces an error and aborts 
 execution.
 
 caveat: it is assumed that x is already is
 converted to s, i.e. the coordinates are
 frational.
 --------------------------------------------*/
void Atoms::xchng_prtl(class VecLst* list)
{
    for(int idim=0;idim<dimension;idim++)
        xchng_prtl(idim,list);

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::xchng_prtl(int idim,class VecLst* list)
{
    if(tot_p_grid[idim]==1)
        return;
    int icomp;
    type0* s;
    vectors[0]->ret(s);
    int s_dim=vectors[0]->dim;
    int snd_buff_0_size=0;
    int snd_buff_1_size=0;
    int rcv_buff_size,rcv_p,snd_p;
    MPI_Request request;
    MPI_Status status;
    int no_atms;
    int iatm=0;
    
    while (iatm<natms)
    {
        icomp=s_dim*iatm;
        if(s[icomp+idim]<s_lo[idim])
        {

            if(my_loc[idim]!=tot_p_grid[idim]-1)
                pack_cmplt_prtl(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
            else
            {
                if(0.5<=s[icomp+idim])
                    pack_cmplt_prtl(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
                else
                    pack_cmplt_prtl(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
            }
            
        }
        else if(s[icomp+idim]>=s_hi[idim])
        {
            if(my_loc[idim]!=0)
                pack_cmplt_prtl(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
            else
            {
                if(s[icomp+idim]<0.5)
                    pack_cmplt_prtl(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
                else
                    pack_cmplt_prtl(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
            }
        }
        else iatm++;
    }
    

    rcv_p=neigh_p[idim][1];
    snd_p=neigh_p[idim][0];
    MPI_Sendrecv(&snd_buff_0_size,1,
                 MPI_INT,snd_p,0,&rcv_buff_size,
                 1,MPI_INT,rcv_p,0,world,
                 &status);
    
    buff_size_management(rcv_buff,rcv_buff_capacity,rcv_buff_size);

    if(rcv_buff_size)
        MPI_Irecv(rcv_buff,rcv_buff_size,
                  MPI_BYTE,rcv_p,0,world,
                  &request);
    
    if(snd_buff_0_size)
        MPI_Send(snd_buff_0,snd_buff_0_size,
                 MPI_BYTE,snd_p,0,world);
    snd_buff_0_size=0;
    if (rcv_buff_size)
    {
        
        MPI_Wait(&request,&status);
        no_atms=rcv_buff_size/list->byte_size;
        
        unpack_prtl(rcv_buff,no_atms,list);
        
        vectors[0]->ret(s);
        icomp=(natms-no_atms)*s_dim+idim;
        for(int i=0;i<no_atms;i++)
        {
            if(s[icomp]<s_lo[idim])
                error->abort("atom moved more than one processor");
            icomp+=s_dim;
        }
        
    }
    
    
    rcv_p=neigh_p[idim][0];
    snd_p=neigh_p[idim][1];
    MPI_Sendrecv(&snd_buff_1_size,1,
                 MPI_INT,snd_p,0,&rcv_buff_size,
                 1,MPI_INT,rcv_p,0,world,
                 &status);
    
    buff_size_management(rcv_buff,rcv_buff_capacity,rcv_buff_size);
    
    if(rcv_buff_size)
        MPI_Irecv(rcv_buff,rcv_buff_size,
                  MPI_BYTE,rcv_p,0,world,
                  &request);
    if(snd_buff_1_size)
        MPI_Send(snd_buff_1,snd_buff_1_size,
                 MPI_BYTE,snd_p,0,world);
    snd_buff_1_size=0;
    if (rcv_buff_size)
    {
        MPI_Wait(&request,&status);
        
        no_atms=rcv_buff_size/list->byte_size;
        
        unpack_prtl(rcv_buff,no_atms,list);
        
        vectors[0]->ret(s);
        icomp=(natms-no_atms)*s_dim+idim;
        for(int i=0;i<no_atms;i++)
        {
            if(s_hi[idim]<=s[icomp])
                error->abort("atom moved more than one processor");
            icomp+=s_dim;
        }

    }
 
}
/*--------------------------------------------
 updates the phantom atoms of a specific
 vector.
 --------------------------------------------*/
void Atoms::update_ph(int ivec)
{
    update_ph(&ivec,1,vectors[ivec]->byte_size);
}
/*--------------------------------------------
 updates the phantom atoms of a group of vectors
 --------------------------------------------*/
void Atoms::update_ph(int* vec_list,int no_vecs,int vec_byte_size)
{
    natms_ph=0;
    for(int iswap=0;iswap<ph_lst->no_swaps;iswap++)
        xchng_ph(iswap,vec_list,no_vecs,vec_byte_size);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::del(int ivec)
{
    if(ivec<0 || no_vecs-1<ivec)
    {
        error->abort("Error: index should be between 0 and no_vecs-1");
    }
    tot_byte_size-=vectors[ivec]->byte_size;
    
    Avec** new_vectors;
    new_vectors = new Avec*[no_vecs-1];
    for(int i=0;i<ivec;i++)
        new_vectors[i]=vectors[i];
    
    for(int i=ivec+1;i<no_vecs;i++)
        new_vectors[i-1]=vectors[i];
    
    delete [] vectors;
    vectors=new_vectors;
    
    no_vecs--;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Atoms::find(const char* name)
{
    for (int i=0;i<no_vecs;i++)
        if(strcmp(vectors[i]->name,name)==0)
            return i;
    
    error->abort("cannot find atomic vector %s",name);
    
    return -1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Atoms::find_exist(const char* name)
{
    for (int i=0;i<no_vecs;i++)
        if(strcmp(vectors[i]->name,name)==0)
            return i;
    return -1;
}
/*----------------------------------------------------------------------------------------------
  _____       ___   _____   _   _        __  _   _   __   _   _____       ___   _____   _   _
 |  _  \     /   | /  ___| | | / /      / / | | | | |  \ | | |  _  \     /   | /  ___| | | / /
 | |_| |    / /| | | |     | |/ /      / /  | | | | |   \| | | |_| |    / /| | | |     | |/ /
 |  ___/   / / | | | |     | |\ \     / /   | | | | | |\   | |  ___/   / / | | | |     | |\ \
 | |      / /  | | | |___  | | \ \   / /    | |_| | | | \  | | |      / /  | | | |___  | | \ \
 |_|     /_/   |_| \_____| |_|  \_\ /_/     \_____/ |_|  \_| |_|     /_/   |_| \_____| |_|  \_\
 
 ----------------------------------------------------------------------------------------------*/
/*--------------------------------------------
 note: I have changed the order of list of
 vecs and atoms; in order to make it more
 consistant with pack which does it atom by
 atom; it might have been ok up to this point
 I have unpacked the atoms one by one.
 
 also the growing seems a bit out of place; it
 seems that it also grows the unneccessary
 vectors; the vectors that are not in the
 VecLst need not to be grown they only need to
 be resized;
 
 I fixed the growth; however, keep an eye on
 that.
 
 these whole changes might cause a shit storm.
 
 used by: 
 xchng_prtl(int,class VecLst*)
 --------------------------------------------*/
inline void Atoms::unpack_prtl(char*& buff,int atm_list_size
,class VecLst* list)
{
    int buff_pos=0;
    
    if(natms+atm_list_size-atm_vec_size>0)
        grow(0,natms+atm_list_size-atm_vec_size,list);
    
    for(int i=0;i<list->no_vecs;i++)
        vectors[list->vec_list[i]]->unpack(buff,buff_pos,list->byte_size,atm_list_size);

    natms+=atm_list_size;
    
}
/*--------------------------------------------
 used by:
 xchng_cmplt(int,class VecLst*)
 --------------------------------------------*/
inline void Atoms::unpack_cmplt(char*& buff,int buff_pos,class VecLst* list)
{
    
    if(natms+1-atm_vec_size>0)
        grow(0,natms+1-atm_vec_size,list);
    
    for(int i=0;i<list->no_vecs;i++)
        vectors[list->vec_list[i]]->unpack(buff,buff_pos);
    
    natms++;
}
/*--------------------------------------------
 this fucntion packs individual atoms from the
 vectors in VecLst; it also replaces the said
 atom with the last one
 
 used by:
 xchng_cmplt(int,class VecLst*)
 xchng_prtl(int,class VecLst*)
 --------------------------------------------*/
inline void Atoms::pack_cmplt_prtl(char*& buff,int& buff_pos,int& buff_capacity
,int iatm,class VecLst* list)
{
    
    if(buff_pos+list->byte_size>buff_capacity)
    {
        GROW(buff,buff_capacity,buff_pos+list->byte_size);
        buff_capacity=buff_pos+list->byte_size;
    }
    
    
    for(int i=0;i<list->no_vecs;i++)
        vectors[list->vec_list[i]]->pack_n_mv(buff,buff_pos,iatm);
    
    natms--;
    
}
/*--------------------------------------------
 this function packs a list of atoms from the
 atomic vectors inside VecLst (from phantom
 vectors).
 
 used by:
 xchng_ph(int,int,int*,int,class VecLst*)
 --------------------------------------------*/
inline void Atoms::pack_ph(char*& buff,class VecLst* list
,int* atm_list,int atm_list_size)
{
    int buff_pos=0;
    
    for(int i=0;i<list->ph_no_vecs;i++)
        vectors[list->ph_vec_list[i]]->pack(buff,buff_pos,atm_list,atm_list_size);
    
}
/*--------------------------------------------
 
 used by:
 update_ph(int*,int,int)
 --------------------------------------------*/
inline void Atoms::pack_ph(char*& buff,int* vec_list,int no_vecs
,int* atm_list,int atm_list_size)
{
    int buff_pos=0;
    
    for(int i=0;i<no_vecs;i++)
        vectors[vec_list[i]]->pack(buff,buff_pos,atm_list,atm_list_size);
    
}
/*--------------------------------------------
 this function unpacks the phantom atoms
 gathered by pack_ph()
 
 used by:
 xchng_ph(int,int,int*,int,class VecLst*)
 --------------------------------------------*/
inline void Atoms::unpack_ph(char*& buff
,class VecLst* list,int atm_list_size)
{
    
    if(natms_ph+natms+atm_list_size-atm_vec_ph_size>0)
    {
        grow(1,natms_ph+natms+atm_list_size-atm_vec_ph_size,list);
        atm_vec_ph_size=natms_ph+natms+atm_list_size;
    }
    
    int buff_pos=0;

    for(int i=0;i<list->ph_no_vecs;i++)
        vectors[list->ph_vec_list[i]]->unpack_ph(buff,buff_pos,atm_list_size);
    
    natms_ph+=atm_list_size;
    
}
/*--------------------------------------------
 
 used by:
 update_ph(int*,int,int)
 --------------------------------------------*/
inline void Atoms::unpack_ph(char*& buff,int* vec_list
,int no_vecs,int atm_list_size)
{
    int buff_pos=0;
    
    for(int i=0;i<no_vecs;i++)
        vectors[vec_list[i]]->unpack_ph(buff,buff_pos,atm_list_size);
    
    natms_ph+=atm_list_size;
}
/*--------------------------------------------
 
 unpack from outside for read
 --------------------------------------------*/
void Atoms::unpack_read(char*& buff,int atm_list_size
,class VecLst* list)
{
    unpack_prtl(buff,atm_list_size,list);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::store_0()
{
    type0* x;
    vectors[0]->ret(x);
    type0* x_0;
    vectors[1]->ret(x_0);
    
    int x_dim=vectors[0]->dim;
    int x_0_dim=vectors[1]->dim;
    int x_comp=0;
    int x_0_comp=0;
    for(int i=0;i<natms;i++)
    {
        for(int idim=0;idim<dimension;idim++)
            x_0[x_0_comp+idim]=x[x_comp+idim];
        x_comp+=x_dim;
        x_0_comp+=x_0_dim;
    }
}
/*--------------------------------------------
 please note that if the box's size changes,
 this function assumes the you have already 
 inversed the H matrix and calculated the B
 matrix; if you have truble inersing your H
 matrix, you can use the invert(double**,
 double**,int). however, please be advised 
 this procedure is numerical and you need to
 adjust TOLERANCE value.
 --------------------------------------------*/
void Atoms::update(int box_change,class VecLst* list)
{
    
    type0* x;
    vectors[0]->ret(x);
    type0* x_0;
    vectors[1]->ret(x_0);
    
    int x_dim=vectors[0]->dim;
    int x_0_dim=vectors[1]->dim;


    type0 dx;
    type0 sq_dx;
    
    int iatm=0;
    int check=0;
    type0 sq_skin=skin*skin;
    
    int x_comp=0;
    int x_0_comp=0;
    while (iatm<natms&& check==0)
    {
        
        sq_dx=0.0;
        for(int idim=0;idim<dimension;idim++)
        {
            dx=x[x_comp+idim]-x_0[x_0_comp+idim];
            sq_dx+=dx*dx;
        }
        if(sq_dx>sq_skin)
            check=1;
        
        iatm++;
        x_comp+=x_dim;
        x_0_comp+=x_0_dim;
    }
    
    int check_all;
    MPI_Allreduce(&check,&check_all,1,
    MPI_INT,MPI_MAX,world);

    
    if(check_all)
    {
        x2s(natms);
        xchng_prtl(list);
        setup_ph(box_change,list);
        neighbor->create_list(box_change,1);
        /*no need to convert s back to x,
         the neighbor list takes care of it*/

        x_comp=0;
        x_0_comp=0;
 
        vectors[0]->ret(x);
        vectors[1]->ret(x_0);
        for(int i=0;i<natms;i++)
        {
            for(int idim=0;idim<dimension;idim++)
                x_0[x_0_comp+idim]=x[x_comp+idim];
            x_comp+=x_dim;
            x_0_comp+=x_0_dim;
        }
    }
    else
    {
        

        update_ph(list->update_every_ph_vec_list
               ,list->update_every_ph_no_vecs
               ,list->update_every_ph_byte_size);
        

    }
}
/*--------------------------------------------
 transform x 2 s
 --------------------------------------------*/
void Atoms::x2s(int no)
{
    type0* x;
    vectors[0]->ret(x);
    int x_dim=vectors[0]->dim;
    int icomp=0;
    for(int i=0;i<no;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            x[icomp+j]=x[icomp+j]*B[j][j];
            for(int k=j+1;k<dimension;k++)
                x[icomp+j]+=B[k][j]*x[icomp+k];
            
        }
        
        for(int j=0;j<dimension;j++)
        {
            while(x[icomp+j]<0.0)
                x[icomp+j]++;
            while(x[icomp+j]>=1.0)
                x[icomp+j]--;
            
        }
               
        icomp+=x_dim;
    }
}
/*--------------------------------------------
 transform s 2 x
 --------------------------------------------*/
void Atoms::s2x(int no)
{
    type0* s;
    vectors[0]->ret(s);
    int s_dim=vectors[0]->dim;
    int icomp=0;
    for(int i=0;i<no;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            s[icomp+j]=s[icomp+j]*H[j][j];
            for(int k=j+1;k<dimension;k++)
                s[icomp+j]+=H[k][j]*s[icomp+k];
        }
        icomp+=s_dim;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::grow(int ph,int no_atms
,class VecLst* list)
{
    int* vec_list;
    int vec_list_size;
    
    if(ph)
    {
        vec_list=list->ph_vec_list;
        vec_list_size=list->ph_no_vecs;
    }
    else
    {
        vec_list=list->vec_list;
        vec_list_size=list->no_vecs;
    }
    
    int ivec=0;
    for(int i=0;i<no_vecs;i++)
    {
        if(i==vec_list[ivec])
        {
            vectors[vec_list[ivec]]->grow(no_atms);
            ivec++;
        }
        else
        {
            if(ph)
            {
                if(vectors[i]->ph)
                    vectors[i]->resize(no_atms);
            }
            else
                vectors[i]->resize(no_atms);
        }
    }
    
    if(ph)
        atm_vec_ph_size+=no_atms;
    else
    {
        atm_vec_size+=no_atms;
        atm_vec_ph_size+=no_atms;
    }
    

}

/*--------------------------------------------
 autogrid the domain
 --------------------------------------------*/
void Atoms::auto_grid_proc()
{
    
    int eq_p_per_n=1;
    for(int i=0;i<tot_n;i++)
        if(p_per_n[i]!=p_per_n[0])
            eq_p_per_n=0;

    type0 vol;
    type0* area;
    CREATE1D(area,dimension);
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
    
    XMath* xmath= new XMath(mapp);
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
            int* tmp_tot_p_grid;
            CREATE1D(tmp_tot_p_grid,dimension-1);
            
            int* dim_indx;
            CREATE1D(dim_indx,dimension-1);
            
            int* my_tmp_loc;
            CREATE1D(my_tmp_loc,dimension-1);
            
            int* my_tmp_denom;
            CREATE1D(my_tmp_denom,dimension-1);
            
            
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
            for(int i=0;i<p_per_n[my_n_no];i++)
                if(n_p_grid[my_n_no][i]==my_p_no)
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
            my_loc[prin_dimension]=my_n_no;
            
            
            
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
                
                neigh_p[dim_indx[i]][1]=n_p_grid[my_n_no][nop];
                neigh_p[dim_indx[i]][0]=n_p_grid[my_n_no][nom];
            }
            
            tmp0=my_n_no-1;
            if(tmp0==-1)
                tmp0=tot_n-1;
            neigh_p[prin_dimension][0]=n_p_grid[tmp0][my_p_in_my_node];
            
            tmp1=my_n_no+1;
            if(tmp1==tot_n)
                tmp1=0;
            neigh_p[prin_dimension][1]=n_p_grid[tmp1][my_p_in_my_node];
            
            
            delete [] dim_indx;
            delete [] tmp_tot_p_grid;
            delete [] my_tmp_loc;
            delete [] my_tmp_denom;
            delete [] area;
            
            
        }
        else
        {
            tot_p_grid[0]=tot_p;
            
            int my_p_in_my_node=-1;
            for(int i=0;i<p_per_n[my_n_no];i++)
                if(n_p_grid[my_n_no][i]==my_p_no)
                    my_p_in_my_node=i;
            
            
            my_loc[0]=p_per_n[my_n_no]*my_n_no;
            my_loc[0]+=my_p_in_my_node;
            
            int tmp0p,tmp1p,tmp0n,tmp1n;
            tmp0p=my_p_in_my_node;
            tmp1p=my_p_in_my_node;
            tmp0n=my_n_no;
            tmp1n=my_n_no;
            
            tmp0p--;
            if(tmp0p==-1)
            {
                tmp0p=p_per_n[my_n_no]-1;
                tmp0n--;
                if(tmp0n==-1)
                    tmp0n=tot_n-1;
            }
            
            tmp1p++;
            if(tmp1p==p_per_n[my_n_no])
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
        
        delete [] area;
        
        for(int i=0;i<fac_list_size;i++)
            delete [] fac_list[i];
        delete [] fac_list;
        
        int* list;
        CREATE1D(list,dimension);
        
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
    

    
    for(int i=0;i<dimension;i++)
    {
        s_lo[i]=
        static_cast<type0>(my_loc[i])
        /static_cast<type0>(tot_p_grid[i]);
        
        s_hi[i]=
        (static_cast<type0>(my_loc[i])+1.0)
        /static_cast<type0>(tot_p_grid[i]);
    }
    if(my_p_no==0)
    {
        fprintf(output,"autogrid performed: ");
        
        for(int i=0;i<dimension-1;i++)
            fprintf(output,"%d\u00D7",tot_p_grid[i]);
        
        fprintf(output,"%d",tot_p_grid[dimension-1]);
        
        fprintf(output,"\n");
    }
    
    
}
/*--------------------------------------------
 manually grid the domain
 --------------------------------------------*/
void Atoms::man_grid_proc(int narg,char** args)
{
    if(narg!=dimension+1)
        error->abort("grid should at least"
        " have %d arguement",dimension);
    
    for(int i=0;i<dimension;i++)
        tot_p_grid[i]=atoi(args[i+1]);

    int tmp=1;
    for(int i=0;i<dimension;i++)
    {
        if(tot_p_grid[i]<1)
            error->abort("number of processors "
            "in dimension %d for grid cannot be "
            "less than 1",i);
        tmp*=tot_p_grid[i];
    }
    if(tmp!=tot_p)
    {
        /*
        char err_msg[MAXCHAR];
        sprintf(err_msg,"%d\u00D7",tot_p_grid[0]);
        for(int i=1;i<dimension-1;i++)
            sprintf(err_msg,"%s%d\u00D7",err_msg,tot_p_grid[i]);
        sprintf(err_msg,"%s%d",err_msg,tot_p_grid[dimension-1]);
        
        error->abort("for grid %d\u2260%s",tot_p,err_msg);
         */
        error->abort("for grid total number of "
        "processors should be equal to the "
        "product of arguments");
    }

    
        
    int* list;
    CREATE1D(list,dimension);
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
    
    if(no_vecs)
    {
        CREATE1D(list,no_vecs);
        for(int i=0;i<no_vecs;i++)
            list[i]=i;
        VecLst* vecs_comm;
        vecs_comm=new VecLst(mapp,list,no_vecs);
        
        xchng_cmplt(vecs_comm);
        
        delete vecs_comm;
        if(no_vecs)
            delete [] list;
    }
    
    
    if(my_p_no==0)
    {

        fprintf(output,"manual grid performed: ");
        
        for(int i=0;i<dimension-1;i++)
            fprintf(output,"%d\u00D7",tot_p_grid[i]);
        
        fprintf(output,"%d",tot_p_grid[dimension-1]);
        
        fprintf(output,"\n");
    }
    
    
}
/*--------------------------------------------
 add skin command
 --------------------------------------------*/
void Atoms::chng_skin(int narg,char** args)
{
    if(narg!=2)
        error->abort("wrong number of inputs for skin command");
    type0 s=atof(args[1]);
    if(s<=0.0)
        error->abort("skin cannot be equal or less than zero");
    skin=s;
}
/*--------------------------------------------
 buffer size management
 --------------------------------------------*/
inline void Atoms::buff_size_management(char*& buff
,int& curr_size,int size_needed)
{
    if(size_needed<curr_size)
        return;
    if(curr_size)
        delete [] buff;
    CREATE1D(buff,size_needed);
    curr_size=size_needed;
}
/*--------------------------------------------
 initialize every necessary thing before a run
 
 caveat: it is assumed that x is already is
 converted to s, i.e. the coordinates are
 frational. neighbor->create_list(0,1) will 
 change fractional coordinates back to normal
 --------------------------------------------*/
void Atoms::init(class VecLst* vecs_comm)
{
    x2s(natms);
    
    /*
     0. reset to make sure that each processor
     has the atoms that really belongs to
     it.
     */
    xchng_cmplt(vecs_comm);
    /*--------------------------------------*/
    
    /*
     1. setup the cutoffs for neighbor lists
     and also the maximum for communication
     in atoms class (this is needed for
     finding the phantom atoms and updating
     them);
     */
    
    /*
     1.1. assign the cut_sk_sq (for neighbor list)
     1.2. find the maximum cuttf off and
     atoms->set_max_cutoff(max_cutoff);
     1.3. set_max_cutoff(max_cutoff) finds boundary of of
     domain for phantom atoms (s)
     */
    forcefield->init();
    
    /*
     additional part to remove the skin from
     ff_*.cpp
     ff_*.cpp must only apply the cut_sq[]
     we will setup cut_sk_sq[] here
     */
    int arr_size=(atom_types->no_types)
    *(atom_types->no_types+1)/2;
    type0* cut_sk_sq=forcefield->cut_sk_sq;
    type0* cut_sq=forcefield->cut_sq;
    type0 tmp;
    type0 max_cut=0.0;
    for(int i=0;i<arr_size;i++)
    {
        tmp=sqrt(cut_sq[i])+skin;
        cut_sk_sq[i]=tmp*tmp;
        max_cut=MAX(max_cut,sqrt(cut_sq[i]));
    }
    set_max_cutoff(max_cut);
    
    
    /*--------------------------------------*/
    
    /*
     2. find the phantom atoms and
     communicate them between the processors
     using the values assigned by the
     previous step;
     */
    setup_ph(1,vecs_comm);
    /*--------------------------------------*/
    
    /*
     3. creates the neighbor list for the
     bins;
     */
    neighbor->init();
    /*--------------------------------------*/
    /*
     4. creates the neighbor list
     */
    neighbor->create_list(0,1);
    /*--------------------------------------*/
    /*
     5. store the first projections
     */
    store_0();
    /*--------------------------------------*/
}
/*---------------------------------------------------------------------------
  _     _   _____   _____   _       _____   _____
 | |   / / | ____| /  ___| | |     /  ___/ |_   _|
 | |  / /  | |__   | |     | |     | |___    | |
 | | / /   |  __|  | |     | |     \___  \   | |
 | |/ /    | |___  | |___  | |___   ___| |   | |
 |___/     |_____| \_____| |_____| /_____/   |_|
 ---------------------------------------------------------------------------*/
/*--------------------------------------------
 
 --------------------------------------------*/
VecLst::VecLst(MAPP* mapp,int no_args,...)
:InitPtrs(mapp)
{
    va_list args;
    va_start (args,no_args);
    
    no_vecs=no_args;
    CREATE1D(vec_list,no_vecs);
    
    for(int i=0;i<no_args;i++)
    {
        vec_list[i]=va_arg(args,int);
        if(vec_list[i]<0 || vec_list[i] >(atoms->no_vecs-1))
        {
            error->abort("wrong index for the vectors: %i",vec_list[i]);
        }
        
    }
    
    va_end (args);
    
    int tmp;
    for(int i=0;i<no_vecs;i++)
    {
        for(int j=i+1;j<no_vecs;j++)
        {
            if(vec_list[i]==vec_list[j])
                error->abort("duplicate atomic vecotr index in vector list");
            
            if(vec_list[i]>vec_list[j])
            {
                tmp=vec_list[i];
                vec_list[i]=vec_list[j];
                vec_list[j]=tmp;
            }
        }
    }
    
    byte_size=0;
    for(int i=0;i<no_vecs;i++)
        byte_size+=atoms->vectors[vec_list[i]]->byte_size;
    
    ph_no_vecs=0;
    for(int i=0;i<no_vecs;i++)
        if(atoms->vectors[vec_list[i]]->ph)
            ph_no_vecs++;
    
    CREATE1D(ph_vec_list,ph_no_vecs);
    
    int pos=0;
    ph_byte_size=0;
    for(int i=0;i<no_vecs;i++)
        if(atoms->vectors[vec_list[i]]->ph)
        {
            ph_vec_list[pos++]=vec_list[i];
            ph_byte_size+=atoms->vectors[vec_list[i]]->byte_size;
        }
    
    update_every_ph_no_vecs=0;
    update_every_ph_byte_size=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
VecLst::VecLst(MAPP* mapp,int* list,int no)
:InitPtrs(mapp)
{
    
    no_vecs=no;
    CREATE1D(vec_list,no_vecs);
    for(int i=0;i<no_vecs;i++)
        vec_list[i]=list[i];
    
    int tmp;
    for(int i=0;i<no_vecs;i++)
    {
        for(int j=i+1;j<no_vecs;j++)
        {
            if(vec_list[i]==vec_list[j])
                error->abort("duplicate atomic vecotr index in vector list");
            
            if(vec_list[i]>vec_list[j])
            {
                tmp=vec_list[i];
                vec_list[i]=vec_list[j];
                vec_list[j]=tmp;
            }
        }
    }
    
    byte_size=0;
    for(int i=0;i<no_vecs;i++)
        byte_size+=atoms->vectors[vec_list[i]]->byte_size;
    
    ph_no_vecs=0;
    for(int i=0;i<no_vecs;i++)
        if(atoms->vectors[vec_list[i]]->ph)
            ph_no_vecs++;
    
    CREATE1D(ph_vec_list,ph_no_vecs);
    
    int pos=0;
    ph_byte_size=0;
    for(int i=0;i<no_vecs;i++)
        if(atoms->vectors[vec_list[i]]->ph)
        {
            ph_vec_list[pos++]=vec_list[i];
            ph_byte_size+=atoms->vectors[vec_list[i]]->byte_size;
        }
    
    update_every_ph_no_vecs=0;
    update_every_ph_byte_size=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
VecLst::~VecLst()
{
    if(no_vecs)
        delete [] vec_list;
    if(ph_no_vecs)
        delete [] ph_vec_list;
    if(update_every_ph_no_vecs)
        delete [] update_every_ph_vec_list;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::add_update(int vec_no)
{
    int chk=0;
    for(int i=0;i<ph_no_vecs;i++)
        if(ph_vec_list[i]==vec_no)
            chk=1;
    
    if(chk==0)
        error->abort("update atomic vector should "
        "be of phantom kind");
    
    
    for(int i=0;i<update_every_ph_no_vecs;i++)
        if(update_every_ph_vec_list[i]==vec_no)
            error->abort("duplicate update atomic vector");
    
    
    GROW(update_every_ph_vec_list
    ,update_every_ph_no_vecs
    ,update_every_ph_no_vecs+1);
    update_every_ph_vec_list[update_every_ph_no_vecs]=vec_no;
    update_every_ph_no_vecs++;
    
    
    int tmp;
    for(int i=0;i<update_every_ph_no_vecs;i++)
    {
        for(int j=i+1;j<update_every_ph_no_vecs;j++)
        {
            if(update_every_ph_vec_list[i]==update_every_ph_vec_list[j])
                error->abort("duplicate atomic vecotr index in vector list");
            
            if(update_every_ph_vec_list[i]>update_every_ph_vec_list[j])
            {
                tmp=update_every_ph_vec_list[i];
                update_every_ph_vec_list[i]=update_every_ph_vec_list[j];
                update_every_ph_vec_list[j]=tmp;
            }
        }
    }
    
    update_every_ph_byte_size+=atoms->vectors[vec_no]->byte_size;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::del_update(int vec_no)
{
    int chk=0;
    for(int i=0;i<ph_no_vecs;i++)
        if(ph_vec_list[i]==vec_no)
            chk=1;
    
    if(chk==0)
        error->abort("update atomic vector not found %i",vec_no);

    int* vec_tmp;
    CREATE1D(vec_tmp,update_every_ph_no_vecs-1);
    
    int k=0;
    for(int i=0;i<update_every_ph_no_vecs;i++)
    {
        if(ph_vec_list[i]!=vec_no)
        {
            vec_tmp[k]=ph_vec_list[i];
            k++;
        }
    }
    
    delete [] ph_vec_list;
    
    ph_vec_list=vec_tmp;
    
    update_every_ph_byte_size-=atoms->vectors[vec_no]->byte_size;
    update_every_ph_no_vecs--;
}
/*---------------------------------------------------------------------------
  _____   _          __      ___   _____   _       _____   _____
 /  ___/ | |        / /     /   | |  _  \ | |     /  ___/ |_   _|
 | |___  | |  __   / /     / /| | | |_| | | |     | |___    | |
 \___  \ | | /  | / /     / / | | |  ___/ | |     \___  \   | |
  ___| | | |/   |/ /     / /  | | | |     | |___   ___| |   | |
 /_____/ |___/|___/     /_/   |_| |_|     |_____| /_____/   |_|
 ---------------------------------------------------------------------------*/

/*--------------------------------------------
constructor
--------------------------------------------*/
SwapLst::SwapLst(MAPP* mapp):InitPtrs(mapp)
{
    grow_size=SWAPGROWTH;
    no_swaps=0;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
SwapLst::SwapLst(MAPP* mapp,int nswaps):InitPtrs(mapp)
{
    no_swaps=nswaps;
    CREATE1D(snd_list,no_swaps);
    
    CREATE1D(snd_list_capacity,no_swaps);
    
    CREATE1D(snd_size,no_swaps);
    CREATE1D(rcv_size,no_swaps);
    
    CREATE1D(snd_p,no_swaps);
    CREATE1D(rcv_p,no_swaps);
    
    CREATE1D(pbc_correction,no_swaps);
    
    grow_size=SWAPGROWTH;
    
    for(int i=0;i<no_swaps;i++)
        snd_list_capacity[i]=snd_size[i]=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
SwapLst::~SwapLst()
{
    if(no_swaps)
    {
        for(int i=0;i<no_swaps;i++)
            if(snd_list_capacity[i])
                delete [] snd_list[i];
        delete [] snd_list;
       
        delete [] snd_list_capacity;
        
        delete [] snd_size;
        delete [] rcv_size;
        
        delete [] snd_p;
        delete [] rcv_p;
        
        delete [] pbc_correction;
    }
    
}
/*--------------------------------------------
 create a new list with new number of swaps
 --------------------------------------------*/
void SwapLst::newlist(int nswaps)
{
    if(nswaps!= no_swaps)
    {
        if(no_swaps)
        {
            for(int i=0;i<no_swaps;i++)
                if(snd_list_capacity[i])
                    delete [] snd_list[i];
            delete [] snd_list;
            
            delete [] snd_list_capacity;
            
            delete [] snd_size;
            delete [] rcv_size;
            
            delete [] snd_p;
            delete [] rcv_p;
            
            delete [] pbc_correction;
        }
        
        no_swaps=nswaps;
        CREATE1D(snd_list,no_swaps);
        
        CREATE1D(snd_list_capacity,no_swaps);
        
        CREATE1D(snd_size,no_swaps);
        CREATE1D(rcv_size,no_swaps);
        
        CREATE1D(snd_p,no_swaps);
        CREATE1D(rcv_p,no_swaps);
        
        CREATE1D(pbc_correction,no_swaps);
    }
    
    grow_size=SWAPGROWTH;
    
    for(int i=0;i<no_swaps;i++)
        snd_list_capacity[i]=snd_size[i]=0;
    
}
/*--------------------------------------------
 add element to a swap
 --------------------------------------------*/
void SwapLst::add(int& iswap,int& val)
{
    if(snd_size[iswap]==snd_list_capacity[iswap])
    {
        GROW(snd_list[iswap],snd_list_capacity[iswap],snd_list_capacity[iswap]+grow_size);
        snd_list_capacity[iswap]+=grow_size;
    }
    snd_list[iswap][snd_size[iswap]]=val;
    snd_size[iswap]++;
}
/*--------------------------------------------
 bring the cursor to the very beginning of all
 of the lists
 --------------------------------------------*/
void SwapLst::reset()
{
    for(int iswap=0;iswap<no_swaps;iswap++)
        snd_size[iswap]=0;
}
/*--------------------------------------------
 bring the cursor to the very beginning of one
 list
 --------------------------------------------*/
void SwapLst::reset(int iswap)
{
    snd_size[iswap]=0;
}


/*----------------------------------------------------------------------------------------------------------
  _____   _____        _____   _____        _____    _____       ___  ___   _____   _     _   _____   _____
 |_   _| /  _  \      |  _  \ | ____|      |  _  \  | ____|     /   |/   | /  _  \ | |   / / | ____| |  _  \
   | |   | | | |      | |_| | | |__        | |_| |  | |__      / /|   /| | | | | | | |  / /  | |__   | | | |
   | |   | | | |      |  _  { |  __|       |  _  /  |  __|    / / |__/ | | | | | | | | / /   |  __|  | | | |
   | |   | |_| |      | |_| | | |___       | | \ \  | |___   / /       | | | |_| | | |/ /    | |___  | |_| |
   |_|   \_____/      |_____/ |_____|      |_|  \_\ |_____| /_/        |_| \_____/ |___/     |_____| |_____/
 ----------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------
 updates the phantom atoms of a specific
 vector.
 --------------------------------------------*/
/*
void Atoms::update_ph_old(int ivec)
{
    update_ph(&ivec,1,vectors[ivec]->byte_size);
}
 */
/*--------------------------------------------
 updates the phantom atoms of a group of vectors
 --------------------------------------------*/
/*
void Atoms::update_ph_old(int* vec_list,int no_vecs,int vec_byte_size)
{
    int snd_proc,rcv_proc;
    int snd_buff_size,rcv_buff_size;
    int last_atm_no,no_new_atms;
    int iswap,idim,idir,icomm,ivec,iatm,jdim;
    
    MPI_Request request[2];
    MPI_Status status[2];
    
    natms_ph=0;
    iswap=0;
    if(no_vecs==1)
    {
        for(idim=0;idim<dimension;idim++)
        {
            if(tot_p_grid[idim]!=1)
            {
                for(idir=0;idir<2;idir++)
                {
                    for(icomm=0;icomm<comm_need[idim][idir];icomm++)
                    {
                        snd_proc=neigh_p[idim][idir];
                        rcv_proc=neigh_p[idim][1-idir];
                        snd_buff_size=vec_byte_size*ph_lst->snd_size[iswap];
                        
                        buff_size_management(snd_ph_buff,snd_ph_buff_capacity,snd_buff_size);
                        
                        pack_ph(snd_ph_buff,vec_list,no_vecs,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
                        
                        
                        
                        MPI_Sendrecv(&snd_buff_size,1,MPI_INT,snd_proc,0,&rcv_buff_size,
                                     1,MPI_INT,rcv_proc,0,world,&status[0]);
                        
                        if (rcv_buff_size)
                        {
                            MPI_Irecv(vectors[vec_list[0]]->ret(natms+natms_ph)
                                      ,rcv_buff_size,MPI_BYTE,rcv_proc,0,world,&request[0]);
                            
                        }
                        if (snd_buff_size)
                        {
                            MPI_Isend(snd_ph_buff,snd_buff_size
                                      ,MPI_BYTE,snd_proc,0,world,&request[1]);
                        }
                        
                        no_new_atms=rcv_buff_size/vec_byte_size;
                        
                        if (rcv_buff_size&&snd_buff_size==0)
                            MPI_Wait(&request[0],&status[0]);
                        else if (rcv_buff_size==0&&snd_buff_size)
                            MPI_Wait(&request[1],&status[1]);
                        else if (rcv_buff_size&&snd_buff_size)
                            MPI_Waitall(2,request,status);
                        
                        natms_ph+=no_new_atms;
                        if(vec_list[0]==0)
                        {
                            if((idir==0&&my_loc[idim]==tot_p_grid[idim]-1) ||
                               (idir==1&&my_loc[idim]==0))
                            {
                                last_atm_no=natms+natms_ph;
                                int strt_atm_no=last_atm_no-no_new_atms;
                                type0* x;
                                vectors[0]->ret(x);
                                int x_dim=vectors[0]->dim;
                                
                                if(idir==0)
                                    for (iatm=strt_atm_no;iatm<last_atm_no;iatm++)
                                        for(jdim=0;jdim<=idim;jdim++)
                                            x[iatm*x_dim+jdim]+=H[idim][jdim];
                                else if(idir==1)
                                    for (iatm=strt_atm_no;iatm<last_atm_no;iatm++)
                                        for(jdim=0;jdim<=idim;jdim++)
                                            x[iatm*x_dim+jdim]-=H[idim][jdim];
                            }
                        }
                        iswap++;
                    }
                }
            }
            else
            {
                for(idir=0;idir<2;idir++)
                {
                    for(icomm=0;icomm<comm_need[idim][idir];icomm++)
                    {
                        last_atm_no=natms+natms_ph;
                        
                        for(ivec=0;ivec<no_vecs;ivec++)
                            for(iatm=0;iatm<ph_lst->snd_size[iswap];iatm++)
                                vectors[vec_list[ivec]]->copy(last_atm_no+iatm,ph_lst->snd_list[iswap][iatm]);
                        
                        if(vec_list[0]==0)
                        {
                            
                            type0* x;
                            vectors[0]->ret(x);
                            int x_dim=vectors[0]->dim;
                            
                            if(idir==0)
                                for (iatm=natms+natms_ph;iatm<natms+natms_ph+ph_lst->snd_size[iswap];iatm++)
                                    for(jdim=0;jdim<=idim;jdim++)
                                        x[iatm*x_dim+jdim]+=H[idim][jdim];
                            
                            else if(idir==1)
                                for (iatm=natms+natms_ph;iatm<natms+natms_ph+ph_lst->snd_size[iswap];iatm++)
                                    for(jdim=0;jdim<=idim;jdim++)
                                        x[iatm*x_dim+jdim]-=H[idim][jdim];
                            
                        }
                        natms_ph+=ph_lst->snd_size[iswap];
                        iswap++;
                    }
                }
            }
            
        }
    }
    else
    {
        for(idim=0;idim<dimension;idim++)
        {
            if(tot_p_grid[idim]!=1)
            {
                for(idir=0;idir<2;idir++)
                {
                    for(icomm=0;icomm<comm_need[idim][idir];icomm++)
                    {
                        snd_proc=neigh_p[idim][idir];
                        rcv_proc=neigh_p[idim][1-idir];
                        snd_buff_size=vec_byte_size*ph_lst->snd_size[iswap];
                        
                        buff_size_management(snd_ph_buff,snd_ph_buff_capacity,snd_buff_size);
                        
                        pack_ph(snd_ph_buff,vec_list,no_vecs,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap]);
                        
                        MPI_Sendrecv(&snd_buff_size,1,
                                     MPI_INT,snd_proc,0,&rcv_buff_size,
                                     1,MPI_INT,rcv_proc,0,world,&status[0]);
                        
                        buff_size_management(rcv_ph_buff,rcv_ph_buff_capacity,rcv_buff_size);
                        
                        if (rcv_buff_size)
                        {
                            MPI_Irecv(rcv_ph_buff,rcv_buff_size,
                                      MPI_BYTE,rcv_proc,0,world,
                                      &request[0]);
                            
                        }
                        if (snd_buff_size)
                        {
                            MPI_Isend(snd_ph_buff,snd_buff_size,
                                      MPI_BYTE,snd_proc,0,world,
                                      &request[1]);
                        }
                        
                        no_new_atms=rcv_buff_size/vec_byte_size;
                        
                        if (rcv_buff_size&&snd_buff_size==0)
                        {
                            MPI_Wait(&request[0],&status[0]);
                            unpack_ph(rcv_ph_buff,vec_list,no_vecs,no_new_atms);
                            
                        }
                        else if (rcv_buff_size==0&&snd_buff_size)
                        {
                            MPI_Wait(&request[1],&status[1]);
                        }
                        else if (rcv_buff_size&&snd_buff_size)
                        {
                            MPI_Waitall(2,request,status);
                            unpack_ph(rcv_ph_buff,vec_list,no_vecs,no_new_atms);
                        }
                        
                        if(vec_list[0]==0)
                        {
                            if((idir==0&&my_loc[idim]==tot_p_grid[idim]-1) ||
                               (idir==1&&my_loc[idim]==0))
                            {
                                last_atm_no=natms+natms_ph;
                                int strt_atm_no=last_atm_no-no_new_atms;
                                type0* x;
                                vectors[0]->ret(x);
                                int x_dim=vectors[0]->dim;
                                
                                if(idir==0)
                                    for (iatm=strt_atm_no;iatm<last_atm_no;iatm++)
                                        for(jdim=0;jdim<=idim;jdim++)
                                            x[iatm*x_dim+jdim]+=H[idim][jdim];
                                else if(idir==1)
                                    for (iatm=strt_atm_no;iatm<last_atm_no;iatm++)
                                        for(jdim=0;jdim<=idim;jdim++)
                                            x[iatm*x_dim+jdim]-=H[idim][jdim];
                            }
                        }
                        iswap++;
                    }
                }
            }
            else
            {
                for(idir=0;idir<2;idir++)
                {
                    for(icomm=0;icomm<comm_need[idim][idir];icomm++)
                    {
                        last_atm_no=natms+natms_ph;
                        
                        for(ivec=0;ivec<no_vecs;ivec++)
                            for(iatm=0;iatm<ph_lst->snd_size[iswap];iatm++)
                                vectors[vec_list[ivec]]->copy(last_atm_no+iatm,ph_lst->snd_list[iswap][iatm]);
                        
                        if(vec_list[0]==0)
                        {
                            type0* x;
                            vectors[0]->ret(x);
                            int x_dim=vectors[0]->dim;
                            
                            if(idir==0)
                                for (iatm=natms+natms_ph;iatm<natms+natms_ph+ph_lst->snd_size[iswap];iatm++)
                                    for(jdim=0;jdim<=idim;jdim++)
                                        x[iatm*x_dim+jdim]+=H[idim][jdim];
                            else if(idir==1)
                                for (iatm=natms+natms_ph;iatm<natms+natms_ph+ph_lst->snd_size[iswap];iatm++)
                                    for(jdim=0;jdim<=idim;jdim++)
                                        x[iatm*x_dim+jdim]-=H[idim][jdim];
                            
                        }
                        natms_ph+=ph_lst->snd_size[iswap];
                        iswap++;
                    }
                }
            }
        }
    }
    
}
 */
/*--------------------------------------------
 this function:
 0. finds the swap lists for each dimension
 & direction
 1. uses xchng_ph(int,int,int*,int
 ,class VecLst*) to transfer said lists to
 neighboring processors
 
 caveat: it is assumed that x is already is
 converted to s, i.e. the coordinates are
 frational.
 --------------------------------------------*/
/*
void Atoms::setup_ph_old(int box_chng,class VecLst* list)
{
 
    if(box_chng)
    {
        type0 tmp;
        for (int i=0;i<dimension;i++)
        {
            tmp=0.0;
            for(int j=i;j<dimension;j++)
                tmp+=B[j][i]*B[j][i];
            
            cut_ph_s[i]=tot_cut_ph*sqrt(tmp);
        }
        
        for (int i=0;i<dimension;i++)
        {
            s_ph_hi[i]=s_hi[i]-cut_ph_s[i];
            s_ph_lo[i]=s_lo[i]+cut_ph_s[i];
        }
        for (int i=0;i<dimension;i++)
            comm_need[i][0]=comm_need[i][1]=static_cast<int>
            (cut_ph_s[i]/(s_hi[i]-s_lo[i]))+1;
        
        int nswap=0;
        for(int i=0;i<dimension;i++)
            nswap+=comm_need[i][0]+comm_need[i][1];
        
        ph_lst->newlist(nswap);
    }
    
    
    
    
    natms_ph=0;
    ph_lst->reset();
    int lo_r_dir,hi_r_dir,lo_r_dim,hi_r_dim,tmp_r;
    int icomm,iatm,idim;
    type0* s;
    int s_dim=vectors[0]->dim;
    int iswap=0;
    
    lo_r_dim=0;
    hi_r_dim=natms+natms_ph;
    for(idim=0;idim<dimension;idim++)
    {
        lo_r_dir=lo_r_dim;
        hi_r_dir=hi_r_dim;
 
        for(icomm=0;icomm<comm_need[idim][0];icomm++)
        {
            tmp_r=natms+natms_ph;
            vectors[0]->ret(s);
            for(iatm=lo_r_dir;iatm<hi_r_dir;iatm++)
                if(s[iatm*s_dim+idim]<=s_ph_lo[idim])
                    ph_lst->add(iswap,iatm);
            
            xchng_ph_old(0,idim,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap],list);
            hi_r_dir=natms+natms_ph;
            lo_r_dir=tmp_r;
            iswap++;
        }
        
        lo_r_dir=lo_r_dim;
        hi_r_dir=hi_r_dim;
        
        for(icomm=0;icomm<comm_need[idim][1];icomm++)
        {
            tmp_r=natms+natms_ph;
            vectors[0]->ret(s);
            for(iatm=lo_r_dir;iatm<hi_r_dir;iatm++)
                if(s_ph_hi[idim]<=s[iatm*s_dim+idim])
                    ph_lst->add(iswap,iatm);
            xchng_ph_old(1,idim,ph_lst->snd_list[iswap],ph_lst->snd_size[iswap],list);
            hi_r_dir=natms+natms_ph;
            lo_r_dir=tmp_r;
            iswap++;
        }
        
        lo_r_dim=0;
        hi_r_dim=natms+natms_ph;
    }
    
}
*/
/*--------------------------------------------
 this function:
 0. sends/recieves phantom atoms to/from
 neighboring processors.
 1. if needed makes correction to positions (s)
 due to periodic boundary conditions
 --------------------------------------------*/
/*
void Atoms::xchng_ph_old(int dir,int dim,
                         int* atm_list,int atm_list_size,class VecLst* list)
{
    
    if(tot_p_grid[dim]!=1)
    {
        int snd_proc=neigh_p[dim][dir];
        int rcv_proc=neigh_p[dim][1-dir];
        int snd_buff_size=list->ph_byte_size
        *atm_list_size;
        int rcv_buff_size;
        buff_size_management(snd_ph_buff,snd_ph_buff_capacity,snd_buff_size);
        
        MPI_Request request[2];
        MPI_Status status[2];
        MPI_Sendrecv(&snd_buff_size,1,
                     MPI_INT,snd_proc,0,&rcv_buff_size,
                     1,MPI_INT,rcv_proc,0,world,
                     &status[0]);
        
        buff_size_management(rcv_ph_buff,rcv_ph_buff_capacity,rcv_buff_size);
        
        pack_ph(snd_ph_buff,list,atm_list,atm_list_size);
        
        
        if (rcv_buff_size)
        {
            MPI_Irecv(rcv_ph_buff,rcv_buff_size,
                      MPI_BYTE,rcv_proc,0,world,
                      &request[0]);
        }
        if (snd_buff_size)
        {
            MPI_Isend(snd_ph_buff,snd_buff_size,
                      MPI_BYTE,snd_proc,0,world,
                      &request[1]);
        }
        
        int no_new_atms=rcv_buff_size/list->ph_byte_size;
        
        if (rcv_buff_size&&snd_buff_size==0)
        {
            MPI_Wait(&request[0],&status[0]);
            unpack_ph(rcv_ph_buff,list,no_new_atms);
        }
        else if (rcv_buff_size==0&&snd_buff_size)
        {
            MPI_Wait(&request[1],&status[1]);
        }
        else if (rcv_buff_size&&snd_buff_size)
        {
            MPI_Waitall(2,request,status);
            unpack_ph(rcv_ph_buff,list, no_new_atms);
        }
        
        //correction due to periodic boundary conditions
        
        if((dir==0&&my_loc[dim]==tot_p_grid[dim]-1) ||
           (dir==1&&my_loc[dim]==0))
        {
            int last_atm_no=natms+natms_ph;
            int strt_atm_no=last_atm_no-no_new_atms;
            type0* s;
            vectors[0]->ret(s);
            int s_dim=vectors[0]->dim;
            if(dir==0)
                for (int iatm=strt_atm_no;iatm<last_atm_no;iatm++)
                    s[iatm*s_dim+dim]++;
            else if(dir==1)
                for (int iatm=strt_atm_no;iatm<last_atm_no;iatm++)
                    s[iatm*s_dim+dim]--;
        }
    }
    else
    {
        int new_size=natms_ph+natms
        +atm_list_size;
        int old_size=atm_vec_ph_size;
        
        if(new_size>old_size)
        {
            grow(1,new_size-old_size,list);
        }
        
        int last_atm_no=natms+natms_ph;
        
        for(int ivec=0;ivec<list->ph_no_vecs;ivec++)
            for(int iatm=0;iatm<atm_list_size;iatm++)
                vectors[list->ph_vec_list[ivec]]->copy(last_atm_no+iatm,atm_list[iatm]);
        
        
        type0* s;
        vectors[0]->ret(s);
        int s_dim=vectors[0]->dim;
        
        if(dir==0)
            for (int iatm=natms+natms_ph;iatm<natms+natms_ph+atm_list_size;iatm++)
                s[iatm*s_dim+dim]++;
        else if(dir==1)
            for (int iatm=natms+natms_ph;iatm<natms+natms_ph+atm_list_size;iatm++)
                s[iatm*s_dim+dim]--;
        
        natms_ph+=atm_list_size;
    }
    
}
*/

