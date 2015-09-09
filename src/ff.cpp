/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "ff.h"
#include "atom_types.h"
#include "timer.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField::
ForceField(MAPP* mapp) : InitPtrs(mapp)
{
    cut_alloc=0;
    ns_alloc=0;
    if(atoms->tot_natms==0 ||
       atom_types->no_types==0)
        error->abort("system configuration "
        "should be loaded before initiating ff");

    if(atoms->dimension!=3)
        error->abort("the dimension of the box for ff");
    
    int no_types=atom_types->no_types;
    if(no_types)
    {
        CREATE1D(cut_sk_sq,no_types*(no_types+1)/2);
        CREATE1D(cut_sq,no_types*(no_types+1)/2);
        cut_alloc=1;
    }
    int dim=atoms->dimension;
    if(dim)
    {
        CREATE1D(nrgy_strss,dim*(dim+1)/2+1);
        ns_alloc=1;
    }

    image_size=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField::~ForceField()
{
    if(cut_alloc)
    {
        delete [] cut_sk_sq;
        delete [] cut_sq;        
    }
    if(ns_alloc)
        delete [] nrgy_strss;
    if(image_size)
        delete [] image;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField::force_calc_timer(int chk,type0* enst)
{
    timer->start(FORCE_TIME_mode);
    if(chk==2) image_calc();
    if(chk==3)
    {
        force_calc(1,enst);
        type0** B=atoms->B;
        enst[1]=-(enst[1]*B[0][0]+enst[6]*B[1][0]+enst[5]*B[2][0]);
        enst[2]=-(enst[2]*B[1][1]+enst[4]*B[2][1]);
        enst[3]=-(enst[3]*B[2][2]);
        enst[4]=-(enst[4]*B[2][2]);
        enst[6]=-(enst[6]*B[1][1]+enst[5]*B[2][1]);
        enst[5]=-(enst[5]*B[2][2]);
    }
    else force_calc(chk,enst);
    timer->stop(FORCE_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceField::energy_calc_timer()
{
    type0 en;
    timer->start(FORCE_TIME_mode);
    en=energy_calc();
    timer->stop(FORCE_TIME_mode);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
type0 ForceField::g_calc_timer(int chk
,type0 alpha,type0* a,type0* g,type0* enst)
{
    type0 en;
    timer->start(FORCE_TIME_mode);
    en=g_calc(chk,alpha,a,g,enst);
    timer->stop(FORCE_TIME_mode);
    return en;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField::c_d_calc_timer(int chk,type0* enst)
{
    timer->start(FORCE_TIME_mode);
    c_d_calc(chk,enst);
    timer->stop(FORCE_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField::create_2nd_neigh_lst_timer()
{
    timer->start(FORCE_TIME_mode);
    create_2nd_neigh_lst();
    timer->stop(FORCE_TIME_mode);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField::image_calc()
{
    int tot_natms=atoms->natms+atoms->natms_ph;
    int dimension=atoms->dimension;
    int x_dim=atoms->vectors[0]->dim;
    int icurs=0;
    int icurs_img=0;
    type0* x;
    type0* xtmp;
    type0** H=atoms->H;
    atoms->vectors[0]->ret(x);
    

    int new_size=dimension*(tot_natms);
    if(new_size>image_size)
    {
        if(image_size)
            delete [] image;
        
        CREATE1D(image,new_size);
        image_size=new_size;
    }
    
    
    CREATE1D(xtmp,dimension);
    
    for(int i=0;i<tot_natms;i++)
    {
        for(int j=0;j<dimension;j++)
            xtmp[j]=x[icurs+j];
        
        for(int j=dimension-1;j>-1;j--)
        {
            image[icurs_img+j]=0.0;
            while(xtmp[j]<0.0)
            {
                for(int k=0;k<j+1;k++)
                    xtmp[k]+=H[j][k];
                
                image[icurs_img+j]--;
            }
            
            while(xtmp[j]>=H[j][j])
            {
                for(int k=0;k<j+1;k++)
                    xtmp[k]-=H[j][k];
                
                image[icurs_img+j]++;
            }
        }
      
        icurs_img+=dimension;
        icurs+=x_dim;
    }
    if(dimension)
        delete [] xtmp;
}


