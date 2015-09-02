#include "ls.h"
#include "ff.h"
#include "clock.h"
#include "atom_types.h"
#define INV_SQ_2 0.7071067811865475
#define TINY 1.0e-15
#define TINY_INV 1.0e15
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearch::LineSearch(MAPP* mapp)
:InitPtrs(mapp)
{
    max_dx=1.0;
    prev_val=0.0;
    golden=0.5+0.5*sqrt(5.0);
    epsilon=numeric_limits<type0>::epsilon();
    epsilon_3_4=pow(epsilon,0.75);
    if(atoms->dimension!=3)
        xmath=new XMath(mapp);
    
    prev_val=-1.0;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
LineSearch::~LineSearch()
{
    if(atoms->dimension!=3)
        delete xmath;
    
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 LineSearch::energy(type0 alpha)
{
    type0* x;
    type0* x_prev;
    type0* h;
    atoms->vectors[0]->ret(x);
    atoms->vectors[x_prev_n]->ret(x_prev);
    atoms->vectors[h_n]->ret(h);
    
    for(int i=0;i<x_dim*atoms->natms;i++)
        x[i]=x_prev[i]+alpha*h[i];
    
    if(chng_box)
    {
        type0** H=atoms->H;
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
            {
                H[i][j]=H_prev[i][j]+alpha*h_H[i][j];
            }
        
        if(dim==3)
            M3INV_TRI_LOWER(atoms->H,atoms->B);
        else
            xmath->invert_lower_triangle(atoms->H,atoms->B,dim);
    }
    
    atoms->update(chng_box,vecs_comm);
    return forcefield->energy_calc_timer();
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
type0 LineSearch::energy(type0 alpha,type0& drev)
{
    type0* x;
    type0* x_prev;
    type0* h;
    type0* f;
    atoms->vectors[0]->ret(x);
    atoms->vectors[x_prev_n]->ret(x_prev);
    atoms->vectors[h_n]->ret(h);
    type0 drev_lcl;
    
    for(int i=0;i<x_dim*atoms->natms;i++)
        x[i]=x_prev[i]+alpha*h[i];
    
    if(chng_box)
    {
        type0** H=atoms->H;
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
            {
                H[i][j]=H_prev[i][j]+alpha*h_H[i][j];
            }
        
        if(dim==3)
            M3INV_TRI_LOWER(atoms->H,atoms->B);
        else
            xmath->invert_lower_triangle(atoms->H,atoms->B,dim);
    }
    
    atoms->update(chng_box,vecs_comm);
    
    atoms->vectors[h_n]->ret(h);
    atoms->vectors[f_n]->ret(f);
    for(int i=0;i<atoms->natms*x_dim;i++)
        f[i]=0.0;
    for(int i=0;i<dim*(dim+1)/2+1;i++)
        nrgy_strss[i]=0.0;
    
    forcefield->force_calc(2,nrgy_strss);
    
    drev_lcl=0.0;
    for(int i=0;i<atoms->natms*x_dim;i++)
        drev_lcl+=f[i]*h[i];
    MPI_Allreduce(&drev_lcl,&drev,1,MPI_TYPE0,MPI_SUM,world);
    
    if(chng_box && dim==3)
    {
        drev+=h_H[0][0]*nrgy_strss[1]
        +h_H[1][1]*nrgy_strss[2]
        +h_H[2][2]*nrgy_strss[3]
        +h_H[2][1]*nrgy_strss[4]
        +h_H[2][0]*nrgy_strss[5]
        +h_H[1][0]*nrgy_strss[5];
    }
    drev*=-1.0;
    
    return nrgy_strss[0];
}
/*--------------------------------------------
 find maximum h
 --------------------------------------------*/
void LineSearch::init_manip(type0& dfa,type0& h_norm,type0& max_a)
{
    type0* h;
    atoms->vectors[h_n]->ret(h);
    type0* f;
    atoms->vectors[f_n]->ret(f);
    type0 max_a_lcl=max_dx;
    type0 h_norm_lcl,dfa_lcl;
    
    h_norm_lcl=0.0;
    dfa_lcl=0.0;
    for(int i=0;i<atoms->natms*x_dim;i++)
    {
        h_norm_lcl+=h[i]*h[i];
        dfa_lcl+=h[i]*f[i];
    }
    
    MPI_Allreduce(&h_norm_lcl,&h_norm,1,MPI_TYPE0,MPI_SUM,world);
    MPI_Allreduce(&dfa_lcl,&dfa,1,MPI_TYPE0,MPI_SUM,world);
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
            {
                h_norm+=h_H[i][j]*h_H[i][j];
                dfa+=h_H[i][j]*f_H[i][j];
            }
    }
    
    dfa*=-1.0;

    if(h_norm==0.0)
    {
        max_a=0.0;
        dfa=0.0;
        return;
    }
    
    if(dfa>=0.0)
    {
        max_a=0.0;
        dfa=1.0;
        return;
    }
    
    h_norm=sqrt(h_norm);
    
    if(mapp->mode==DMD_mode)
    {
        type0* x;
        int icmp=0;
        type0 a_max=forcefield->alpha_max;
        type0 a_min=forcefield->alpha_min;
        max_a_lcl=a_max;
        atoms->vectors[0]->ret(x);
        
        for(int i=0;i<atoms->natms;i++)
        {
            for(int j=0;j<atoms->dimension;j++)
            {
                max_a_lcl=MIN(max_a_lcl,0.999*fabs(max_dx/h[icmp]));
                icmp++;
            }
            for(int j=0;j<atom_types->no_types;j++)
            {
                if(h[icmp]>0.0)
                    max_a_lcl=MIN(max_a_lcl,0.999*(a_max-x[icmp])/h[icmp]);
                if(h[icmp]<0.0)
                    max_a_lcl=MIN(max_a_lcl,0.999*(a_min-x[icmp])/h[icmp]);
                icmp++;
            }
        }
    }
    else
    {
        for(int i=0;i<atoms->natms*x_dim;i++)
            max_a_lcl=MIN(max_a_lcl,0.999*fabs(max_dx/h[i]));
    }
    
    MPI_Allreduce(&max_a_lcl,&max_a,1,MPI_TYPE0,MPI_MAX,world);
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                if(h_H[i][j]!=0.0)
                    max_a=MIN(max_a,0.999*fabs(max_dx/h_H[i][j]));
    }

}
/*--------------------------------------------
 bracketing routine
 --------------------------------------------*/
int LineSearch::bracket(type0 dfa,type0 max_a,type0& a,
type0& b,type0& c,type0& fa,type0& fb,type0& fc)
{
    
    type0 u,fu,r,q,ulim;
    int uphill_iter=4,iter;

    b=MIN(1.e-2*max_a,-epsilon_3_4/dfa);

    //test(fa,dfa,max_a);
    
    r=u=b;
    b=0.0;
    fb=fa;
    iter=uphill_iter;
    

    while(fb>=fa && u<max_a && iter)
    {
        fu=energy(u);
        if(fa<=fu)
            iter--;
        else
            iter=uphill_iter;
        b=u;
        fb=fu;
        u+=r;
    }
    
    if(u>max_a && fb>=fa)
        return B_F_MAX_ALPHA;
    
    if(fb>fa)
    {
        return B_F_DOWNHILL;
    }
    
    fc=fb;
    while (fb>=fc)
    {
        
        c=b+golden*(b-a);
        if(c>=max_a)
        {
            c=max_a;
            fc=energy(c);
            return B_S;
        }
        
        fc=energy(c);
        
        if(fc>fb)
            continue;
        
        ulim=MIN(b+(golden+2.0)*(b-a),max_a);
        
        r=(b-a)*(fb-fc);
        q=(b-c)*(fb-fa);
        
        u=0.5*b+(c*q-a*r)/(2.0*(q-r));
        
        if(b<u && u<c)
        {
            fu=energy(u);
            if(fu<fc)
            {
                a=b;
                b=u;
                fa=fb;
                fb=fu;
                return B_S;
            }
            else if (fu>fb)
            {
                c=u;
                fc=fu;
                return B_S;
            }
            
            a=b;
            b=c;
            
            fa=fb;
            fb=fc;
        }
        else if(u>c)
        {
            u=MIN(u,ulim);
            fu=energy(u);
            
            a=b;
            b=c;
            c=u;
            
            fa=fb;
            fb=fc;
            fc=fu;
            
            if(fu>fc)
            {
                return B_S;
            }
        }
        else
        {
            a=b;
            b=c;
            
            fa=fb;
            fb=fc;
        }
    }
    
    return B_S;
}
/*--------------------------------------------
 reset to initial position
 --------------------------------------------*/
void LineSearch::reset()
{
    type0* x;
    type0* x_prev;
    atoms->vectors[0]->ret(x);
    atoms->vectors[x_prev_n]->ret(x_prev);
    
    memcpy(x,x_prev,atoms->natms*x_dim*sizeof(type0));
    
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
            {
                atoms->H[i][j]=H_prev[i][j];
                atoms->B[i][j]=B_prev[i][j];
            }
    }
    
    atoms->update(chng_box,vecs_comm);
}

/*--------------------------------------------
 reset to initial position
 --------------------------------------------*/
void LineSearch::test(type0 fa,type0 dfa,type0 max_a)
{
    int no=100;
    type0 frac=1.0e-2*max_a;
    type0 dfu,fu,u=0.0;

    printf("dfa %e \n",dfa);
    printf("u fu u*dfa\n");
    
    
    for(int i=0;i<no;i++)
    {
        fu=energy(u,dfu);
        printf("%22.20lf %22.20lf %22.20lf  %22.20lf \n",u,fu-fa,u*dfa,dfu);
        u+=frac;
    }
    //error->abort("");
    
    
}
