//
//  ff_meam.cpp
//  MAPP
//
//  Created by Sina on 7/23/14.
//  Copyright (c) 2014 Li Group/Sina. All rights reserved.
//
#include "neighbor.h"
#include "ff_meam.h"
#include "atom_types.h"
#include "memory.h"
#include "error.h"
#include <iostream>
#include <fstream>
using namespace MAPP_NS;
enum{FCC,BCC,HCP,DIM,DIAMOND,B1,C11,L12,B2};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_meam::
ForceField_meam(MAPP* mapp):ForceFieldMD(mapp)
{
    max_pairs=0;
    no_types=0;
    rho_dim=38;
    nr=1000;
    third=1.0/3.0;
    sixth=1.0/6.0;
    
    v2d[0]=1;
    v2d[1]=2;
    v2d[2]=2;
    v2d[3]=1;
    v2d[4]=2;
    v2d[5]=1;
    
    v3d[0]=1;
    v3d[1]=3;
    v3d[2]=3;
    v3d[3]=3;
    v3d[4]=6;
    v3d[5]=3;
    v3d[6]=1;
    v3d[7]=3;
    v3d[8]=3;
    v3d[9]=1;

    int comp0=0;
    int comp1=0;
    for(int i=0;i<3;i++)
    {
        for(int j=i;j<3;j++)
        {
            vind2d[i][j]=comp0;
            vind2d[j][i]=comp0;
            comp0++;
            for(int k=j;k<3;k++)
            {
                vind3d[i][j][k]=comp1;
                vind3d[i][k][j]=comp1;
                vind3d[j][i][k]=comp1;
                vind3d[j][k][i]=comp1;
                vind3d[k][i][j]=comp1;
                vind3d[k][j][i]=comp1;
                comp1++;
            }
        }
    }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_meam::~ForceField_meam()
{
    deallocate();
}
/*--------------------------------------------
 allocation
 --------------------------------------------*/
void ForceField_meam::allocate()
{
    if(no_types==atom_types->no_types)
        return;
    
    deallocate();
    
    no_types=atom_types->no_types;
    
    CREATE1D(c_min,no_types);
    CREATE1D(c_max,no_types);
    for(int i=0;i<no_types;i++)
    {
        CREATE1D(c_min[i],no_types);
        CREATE1D(c_max[i],no_types);
    }
    
    for(int j=0;j<no_types;j++)
    {
        for(int i=0;i<no_types;i++)
        {
            CREATE1D(c_min[i][j],no_types);
            CREATE1D(c_max[i][j],no_types);
        }
    }

    int tot_types=(no_types+1)*no_types/2;
    CREATE1D(phirar,tot_types);
    for(int i=0;i<tot_types;i++)
        CREATE1D(phirar[i],nr);
    
    for(int i=0;i<tot_types;i++)
        for(int j=0;j<nr;j++)
            CREATE1D(phirar[i][j],7);
    
    CREATE_2D(re_meam,no_types,no_types);
    CREATE_2D(ebound_meam,no_types,no_types);
    CREATE_2D(Ec_meam,no_types,no_types);
    CREATE_2D(alpha_meam,no_types,no_types);
    CREATE_2D(delta_meam,no_types,no_types);
    CREATE_2D(attrac_meam,no_types,no_types);
    CREATE_2D(repuls_meam,no_types,no_types);
    CREATE_2D(lattice,no_types,no_types);
    CREATE_2D(nn2_meam,no_types,no_types);
    CREATE_2D(zbl_meam,no_types,no_types);
    
    
    CREATE1D(rho0_meam,no_types);
    CREATE1D(beta0_meam,no_types);
    CREATE1D(beta1_meam,no_types);
    CREATE1D(beta2_meam,no_types);
    CREATE1D(beta3_meam,no_types);
    CREATE1D(t0_meam,no_types);
    CREATE1D(t1_meam,no_types);
    CREATE1D(t2_meam,no_types);
    CREATE1D(t3_meam,no_types);
    CREATE1D(Z_meam,no_types);
    CREATE1D(rho_ref_meam,no_types);
    CREATE1D(A_meam,no_types);
    CREATE1D(ibar_meam,no_types);
    CREATE1D(ielt_meam,no_types);
    CREATE1D(type_ref,no_types);
    
    
    for(int i=0;i<no_types;i++)
    {
        for(int j=0;j<no_types;j++)
        {
            Ec_meam[i][j]=0.0;
            alpha_meam[i][j]=0.0;
            lattice[i][j]=-1;
            re_meam[i][j]=0.0;
        }
    }
}
/*--------------------------------------------
 deallocation
 --------------------------------------------*/
void ForceField_meam::deallocate()
{
    if(!no_types)
        return;
    
    for(int j=0;j<no_types;j++)
    {
        for(int i=0;i<no_types;i++)
        {
            delete [] c_min[i][j];
            delete [] c_max[i][j];
        }
    }
    
    for(int i=0;i<no_types;i++)
    {
        delete [] c_min[i];
        delete [] c_max[i];
    }
    delete [] c_min;
    delete [] c_max;
    
    int tot_types=(no_types+1)*no_types/2;
    for(int i=0;i<tot_types;i++)
        for(int j=0;j<nr;j++)
            delete [] phirar[i][j];
    
    for(int i=0;i<tot_types;i++)
        delete [] phirar[i];
    delete [] phirar;
    
    DEL_2D(re_meam);
    DEL_2D(ebound_meam);
    DEL_2D(Ec_meam);
    DEL_2D(alpha_meam);
    DEL_2D(delta_meam);
    DEL_2D(attrac_meam);
    DEL_2D(repuls_meam);
    DEL_2D(lattice);
    DEL_2D(nn2_meam);
    DEL_2D(zbl_meam);
    
    delete [] rho0_meam;
    delete [] beta0_meam;
    delete [] beta1_meam;
    delete [] beta2_meam;
    delete [] beta3_meam;
    delete [] t0_meam;
    delete [] t1_meam;
    delete [] t2_meam;
    delete [] t3_meam;
    delete [] Z_meam;
    delete [] rho_ref_meam;
    delete [] A_meam;
    delete [] ibar_meam;
    delete [] ielt_meam;
    delete [] type_ref;
    
    no_types=0;
}
/*--------------------------------------------
 G_gam
 --------------------------------------------*/
void ForceField_meam::G_gam(type0 Gamma,int ibar
,type0 gsmooth_factor_0,type0& G)
{
    type0 gsmooth_switchpoint;
    if(ibar==0||ibar==4)
    {
        gsmooth_switchpoint=-gsmooth_factor_0/(gsmooth_factor_0+1.0);
        if(Gamma<gsmooth_switchpoint)
        {
            G=pow(gsmooth_switchpoint/Gamma,gsmooth_factor_0)
            /(gsmooth_factor+1.0);
            G=sqrt(G);
        }
        else
            G=sqrt(1.0+Gamma);
    }
    else if(ibar==1)
    {
        G=exp(Gamma*0.5);
    }
    else if(ibar==3)
    {
        G=2.0/(1.0+exp(-Gamma));
    }
    else if(ibar==-5)
    {
        if((1.0+Gamma)>=0.0)
        {
            G = sqrt(1.0+Gamma);
        }
        else
        {
            G=-sqrt(-1.0-Gamma);
        }
    }
    else
        error->abort("ibar error in MEAM");
}
/*--------------------------------------------
 dG_gam
 --------------------------------------------*/
void ForceField_meam::dG_gam(type0 Gamma,int ibar
,type0 gsmooth_factor,type0& G,type0& dG)
{
    type0 gsmooth_switchpoint;
    if(ibar==0||ibar==4)
    {
        gsmooth_switchpoint=-gsmooth_factor/(gsmooth_factor+1.0);
        if(Gamma<gsmooth_switchpoint)
        {
            G=pow(gsmooth_switchpoint/Gamma,gsmooth_factor)
            /(gsmooth_factor+1.0);
            G=sqrt(G);
            dG=0.5/G;
        }
        else
        {
            G=sqrt(1.0+Gamma);
            dG=0.5/G;
        }
    }
    else if(ibar==1)
    {
        G=exp(Gamma*0.5);
        dG=0.5*G;
    }
    else if(ibar==3)
    {
        G=2.0/(1.0+exp(-Gamma));
        dG=G*(1.0-0.5*G);
    }
    else if(ibar==-5)
    {
        if((1.0+Gamma)>=0.0)
        {
            G=sqrt(1.0+Gamma);
            dG=0.5/G;
        }
        else
        {
            G=-sqrt(-1.0-Gamma);
            dG=-0.5/G;
        }
    }
    else
        error->abort("ibar error in MEAM");
}
/*--------------------------------------------
 fcut
 --------------------------------------------*/
void ForceField_meam::fcut(type0 xi,type0& fc)
{
    if(xi>=1.0)
    {
        fc=1.0;
    }
    else if(xi<=0.0)
    {
        fc=0.0;
    }
    else
    {
        type0 a=1.0-xi;
        a*=a;
        a*=a;
        a=1.0-a;
        fc=a*a;
    }
}
/*--------------------------------------------
 dfcut
 --------------------------------------------*/
void ForceField_meam::dfcut(type0 xi,type0& fc
,type0& dfc)
{
    if(xi>1.0)
    {
        fc=1.0;
        dfc=0.0;
    }
    else if(xi<0.0)
    {
        fc=0.0;
        dfc=0.0;
    }
    else
    {
        type0 a=1.0-xi;
        type0 a3=a*a*a;
        type0 a4=a*a*a*a;
        fc=(1.0-a4)*(1.0-a4);
        dfc=8*(1.0-a4)*a3;
    }
}
/*--------------------------------------------
 dCfunc
 --------------------------------------------*/
void ForceField_meam::dCfunc(type0 rij2
,type0 rjk2,type0 rki2,type0& dCikj)
{
    type0 rij4=rij2*rij2;
    type0 a=rki2-rjk2;
    type0 b=rki2+rjk2;
    type0 denom=rij4-a*a;
    denom *= denom;
    dCikj=-4.0*(-2.0*rij2*a*a +rij4*b+a*a*b)/denom;
}
/*--------------------------------------------
 dCfunc2
 --------------------------------------------*/
void ForceField_meam::dCfunc2(type0 rij2
,type0 rjk2,type0 rki2,type0& dCikj1,type0& dCikj2)
{
    type0 rij4=rij2*rij2;
    type0 rjk4=rjk2*rjk2;
    type0 rki4=rki2*rki2;
    type0 a=rki2-rjk2;
    type0 denom=rij4-a*a;
    denom *= denom;
    
    dCikj1=4.0*rij2*(rij4+rki4
                     +2.0*rki2*rjk2-3.0*rjk4
                     -2.0*rij2*a)/denom;
    
    dCikj2=4.0*rij2*(rij4-3.0*rki4
                     +2.0*rki2*rjk2+rjk4
                     +2.0*rij2*a)/denom;
}
/*--------------------------------------------
 dCfunc2
 checked
 --------------------------------------------*/
void ForceField_meam::
get_shpfcn(type0* s,int latt)
{
    if(latt==FCC ||
       latt==BCC ||
       latt==B1  ||
       latt==B2)
    {
        s[0]=0.0;
        s[1]=0.0;
        s[2]=0.0;
    }
    else if(latt==HCP)
    {
        s[0]=0.0;
        s[1]=0.0;
        s[2]=1.0/3.0;
    }
    else if(latt==DIAMOND)
    {
        s[0]=0.0;
        s[1]=0.0;
        s[2]=32.0/9.0;
    }
    else if(latt==DIM)
    {
        s[0]=1.0;
        s[1]=2.0/3.0;
        s[2]=0.4;
    }
}
/*--------------------------------------------
 Zij
 --------------------------------------------*/
void ForceField_meam::get_Zij(type0& Zij,int latt)
{
    if(latt==FCC)
        Zij=12.0;
    else if(latt==BCC)
        Zij=8.0;
    else if(latt==HCP)
        Zij=12.0;
    else if(latt==DIM)
        Zij=1.0;
    else if(latt==DIAMOND)
        Zij=4.0;
    else if(latt==B1)
        Zij=6.0;
    else if(latt==C11)
        Zij=10.0;
    else if(latt==L12)
        Zij=12.0;
    else if(latt==B2)
        Zij=8.0;
    
}
/*--------------------------------------------
 Zij
 --------------------------------------------*/
void ForceField_meam::get_Zij2(type0& Zij2,type0& a,
type0& S,int latt,type0 cmin,type0 cmax)
{
    int num;
    if(latt==FCC)
    {
        Zij2=6.0;
        a=sqrt(2.0);
        num=4;
    }
    else if(latt==BCC)
    {
        Zij2=6.0;
        a=2.0/sqrt(3.0);
        num=4;
    }
    else if(latt==HCP)
    {
        Zij2=6.0;
        a=sqrt(2.0);
        num=4;
    }
    else if(latt==DIM)
    {
        Zij2=0.0;
        a=1.0;
        S=0.0;
        return;
    }
    else if(latt==DIAMOND)
    {
        Zij2=0.0;
        a=sqrt(8.0/3.0);
        num=4;
    }
    else if(latt==B1)
    {
        Zij2=12.0;
        a=sqrt(2.0);
        num=2;
    }
    else if(latt==L12)
    {
        Zij2=6.0;
        a=sqrt(2.0);
        num=4;
    }
    else if(latt==B2)
    {
        Zij2=6.0;
        a=2.0/sqrt(3.0);
        num=4;
    }
    else error->abort("wrong type");
    
    type0 C=4.0/(a*a)-1.0;
    type0 x=(C-cmin)/(cmax-cmin);
    type0 sijk;
    fcut(x,sijk);
    S=pow(sijk,num);
    
}
/*--------------------------------------------
 zbl
 --------------------------------------------*/
type0 ForceField_meam::zbl(type0 r,type0 z1
,type0 z2)
{
    type0 a=0.4685/(pow(z1,0.23)+pow(z2,0.23));
    type0 zbl=0.0;
    type0 x=r/a;
    
    zbl=0.028171*exp(-0.20162*x)
    +0.28022*exp(-0.40290*x)
    +0.50986*exp(-0.94229*x)
    +0.18175*exp(-3.1998*x);
    if(r>0.0) zbl*=z1*z2*14.3997/r;
    
    return zbl;
}
/*--------------------------------------------
 Rose energy function
 --------------------------------------------*/
type0 ForceField_meam::erose(type0 r,type0 re,
type0 alpha,type0 Ec,type0 repuls,type0 attrac,
int form)
{
    type0 a3,astar,erose=0.0;
    
    if(r>0.0)
    {
        astar=alpha*(r/re-1.0);
        a3=0.0;
        if(astar>=0.0)
            a3=attrac;
        else
            a3=repuls;
        
        if(form==1)
            erose=-Ec*(1.0+astar+(repuls/r-attrac)*
            astar*astar*astar)*exp(-astar);
        else if(form==2)
            erose=-Ec*(1.0+astar+
            a3*astar*astar*astar)*exp(-astar);
        else
            erose=-Ec*(1.0+astar+
            a3*(astar*astar*astar)/(r/re))*exp(-astar);
    }
    
    return erose;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void ForceField_meam::init()
{
    neighbor->pair_wise=false;
    rho_vec_ptr=new Vec<type0>(atoms,1);
    rho_ptr=new Vec<type0>(atoms,rho_dim);
}
/*--------------------------------------------
 fin
 --------------------------------------------*/
void ForceField_meam::fin()
{
    delete rho_ptr;
    delete rho_vec_ptr;
    
    if(max_pairs)
    {
        delete [] scrfcn;
        delete [] dscrfcn;
        delete [] fcpair;
        max_pairs=0;
    }
}
/*--------------------------------------------
 init xchng
 --------------------------------------------*/
void ForceField_meam::init_xchng()
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 fin xchng
 --------------------------------------------*/
void ForceField_meam::fin_xchng()
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 pre xchng energy
 --------------------------------------------*/
void ForceField_meam::pre_xchng_energy(GCMC*)
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 xchng energy
 --------------------------------------------*/
type0 ForceField_meam::xchng_energy(GCMC*)
{
    error->abort("exchange has not been set for this forcefield");
    return 0.0;
}
/*--------------------------------------------
 post xchng energy
 --------------------------------------------*/
void ForceField_meam::post_xchng_energy(GCMC*)
{
    error->abort("exchange has not been set for this forcefield");
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void ForceField_meam::force_calc
(bool st_clc)
{
    nrgy_strss_lcl[0]=0.0;
    if(st_clc)
        for(int i=1;i<7;i++)
            nrgy_strss_lcl[i]=0.0;
    
    reset();
    
    int* my_list;
    int my_list_size;
    
    int jatm,katm;
    
    type0 xij,xjk,xki,yij,yjk,yki,zij,zjk,zki;
    type0 rij2,rjk2,rki2;
    type0 rij,rij3;
    type0 drhoa0i,drhoa0j,drhoa1i,drhoa1j,drhoa2i,drhoa2j,drhoa3i,drhoa3j,recip;
    type0 arg1i1,arg1j1,arg1i2,arg1j2,arg1i3,arg1j3,arg3i3,arg3j3,arg;
    int nv2_comp,nv3_comp;
    type0 drho0dr1,drho0dr2,a1,a2,a3,a3a,drho1dr1,drho1dr2,drho2dr1,drho2dr2,drho3dr1,drho3dr2;
    type0 drho0ds1,drho0ds2,drho1ds1,drho1ds2,drho2ds1,drho2ds2,drho3ds1,drho3ds2;
    type0 dt1ds1,dt1ds2,dt2ds1,dt2ds2,dt3ds1,dt3ds2;
    type0 dsij1,dsij2,dcikj1,dcikj2,dfc,force1,force2;
    type0 drhods1=0.0,drhods2=0.0;
    type0 re_meam_inv_i,re_meam_inv_j;
    type0 dUdrij,dUdsij,phi,phip,force,forcem;
    
    int icomp,jcomp,kcomp,itype,jtype,ktype;
    type0 coef1,coef2;
    int icomp_rho,jcomp_rho;
    
    
    type0 fcij,dfcij,sij,dsij;
    type0 rnorm,rbound,x_ki,x_jk,a;
    type0 cikj,sikj=0.0,cmax,cmin,dfikj,dcikj;
    type0 ai,aj,ro0i,ro0j,rhoa0i,rhoa0j,rhoa1i,rhoa1j
    ,rhoa2i,rhoa2j,rhoa3i,rhoa3j;
    type0 Z;
    type0 G,dG,Gbar,dGbar=0.0,gam,rho_bkgd,rhob,denom,B;
    int iArho1_comp,jArho1_comp,iArho2_comp
    ,jArho2_comp,iArho3_comp,jArho3_comp
    ,iArhob3_comp,jArhob3_comp;
    type0 A1i,A2i,A3i,A1j,A2j,A3j,B1i,B1j;
    type0 drhodr1,drhodr2;
    
    type0 t1i,t2i,t3i,t1j,t2j,t3j;
    type0 a1i,a2i,a3i,a1j,a2j,a3j,dt1dr1,dt1dr2,dt2dr1,dt2dr2,dt3dr1,dt3dr2;
    
    type0 delij[3],delik[3],deljk[3];
    type0 drhodrm1[3],drhodrm2[3],drho1drm1[3],drho1drm2[3],drho2drm1[3],drho2drm2[3],drho3drm1[3],drho3drm2[3];
    type0 s[3],si[3],sj[3];
    type0 dUdrijm[3];
    type0 fi[3],fj[3];
    type0 v[6];
    
    type0* x=mapp->x->begin();
    type0* fvec=f->begin();
    type0* rho=rho_ptr->begin();
    type0* rho_vec=rho_vec_ptr->begin();
    md_type* type=mapp->type->begin();
    
    int kn,kk,curs;
    type0 pp;
    type0* coef;
    
    
    int istart=0;
    
    type0 a4,rij4,b;
    
    for(int iatm=0;iatm<atoms->natms;iatm++)
    {
        my_list_size=neighbor->neighbor_list_size[iatm];
        my_list=neighbor->neighbor_list[iatm];
        icomp=iatm*3;
        icomp_rho=iatm*rho_dim;
        itype=type[iatm];
        for(int jn=0;jn<my_list_size;jn++)
        {
            jatm=my_list[jn];
            if(jatm>iatm)
            {
                fcpair[istart]=scrfcn[istart]=dscrfcn[istart]=0.0;
                jtype=type[jatm];
                jcomp=jatm*3;
                jcomp_rho=jatm*rho_dim;
                
                xij=x[icomp]-x[jcomp];
                yij=x[icomp+1]-x[jcomp+1];
                zij=x[icomp+2]-x[jcomp+2];
                delij[0]=-xij;
                delij[1]=-yij;
                delij[2]=-zij;
                rij2=xij*xij+yij*yij+zij*zij;
                rij=sqrt(rij2);
                
                if(rij<rc_meam)
                {
                    rnorm=(rc_meam-rij)*delr_meam_inv;
                    sij=1.0;
                    dsij=0.0;
                    rbound = ebound_meam[itype][jtype]*rij2;
                    rij4=rij2*rij2;
                    kn=0;
                    dfcut(rnorm,fcij,dfcij);
                    dfcij*=delr_meam_inv;
                    while(kn<my_list_size && sij!=0.0)
                    {
                        katm=my_list[kn];
                        
                        if(katm!=jatm)
                        {
                            ktype=type[katm];
                            kcomp=3*katm;
                            
                            xjk=x[jcomp]-x[kcomp];
                            yjk=x[jcomp+1]-x[kcomp+1];
                            zjk=x[jcomp+2]-x[kcomp+2];
                            rjk2=xjk*xjk+yjk*yjk+zjk*zjk;
                            
                            if(rjk2<rbound)
                            {
                                xki=x[kcomp]-x[icomp];
                                yki=x[kcomp+1]-x[icomp+1];
                                zki=x[kcomp+2]-x[icomp+2];
                                rki2=xki*xki+yki*yki+zki*zki;
                                
                                if(rki2<rbound)
                                {
                                    x_ki=rki2/rij2;
                                    x_jk=rjk2/rij2;
                                    a=1.0-(x_ki-x_jk)*(x_ki-x_jk);
                                    if(a>0)
                                    {
                                        cikj=(2.0*(x_ki+x_jk)+a-2.0)/a;
                                        cmax=c_max[itype][jtype][ktype];
                                        cmin=c_min[itype][jtype][ktype];
                                        if(cikj<cmax&&cikj>cmin)
                                        {
                                            cikj=(cikj-cmin)/(cmax-cmin);
                                            dfcut(cikj,sikj,dfikj);
                                            
                                            if(cikj>1.0)
                                            {
                                                sikj=1.0;
                                                dfikj=0.0;
                                            }
                                            else if(cikj<0.0)
                                            {
                                                sikj=0.0;
                                                dfikj=0.0;
                                            }
                                            else
                                            {
                                                a=1.0-cikj;
                                                a3=a*a*a;
                                                a4=a*a*a*a;
                                                sikj=(1.0-a4)*(1.0-a4);
                                                dfikj=8*(1.0-a4)*a3;
                                            }
                                            
                                            sij*=sikj;
                                            coef1=dfikj/((cmax-cmin)*sikj);
                                            dCfunc(rij2,rki2,rjk2,dcikj);
                                            
                                            a=rjk2-rki2;
                                            b=rki2+rjk2;
                                            denom=rij4-a*a;
                                            denom *= denom;
                                            dcikj=-4.0*(-2.0*rij2*a*a +rij4*b+a*a*b)/denom;
                                            
                                            dsij+=coef1*dcikj;
                                        }
                                        else if(cikj<cmin)
                                        {
                                            sij=0.0;
                                            dsij=0.0;
                                        }
                                    }
                                }
                            }
                        }
                        kn++;
                    }
                    coef1=sij*fcij;
                    coef2=sij*dfcij/rij;
                    dsij*=coef1;
                    dsij-=coef2;
                    if(sij==1.0 || sij==0.0)
                        dsij=0.0;
                    if(sij*fcij==1.0 || sij*fcij==0.0)
                        dsij=0.0;
                }
                else
                {
                    fcij=0.0;
                    dfcij=0.0;
                    sij=0.0;
                    dsij=0.0;
                }
                
                fcpair[istart]=fcij;
                scrfcn[istart]=sij;
                dscrfcn[istart]=dsij;
                
                if(sij*fcij!=0.0)
                {
                    sij*=fcij;
                    ai=rij/re_meam[itype][itype]-1.0;
                    aj=rij/re_meam[jtype][jtype]-1.0;
                    ro0i=rho0_meam[itype];
                    ro0j=rho0_meam[jtype];
                    
                    rhoa0i=sij*ro0i*exp(-beta0_meam[itype]*ai);
                    rhoa1i=sij*ro0i*exp(-beta1_meam[itype]*ai);
                    rhoa2i=sij*ro0i*exp(-beta2_meam[itype]*ai);
                    rhoa3i=sij*ro0i*exp(-beta3_meam[itype]*ai);
                    
                    rhoa0j=sij*ro0j*exp(-beta0_meam[jtype]*aj);
                    rhoa1j=sij*ro0j*exp(-beta1_meam[jtype]*aj);
                    rhoa2j=sij*ro0j*exp(-beta2_meam[jtype]*aj);
                    rhoa3j=sij*ro0j*exp(-beta3_meam[jtype]*aj);
                    
                    
                    if(ialloy==1)
                    {
                        rhoa1i*=t1_meam[itype];
                        rhoa2i*=t2_meam[itype];
                        rhoa3i*=t3_meam[itype];
                        
                        rhoa1j*=t1_meam[jtype];
                        rhoa2j*=t2_meam[jtype];
                        rhoa3j*=t3_meam[jtype];
                    }
                    
                    rho[icomp_rho]+=rhoa0j;
                    if(jatm<atoms->natms)
                        rho[jcomp_rho]+=rhoa0i;
                    
                    
                    if(ialloy!=2)
                    {
                        
                        rho[icomp_rho+27]+=t1_meam[jtype]*rhoa0j;
                        rho[icomp_rho+28]+=t2_meam[jtype]*rhoa0j;
                        rho[icomp_rho+29]+=t3_meam[jtype]*rhoa0j;
                        if(jatm<atoms->natms)
                        {
                            rho[jcomp_rho+27]+=t1_meam[itype]*rhoa0i;
                            rho[jcomp_rho+28]+=t2_meam[itype]*rhoa0i;
                            rho[jcomp_rho+29]+=t3_meam[itype]*rhoa0i;
                        }
                        
                    }
                    if(ialloy==1)
                    {
                        
                        rho[icomp_rho+30]+=t1_meam[jtype]*t1_meam[jtype]*rhoa0j;
                        rho[icomp_rho+31]+=t2_meam[jtype]*t2_meam[jtype]*rhoa0j;
                        rho[icomp_rho+32]+=t3_meam[jtype]*t3_meam[jtype]*rhoa0j;
                        if(jatm<atoms->natms)
                        {
                            rho[jcomp_rho+30]+=t1_meam[itype]*t1_meam[itype]*rhoa0i;
                            rho[jcomp_rho+31]+=t2_meam[itype]*t2_meam[itype]*rhoa0i;
                            rho[jcomp_rho+32]+=t3_meam[itype]*t3_meam[itype]*rhoa0i;
                        }
                    }
                    
                    rho[icomp_rho+13]+=rhoa2j;
                    if(jatm<atoms->natms)
                        rho[jcomp_rho+13]+=rhoa2i;
                    
                    
                    iArho1_comp=icomp_rho+4;
                    iArho2_comp=icomp_rho+7;
                    iArho3_comp=icomp_rho+14;
                    iArhob3_comp=icomp_rho+24;
                    
                    jArho1_comp=jcomp_rho+4;
                    jArho2_comp=jcomp_rho+7;
                    jArho3_comp=jcomp_rho+14;
                    jArhob3_comp=jcomp_rho+24;
                    
                    A1j=rhoa1j/rij;
                    B1j=rhoa3j/rij;
                    A2j=rhoa2j/rij2;
                    A3j=rhoa3j/(rij2*rij);
                    
                    for(int i=0;i<3;i++)
                    {
                        rho[iArho1_comp++]+=delij[i]*A1j;
                        rho[iArhob3_comp++]+=delij[i]*B1j;
                        for(int j=i;j<3;j++)
                        {
                            rho[iArho2_comp++]+=delij[i]*delij[j]*A2j;
                            for(int k=j;k<3;k++)
                            {
                                rho[iArho3_comp++]+=delij[i]*delij[j]*delij[k]*A3j;
                            }
                        }
                    }
                    if(jatm<atoms->natms)
                    {
                        
                        A1i=rhoa1i/rij;
                        B1i=rhoa3i/rij;
                        A2i=rhoa2i/rij2;
                        A3i=rhoa3i/(rij2*rij);
                        
                        for(int i=0;i<3;i++)
                        {
                            rho[jArho1_comp++]-=delij[i]*A1i;
                            rho[jArhob3_comp++]-=delij[i]*B1i;
                            for(int j=i;j<3;j++)
                            {
                                rho[jArho2_comp++]+=delij[i]*delij[j]*A2i;
                                for(int k=j;k<3;k++)
                                {
                                    rho[jArho3_comp++]-=delij[i]*delij[j]*delij[k]*A3i;
                                }
                            }
                        }
                    }
                    
                }
                istart++;
            }
        }
        
        
        
        rho[icomp_rho+1]=0.0;
        rho[icomp_rho+2]=-rho[icomp_rho+13]*rho[icomp_rho+13]/3.0;
        rho[icomp_rho+3]=0.0;
        for(int i=0;i<3;i++)
        {
            rho[icomp_rho+1]+=rho[icomp_rho+4+i]*rho[icomp_rho+4+i];
            rho[icomp_rho+3]-=0.6*rho[icomp_rho+24+i]*rho[icomp_rho+24+i];
        }
        
        for(int i=0;i<6;i++)
        {
            rho[icomp_rho+2]+=v2d[i]*rho[icomp_rho+7+i]*rho[icomp_rho+7+i];
        }
        
        for(int i=0;i<10;i++)
        {
            rho[icomp_rho+3]+=v3d[i]*rho[icomp_rho+14+i]*rho[icomp_rho+14+i];
        }
        
        if(rho[icomp_rho]>0.0)
        {
            if(ialloy==1)
            {
                rho[icomp_rho+27]=rho[icomp_rho+27]/rho[icomp_rho+30];
                rho[icomp_rho+28]=rho[icomp_rho+28]/rho[icomp_rho+31];
                rho[icomp_rho+29]=rho[icomp_rho+29]/rho[icomp_rho+32];
            }
            else if(ialloy==2)
            {
                rho[icomp_rho+27]=t1_meam[itype];
                rho[icomp_rho+28]=t2_meam[itype];
                rho[icomp_rho+29]=t3_meam[itype];
            }
            else
            {
                rho[icomp_rho+27]=rho[icomp_rho+27]/rho[icomp_rho];
                rho[icomp_rho+28]=rho[icomp_rho+28]/rho[icomp_rho];
                rho[icomp_rho+29]=rho[icomp_rho+29]/rho[icomp_rho];
            }
        }
        
        
        rho[icomp_rho+33]=rho[icomp_rho+27]*rho[icomp_rho+1]
        +rho[icomp_rho+28]*rho[icomp_rho+2]+rho[icomp_rho+29]*rho[icomp_rho+3];
        
        if(rho[icomp_rho]>0.0)
            rho[icomp_rho+33]=rho[icomp_rho+33]/(rho[icomp_rho]*rho[icomp_rho]);
        
        
        
        Z=Z_meam[itype];
        
        G_gam(rho[icomp_rho+33],ibar_meam[itype],gsmooth_factor,G);
        
        if(ibar_meam[itype]<=0)
            Gbar=1.0;
        else
        {
            get_shpfcn(s,lattice[itype][itype]);
            if(mix_ref_t==1)
                gam=(rho[icomp_rho+27]*s[0]+rho[icomp_rho+28]*s[1]+rho[icomp_rho+29]*s[2])/(Z*Z);
            else
                gam=(t1_meam[itype]*s[0]+t2_meam[itype]*s[1]+t3_meam[itype]*s[2])/(Z*Z);
            
            G_gam(gam,ibar_meam[itype],gsmooth_factor,Gbar);
        }
        
        rho_vec[iatm]=rho[icomp_rho]*G;
        
        if(mix_ref_t==1)
        {
            if(ibar_meam[itype]<=0.0)
            {
                Gbar=1.0;
                dGbar=0.0;
            }
            else
            {
                get_shpfcn(s,lattice[itype][itype]);
                gam=(rho[icomp_rho+27]*s[0]+rho[icomp_rho+28]*s[1]+rho[icomp_rho+29]*s[2])/(Z*Z);
                dG_gam(gam,ibar_meam[itype],gsmooth_factor,Gbar,dGbar);
            }
            rho_bkgd=rho0_meam[itype]*Z*Gbar;
        }
        else
        {
            if(bkgd_dyn==1)
            {
                rho_bkgd=rho0_meam[itype]*Z;
            }
            else
                rho_bkgd=rho_ref_meam[itype];
        }
        
        rhob=rho_vec[iatm]/rho_bkgd;
        denom=1.0/rho_bkgd;
        dG_gam(rho[icomp_rho+33],ibar_meam[itype],gsmooth_factor,G,dG);
        rho[icomp_rho+34]=(G-2.0*dG*rho[icomp_rho+33])*denom;
        
        if(rho[icomp_rho]!=0)
            rho[icomp_rho+35]=(dG/rho[icomp_rho])*denom;
        else
            rho[icomp_rho+35]=0.0;
        
        if(mix_ref_t==1)
            rho[icomp_rho+36]=rho[icomp_rho]*G*dGbar*denom/(Gbar*Z*Z);
        else
            rho[icomp_rho+36]=0.0;
        
        B=A_meam[itype]*Ec_meam[itype][itype];
        
        if(rhob!=0.0)
        {
            if(emb_lin_neg==1 && rhob<=0.0)
            {
                rho[icomp_rho+37]=-B;
                nrgy_strss_lcl[0]-=B*rhob;
            }
            else
            {
                rho[icomp_rho+37]=B*(log(rhob)+1.0);
                nrgy_strss_lcl[0]+=B*rhob*log(rhob);
            }
        }
        else
        {
            if(emb_lin_neg==1)
                rho[icomp_rho+37]=-B;
            else
                rho[icomp_rho+37]=B;
        }
        
    }
    
    /*----------------------------------------------------------------------------*/
    
    
    atoms->update(rho_ptr);
    
    
    /*----------------------------------------------------------------------------*/
    
    istart=0;
    for(int iatm=0;iatm<atoms->natms;iatm++)
    {
        my_list=neighbor->neighbor_list[iatm];
        my_list_size=neighbor->neighbor_list_size[iatm];
        icomp=iatm*3;
        icomp_rho=iatm*rho_dim;
        itype=type[iatm];
        
        for(int jn=0;jn<my_list_size;jn++)
        {
            jatm=my_list[jn];
            if(jatm>iatm)
            {
                if(scrfcn[istart]!=0.0)
                {
                    jtype=type[jatm];
                    jcomp=jatm*3;
                    jcomp_rho=jatm*rho_dim;
                    
                    xij=x[icomp]-x[jcomp];
                    yij=x[icomp+1]-x[jcomp+1];
                    zij=x[icomp+2]-x[jcomp+2];
                    delij[0]=-xij;
                    delij[1]=-yij;
                    delij[2]=-zij;
                    rij2=xij*xij+yij*yij+zij*zij;
                    rij=sqrt(rij2);
                    if(rij<rc_meam)
                    {
                        sij=scrfcn[istart]*fcpair[istart];
                        curs=COMP(itype,jtype);
                        pp=rij*dr_inv;
                        kk=static_cast<int>(pp);
                        kk=MIN(kk,nr-2);
                        pp-=static_cast<type0>(kk);
                        pp=MIN(pp,1.0);
                        coef=phirar[curs][kk];
                        phi=((coef[3]*pp+coef[2])*pp+coef[1])*pp+coef[0];
                        phip=(coef[6]*pp+coef[5])*pp+coef[4];
                        
                        if(jatm<atoms->natms)
                        {
                            nrgy_strss_lcl[0]+=phi*sij;
                        }
                        else
                        {
                            nrgy_strss_lcl[0]+=0.5*phi*sij;
                        }
                        
                        /*
                         energy sij*phi
                         */
                        //printf("%20.15lf\n",phi);
                        recip=1.0/rij;
                        re_meam_inv_i=1.0/re_meam[itype][itype];
                        ai=rij*re_meam_inv_i-1.0;
                        ro0i=rho0_meam[itype];
                        rhoa0i=ro0i*exp(-beta0_meam[itype]*ai);
                        drhoa0i=-beta0_meam[itype]*re_meam_inv_i*rhoa0i;
                        
                        rhoa1i=ro0i*exp(-beta1_meam[itype]*ai);
                        drhoa1i=-beta1_meam[itype]*re_meam_inv_i*rhoa1i;
                        
                        rhoa2i=ro0i*exp(-beta2_meam[itype]*ai);
                        drhoa2i=-beta2_meam[itype]*re_meam_inv_i*rhoa2i;
                        
                        rhoa3i=ro0i*exp(-beta3_meam[itype]*ai);
                        drhoa3i=-beta3_meam[itype]*re_meam_inv_i*rhoa3i;
                        
                        if(itype!=jtype)
                        {
                            
                            re_meam_inv_j=1.0/re_meam[jtype][jtype];
                            aj=rij*re_meam_inv_j-1.0;
                            ro0j=rho0_meam[jtype];
                            rhoa0j=ro0j*exp(-beta0_meam[jtype]*aj);
                            drhoa0j=-beta0_meam[jtype]*re_meam_inv_j*rhoa0j;
                            
                            rhoa1j=ro0j*exp(-beta1_meam[jtype]*aj);
                            drhoa1j=-beta1_meam[jtype]*re_meam_inv_j*rhoa1j;
                            
                            rhoa2j=ro0j*exp(-beta2_meam[jtype]*aj);
                            drhoa2j=-beta2_meam[jtype]*re_meam_inv_j*rhoa2j;
                            
                            rhoa3j=ro0j*exp(-beta3_meam[jtype]*aj);
                            drhoa3j=-beta3_meam[jtype]*re_meam_inv_j*rhoa3j;
                        }
                        else
                        {
                            rhoa0j=rhoa0i;
                            drhoa0j=drhoa0i;
                            rhoa1j=rhoa1i;
                            drhoa1j=drhoa1i;
                            rhoa2j=rhoa2i;
                            drhoa2j=drhoa2i;
                            rhoa3j=rhoa3i;
                            drhoa3j=drhoa3i;
                        }
                        
                        if(ialloy==1)
                        {
                            rhoa1j*=t1_meam[jtype];
                            rhoa2j*=t2_meam[jtype];
                            rhoa3j*=t3_meam[jtype];
                            rhoa1i*=t1_meam[itype];
                            rhoa2i*=t2_meam[itype];
                            rhoa3i*=t3_meam[itype];
                            drhoa1j*=t1_meam[jtype];
                            drhoa2j*=t2_meam[jtype];
                            drhoa3j*=t3_meam[jtype];
                            drhoa1i*=t1_meam[itype];
                            drhoa2i*=t2_meam[itype];
                            drhoa3i*=t3_meam[itype];
                        }
                        
                        
                        
                        
                        arg1i1=0.0;
                        arg1j1=0.0;
                        arg1i2=0.0;
                        arg1j2=0.0;
                        arg1i3=0.0;
                        arg1j3=0.0;
                        arg3i3=0.0;
                        arg3j3=0.0;
                        
                        nv2_comp=0;
                        nv3_comp=0;
                        for(int i=0;i<3;i++)
                        {
                            for(int j=i;j<3;j++)
                            {
                                for(int k=j;k<3;k++)
                                {
                                    arg=delij[i]*delij[j]*delij[k]*v3d[nv3_comp];
                                    arg1i3+=rho[icomp_rho+14+nv3_comp]*arg;
                                    arg1j3-=rho[jcomp_rho+14+nv3_comp]*arg;
                                    nv3_comp++;
                                }
                                arg=delij[i]*delij[j]*v2d[nv2_comp];
                                arg1i2+=rho[icomp_rho+7+nv2_comp];
                                arg1j2+=rho[jcomp_rho+7+nv2_comp];
                                nv2_comp++;
                            }
                            arg1i1+=rho[icomp_rho+4+i]*delij[i];
                            arg1j1-=rho[jcomp_rho+4+i]*delij[i];
                            arg3i3+=rho[icomp_rho+24+i]*delij[i];
                            arg3j3-=rho[jcomp_rho+24+i]*delij[i];
                        }
                        
                        
                        
                        
                        
                        
                        //rho0 terms
                        drho0dr1=drhoa0j*sij;
                        drho0dr2=drhoa0i*sij;
                        
                        //rho1 terms
                        a1=2.0*sij/rij;
                        drho1dr1=a1*(drhoa1j-rhoa1j/rij)*arg1i1;
                        drho1dr2=a1*(drhoa1i-rhoa1i/rij)*arg1j1;
                        
                        
                        for(int i=0;i<3;i++)
                        {
                            drho1drm1[i]=a1*rhoa1j*rho[icomp_rho+4+i];
                            drho1drm2[i]=-a1*rhoa1i*rho[jcomp_rho+4+i];
                        }
                        
                        //rho2 terms
                        a2=2.0*sij/rij2;
                        drho2dr1=a2*(drhoa2j-2.0*rhoa2j/rij)*arg1i2-2.0/3.0*rho[icomp_rho+13]*drhoa2j*sij;
                        drho2dr2=a2*(drhoa2i-2.0*rhoa2i/rij)*arg1j2-2.0/3.0*rho[jcomp_rho+13]*drhoa2i*sij;
                        a2=4.0*sij/rij2;
                        
                        
                        
                        for(int i=0;i<3;i++)
                        {
                            drho2drm1[i]=0.0;
                            drho2drm2[i]=0.0;
                            for(int j=0;j<3;j++)
                            {
                                drho2drm1[i]+=rho[icomp_rho+7+vind2d[i][j]]*delij[j];
                                drho2drm2[i]-=rho[jcomp_rho+7+vind2d[i][j]]*delij[j];
                            }
                            drho2drm1[i]*=a2*rhoa2j;
                            drho2drm2[i]*=-a2*rhoa2i;
                        }
                        
                        rij3=rij*rij2;
                        a3=2.0*sij/rij3;
                        a3a=1.2*sij/rij;
                        drho3dr1 = a3*(drhoa3j-3.0*rhoa3j/rij)*arg1i3- a3a*(drhoa3j - rhoa3j/rij)*arg3i3;
                        drho3dr2 = a3*(drhoa3i-3.0*rhoa3i/rij)*arg1j3- a3a*(drhoa3i - rhoa3i/rij)*arg3j3;
                        a3 = 6.0*sij/rij3;
                        a3a = 1.2*sij/rij;
                        
                        
                        for(int i=0;i<3;i++)
                        {
                            drho3drm1[i]=0.0;
                            drho3drm2[i]=0.0;
                            nv2_comp=0;
                            for(int j=0;j<3;j++)
                            {
                                for(int k=j;k<3;k++)
                                {
                                    arg = delij[j]*delij[k]*v2d[nv2_comp];
                                    drho3drm1[i]+=rho[icomp_rho+14+vind3d[i][j][k]]*arg;
                                    drho3drm2[i]+=rho[jcomp_rho+14+vind3d[i][j][k]]*arg;
                                    nv2_comp++;
                                }
                            }

                            
                            drho3drm1[i]=(a3*drho3drm1[i]-a3a*rho[icomp_rho+24+i])*rhoa3j;
                            drho3drm2[i]=-(a3*drho3drm2[i]-a3a*rho[jcomp_rho+24+i])*rhoa3i;
                        }
                        
                        
                        t1i=rho[icomp_rho+27];
                        t2i=rho[icomp_rho+28];
                        t3i=rho[icomp_rho+29];
                        t1j=rho[jcomp_rho+27];
                        t2j=rho[jcomp_rho+28];
                        t3j=rho[jcomp_rho+29];
                        
                        
                        
                        if(ialloy==1)
                        {
                            a1i=0.0;
                            a2i=0.0;
                            a3i=0.0;
                            a1j=0.0;
                            a2j=0.0;
                            a3j=0.0;
                            
                            if(rho[icomp_rho+30]!=0.0)
                                a1i=drhoa0j*sij/rho[icomp_rho+30];
                            if(rho[jcomp_rho+30]!=0.0)
                                a1j=drhoa0i*sij/rho[jcomp_rho+30];
                            if(rho[icomp_rho+31]!=0.0)
                                a2i=drhoa0j*sij/rho[icomp_rho+31];
                            if(rho[jcomp_rho+31]!=0.0)
                                a2j=drhoa0i*sij/rho[jcomp_rho+31];
                            if(rho[icomp_rho+32]!=0.0)
                                a3i=drhoa0j*sij/rho[icomp_rho+32];
                            if(rho[jcomp_rho+32]!=0.0)
                                a3j=drhoa0i*sij/rho[jcomp_rho+32];
                            
                            dt1dr1=a1i*(t1_meam[jtype]-t1i*t1_meam[jtype]*t1_meam[jtype]);
                            dt1dr2=a1j*(t1_meam[itype]-t1j*t1_meam[itype]*t1_meam[itype]);
                            dt2dr1=a2i*(t2_meam[jtype]-t2i*t2_meam[jtype]*t2_meam[jtype]);
                            dt2dr2=a2j*(t2_meam[itype]-t2j*t2_meam[itype]*t2_meam[itype]);
                            dt3dr1=a3i*(t3_meam[jtype]-t3i*t3_meam[jtype]*t3_meam[jtype]);
                            dt3dr2=a3j*(t3_meam[itype]-t3j*t3_meam[itype]*t3_meam[itype]);
                        }
                        else if(ialloy==2)
                        {
                            dt1dr1=0.0;
                            dt1dr2=0.0;
                            dt2dr1=0.0;
                            dt2dr2=0.0;
                            dt3dr1=0.0;
                            dt3dr2=0.0;
                        }
                        else
                        {
                            ai=0.0;
                            if(rho[icomp_rho]!=0.0)
                                ai=drhoa0j*sij/rho[icomp_rho];
                            aj=0.0;
                            if(rho[jcomp_rho]!=0.0)
                                aj=drhoa0i*sij/rho[jcomp_rho];
                            
                            dt1dr1=ai*(t1_meam[jtype]-t1i);
                            dt1dr2=aj*(t1_meam[itype]-t1j);
                            dt2dr1=ai*(t2_meam[jtype]-t2i);
                            dt2dr2=aj*(t2_meam[itype]-t2j);
                            dt3dr1=ai*(t3_meam[jtype]-t3i);
                            dt3dr2=aj*(t3_meam[itype]-t3j);
                            
                        }
                        
                        get_shpfcn(si,lattice[itype][itype]);
                        get_shpfcn(sj,lattice[jtype][jtype]);
                        
                        
                        drhodr1=rho[icomp_rho+34]*drho0dr1+rho[icomp_rho+35]
                        *(dt1dr1*rho[icomp_rho+1]+t1i*drho1dr1
                          +dt2dr1*rho[icomp_rho+2]+t2i*drho2dr1
                          +dt3dr1*rho[icomp_rho+3]+t3i*drho3dr1)
                        -rho[icomp_rho+36]*
                        (si[0]*dt1dr1+si[1]*dt2dr1+si[2]*dt3dr1);
                        
                        drhodr2=rho[jcomp_rho+34]*drho0dr2+rho[jcomp_rho+35]
                        *(dt1dr2*rho[jcomp_rho+1]+t1j*drho1dr2
                          +dt2dr2*rho[jcomp_rho+2]+t2j*drho2dr2
                          +dt3dr2*rho[jcomp_rho+3]+t3j*drho3dr2)
                        -rho[jcomp_rho+36]*
                        (sj[0]*dt1dr2+sj[1]*dt2dr2+sj[2]*dt3dr2);
                        
                        for(int i=0;i<3;i++)
                        {
                            drhodrm1[i]=rho[icomp_rho+35]*(t1i*drho1drm1[i]
                                                           +t2i*drho2drm1[i]+t3i*drho3drm1[i]);
                            drhodrm2[i]=rho[jcomp_rho+35]*(t1j*drho1drm2[i]
                                                           +t2j*drho2drm2[i]+t3j*drho3drm2[i]);
                        }
                        
                        
                        if(dscrfcn[istart]!=0.0)
                        {
                            drho0ds1=rhoa0j;
                            drho0ds2=rhoa0i;
                            a1=2.0/rij;
                            drho1ds1=a1*rhoa1j*arg1i1;
                            drho1ds2=a1*rhoa1i*arg1j1;
                            a2=2.0/rij2;
                            
                            drho2ds1=a2*rhoa2j*arg1i2
                            -2.0/3.0*rho[icomp_rho+13]*rhoa2j;
                            drho2ds2=a2*rhoa2i*arg1j2
                            -2.0/3.0*rho[jcomp_rho+13]*rhoa2i;
                            
                            a3=2.0/rij3;
                            a3a=1.2/rij;
                            
                            drho3ds1=a3*rhoa3j*arg1i3-a3a*rhoa3j*arg3i3;
                            drho3ds2=a3*rhoa3i*arg1j3-a3a*rhoa3i*arg3j3;
                            
                            if(ialloy==1)
                            {
                                a1i=0.0;
                                a1j=0.0;
                                a2i=0.0;
                                a2j=0.0;
                                a3i=0.0;
                                a3j=0.0;
                                if(rho[icomp_rho+30]!=0.0)
                                    a1i=drhoa0j/rho[icomp_rho+30];
                                if(rho[jcomp_rho+30]!=0.0)
                                    a1j=drhoa0i/rho[jcomp_rho+30];
                                if(rho[icomp_rho+31]!=0.0)
                                    a2i=drhoa0j/rho[icomp_rho+31];
                                if(rho[jcomp_rho+31]!=0.0)
                                    a2j=drhoa0i/rho[jcomp_rho+31];
                                if(rho[icomp_rho+32]!=0.0)
                                    a3i=drhoa0j/rho[icomp_rho+32];
                                if(rho[jcomp_rho+32]!=0.0)
                                    a3j=drhoa0i/rho[jcomp_rho+32];
                                
                                dt1ds1=a1i*(t1_meam[jtype]-t1i*t1_meam[jtype]*t1_meam[jtype]);
                                dt1ds2=a1j*(t1_meam[itype]-t1j*t1_meam[itype]*t1_meam[itype]);
                                dt2ds1=a2i*(t2_meam[jtype]-t2i*t2_meam[jtype]*t1_meam[jtype]);
                                dt2ds2=a2j*(t2_meam[itype]-t2j*t2_meam[itype]*t1_meam[itype]);
                                dt3ds1=a3i*(t3_meam[jtype]-t3i*t3_meam[jtype]*t1_meam[jtype]);
                                dt3ds2=a3j*(t3_meam[itype]-t3j*t3_meam[itype]*t1_meam[itype]);
                                
                            }
                            else if(ialloy==2)
                            {
                                dt1ds1=0.0;
                                dt1ds2=0.0;
                                dt2ds1=0.0;
                                dt2ds2=0.0;
                                dt3ds1=0.0;
                                dt3ds2=0.0;
                            }
                            else
                            {
                                ai=0.0;
                                if(rho[icomp_rho]!=0.0)
                                    ai=rhoa0j/rho[icomp_rho];
                                aj=0.0;
                                if(rho[jcomp_rho]!=0.0)
                                    aj=rhoa0i/rho[jcomp_rho];
                                
                                dt1ds1=ai*(t1_meam[jtype]-t1i);
                                dt1ds2=aj*(t1_meam[itype]-t1j);
                                dt2ds1=ai*(t2_meam[jtype]-t2i);
                                dt2ds2=aj*(t2_meam[itype]-t2j);
                                dt3ds1=ai*(t3_meam[jtype]-t3i);
                                dt3ds2=aj*(t3_meam[itype]-t3j);
                            }
                            
                            
                            drhods1=rho[icomp_rho+34]*drho0ds1+rho[icomp_rho+35]
                            *(dt1ds1*rho[icomp_rho+1]+t1i*drho1ds1
                              +dt2ds1*rho[icomp_rho+2]+t2i*drho2ds1
                              +dt3ds1*rho[icomp_rho+3]+t3i*drho3ds1)
                            -rho[icomp_rho+36]*
                            (si[0]*dt1ds1+si[1]*dt2ds1+si[2]*dt3ds1);
                            
                            drhods2=rho[jcomp_rho+34]*drho0ds2+rho[jcomp_rho+35]
                            *(dt1ds2*rho[jcomp_rho+1]+t1j*drho1ds2
                              +dt2ds2*rho[jcomp_rho+2]+t2j*drho2ds2
                              +dt3ds2*rho[jcomp_rho+3]+t3j*drho3ds2)
                            -rho[jcomp_rho+36]*
                            (sj[0]*dt1ds2+sj[1]*dt2ds2+sj[2]*dt3ds2);
                        }
                        
                        
                        dUdrij=phip*sij+rho[icomp_rho+37]*drhodr1
                        +rho[jcomp_rho+37]*drhodr2;
                        
                        dUdsij=0.0;
                        if(dscrfcn[istart]!=0.0)
                            dUdsij=phi+rho[icomp_rho+37]*drhods1
                            +rho[jcomp_rho+37]*drhods2;
                        
                        for(int i=0;i<3;i++)
                            dUdrijm[i]=rho[icomp_rho+37]*drhodrm1[i]
                            +rho[jcomp_rho+37]*drhodrm2[i];
                        
                        force=dUdrij*recip+dUdsij*dscrfcn[istart];
                        
                        for(int i=0;i<3;i++)
                        {
                            forcem=delij[i]*force+dUdrijm[i];
                            fvec[icomp+i]+=forcem;
                            if(jatm<atoms->natms)
                                fvec[jcomp+i]-=forcem;
                        }
                        
                        /*
                         stress and all the other shit
                         
                         remember to fix it
                         */
                        if(st_clc)
                        {
                            fi[0]=delij[0]*force+dUdrijm[0];
                            fi[1]=delij[1]*force+dUdrijm[1];
                            fi[2]=delij[2]*force+dUdrijm[2];
                            
                            if(jatm<atoms->natms)
                            {
                                nrgy_strss_lcl[1]+=delij[0]*fi[0];
                                nrgy_strss_lcl[2]+=delij[1]*fi[1];
                                nrgy_strss_lcl[3]+=delij[2]*fi[2];
                                nrgy_strss_lcl[4]+=0.5*(delij[1]*fi[2]+delij[2]*fi[1]);
                                nrgy_strss_lcl[5]+=0.5*(delij[0]*fi[2]+delij[2]*fi[0]);
                                nrgy_strss_lcl[6]+=0.5*(delij[0]*fi[1]+delij[1]*fi[0]);
                            }
                            else
                            {
                                nrgy_strss_lcl[1]+=0.5*delij[0]*fi[0];
                                nrgy_strss_lcl[2]+=0.5*delij[1]*fi[1];
                                nrgy_strss_lcl[3]+=0.5*delij[2]*fi[2];
                                nrgy_strss_lcl[4]+=0.25*(delij[1]*fi[2]+delij[2]*fi[1]);
                                nrgy_strss_lcl[5]+=0.25*(delij[0]*fi[2]+delij[2]*fi[0]);
                                nrgy_strss_lcl[6]+=0.25*(delij[0]*fi[1]+delij[1]*fi[0]);
                            }
                        }
                        
                        if(sij!=0.0 && sij!=1.0)
                        {
                            rbound=rij2*ebound_meam[itype][jtype];
                            
                            for(int kn=0;kn<my_list_size;kn++)
                            {
                                katm=my_list[kn];
                                
                                if(katm!=jatm)
                                {
                                    dsij1=dsij2=0.0;
                                    ktype=type[katm];
                                    kcomp=3*katm;
                                    
                                    xjk=x[jcomp]-x[kcomp];
                                    yjk=x[jcomp+1]-x[kcomp+1];
                                    zjk=x[jcomp+2]-x[kcomp+2];
                                    rjk2=xjk*xjk+yjk*yjk+zjk*zjk;
                                    
                                    if(rjk2<rbound)
                                    {
                                        xki=x[kcomp]-x[icomp];
                                        yki=x[kcomp+1]-x[icomp+1];
                                        zki=x[kcomp+2]-x[icomp+2];
                                        rki2=xki*xki+yki*yki+zki*zki;
                                        
                                        if(rki2<rbound)
                                        {
                                            x_ki=rki2/rij2;
                                            x_jk=rjk2/rij2;
                                            a=1.0-(x_ki-x_jk)*(x_ki-x_jk);
                                            if(a!=0.0)
                                            {
                                                cikj=(2.0*(x_ki+x_jk)+a-2.0)/a;
                                                cmax=c_max[itype][jtype][ktype];
                                                cmin=c_min[itype][jtype][ktype];
                                                if(cikj<cmax&&cikj>cmin)
                                                {
                                                    cikj=(cikj-cmin)/(cmax-cmin);
                                                    dfcut(cikj,sikj,dfc);
                                                    dCfunc2(rij2,rki2,rjk2,dcikj1,dcikj2);
                                                    a=(sij/(cmax-cmin))*(dfc/sikj);
                                                    dsij1=a*dcikj1;
                                                    dsij2=a*dcikj2;
                                                }
                                            }
                                        }
                                        
                                        if(dsij1!=0.0 || dsij2!=0.0)
                                        {
                                            force1=dUdsij*dsij1;
                                            force2=dUdsij*dsij2;
                                            delik[0]=xki;
                                            deljk[0]=-xjk;
                                            delik[1]=yki;
                                            deljk[1]=-yjk;
                                            delik[2]=zki;
                                            deljk[2]=-zjk;
                                            for(int i=0;i<3;i++)
                                                fvec[icomp+i]+=force1*delik[i];
                                            if(jatm<atoms->natms)
                                                for(int i=0;i<3;i++)
                                                    fvec[jcomp+i]+=force2*deljk[i];
                                            
                                            if(katm<atoms->natms)
                                                for(int i=0;i<3;i++)
                                                    fvec[kcomp+i]-=force1*delik[i]+force2*deljk[i];
                                            
                                            /*
                                             the stress contrbution crap
                                             remember to set it
                                             */
                                            if(st_clc)
                                            {
                                                fi[0]=force1*delik[0];
                                                fi[1]=force1*delik[1];
                                                fi[2]=force1*delik[2];
                                                fj[0]=force2*deljk[0];
                                                fj[1]=force2*deljk[1];
                                                fj[2]=force2*deljk[2];
                                                
                                                v[0]=-third*(delik[0]*fi[0]+deljk[0]*fj[0]);
                                                v[1]=-third*(delik[1]*fi[1]+deljk[1]*fj[1]);
                                                v[2]=-third*(delik[2]*fi[2]+deljk[2]*fj[2]);
                                                v[3]=-sixth*(delik[1]*fi[2]+deljk[1]*fj[2]
                                                +delik[2]*fi[1]+deljk[2]*fj[1]);
                                                v[4]=-sixth*(delik[0]*fi[2]+deljk[0]*fj[2]
                                                +delik[2]*fi[0]+deljk[2]*fj[0]);
                                                v[5]=-sixth*(delik[0]*fi[1]+deljk[0]*fj[1]
                                                +delik[1]*fi[0]+deljk[1]*fj[0]);
                                                
                                                for(int i=0;i<6;i++)
                                                    nrgy_strss_lcl[i+1]-=v[i];
                                                if(jatm<atoms->natms)
                                                    for(int i=0;i<6;i++)
                                                        nrgy_strss_lcl[i+1]-=v[i];
                                                if(katm<atoms->natms)
                                                    for(int i=0;i<6;i++)
                                                        nrgy_strss_lcl[i+1]-=v[i];
                                                
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                istart++;
            }
        }
    }
    
    if(st_clc)
    {
        for(int i=0;i<7;i++)
            nrgy_strss[i]=0.0;
        
        MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,7,MPI_TYPE0,MPI_SUM,world);
    }
    else
    {
        nrgy_strss[0]=0.0;
        MPI_Allreduce(nrgy_strss_lcl,nrgy_strss,1,MPI_TYPE0,MPI_SUM,world);
    }
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
type0 ForceField_meam::energy_calc()
{
    
    type0 en=0.0;
    type0 en_tot=0.0;
    
    int* my_list;
    int my_list_size;
    
    int jatm,katm;
    
    type0 xij,xjk,xki,yij,yjk,yki,zij,zjk,zki;
    type0 rij2,rjk2,rki2;
    type0 rij;
    type0 delij[3],s[3];
    
    int icomp,jcomp,kcomp,itype,jtype,ktype;
    type0 coef1,coef2;
    int icomp_rho,jcomp_rho;
    
    
    type0 fcij,dfcij,sij,dsij;
    type0 rnorm,rbound,x_ki,x_jk,a;
    type0 cikj,sikj=0.0,cmax,cmin,dfikj,dcikj;
    type0 ai,aj,ro0i,ro0j,rhoa0i,rhoa0j,rhoa1i,rhoa1j
    ,rhoa2i,rhoa2j,rhoa3i,rhoa3j;
    type0 Z;
    type0 G,dG,Gbar,dGbar=0.0,gam,rho_bkgd,rhob,denom,B;
    int iArho1_comp,jArho1_comp,iArho2_comp
    ,jArho2_comp,iArho3_comp,jArho3_comp
    ,iArhob3_comp,jArhob3_comp;
    type0 A1i,A2i,A3i,A1j,A2j,A3j,B1i,B1j;
    
    
    type0* x=mapp->x->begin();
    type0* rho=rho_ptr->begin();
    type0* rho_vec=rho_vec_ptr->begin();
    md_type* type=mapp->type->begin();
    
    
    for(int i=0;i<atoms->natms;i++)
        rho_vec[i]=0.0;
    for(int i=0;i<(atoms->natms+atoms->natms_ph)*rho_dim;i++)
        rho[i]=0.0;
    
    type0 phi;
    int kn,kk,curs;
    type0 pp;
    type0* coef;
    
    for(int iatm=0;iatm<atoms->natms;iatm++)
    {
        my_list_size=neighbor->neighbor_list_size[iatm];
        my_list=neighbor->neighbor_list[iatm];
        icomp=iatm*3;
        icomp_rho=iatm*rho_dim;
        itype=type[iatm];
        for(int jn=0;jn<my_list_size;jn++)
        {
            jatm=my_list[jn];
            if(jatm>iatm)
            {
                fcij=sij=dsij=0.0;
                jtype=type[jatm];
                jcomp=jatm*3;
                jcomp_rho=jatm*rho_dim;
                
                xij=x[icomp]-x[jcomp];
                yij=x[icomp+1]-x[jcomp+1];
                zij=x[icomp+2]-x[jcomp+2];
                delij[0]=-xij;
                delij[1]=-yij;
                delij[2]=-zij;
                rij2=xij*xij+yij*yij+zij*zij;
                rij=sqrt(rij2);
                
                if(rij<rc_meam)
                {
                    rnorm=(rc_meam-rij)*delr_meam_inv;
                    sij=1.0;
                    dsij=0.0;
                    rbound = ebound_meam[itype][jtype]*rij2;
                    
                    kn=0;
                    dfcut(rnorm,fcij,dfcij);
                    dfcij*=delr_meam_inv;
                    while(kn<my_list_size && sij!=0.0)
                    {
                        katm=my_list[kn];
                        
                        if(katm!=jatm)
                        {
                            ktype=type[katm];
                            kcomp=3*katm;
                            
                            xjk=x[jcomp]-x[kcomp];
                            yjk=x[jcomp+1]-x[kcomp+1];
                            zjk=x[jcomp+2]-x[kcomp+2];
                            rjk2=xjk*xjk+yjk*yjk+zjk*zjk;
                            
                            if(rjk2<rbound)
                            {
                                xki=x[kcomp]-x[icomp];
                                yki=x[kcomp+1]-x[icomp+1];
                                zki=x[kcomp+2]-x[icomp+2];
                                rki2=xki*xki+yki*yki+zki*zki;
                                
                                if(rki2<rbound)
                                {
                                    x_ki=rki2/rij2;
                                    x_jk=rjk2/rij2;
                                    a=1.0-(x_ki-x_jk)*(x_ki-x_jk);
                                    if(a>0)
                                    {
                                        cikj=(2.0*(x_ki+x_jk)+a-2.0)/a;
                                        cmax=c_max[itype][jtype][ktype];
                                        cmin=c_min[itype][jtype][ktype];
                                        if(cikj<cmax&&cikj>cmin)
                                        {
                                            cikj=(cikj-cmin)/(cmax-cmin);
                                            dfcut(cikj,sikj,dfikj);
                                            sij*=sikj;
                                            coef1=dfikj/((cmax-cmin)*sikj);
                                            dCfunc(rij2,rki2,rjk2,dcikj);
                                            dsij+=coef1*dcikj;
                                        }
                                        else if(cikj<cmin)
                                        {
                                            sij=0.0;
                                            dsij=0.0;
                                        }
                                    }
                                }
                            }
                        }
                        kn++;
                    }
                    coef1=sij*fcij;
                    coef2=sij*dfcij/rij;
                    dsij*=coef1;
                    dsij-=coef2;
                    if(sij==1.0 || sij==0.0)
                        dsij=0.0;
                    if(sij*fcij==1.0 || sij*fcij==0.0)
                        dsij=0.0;
                }
                else
                {
                    fcij=0.0;
                    dfcij=0.0;
                    sij=0.0;
                    dsij=0.0;
                }
                
                
                
                if(sij*fcij!=0.0)
                {
                    sij*=fcij;
                    
                    curs=COMP(itype,jtype);
                    pp=rij*dr_inv;
                    kk=static_cast<int>(pp);
                    kk=MIN(kk,nr-2);
                    pp-=static_cast<type0>(kk);
                    pp=MIN(pp,1.0);
                    coef=phirar[curs][kk];
                    phi=((coef[3]*pp+coef[2])*pp+coef[1])*pp+coef[0];
                    
                    
                    if(jatm<atoms->natms)
                        en+=phi*sij;
                    else
                        en+=0.5*phi*sij;
                    
                    ai=rij/re_meam[itype][itype]-1.0;
                    aj=rij/re_meam[jtype][jtype]-1.0;
                    ro0i=rho0_meam[itype];
                    ro0j=rho0_meam[jtype];
                    
                    rhoa0i=sij*ro0i*exp(-beta0_meam[itype]*ai);
                    rhoa1i=sij*ro0i*exp(-beta1_meam[itype]*ai);
                    rhoa2i=sij*ro0i*exp(-beta2_meam[itype]*ai);
                    rhoa3i=sij*ro0i*exp(-beta3_meam[itype]*ai);
                    
                    rhoa0j=sij*ro0j*exp(-beta0_meam[jtype]*aj);
                    rhoa1j=sij*ro0j*exp(-beta1_meam[jtype]*aj);
                    rhoa2j=sij*ro0j*exp(-beta2_meam[jtype]*aj);
                    rhoa3j=sij*ro0j*exp(-beta3_meam[jtype]*aj);
                    
                    
                    if(ialloy==1)
                    {
                        rhoa1i*=t1_meam[itype];
                        rhoa2i*=t2_meam[itype];
                        rhoa3i*=t3_meam[itype];
                        
                        rhoa1j*=t1_meam[jtype];
                        rhoa2j*=t2_meam[jtype];
                        rhoa3j*=t3_meam[jtype];
                    }
                    
                    rho[icomp_rho]+=rhoa0j;
                    if(jatm<atoms->natms)
                        rho[jcomp_rho]+=rhoa0i;
                    
                    
                    if(ialloy!=2)
                    {
                        
                        rho[icomp_rho+27]+=t1_meam[jtype]*rhoa0j;
                        rho[icomp_rho+28]+=t2_meam[jtype]*rhoa0j;
                        rho[icomp_rho+29]+=t3_meam[jtype]*rhoa0j;
                        if(jatm<atoms->natms)
                        {
                            rho[jcomp_rho+27]+=t1_meam[itype]*rhoa0i;
                            rho[jcomp_rho+28]+=t2_meam[itype]*rhoa0i;
                            rho[jcomp_rho+29]+=t3_meam[itype]*rhoa0i;
                        }
                        
                    }
                    if(ialloy==1)
                    {
                        
                        rho[icomp_rho+30]+=t1_meam[jtype]*t1_meam[jtype]*rhoa0j;
                        rho[icomp_rho+31]+=t2_meam[jtype]*t2_meam[jtype]*rhoa0j;
                        rho[icomp_rho+32]+=t3_meam[jtype]*t3_meam[jtype]*rhoa0j;
                        if(jatm<atoms->natms)
                        {
                            rho[jcomp_rho+30]+=t1_meam[itype]*t1_meam[itype]*rhoa0i;
                            rho[jcomp_rho+31]+=t2_meam[itype]*t2_meam[itype]*rhoa0i;
                            rho[jcomp_rho+32]+=t3_meam[itype]*t3_meam[itype]*rhoa0i;
                        }
                    }
                    
                    rho[icomp_rho+13]+=rhoa2j;
                    if(jatm<atoms->natms)
                        rho[jcomp_rho+13]+=rhoa2i;
                    
                    
                    iArho1_comp=icomp_rho+4;
                    iArho2_comp=icomp_rho+7;
                    iArho3_comp=icomp_rho+14;
                    iArhob3_comp=icomp_rho+24;
                    
                    jArho1_comp=jcomp_rho+4;
                    jArho2_comp=jcomp_rho+7;
                    jArho3_comp=jcomp_rho+14;
                    jArhob3_comp=jcomp_rho+24;
                    
                    A1j=rhoa1j/rij;
                    B1j=rhoa3j/rij;
                    A2j=rhoa2j/rij2;
                    A3j=rhoa3j/(rij2*rij);
                    
                    for(int i=0;i<3;i++)
                    {
                        rho[iArho1_comp++]+=delij[i]*A1j;
                        rho[iArhob3_comp++]+=delij[i]*B1j;
                        for(int j=i;j<3;j++)
                        {
                            rho[iArho2_comp++]+=delij[i]*delij[j]*A2j;
                            for(int k=j;k<3;k++)
                            {
                                rho[iArho3_comp++]+=delij[i]*delij[j]*delij[k]*A3j;
                            }
                        }
                    }
                    
                    
                    if(jatm<atoms->natms)
                    {
                        
                        A1i=rhoa1i/rij;
                        B1i=rhoa3i/rij;
                        A2i=rhoa2i/rij2;
                        A3i=rhoa3i/(rij2*rij);
                        
                        for(int i=0;i<3;i++)
                        {
                            rho[jArho1_comp++]-=delij[i]*A1i;
                            rho[jArhob3_comp++]-=delij[i]*B1i;
                            for(int j=i;j<3;j++)
                            {
                                rho[jArho2_comp++]+=delij[i]*delij[j]*A2i;
                                for(int k=j;k<3;k++)
                                {
                                    rho[jArho3_comp++]-=delij[i]*delij[j]*delij[k]*A3i;
                                }
                            }
                        }
                    }
                    
                }
                
            }
        }
        
        
        
        rho[icomp_rho+1]=0.0;
        rho[icomp_rho+2]=-rho[icomp_rho+13]*rho[icomp_rho+13]/3.0;
        rho[icomp_rho+3]=0.0;
        for(int i=0;i<3;i++)
        {
            rho[icomp_rho+1]+=rho[icomp_rho+4+i]*rho[icomp_rho+4+i];
            rho[icomp_rho+3]-=0.6*rho[icomp_rho+24+i]*rho[icomp_rho+24+i];
        }
        
        for(int i=0;i<6;i++)
        {
            rho[icomp_rho+2]+=v2d[i]*rho[icomp_rho+7+i]*rho[icomp_rho+7+i];
        }
        
        for(int i=0;i<10;i++)
        {
            rho[icomp_rho+3]+=v3d[i]*rho[icomp_rho+14+i]*rho[icomp_rho+14+i];
        }
        
        if(rho[icomp_rho]>0.0)
        {
            if(ialloy==1)
            {
                rho[icomp_rho+27]=rho[icomp_rho+27]/rho[icomp_rho+30];
                rho[icomp_rho+28]=rho[icomp_rho+28]/rho[icomp_rho+31];
                rho[icomp_rho+29]=rho[icomp_rho+29]/rho[icomp_rho+32];
            }
            else if(ialloy==2)
            {
                rho[icomp_rho+27]=t1_meam[itype];
                rho[icomp_rho+28]=t2_meam[itype];
                rho[icomp_rho+29]=t3_meam[itype];
            }
            else
            {
                rho[icomp_rho+27]=rho[icomp_rho+27]/rho[icomp_rho];
                rho[icomp_rho+28]=rho[icomp_rho+28]/rho[icomp_rho];
                rho[icomp_rho+29]=rho[icomp_rho+29]/rho[icomp_rho];
            }
        }
        
        
        rho[icomp_rho+33]=rho[icomp_rho+27]*rho[icomp_rho+1]
        +rho[icomp_rho+28]*rho[icomp_rho+2]+rho[icomp_rho+29]*rho[icomp_rho+3];
        
        if(rho[icomp_rho]>0.0)
            rho[icomp_rho+33]=rho[icomp_rho+33]/(rho[icomp_rho]*rho[icomp_rho]);
        
        
        
        Z=Z_meam[itype];
        
        G_gam(rho[icomp_rho+33],ibar_meam[itype],gsmooth_factor,G);
        
        if(ibar_meam[itype]<=0)
            Gbar=1.0;
        else
        {
            get_shpfcn(s,lattice[itype][itype]);
            if(mix_ref_t==1)
                gam=(rho[icomp_rho+27]*s[0]+rho[icomp_rho+28]*s[1]+rho[icomp_rho+29]*s[2])/(Z*Z);
            else
                gam=(t1_meam[itype]*s[0]+t2_meam[itype]*s[1]+t3_meam[itype]*s[2])/(Z*Z);
            
            G_gam(gam,ibar_meam[itype],gsmooth_factor,Gbar);
        }
        
        rho_vec[iatm]=rho[icomp_rho]*G;
        
        if(mix_ref_t==1)
        {
            if(ibar_meam[itype]<=0.0)
            {
                Gbar=1.0;
                dGbar=0.0;
            }
            else
            {
                get_shpfcn(s,lattice[itype][itype]);
                gam=(rho[icomp_rho+27]*s[0]+rho[icomp_rho+28]*s[1]+rho[icomp_rho+29]*s[2])/(Z*Z);
                dG_gam(gam,ibar_meam[itype],gsmooth_factor,Gbar,dGbar);
            }
            rho_bkgd=rho0_meam[itype]*Z*Gbar;
        }
        else
        {
            if(bkgd_dyn==1)
            {
                rho_bkgd=rho0_meam[itype]*Z;
            }
            else
                rho_bkgd=rho_ref_meam[itype];
        }
        
        rhob=rho_vec[iatm]/rho_bkgd;
        denom=1.0/rho_bkgd;
        dG_gam(rho[icomp_rho+33],ibar_meam[itype],gsmooth_factor,G,dG);
        rho[icomp_rho+34]=(G-2.0*dG*rho[icomp_rho+33])*denom;
        
        if(rho[icomp_rho]!=0)
            rho[icomp_rho+35]=(dG/rho[icomp_rho])*denom;
        else
            rho[icomp_rho+35]=0.0;
        
        if(mix_ref_t==1)
            rho[icomp_rho+36]=rho[icomp_rho]*G*dGbar*denom/(Gbar*Z*Z);
        else
            rho[icomp_rho+36]=0.0;
        
        B=A_meam[itype]*Ec_meam[itype][itype];
        
        if(rhob!=0.0)
        {
            if(emb_lin_neg==1 && rhob<=0.0)
            {
                rho[icomp_rho+37]=-B;
                en-=B*rhob;
            }
            else
            {
                rho[icomp_rho+37]=B*(log(rhob)+1.0);
                en+=B*rhob*log(rhob);
            }
        }
        else
        {
            if(emb_lin_neg==1)
                rho[icomp_rho+37]=-B;
            else
                rho[icomp_rho+37]=B;
        }
        
    }
    
    
    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    return en_tot;
}
/*--------------------------------------------
 coefficients, file readings and stuff
 --------------------------------------------*/
void ForceField_meam::coef(int narg,char** args)
{
    if(narg!=3)
        error->abort("wrong coeff command for MEAM FF");
    cut_off_alloc();
    allocate();
    read_global(args[1]);
    read_local(args[2]);
    setup();
}
/*--------------------------------------------
 coefficients, file readings and stuff
 --------------------------------------------*/
void ForceField_meam::setup()
{
    for(int itype=0;itype<no_types;itype++)
        t1_meam[itype]+=static_cast<type0>(augt1)*0.6*t3_meam[itype];
    
    alloy_params();
    compute_reference_density();
    dr=1.1*rc_meam/static_cast<type0>(nr);
    dr_inv=1.0/dr;
    compute_pair_meam();
    for(int i=0;i<no_types;i++)
        for(int j=0;j<no_types;j++)
        {
            cut[i][j]=rc_meam;
            cut_sq[i][j]=rc_meam*rc_meam;
        }
}
/*--------------------------------------------
 reset the vectors
 --------------------------------------------*/
void ForceField_meam::reset()
{
    type0* rho=rho_ptr->begin();
    type0* rho_vec=rho_vec_ptr->begin();
    
    for(int i=0;i<atoms->natms;i++)
        rho_vec[i]=0.0;
    for(int i=0;i<(atoms->natms+atoms->natms_ph)*rho_dim;i++)
        rho[i]=0.0;
    
    int no_pairs=neighbor->no_pairs;
    if(no_pairs>max_pairs)
    {
        if(max_pairs)
        {
            delete [] scrfcn;
            delete [] dscrfcn;
            delete [] fcpair;
        }
        CREATE1D(scrfcn,no_pairs);
        CREATE1D(dscrfcn,no_pairs);
        CREATE1D(fcpair,no_pairs);
        max_pairs=no_pairs;
    }
}
/*--------------------------------------------
 read the global file and broad cast the
 parameters
 --------------------------------------------*/
void ForceField_meam::read_global(char* file_name)
{
    long byte_sz;
    char* buff;
    if(atoms->my_p==0)
    {
        std::fstream file(file_name,std::ios::in);
        file.seekg(0,std::ios::end);
        byte_sz=file.tellg();
    }
    
    MPI_Bcast(&byte_sz,1,MPI_LONG,0,world);
    byte_sz++;
    CREATE1D(buff,byte_sz);
    
    FILE* fp;
    mapp->open_file(fp,file_name,"r");
    
    if(atoms->my_p==0)
    {
        fread (buff,1,byte_sz-1,fp);
        buff[byte_sz-1]='\0';
        fclose(fp);
    }
    
    int sz_int=static_cast<int>(byte_sz);
    long sz_long=static_cast<long>(sz_int);
    long byte_sz_=byte_sz;
    char* buff_=buff;
    while(byte_sz_)
    {
        MPI_Bcast(buff_,sz_int,MPI_CHAR,0, world);
        byte_sz_-=sz_long;
        buff_+=sz_int;
        if(byte_sz_<sz_long)
            sz_int=static_cast<int>(byte_sz_);
    }
    
    char* pch;
    char* new_buff;
    char* new_buff_;
    char* st;
    CREATE1D(new_buff,byte_sz);
    
    buff_=buff;
    new_buff_=new_buff;
    st=buff;
    pch=strchr(buff_,'#');
    while (pch!=NULL)
    {
        memcpy(new_buff_,st,pch-st);
        new_buff_+=pch-st;
        pch=strchr(pch+1,'\n');
        st=pch;
        
        pch=strchr(pch+1,'#');
    }
    memcpy(new_buff_,st,buff+byte_sz-st);
    std::swap(new_buff,buff);
    
    buff_=buff;
    new_buff_=new_buff;
    st=buff;
    
    pch=strchr(buff_,'\\');
    while (pch!=NULL)
    {
        memcpy(new_buff_,st,pch-st);
        new_buff_+=pch-st;
        pch=strchr(pch+1,'\n');
        st=pch+1;
        
        pch=strchr(pch+1,'\\');
    }
    memcpy(new_buff_,st,buff+byte_sz-st);
    
    
    int buff_length=static_cast<int>(strlen(new_buff)+1);
    if(byte_sz)
        delete [] buff;
    CREATE1D(buff,buff_length);
    memcpy(buff,new_buff,buff_length*sizeof(char));
    if(byte_sz)
        delete [] new_buff;
    
    char** args=NULL;
    int nargs=0;
    pch=strtok(buff,"' \n\t\r");
    while (pch!=NULL)
    {
        GROW(args,nargs,nargs+1);
        args[nargs]=pch;
        nargs++;
        pch=strtok(NULL,"' \n\t\r");
    }
    
    if(nargs%19!=0)
        error->abort("number of enteries in file %s must be a multiple of 19",file_name);
    int nelem=nargs/19;
    
    for(int i=0;i<no_types;i++)
        type_ref[i]=-1;
    int no_types_=0;
    
    for(int ielem=0;ielem<nelem;ielem++)
    {
        int itype=atom_types->find_type(args[ielem*19]);
        if(itype!=-1)
        {
            int itype_=0;
            while (itype_<no_types_&& type_ref[itype_]!=itype)
                itype_++;

            if(itype_!=no_types_)
                continue;

            type_ref[no_types_++]=itype;
            
            char** args_=args+ielem*19;
            
            if(strcmp(args_[1],"fcc")==0) lattice[itype][itype]=FCC;
            else if(strcmp(args_[1],"bcc")==0) lattice[itype][itype]=BCC;
            else if(strcmp(args_[1],"hcp")==0) lattice[itype][itype]=HCP;
            else if(strcmp(args_[1],"dim")==0) lattice[itype][itype]=DIM;
            else if(strcmp(args_[1],"dia")==0) lattice[itype][itype]=DIAMOND;
            else error->abort("unknown lattice");
            
            Z_meam[itype]=atof(args_[2]);
            ielt_meam[itype]=atoi(args_[3]);
            
            alpha_meam[itype][itype]=atof(args_[5]);
            beta0_meam[itype]=atof(args_[6]);
            beta1_meam[itype]=atof(args_[7]);
            beta2_meam[itype]=atof(args_[8]);
            beta3_meam[itype]=atof(args_[9]);
            
            if(lattice[itype][itype]==FCC)
                re_meam[itype][itype]=atof(args_[10])/sqrt(2.0);
            else if(lattice[itype][itype]==BCC)
                re_meam[itype][itype]=atof(args_[10])*0.5*sqrt(3.0);
            else if(lattice[itype][itype]==HCP)
                re_meam[itype][itype]=atof(args_[10]);
            else if(lattice[itype][itype]==DIM)
                re_meam[itype][itype]=atof(args_[10]);
            else if(lattice[itype][itype]==DIAMOND)
                re_meam[itype][itype]=atof(args_[10])*0.25*sqrt(3.0);
            
            Ec_meam[itype][itype]=atof(args_[11]);
            A_meam[itype]=atof(args_[12]);
            t0_meam[itype]=atof(args_[13]);
            t1_meam[itype]=atof(args_[14]);
            t2_meam[itype]=atof(args_[15]);
            t3_meam[itype]=atof(args_[16]);
            rho0_meam[itype]=atof(args_[17]);
            ibar_meam[itype]=atoi(args_[18]);

        }
    }
    if(buff_length)
        delete [] buff;
    if(nargs)
        delete [] args;
    
    
    for(int itype=0;itype<no_types;itype++)
        if(type_ref[itype]==-1)
            error->abort("properties of "
            "element %s was not found in "
            "%s file",atom_types->atom_names[itype],file_name);


    for(int i=0;i<no_types;i++)
    {
        for(int j=0;j<no_types;j++)
        {
            for(int k=0;k<no_types;k++)
            {
                c_min[i][j][k]=2.0;
                c_max[i][j][k]=2.8;
            }
            ebound_meam[i][j]=1.96/1.8;
            
            delta_meam[i][j]=0.0;
            attrac_meam[i][j]=0.0;
            repuls_meam[i][j]=0.0;
            
            nn2_meam[i][j]=0;
            zbl_meam[i][j]=1;
        }
    }
    
    emb_lin_neg=0;
    ialloy=0;
    mix_ref_t=0;
    bkgd_dyn=0;
    augt1=1;
    erose_form=0;
    delr_meam_inv=10.0;
    delr_meam=0.1;
    rc_meam=4.0;
    gsmooth_factor=99.0;
}
/*--------------------------------------------
 read the global file and broad cast the
 parameters
 --------------------------------------------*/
void ForceField_meam::read_local(char* file_name)
{
    if(strcmp(file_name,"NULL")==0)
        return;
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    
    mapp->open_file(fp,file_name,"r");
    
    
    int icmp,jcmp,kcmp;
    type0 tmp;
    int tmp_;
    char tmp_str[20];
    while(mapp->read_line(fp,line)!=-1)
    {
        if(mapp->hash_remover(line)==0)
            continue;
        if(0);
        else if(sscanf(line,"gsmooth_factor = %lf",&tmp)==1)
        {
            gsmooth_factor=tmp;
        }
        else if(sscanf(line,"augt1 = %lf",&tmp)==1)
        {
            augt1=tmp;
        }
        else if(sscanf(line,"delr = %lf",&tmp)==1)
        {
            delr_meam=tmp;
        }
        else if(sscanf(line,"rc = %lf",&tmp)==1)
        {
            rc_meam=tmp;
        }
        else if(sscanf(line,"bkgd_dyn = %d",&tmp_)==1)
        {
            bkgd_dyn=tmp_;
        }
        else if(sscanf(line,"emb_lin_neg = %d",&tmp_)==1)
        {
            emb_lin_neg=tmp_;
        }
        else if(sscanf(line,"mixture_ref_t = %d",&tmp_)==1)
        {
            mix_ref_t=tmp_;
        }
        else if(sscanf(line,"erose_form = %d",&tmp_)==1)
        {
            erose_form=tmp_;
        }
        else if(sscanf(line,"ialloy = %d",&tmp_)==1)
        {
            ialloy=tmp_;
        }
        else if(sscanf(line,"rho0(%d) = %lf",&icmp,&tmp)==2)
        {
            icmp--;
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for rho0(%d)",file_name,icmp+1);

            rho0_meam[type_ref[icmp]]=tmp;
        }
        else if(sscanf(line,"Ec(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for Ec(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp<0 || jcmp>no_types-1)
                error->abort("wrong component in %s file for Ec(%d,%d)",file_name,icmp+1,jcmp+1);
            Ec_meam[type_ref[icmp]][type_ref[jcmp]]
            =Ec_meam[type_ref[jcmp]][type_ref[icmp]]=tmp;
        }
        else if(sscanf(line,"alpha(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for alpha(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp<0 || jcmp>no_types-1)
                error->abort("wrong component in %s file for alpha(%d,%d)",file_name,icmp+1,jcmp+1);
            alpha_meam[type_ref[icmp]][type_ref[jcmp]]
            =alpha_meam[type_ref[jcmp]][type_ref[icmp]]=tmp;
        }
        else if(sscanf(line,"re(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for re(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp<0 || jcmp>no_types-1)
                error->abort("wrong component in %s file for re(%d,%d)",file_name,icmp+1,jcmp+1);
            re_meam[type_ref[icmp]][type_ref[jcmp]]
            =re_meam[type_ref[jcmp]][type_ref[icmp]]=tmp;
        }
        else if(sscanf(line,"delta(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for delta(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp<0 || jcmp>no_types-1)
                error->abort("wrong component in %s file for delta(%d,%d)",file_name,icmp+1,jcmp+1);
            delta_meam[type_ref[icmp]][type_ref[jcmp]]=tmp;
        }
        else if(sscanf(line,"attrac(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for attrac(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp<0 || jcmp>no_types-1)
                error->abort("wrong component in %s file for attrac(%d,%d)",file_name,icmp+1,jcmp+1);
            attrac_meam[type_ref[icmp]][type_ref[jcmp]]=tmp;
        }
        else if(sscanf(line,"repuls(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for repuls(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp<0 || jcmp>no_types-1)
                error->abort("wrong component in %s file for repuls(%d,%d)",file_name,icmp+1,jcmp+1);
            repuls_meam[type_ref[icmp]][type_ref[jcmp]]=tmp;
        }
        else if(sscanf(line,"zbl(%d,%d) = %d",&icmp,&jcmp,&tmp_)==3)
        {
            icmp--;
            jcmp--;
            kcmp=MIN(icmp,jcmp);
            jcmp=MAX(icmp,jcmp);
            icmp=kcmp;
            
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for zbl(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp<0 || jcmp>no_types-1)
                error->abort("wrong component in %s file for zbl(%d,%d)",file_name,icmp+1,jcmp+1);
            zbl_meam[type_ref[icmp]][type_ref[jcmp]]=tmp_;
        }
        else if(sscanf(line,"nn2(%d,%d) = %d",&icmp,&jcmp,&tmp_)==3)
        {
            icmp--;
            jcmp--;
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for nn2(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp<0 || jcmp>no_types-1)
                error->abort("wrong component in %s file for nn2(%d,%d)",file_name,icmp+1,jcmp+1);
            nn2_meam[type_ref[icmp]][type_ref[jcmp]]=nn2_meam[type_ref[jcmp]][type_ref[icmp]]=tmp_;
        }
        else if(sscanf(line,"lattce(%d,%d) = %s",&icmp,&jcmp,tmp_str)==3)
        {
            icmp--;
            jcmp--;
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for lattice(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp<0 || jcmp>no_types-1)
                error->abort("wrong component in %s file for lattice(%d,%d)",file_name,icmp+1,jcmp+1);
            
            if(!strcmp(tmp_str,"'fcc'")) lattice[type_ref[icmp]][type_ref[jcmp]]=FCC;
            else if(!strcmp(tmp_str,"'bcc'")) lattice[type_ref[icmp]][type_ref[jcmp]]=BCC;
            else if(!strcmp(tmp_str,"'hcp'")) lattice[type_ref[icmp]][type_ref[jcmp]]=HCP;
            else if(!strcmp(tmp_str,"'dim'")) lattice[type_ref[icmp]][type_ref[jcmp]]=DIM;
            else if(!strcmp(tmp_str,"'dia'")) lattice[type_ref[icmp]][type_ref[jcmp]]=DIAMOND;
            else if(!strcmp(tmp_str,"'b1'")) lattice[type_ref[icmp]][type_ref[jcmp]]=B1;
            else if(!strcmp(tmp_str,"'c11'")) lattice[type_ref[icmp]][type_ref[jcmp]]=C11;
            else if(!strcmp(tmp_str,"'l12'")) lattice[type_ref[icmp]][type_ref[jcmp]]=L12;
            else if(!strcmp(tmp_str,"'b2'")) lattice[type_ref[icmp]][type_ref[jcmp]]=B2;
            else error->abort("wrong component in %s file for lattce(%i,%i) = %s",file_name,icmp+1,jcmp+1,tmp_str);
            
        }
        else if(sscanf(line,"Cmin(%d,%d,%d) = %lf",&icmp,&jcmp,&kcmp,&tmp)==4)
        {
            icmp--;
            jcmp--;
            kcmp--;
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for Cmin(%d,%d,%d)",file_name,icmp+1,jcmp+1,kcmp+1);
            if(jcmp<0 || jcmp>no_types-1)
                error->abort("wrong component in %s file for Cmin(%d,%d,%d)",file_name,icmp+1,jcmp+1,kcmp+1);
            if(kcmp<0 || kcmp>no_types-1)
                error->abort("wrong component in %s file for Cmin(%d,%d,%d)",file_name,icmp+1,jcmp+1,kcmp+1);
            c_min[type_ref[icmp]][type_ref[jcmp]][type_ref[kcmp]]
            =c_min[type_ref[jcmp]][type_ref[icmp]][type_ref[kcmp]]=tmp;
        }
        else if(sscanf(line,"Cmax(%d,%d,%d) = %lf",&icmp,&jcmp,&kcmp,&tmp)==4)
        {
            icmp--;
            jcmp--;
            kcmp--;
            if(icmp<0 || icmp>no_types-1)
                error->abort("wrong component in %s file for Cmax(%d,%d,%d)",file_name,icmp+1,jcmp+1,kcmp+1);
            if(jcmp<0 || jcmp>no_types-1)
                error->abort("wrong component in %s file for Cmax(%d,%d,%d)",file_name,icmp+1,jcmp+1,kcmp+1);
            if(kcmp<0 || kcmp>no_types-1)
                error->abort("wrong component in %s file for Cmax(%d,%d,%d)",file_name,icmp+1,jcmp+1,kcmp+1);
            c_max[type_ref[icmp]][type_ref[jcmp]][type_ref[kcmp]]
            =c_max[type_ref[jcmp]][type_ref[icmp]][type_ref[kcmp]]=tmp;
        }
        else
            error->abort("invalid line in %s file: %s",file_name,line);
    }
    
    if(atoms->my_p==0)
        fclose(fp);

    delete [] line;
    
    if(emb_lin_neg<0 || emb_lin_neg>1)
        error->abort("emb_lin_neg in %s file can only be 0 or 1",file_name);
    if(bkgd_dyn<0 || bkgd_dyn>1)
        error->abort("bkgd_dyn in %s file can only be 0 or 1",file_name);
    if(ialloy<0 || ialloy>2)
        error->abort("ialloy in %s file can only be 0 or 1 or 2",file_name);
    if(mix_ref_t<0 || mix_ref_t>1)
        error->abort("mixture_ref_t in %s file can only be 0 or 1",file_name);
    if(erose_form<0 || erose_form>2)
        error->abort("erose_form in %s file can only be 0 or 1 or 2",file_name);
    if(augt1<0 || augt1>1)
        error->abort("augt1 in %s file can only be 0 or 1",file_name);
    if(rc_meam<=0.0)
        error->abort("rc in %s file should be greater than 0.0",file_name);
    if(delr_meam<=0.0)
        error->abort("delr in %s file should be greater than 0.0",file_name);
    if(gsmooth_factor<0.5 || gsmooth_factor>99.0)
        error->abort("gsmooth_factor in %s file should be between 0.5 & 99.0",file_name);
    for(int itype=0;itype<no_types;itype++)
        for(int jtype=0;jtype<no_types;jtype++)
        {
            char* iname=atom_types->atom_names[itype];
            char* jname=atom_types->atom_names[jtype];
            if(nn2_meam[itype][jtype]<0 || nn2_meam[itype][jtype]>1)
                error->abort("nn2(%s,%s) in %s file should be 0 or 1",iname,jname,file_name);
            if(re_meam[itype][jtype]<=0.0)
                error->abort("re(%s,%s) in %s file should be greater than 0.0",iname,jname,file_name);
        }
}
/*--------------------------------------------
 checked but check again
 --------------------------------------------*/
void ForceField_meam::alloy_params()
{
    for(int i=0;i<no_types;i++)
    {
        for(int j=0;j<no_types;j++)
        {
            
            if(Ec_meam[i][j]==0.0)
            {
                if(lattice[i][j]==L12)
                {
                    Ec_meam[i][j]=0.25*(3.0*Ec_meam[i][i]+Ec_meam[j][j])-delta_meam[i][j];
                    Ec_meam[j][i]=Ec_meam[i][j];
                }
                else if(lattice[i][j]==C11)
                {
                    if(lattice[i][i]==DIAMOND)
                    {
                        Ec_meam[i][j]=(2.0*Ec_meam[i][i]+Ec_meam[j][j])/3.0-delta_meam[i][j];
                        Ec_meam[j][i]=Ec_meam[i][j];
                    }
                    else
                    {
                        Ec_meam[i][j]=(Ec_meam[i][i]+2.0*Ec_meam[j][j])/3.0-delta_meam[i][j];
                        Ec_meam[j][i]=Ec_meam[i][j];
                    }
                    
                }
                else
                {
                    Ec_meam[i][j]=0.5*(Ec_meam[i][i]+Ec_meam[j][j]);
                    Ec_meam[j][i]=Ec_meam[i][j];
                }
            }
            
            if(alpha_meam[i][j]==0.0)
            {
                alpha_meam[i][j]=0.5*(alpha_meam[i][i]+alpha_meam[j][j]);
                alpha_meam[j][i]=alpha_meam[i][j];
            }
            if(re_meam[i][j]==0.0)
            {
                re_meam[i][j]=0.5*(re_meam[i][i]+re_meam[j][j]);
                re_meam[j][i]=re_meam[i][j];
            }
        }
    }
    
    for(int i=0;i<no_types;i++)
        for(int j=0;j<no_types;j++)
            for(int k=0;k<no_types;k++)
                ebound_meam[i][j]=MAX(ebound_meam[i][j],0.25*(c_max[i][j][k]*c_max[i][j][k])/(c_max[i][j][k]-1.0));
    
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField_meam::compute_pair_meam()
{
    int curs;
    type0 Z1,Z2,r,rarat,C,s111,s112,s221,S11,S22;
    type0 phiaa=0.0,phibb=0.0,arat,scrn,scrn2;
    int nmax=10;
    type0 tmp0,tmp1,tmp2,tmp3,astar,phizbl,frac;
    
    for(int itype=0;itype<no_types;itype++)
    {
        for(int jtype=itype;jtype<no_types;jtype++)
        {
            curs=COMP(itype,jtype);
            
            for(int in=0;in<nr;in++)
            {
                r=static_cast<type0>(in)*dr;
                phirar[curs][in][0]=phi_meam(r,itype,jtype);
                //cout<<in<<" "<<phirar[curs][in][0] <<endl;
                if(nn2_meam[itype][jtype]==1)
                {
                    get_Zij(Z1,lattice[itype][jtype]);
                    get_Zij2(Z2,arat,scrn,lattice[itype][jtype],c_min[itype][itype][jtype],c_max[itype][itype][jtype]);
                    
                    if(lattice[itype][jtype]==B1 ||
                       lattice[itype][jtype]==B2 ||
                       lattice[itype][jtype]==L12)
                    {
                        rarat=r*arat;
                        
                        phiaa=phi_meam(rarat,itype,itype);
                        get_Zij(Z1,lattice[itype][itype]);
                        get_Zij2(Z2,arat,scrn,lattice[itype][itype],c_min[itype][itype][itype],c_max[itype][itype][itype]);
                        
                        tmp0=1.0;
                        tmp1=rarat;
                        tmp2=-Z2*scrn/Z1;
                        tmp3=arat;
                        if(scrn>0.0)
                        {
                            for(int i=0;i<nmax;i++)
                            {
                                tmp0*=tmp2;
                                tmp1*=tmp3;
                                phiaa+=tmp0*phi_meam(tmp1,itype,itype);
                            }
                        }
                        
                        phibb=phi_meam(rarat,jtype,jtype);
                        get_Zij(Z1,lattice[jtype][jtype]);
                        get_Zij2(Z2,arat,scrn,lattice[jtype][jtype],c_min[jtype][jtype][jtype],c_max[jtype][jtype][jtype]);
                        
                        tmp0=1.0;
                        tmp1=rarat;
                        tmp2=-Z2*scrn/Z1;
                        tmp3=arat;
                        if(scrn>0.0)
                        {
                            for(int i=0;i<nmax;i++)
                            {
                                tmp0*=tmp2;
                                tmp1*=tmp3;
                                phibb+=tmp0*phi_meam(tmp1,jtype,jtype);
                            }
                        }
                        
                        if(lattice[itype][jtype]==B1 ||
                           lattice[itype][jtype]==B2)
                        {
                            get_Zij(Z1,lattice[itype][jtype]);
                            get_Zij2(Z2,arat,scrn,lattice[itype][jtype],c_min[itype][itype][jtype],c_max[itype][itype][jtype]);
                            phirar[curs][in][0]-=0.5*Z2*scrn/Z1*phiaa;
                            
                            get_Zij2(Z2,arat,scrn2,lattice[itype][jtype],c_min[jtype][jtype][itype],c_max[jtype][jtype][itype]);
                            phirar[curs][in][0]-=0.5*Z2*scrn2/Z1*phibb;
                            
                        }
                        else if(lattice[itype][jtype]==L12)
                        {
                            C=1.0;
                            get_sijk(C,itype,itype,itype,s111);
                            get_sijk(C,itype,itype,jtype,s112);
                            get_sijk(C,jtype,jtype,itype,s221);
                            S11=s111*s111*s112*s112;
                            S22=s221*s221*s221*s221;
                            phirar[curs][in][0]+=-0.75*S11*phiaa-0.25*S22*phibb;
                        }
                        
                    }
                    else
                    {
                        tmp0=1.0;
                        tmp1=r;
                        tmp2=-Z2*scrn/Z1;
                        tmp3=arat;
                        
                        for(int i=0;i<nmax;i++)
                        {
                            tmp0*=tmp2;
                            tmp1*=tmp3;
                            phirar[curs][in][0]+=tmp0*phi_meam(tmp1,itype,jtype);
                            
                        }
                    }
                    
                }
                if(zbl_meam[itype][jtype]==1)
                {
                    astar=alpha_meam[itype][jtype]*(r/re_meam[itype][jtype]-1.0);
                    
                    if(astar<=-3.0)
                    {
                        phirar[curs][in][0]=zbl(r,static_cast<type0>(ielt_meam[itype]),static_cast<type0>(ielt_meam[jtype]));
                    }
                    else if(astar>-3.0 && astar<-1.0)
                    {
                        fcut(1.0+0.5*(astar+1.0),frac);
                        phizbl=zbl(r,static_cast<type0>(ielt_meam[itype]),static_cast<type0>(ielt_meam[jtype]));
                        phirar[curs][in][0]=frac*phirar[curs][in][0]+(1.0-frac)*phizbl;
                    }
                }
            }
            
            phirar[curs][0][1]=phirar[curs][1][0]-phirar[curs][0][0];
            phirar[curs][1][1]=0.5*(phirar[curs][2][0]-phirar[curs][0][0]);
            phirar[curs][nr-2][1]=0.5*(phirar[curs][nr-1][0]-phirar[curs][nr-3][0]);
            phirar[curs][nr-1][1]=0.0;
            
            for(int in=2;in<nr-2;in++)
            {
                phirar[curs][in][1]=((phirar[curs][in-2][0]-phirar[curs][in+2][0])+
                                     8.0*(phirar[curs][in+1][0]-phirar[curs][in-1][0]))/12.0;
            }
            
            phirar[curs][nr-1][2]=0.0;
            phirar[curs][nr-1][3]=0.0;
            for(int in=0;in<nr-1;in++)
            {
                phirar[curs][in][2]=3.0*(phirar[curs][in+1][0]-phirar[curs][in][0])
                -2.0*phirar[curs][in][1]-phirar[curs][in+1][1];
                phirar[curs][in][3]=phirar[curs][in][1]+phirar[curs][in+1][1]
                -2.0*(phirar[curs][in+1][0]-phirar[curs][in][0]);
            }
            
            for(int in=0;in<nr;in++)
            {
                phirar[curs][in][4]=phirar[curs][in][1]*dr_inv;
                phirar[curs][in][5]=2.0*phirar[curs][in][2]*dr_inv;
                phirar[curs][in][6]=3.0*phirar[curs][in][3]*dr_inv;
            }
            
            
        }
    }
    
}
/*--------------------------------------------
 checked
 --------------------------------------------*/
void ForceField_meam::get_sijk(type0 C,int i,
int j,int k,type0& sijk)
{
    type0 x;
    x=(C-c_min[i][j][k])/(c_max[i][j][k]-c_min[i][j][k]);
    fcut(x,sijk);
}
/*--------------------------------------------
 checked
 --------------------------------------------*/
void ForceField_meam::get_dens_ref(type0 r,
int a,int b,type0& rho01,type0& rho11,
type0& rho21,type0& rho31,type0& rho02,
type0& rho12,type0& rho22,type0& rho32)
{
    type0* s;
    CREATE1D(s,3);
    type0 a1=r/re_meam[a][a]-1.0;
    type0 a2=r/re_meam[b][b]-1.0;
    
    type0 rhoa01,rhoa11,rhoa21,rhoa31,
    rhoa02,rhoa12,rhoa22,rhoa32;
    
    rhoa01=rho0_meam[a]*exp(-beta0_meam[a]*a1);
    rhoa11=rho0_meam[a]*exp(-beta1_meam[a]*a1);
    rhoa21=rho0_meam[a]*exp(-beta2_meam[a]*a1);
    rhoa31=rho0_meam[a]*exp(-beta3_meam[a]*a1);
    rhoa02=rho0_meam[b]*exp(-beta0_meam[b]*a2);
    rhoa12=rho0_meam[b]*exp(-beta1_meam[b]*a2);
    rhoa22=rho0_meam[b]*exp(-beta2_meam[b]*a2);
    rhoa32=rho0_meam[b]*exp(-beta3_meam[b]*a2);
    
    int lat=lattice[a][b];
    
    rho11=0.0;
    rho21=0.0;
    rho31=0.0;
    rho12=0.0;
    rho22=0.0;
    rho32=0.0;
    type0 Zij1nn,Zij2nn;
    type0 denom;
    get_Zij(Zij1nn,lat);
    
    
    if(lat==FCC)
    {
        rho01=12.0*rhoa02;
        rho02=12.0*rhoa01;
    }
    else if(lat==BCC)
    {
        rho01=8.0*rhoa02;
        rho02=8.0*rhoa01;
    }
    else if(lat==B1)
    {
        rho01=6.0*rhoa02;
        rho02=6.0*rhoa01;
    }
    else if(lat==DIAMOND)
    {
        rho01=4.0*rhoa02;
        rho02=4.0*rhoa01;
        rho31=32.0/9.0*rhoa32*rhoa32;
        rho32=32.0/9.0*rhoa31*rhoa31;
    }
    else if(lat==HCP)
    {
        rho01=12.0*rhoa02;
        rho02=12.0*rhoa01;
        rho31=1.0/3.0*rhoa32*rhoa32;
        rho32=1.0/3.0*rhoa31*rhoa31;
    }
    else if(lat==DIM)
    {
        get_shpfcn(s,lat);
        rho01=rhoa02;
        rho02=rhoa01;
        rho11=s[0]*rhoa12*rhoa12;
        rho12=s[0]*rhoa11*rhoa11;
        rho21=s[1]*rhoa22*rhoa22;
        rho22=s[1]*rhoa21*rhoa21;
        rho31=s[2]*rhoa32*rhoa32;
        rho32=s[2]*rhoa31*rhoa31;
    }
    else if(lat==C11)
    {
        rho01=rhoa01;
        rho02=rhoa02;
        rho11=rhoa11;
        rho12=rhoa12;
        rho21=rhoa21;
        rho22=rhoa22;
        rho31=rhoa31;
        rho32=rhoa32;
    }
    else if(lat==L12)
    {
        rho01=8.0*rhoa01+4.0*rhoa02;
        rho02=12.0*rhoa01;
        if(ialloy==1)
        {
            rho21=(8.0/3.0)*(rhoa21*t2_meam[a]-rhoa22*t2_meam[b])*(rhoa21*t2_meam[a]-rhoa22*t2_meam[b]);
            denom=8.0*rhoa01*t2_meam[a]*t2_meam[a]+4.0*rhoa02*t2_meam[b]*t2_meam[b];
            if(denom>0.0)
                rho21*=rho01/denom;
        }
        else
        {
            rho21=8.0/3.0*(rhoa21-rhoa22)*(rhoa21-rhoa22);
        }
        
    }
    else if(lat==B2)
    {
        rho01=8.0*rhoa02;
        rho02=8.0*rhoa01;
    }
    else
    {
        error->abort("unknown lattice %d",lat);
    }
    
    
    
    type0 arat,scrn,rhoa01nn,rhoa02nn,C,s111,s112,s221,S11,S22;
    if(nn2_meam[a][b]==1)
    {
        get_Zij2(Zij2nn,arat,scrn,lat,c_min[a][a][b],c_max[a][a][b]);
        a1=arat*r/re_meam[a][a]-1.0;
        a2=arat*r/re_meam[b][b]-1.0;
        rhoa01nn=rho0_meam[a]*exp(-beta0_meam[a]*a1);
        rhoa02nn=rho0_meam[b]*exp(-beta0_meam[b]*a2);
        if(lat==L12)
        {
            C=1.0;
            get_sijk(C,a,a,a,s111);
            get_sijk(C,a,a,b,s112);
            get_sijk(C,b,b,a,s221);
            S11=s111*s111*s112*s112;
            S22=s221*s221*s221*s221;
            rho01+=6.0*S11*rhoa01nn;
            rho02+=6.0*S22*rhoa02nn;
        }
        else
        {
            rho01+=Zij2nn*scrn*rhoa01nn;
            get_Zij2(Zij2nn,arat,scrn,lat,c_min[b][b][a],c_max[b][b][a]);
            rho02+=Zij2nn*scrn*rhoa02nn;
        }
    }
    
    delete [] s;
}
/*--------------------------------------------
 checked
 --------------------------------------------*/
void ForceField_meam::get_tavref(type0& t11av,
type0& t21av,type0& t31av,type0& t12av,
type0& t22av,type0& t32av,type0 t11,type0 t21,
type0 t31,type0 t12,type0 t22,type0 t32,type0 r
,int a,int b,int latt)
{
    if(ialloy==2)
    {
        t11av=t11;
        t21av=t21;
        t31av=t31;
        t12av=t12;
        t22av=t22;
        t32av=t32;
    }
    else
    {
        if(latt==FCC ||
           latt==BCC ||
           latt==DIAMOND ||
           latt==HCP ||
           latt==B1 ||
           latt==DIM ||
           latt==B2
           )
        {
            t11av=t12;
            t21av=t22;
            t31av=t32;
            t12av=t11;
            t22av=t21;
            t32av=t31;
        }
        else
        {
            type0 a1,a2,rhoa01,rhoa02,rho01;
            a1=r/re_meam[a][a]-1.0;
            a2=r/re_meam[b][b]-1.0;
            rhoa01=rho0_meam[a]*exp(-beta0_meam[a]*a1);
            rhoa02=rho0_meam[b]*exp(-beta0_meam[b]*a2);
            if(latt==L12)
            {
                rho01=8.0*rhoa01+4.0*rhoa02;
                t11av=(8.0*t11*rhoa01+4.0*t12*rhoa02)/rho01;
                t12av=t11;
                t21av=(8.0*t21*rhoa01+4.0*t22*rhoa02)/rho01;
                t22av=t21;
                t31av=(8.0*t31*rhoa01+4.0*t32*rhoa02)/rho01;
                t32av=t31;
            }
            else
                error->abort("unknown lattice type");
        }
    }
}
/*--------------------------------------------
 checked
 --------------------------------------------*/
type0 ForceField_meam::phi_meam(type0 r,int a
,int b)
{
    type0 phi_m;
    type0 Z12,rho01,rho11,rho21,rho31
    ,rho02,rho12,rho22,rho32;
    type0 t11av,t21av,t31av,t12av,t22av,t32av;
    type0 scalfac,arat,scrn,Z1nn,Z2nn;
    type0 rhobar1,rhobar2,Z1,Z2,G1,Gam1,G2,Gam2
    ,rho0_1=0.0,rho0_2=0.0,rho_bkgd1,rho_bkgd2,Eu;
    
    get_Zij(Z12,lattice[a][b]);
    get_dens_ref(r,a,b,rho01,rho11,rho21,rho31
                 ,rho02,rho12,rho22,rho32);
    
    
    if(rho01<=1.0e-14 && rho02<=1.0e-14)
    {
        return 0.0;
    }
    
    if(lattice[a][b]==C11)
    {
        if(ialloy==2)
        {
            t11av=t1_meam[a];
            t21av=t2_meam[a];
            t31av=t3_meam[a];
            t12av=t1_meam[b];
            t22av=t2_meam[b];
            t32av=t3_meam[b];
        }
        else
        {
            scalfac=1.0/(rho01+rho02);
            t11av=scalfac*(t1_meam[a]*rho01+t1_meam[b]*rho02);
            t12av=t11av;
            t21av=scalfac*(t2_meam[a]*rho01+t2_meam[b]*rho02);
            t22av=t21av;
            t31av=scalfac*(t3_meam[a]*rho01+t3_meam[b]*rho02);
            t32av=t31av;
        }
    }
    else
    {
        get_tavref(t11av,t21av,t31av,t12av,t22av,t32av
                   ,t1_meam[a],t2_meam[a],t3_meam[a],t1_meam[b]
                   ,t2_meam[b],t3_meam[b],r,a,b,lattice[a][b]);
        
    }
    
    if(lattice[a][b]==C11)
    {
        if(lattice[a][a]==DIAMOND)
        {
            rhobar1=0.25*Z12*Z12*(rho02+rho01)*(rho02+rho01)+
            t11av*(rho12-rho11)*(rho12-rho11)+
            t21av/6.0*(rho22+rho21)*(rho22+rho21)+
            3.025*t31av*(rho32-rho31)*(rho32-rho31);
            rhobar1=sqrt(rhobar1);
            rhobar2=Z12*rho01*Z12*rho01+2.0/3.0*t21av*rho21*rho21;
            rhobar2=sqrt(rhobar2);
        }
        else
        {
            rhobar2=0.25*Z12*Z12*(rho01+rho02)*(rho01+rho02)+
            t12av*(rho12-rho11)*(rho12-rho11)+
            t22av/6.0*(rho22+rho21)*(rho22+rho21)+
            3.025*t32av*(rho32-rho31)*(rho32-rho31);
            rhobar2=sqrt(rhobar2);
            rhobar1=Z12*rho02*Z12*rho02+2.0/3.0*t22av*rho22*rho22;
            rhobar1=sqrt(rhobar1);
        }
    }
    else
    {
        if(mix_ref_t==1)
        {
            Z1=Z_meam[a];
            Z2=Z_meam[b];
            if(ibar_meam[a]<=0)
            {
                G1=1.0;
            }
            else
            {
                type0* s1;
                CREATE1D(s1,3);
                get_shpfcn(s1,lattice[a][a]);
                Gam1=(s1[0]*t11av+s1[1]*t21av+s1[2]*t31av)/(Z1*Z1);
                G_gam(Gam1,ibar_meam[a],gsmooth_factor,G1);
                delete [] s1;
            }
            
            if(ibar_meam[b]<=0)
            {
                G2=1.0;
            }
            else
            {
                type0* s2;
                CREATE1D(s2,3);
                get_shpfcn(s2,lattice[b][b]);
                Gam2=(s2[0]*t12av+s2[1]*t22av+s2[2]*t32av)/(Z2*Z2);
                G_gam(Gam2,ibar_meam[b],gsmooth_factor,G2);
                delete [] s2;
            }
            rho0_1=rho0_meam[a]*Z1*G1;
            rho0_2=rho0_meam[b]*Z2*G2;
        }
        Gam1=t11av*rho11+t21av*rho21+t31av*rho31;
        Gam1=Gam1/(rho01*rho01);
        Gam2=t12av*rho12+t22av*rho22+t32av*rho32;
        Gam2=Gam2/(rho02*rho02);
        G_gam(Gam1,ibar_meam[a],gsmooth_factor,G1);
        G_gam(Gam2,ibar_meam[b],gsmooth_factor,G2);
        if(mix_ref_t==1)
        {
            rho_bkgd1=rho0_1;
            rho_bkgd1=rho0_2;
        }
        else
        {
            if(bkgd_dyn==1)
            {
                rho_bkgd1=rho0_meam[a]*Z_meam[a];
                rho_bkgd2=rho0_meam[b]*Z_meam[b];
            }
            else
            {
                rho_bkgd1=rho_ref_meam[a];
                rho_bkgd2=rho_ref_meam[b];
            }
            
            
        }
        rhobar1=rho01/rho_bkgd1*G1;
        rhobar2=rho02/rho_bkgd2*G2;
        
    }
    
    type0 F1,F2,phiaa=0.0,phibb=0.0;
    
    if(rhobar1==0.0)
        F1=0.0;
    else
    {
        if(emb_lin_neg==1 && rhobar1<=0.0)
            F1=-A_meam[a]*Ec_meam[a][a]*rhobar1;
        else
            F1=A_meam[a]*Ec_meam[a][a]*rhobar1*log(rhobar1);
    }
    if(rhobar2==0.0)
        F2=0.0;
    else
    {
        if(emb_lin_neg==1 && rhobar2<=0.0)
            F2=-A_meam[b]*Ec_meam[b][b]*rhobar2;
        else
            F2=A_meam[b]*Ec_meam[b][b]*rhobar2*log(rhobar2);
    }
    
    Eu=erose(r,re_meam[a][b],alpha_meam[a][b],Ec_meam[a][b]
             ,repuls_meam[a][b],attrac_meam[a][b],erose_form);
    
    
    
    if(lattice[a][b]==C11)
    {
        
        if(lattice[a][a]==DIAMOND)
        {
            phiaa=phi_meam(r,a,a);
            phi_m=(3.0*Eu-F1-2.0*F2-5.0*phiaa)/Z12;
        }
        else
        {
            phibb=phi_meam(r,b,b);
            phi_m=(3.0*Eu-F1-2.0*F2-5.0*phibb)/Z12;
        }
        
    }
    else if(lattice[a][b]==L12)
    {
        phiaa=phi_meam(r,a,a);
        get_Zij(Z1nn,lattice[a][a]);
        get_Zij2(Z2nn,arat,scrn,lattice[a][a],c_min[a][a][a],c_max[a][a][a]);
        int nmax=10;
        if(scrn>0.0)
        {
            type0 tmp0=1.0;
            type0 tmp1=r;
            type0 tmp2=-Z2nn*scrn/Z1nn;
            type0 tmp3=arat;
            for(int i=0;i<nmax;i++)
            {
                tmp0*=tmp2;
                tmp1*=tmp3;
                phiaa+= tmp2*phi_meam(tmp1,a,a);
            }
        }
        phi_m=Eu/3.0-F1/4.0-F2/12.0-phiaa;
    }
    else
    {
        phi_m=(2.0*Eu-F1-F2)/Z12;
    }
    if(r==0.0)
        phi_m=0.0;
    
    return phi_m;
}
/*--------------------------------------------
 checked
 --------------------------------------------*/
void ForceField_meam::compute_reference_density()
{
    type0 Z,Gbar,gam,rho0,Z2,arat,scrn,rho0_2nn;
    type0 s[3];
    
    for(int itype=0;itype<no_types;itype++)
    {
        Z=Z_meam[itype];
        if(ibar_meam[itype]<=0.0)
        {
            Gbar=1.0;
        }
        else
        {
            get_shpfcn(s,lattice[itype][itype]);
            gam=(t1_meam[itype]*s[0]+
                 t2_meam[itype]*s[1]+
                 t3_meam[itype]*s[2])/(Z*Z);
            G_gam(gam,ibar_meam[itype],gsmooth_factor,Gbar);
            
        }
        rho0=rho0_meam[itype]*Z;
        if(nn2_meam[itype][itype]==1)
        {
            get_Zij2(Z2,arat,scrn,lattice[itype][itype],
                     c_min[itype][itype][itype],c_max[itype][itype][itype]);
            rho0_2nn=rho0_meam[itype]*exp(-beta0_meam[itype]*(arat-1.0));
            rho0+=Z2*rho0_2nn*scrn;
        }
        rho_ref_meam[itype]=rho0*Gbar;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField_meam::compute_phi(type0 r
,int itype,int jtype,type0& phi)
{
    int kk,curs=COMP(itype,jtype);
    type0 pp=r*dr_inv;
    kk=static_cast<int>(pp);
    kk=MIN(kk,nr-2);
    pp-=static_cast<type0>(kk);
    pp=MIN(pp,1.0);
    
    phi=((phirar[curs][kk][3]*pp+phirar[curs][kk][2])*pp
         +phirar[curs][kk][1])*pp+phirar[curs][kk][0];
}
/*--------------------------------------------
 
 --------------------------------------------*/
void ForceField_meam::compute_phi_dphi(type0 r
,int itype,int jtype,type0& phi,type0& dphi)
{
    int kk,curs=COMP(itype,jtype);
    type0 pp=r*dr_inv;
    kk=static_cast<int>(pp);
    kk=MIN(kk,nr-2);
    pp-=static_cast<type0>(kk);
    pp=MIN(pp,1.0);
    
    phi=((phirar[curs][kk][3]*pp+phirar[curs][kk][2])*pp
         +phirar[curs][kk][1])*pp+phirar[curs][kk][0];
    dphi=(phirar[curs][kk][6]*pp+phirar[curs][kk][5])*pp
    +phirar[curs][kk][4];
}
