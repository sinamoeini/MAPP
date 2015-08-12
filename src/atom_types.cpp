/*--------------------------------------------
 Created by Sina on 06/27/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#include "atom_types.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
AtomTypes::AtomTypes(MAPP* mapp):InitPtrs(mapp)
{
    no_types=0;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
AtomTypes::~AtomTypes()
{

    for(int i=0;i<no_types;i++)
    {
        delete [] atom_names[i];
        delete [] clr_rad[i];
    }
    if(no_types)
    {
        delete [] atom_names;
        delete [] clr_rad;
        delete [] mass;
    }
    
}
/*--------------------------------------------
 add a type
 --------------------------------------------*/
int AtomTypes::add_type(type0 m,char* name)
{
    
    
    int type=0;
    for (int i=0;i<no_types;i++)
        if(!strcmp(name,atom_names[i]))
            type=i+1;
    
    if (type)
        return (type-1);
    
    GROW(atom_names,no_types,no_types+1);
    GROW(clr_rad,no_types,no_types+1);
    CREATE1D(clr_rad[no_types],4);
    GROW(mass,no_types,no_types+1);
    int lngth= static_cast<int>(strlen(name))+1;
    CREATE1D(atom_names[no_types],lngth);
    for(int i=0;i<lngth;i++)
        atom_names[no_types][i]=name[i];
    mass[no_types]=m;
    assign_color_rad(name,clr_rad[no_types]);
    
    no_types++;
    return (no_types-1);
   
}
/*--------------------------------------------
 find a type
 --------------------------------------------*/
int AtomTypes::find_type(char* name)
{
    int type=0;
    for (int i=0;i<no_types;i++)
        if(!strcmp(name,atom_names[i]))
            type=i+1;
    
    if(type)
        return (type-1);
    else
    {
        error->abort("atom type %s not found",name);
        return -1;
    }
}
/*--------------------------------------------
 find a type without error
 --------------------------------------------*/
int AtomTypes::find_type_exist(char* name)
{

    for (int i=0;i<no_types;i++)
        if(!strcmp(name,atom_names[i]))
            return i;
    
    return -1;
}
/*--------------------------------------------
 find a type without error
 --------------------------------------------*/
void AtomTypes::assign_color_rad(char* name,type0* mat)
{
    if(0)
    {}
    else if(strcmp(name,"H")==0)
    {
        mat[0]=0.800000;
        mat[1]=0.800000;
        mat[2]=0.800000;
        mat[3]=0.435000;
    }
    else if(strcmp(name,"He")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=0.700000;
    }
    else if(strcmp(name,"Li")==0)
    {
        mat[0]=0.700000;
        mat[1]=0.700000;
        mat[2]=0.700000;
        mat[3]=1.519900;
    }
    else if(strcmp(name,"Be")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.143000;
    }
    else if(strcmp(name,"B")==0)
    {
        mat[0]=0.900000;
        mat[1]=0.400000;
        mat[2]=0.000000;
        mat[3]=0.975000;
    }
    else if(strcmp(name,"C")==0)
    {
        mat[0]=0.350000;
        mat[1]=0.350000;
        mat[2]=0.350000;
        mat[3]=0.655000;
    }
    else if(strcmp(name,"N")==0)
    {
        mat[0]=0.200000;
        mat[1]=0.200000;
        mat[2]=0.800000;
        mat[3]=0.750000;
    }
    else if(strcmp(name,"O")==0)
    {
        mat[0]=0.800000;
        mat[1]=0.200000;
        mat[2]=0.200000;
        mat[3]=0.730000;
    }
    else if(strcmp(name,"F")==0)
    {
        mat[0]=0.700000;
        mat[1]=0.850000;
        mat[2]=0.450000;
        mat[3]=0.720000;
    }
    else if(strcmp(name,"Ne")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.600000;
    }
    else if(strcmp(name,"Na")==0)
    {
        mat[0]=0.600000;
        mat[1]=0.600000;
        mat[2]=0.600000;
        mat[3]=1.857900;
    }
    else if(strcmp(name,"Mg")==0)
    {
        mat[0]=0.600000;
        mat[1]=0.600000;
        mat[2]=0.700000;
        mat[3]=1.604700;
    }
    else if(strcmp(name,"Al")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.431800;
    }
    else if(strcmp(name,"Si")==0)
    {
        mat[0]=0.690196;
        mat[1]=0.768627;
        mat[2]=0.870588;
        mat[3]=1.175800;
    }
    else if(strcmp(name,"P")==0)
    {
        mat[0]=0.100000;
        mat[1]=0.700000;
        mat[2]=0.300000;
        mat[3]=1.060000;
    }
    else if(strcmp(name,"S")==0)
    {
        mat[0]=0.950000;
        mat[1]=0.900000;
        mat[2]=0.200000;
        mat[3]=1.020000;
    }
    else if(strcmp(name,"Cl")==0)
    {
        mat[0]=0.150000;
        mat[1]=0.500000;
        mat[2]=0.100000;
        mat[3]=0.990000;
    }
    else if(strcmp(name,"Ar")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.900000;
    }
    else if(strcmp(name,"K")==0)
    {
        mat[0]=0.576471;
        mat[1]=0.439216;
        mat[2]=0.858824;
        mat[3]=2.262000;
    }
    else if(strcmp(name,"Ca")==0)
    {
        mat[0]=0.800000;
        mat[1]=0.800000;
        mat[2]=0.700000;
        mat[3]=1.975800;
    }
    else if(strcmp(name,"Sc")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.654500;
    }
    else if(strcmp(name,"Ti")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.475500;
    }
    else if(strcmp(name,"V")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.309000;
    }
    else if(strcmp(name,"Cr")==0)
    {
        mat[0]=0.000000;
        mat[1]=0.800000;
        mat[2]=0.000000;
        mat[3]=1.249000;
    }
    else if(strcmp(name,"Mn")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.350000;
    }
    else if(strcmp(name,"Fe")==0)
    {
        mat[0]=0.517647;
        mat[1]=0.576471;
        mat[2]=0.652941;
        mat[3]=1.241100;
    }
    else if(strcmp(name,"Co")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.253500;
    }
    else if(strcmp(name,"Ni")==0)
    {
        mat[0]=0.257255;
        mat[1]=0.266667;
        mat[2]=0.271373;
        mat[3]=1.246000;
    }
    else if(strcmp(name,"Cu")==0)
    {
        mat[0]=0.950000;
        mat[1]=0.790074;
        mat[2]=0.013859;
        mat[3]=1.278000;
    }
    else if(strcmp(name,"Zn")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.332500;
    }
    else if(strcmp(name,"Ga")==0)
    {
        mat[0]=0.900000;
        mat[1]=0.000000;
        mat[2]=1.000000;
        mat[3]=1.350100;
    }
    else if(strcmp(name,"Ge")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.224800;
    }
    else if(strcmp(name,"As")==0)
    {
        mat[0]=1.000000;
        mat[1]=1.000000;
        mat[2]=0.300000;
        mat[3]=1.200000;
    }
    else if(strcmp(name,"Se")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.160000;
    }
    else if(strcmp(name,"Br")==0)
    {
        mat[0]=0.500000;
        mat[1]=0.080000;
        mat[2]=0.120000;
        mat[3]=1.140000;
    }
    else if(strcmp(name,"Kr")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.000000;
    }
    else if(strcmp(name,"Rb")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.470000;
    }
    else if(strcmp(name,"Sr")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.151300;
    }
    else if(strcmp(name,"Y")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.823700;
    }
    else if(strcmp(name,"Zr")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.615600;
    }
    else if(strcmp(name,"Nb")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.431800;
    }
    else if(strcmp(name,"Mo")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.362600;
    }
    else if(strcmp(name,"Tc")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.367500;
    }
    else if(strcmp(name,"Ru")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.352900;
    }
    else if(strcmp(name,"Rh")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.345000;
    }
    else if(strcmp(name,"Pd")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.375500;
    }
    else if(strcmp(name,"Ag")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.444700;
    }
    else if(strcmp(name,"Cd")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.489400;
    }
    else if(strcmp(name,"In")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.666200;
    }
    else if(strcmp(name,"Sn")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.537500;
    }
    else if(strcmp(name,"Sb")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.400000;
    }
    else if(strcmp(name,"Te")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.360000;
    }
    else if(strcmp(name,"I")==0)
    {
        mat[0]=0.500000;
        mat[1]=0.100000;
        mat[2]=0.500000;
        mat[3]=1.330000;
    }
    else if(strcmp(name,"Xe")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.200000;
    }
    else if(strcmp(name,"Cs")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.632500;
    }
    else if(strcmp(name,"Ba")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.170500;
    }
    else if(strcmp(name,"La")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.872500;
    }
    else if(strcmp(name,"Ce")==0)
    {
        mat[0]=0.800000;
        mat[1]=0.800000;
        mat[2]=0.000000;
        mat[3]=1.824300;
    }
    else if(strcmp(name,"Pr")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.836200;
    }
    else if(strcmp(name,"Nd")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.829500;
    }
    else if(strcmp(name,"Pm")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.809000;
    }
    else if(strcmp(name,"Sm")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.804000;
    }
    else if(strcmp(name,"Eu")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.984000;
    }
    else if(strcmp(name,"Gd")==0)
    {
        mat[0]=1.000000;
        mat[1]=0.843137;
        mat[2]=0.000000;
        mat[3]=1.818000;
    }
    else if(strcmp(name,"Tb")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.800500;
    }
    else if(strcmp(name,"Dy")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.795100;
    }
    else if(strcmp(name,"Ho")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.788600;
    }
    else if(strcmp(name,"Er")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.779400;
    }
    else if(strcmp(name,"Tm")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.768700;
    }
    else if(strcmp(name,"Yb")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.939600;
    }
    else if(strcmp(name,"Lu")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.751500;
    }
    else if(strcmp(name,"Hf")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.597300;
    }
    else if(strcmp(name,"Ta")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.428000;
    }
    else if(strcmp(name,"W")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.370500;
    }
    else if(strcmp(name,"Re")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.380000;
    }
    else if(strcmp(name,"Os")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.367600;
    }
    else if(strcmp(name,"Ir")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.357300;
    }
    else if(strcmp(name,"Pt")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.387300;
    }
    else if(strcmp(name,"Au")==0)
    {
        mat[0]=0.900000;
        mat[1]=0.800000;
        mat[2]=0.000000;
        mat[3]=1.441900;
    }
    else if(strcmp(name,"Hg")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.502500;
    }
    else if(strcmp(name,"Tl")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.728300;
    }
    else if(strcmp(name,"Pb")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.750100;
    }
    else if(strcmp(name,"Bi")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.460000;
    }
    else if(strcmp(name,"Po")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.460000;
    }
    else if(strcmp(name,"At")==0)
    {
        mat[0]=0.800000;
        mat[1]=0.200000;
        mat[2]=0.200000;
        mat[3]=4.350000;
    }
    else if(strcmp(name,"Rn")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.430000;
    }
    else if(strcmp(name,"Fr")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.500000;
    }
    else if(strcmp(name,"Ra")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=2.140000;
    }
    else if(strcmp(name,"Ac")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.877500;
    }
    else if(strcmp(name,"Th")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.797500;
    }
    else if(strcmp(name,"Pa")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.608600;
    }
    else if(strcmp(name,"U")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.568300;
    }
    else if(strcmp(name,"Np")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"Pu")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"Am")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"Cm")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"Bk")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"Cf")==0)
    {
        mat[0]=0.100000;
        mat[1]=0.700000;
        mat[2]=0.300000;
        mat[3]=1.460000;
    }
    else if(strcmp(name,"Es")==0)
    {
        mat[0]=0.100000;
        mat[1]=0.300000;
        mat[2]=0.700000;
        mat[3]=1.752000;
    }
    else if(strcmp(name,"Fm")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"Md")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"No")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"Lr")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"Rf")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"Db")==0)
    {
        mat[0]=0.900000;
        mat[1]=0.800000;
        mat[2]=0.000000;
        mat[3]=1.460000;
    }
    else if(strcmp(name,"Sg")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"Bh")==0)
    {
        mat[0]=1.000000;
        mat[1]=1.000000;
        mat[2]=0.000000;
        mat[3]=1.900000;
    }
    else if(strcmp(name,"Hs")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else if(strcmp(name,"Mt")==0)
    {
        mat[0]=0.643137;
        mat[1]=0.666667;
        mat[2]=0.678431;
        mat[3]=1.000000;
    }
    else
    {
        mat[0]=0.500000;
        mat[1]=0.500000;
        mat[2]=0.500000;
        mat[3]=0.435000;
    }
}


