#include "error.h"
#include "md.h"
#include <stdlib.h>
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
MD::MD(MAPP* mapp):InitPtrs(mapp)
{
    ns_alloc=0;
    
    if(forcefield==NULL)
        error->abort("ff should be "
        "initiated before md");
    
    if(mapp->mode!=MD_mode)
        error->abort("md works only "
        "for md mode");
    
    char** args;
    int narg=mapp->parse_line((char*)"KE Temp. "
    "PE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    ke_idx=0;
    temp_idx=1;
    pe_idx=2;
    stress_idx=3;
    
    thermo=new ThermoDynamics(mapp,narg,args);
    for(int i=0;i<narg;i++)
        delete [] args[i];
    delete [] args;
    
    int dim=atoms->dimension;
    if(dim)
    {
        CREATE1D(nrgy_strss,dim*(dim+1)/2+1);
        ns_alloc=1;
    }

}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MD::~MD()
{
    delete thermo;
    if(ns_alloc)
        delete [] nrgy_strss;
}
/*--------------------------------------------
 add time step
 --------------------------------------------*/
void MD::add_dt(int narg,char** args)
{
    if(narg!=2)
        error->abort("wrong command");
    dt=atof(args[1]);
    if(dt<=0.0)
        error->abort("time step cannot be"
        " equal or less than zero");
}
/*--------------------------------------------
 add time step
 --------------------------------------------*/
void MD::add_boltzmann(int narg,char** args)
{
    if(narg!=2)
        error->abort("wrong command");
    boltz=atof(args[1]);
    if(boltz<=0.0)
        error->abort("boltzmann constant cannot"
                     " be equal or less than zero");
}
/*--------------------------------------------
 add time step
 --------------------------------------------*/
void MD::run(int narg,char** args)
{
    if(narg!=2)
        error->abort("wrong command");
    int steps=atoi(args[1]);
    
    if(steps<=0)
        error->abort("number of steps cannot "
        "be equal or less than zero");
    
    init();
    run(steps);
    fin();
    
}
