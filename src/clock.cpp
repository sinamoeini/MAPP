#include "clock.h"
#include "error.h"
#include "thermo_dynamics.h"
#include "atom_types.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock::Clock(MAPP* mapp):InitPtrs(mapp)
{
    ns_alloc=0;
    if(forcefield==NULL)
        error->abort("ff should be "
        "initiated before clock");
    
    if(mapp->mode!=DMD_mode)
        error->abort("clock works only "
        "for md mode");
    
    char** args;
    int narg=mapp->parse_line((char*)
    "Time FE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    
    fe_idx=1;
    stress_idx=2;
    time_idx=0;
    thermo=new ThermoDynamics(mapp,narg,args);
    for(int i=0;i<narg;i++)
        delete [] args[i];
    delete [] args;
    
    cdof_n=atoms->find_exist("cdof");
    
    int dim=atoms->dimension;
    if(dim)
    {
        CREATE1D(nrgy_strss,dim*(dim+1)/2+1);
        ns_alloc=1;
    }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock::~Clock()
{
    delete thermo;
    if(ns_alloc)
        delete [] nrgy_strss;
}

/*--------------------------------------------
 rectify
 --------------------------------------------*/
void Clock::rectify(type0* f)
{
    if(cdof_n==-1)
        return;
    atoms->vectors[cdof_n].ret(cdof);
    
    int tot=(atoms->natms)*(atom_types->no_types);
    for(int i=0;i<tot;i++) if(cdof[i]==1) f[i]=0.0;
}
