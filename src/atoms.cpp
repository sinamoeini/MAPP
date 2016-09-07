#include "atoms.h"
using namespace MAPP_NS;
/*----------------------------------------------
  _     _   _____   _____   _       _____   _____
 | |   / / | ____| /  ___| | |     /  ___/ |_   _|
 | |  / /  | |__   | |     | |     | |___    | |
 | | / /   |  __|  | |     | |     \___  \   | |
 | |/ /    | |___  | |___  | |___   ___| |   | |
 |___/     |_____| \_____| |_____| /_____/   |_|
 
 ----------------------------------------------*/
VecLst::VecLst(Atoms<MAPP_NS::dimension>* atoms)
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

