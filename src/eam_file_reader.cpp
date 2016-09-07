#include <stdlib.h>
#include "eam_file_reader.h"
#include "atom_types.h"
#include "memory.h"
#include "error.h"
#include "atoms.h"
#include "script_reader.h"
enum{NOT_SET,FUNC_FL,SET_FL,FINNIS_FL};
using namespace MAPP_NS;
/*--------------------------------------------
 
 --------------------------------------------*/
EAMFileReader::EAMFileReader()

{
    allocated=0;
    no_types=0;
    nfiles=0;
    
    type2phi=NULL;
    type2rho=NULL;
    type2rho_pair_ij=NULL;
    type2rho_pair_ji=NULL;
    type2phi_pair_ij=NULL;
    type2phi_pair_ji=NULL;
    cut=NULL;
    cut_sq=NULL;
}
/*--------------------------------------------
 
 --------------------------------------------*/
EAMFileReader::~EAMFileReader()
{
    if(allocated)
        deallocate();
    for(int ifile=0;ifile<nfiles;ifile++)
        delete [] files[ifile];
    if(nfiles)
    {
        delete [] types;
        delete [] files;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void EAMFileReader::setup()
{
    deallocate();
    no_types=atom_types->no_types;
    
    if(eam_mode==FINNIS_FL)
    {
        if(nfiles!=1)
            error->abort("number of files for "
            "FS should be one");
        set_fs();
    }
    else if(eam_mode==SET_FL)
    {
        if(nfiles!=1)
            error->abort("number of files for "
            "SetFL should be one");
        
        set_setfl();
    }
    else if(eam_mode==FUNC_FL)
    {
        if(nfiles!=no_types)
            error->abort("number of files for "
            "FuncFL should be euqal to number of "
            "atomic types");
        
        char** files_;
        CREATE1D(files_,no_types);
        for(int itype=0;itype<no_types;itype++)
        {
            int ifile=0;
            while(types[ifile]==itype && ifile<nfiles)
                ifile++;
            if(ifile==nfiles)
                error->abort("file FuncFL for %s "
                "was not found",atom_types->atom_names[itype]);

            files_[itype]=files[itype];
        }
        
        if(nfiles)
            delete files;
        files=files_;
        
        set_funcfl();
    }
    
    set_arrays();
}
/*--------------------------------------------
 
 --------------------------------------------*/
void EAMFileReader::add_file(const char* file,int itype)
{
    int lngth=static_cast<int>(strlen(file))+1;
    GROW(files,nfiles,nfiles+1);
    GROW(types,nfiles,nfiles+1);
    CREATE1D(files[nfiles],lngth);
    memcpy(files[nfiles],file,lngth*sizeof(char));
    types[nfiles]=itype;
    nfiles++;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void EAMFileReader::file_format(const char* format)
{
    if(strcmp(format,"FS")==0)
        eam_mode=FINNIS_FL;
    else if(strcmp(format,"SetFL")==0)
        eam_mode=SET_FL;
    else if(strcmp(format,"FuncFL")==0)
        eam_mode=FUNC_FL;
    else
        error->abort("unknown file format "
        "for ff eam_dmd: %s",format);
}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void EAMFileReader::set_funcfl()
{
    if(nfiles!=no_types)
        error->abort("for FuncFL mode number of"
        " ff eam_dmd files should be equal to the number"
        " of atom types present in the system");
    
    type0* drs;
    type0* drhos;
    int* nrs;
    int* nrhos;
    CREATE1D(drs,no_types);
    CREATE1D(drhos,no_types);
    CREATE1D(nrs,no_types);
    CREATE1D(nrhos,no_types);
    
    type0** tmp_F;
    type0** tmp_rho;
    type0** tmp_zi;
    CREATE1D(tmp_F,no_types);
    CREATE1D(tmp_rho,no_types);
    CREATE1D(tmp_zi,no_types);
    type0* tmp;
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    char** args=NULL;
    int args_cpcty=0;
    int nargs;
    
    type0 mass;
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        ScriptReader::open_file(fp,files[ityp],"r");
        
        for(int i=0;i<2;i++)
            if(ScriptReader::read_line(fp,line)==-1)
                error->abort("%s file ended immaturely",files[ityp]);
        
        nargs=ScriptReader::parse_line(line,args,args_cpcty);
        if(nargs!=4)
            error->abort("invalid line in %s file: %s",files[ityp],line);
        
        mass=atof(args[1]);
        
        if(mass!=atom_types->mass[ityp])
            error->abort("mass of element %s in %s file (%lf) does not match "
            "the mass that is already assigned to the element (%lf), make "
            "sure that the right sequence of ff eam_dmd files are "
            "used",atom_types->atom_names[ityp],files[ityp],mass
            ,atom_types->mass[ityp]);
  
        if(ScriptReader::read_line(fp,line)==-1)
            error->abort("%s file ended immaturely",files[ityp]);
        
        nargs=ScriptReader::parse_line(line,args,args_cpcty);
        if(nargs!=5)
            error->abort("invalid line in %s file: %s",files[ityp],line);
        
        nrhos[ityp]=atoi(args[0]);
        nrs[ityp]=atoi(args[2]);
        drhos[ityp]=atof(args[1]);
        drs[ityp]=atof(args[3]);
        
        if(nrhos[ityp]<5)
            error->abort("nrho in %s file should be larger than 5",files[ityp]);
        if(nrs[ityp]<5)
            error->abort("nr in %s file should be larger than 5",files[ityp]);
        if(drhos[ityp]<=0.0)
            error->abort("drho in %s file should be larger than 0.0",files[ityp]);
        if(drs[ityp]<=0.0)
            error->abort("dr in %s file should be larger than 0.0",files[ityp]);
        
        
        
        int tot=nrhos[ityp]+2*nrs[ityp];
        
        CREATE1D(tmp,tot);
        
        int ipos=0;
        while (ipos<tot)
        {
            if(ScriptReader::read_line(fp,line)==-1)
                error->abort("%s file ended immaturely",files[ityp]);
            
            nargs=ScriptReader::parse_line(line,args,args_cpcty);
            
            if(ipos+nargs>tot)
                error->abort("%s file ended immaturely",files[ityp]);
            
            for(int i=0;i<nargs;i++)
            {
                tmp[ipos]=atof(args[i]);
                ipos++;
            }
        }
        
        if(atoms->my_p==0)
            fclose(fp);
        
        CREATE1D(tmp_F[ityp],nrhos[ityp]);
        CREATE1D(tmp_zi[ityp],nrs[ityp]);
        CREATE1D(tmp_rho[ityp],nrs[ityp]);
        
        
        memcpy(tmp_F[ityp],&tmp[0],nrhos[ityp]*sizeof(type0));
        memcpy(tmp_zi[ityp],&tmp[nrhos[ityp]],nrs[ityp]*sizeof(type0));
        memcpy(tmp_rho[ityp],&tmp[nrhos[ityp]+nrs[ityp]],nrs[ityp]*sizeof(type0));
        
        delete [] tmp;
        
    }
    delete [] line;
    if(args_cpcty)
        delete [] args;
    
    
    type0 maxr=0.0;
    type0 maxrho=0.0;
    type0 maxdr=0.0;
    type0 maxdrho=0.0;
    
    for(int i=0;i<no_types;i++)
    {
        maxr=MAX(maxr,static_cast<type0>(nrs[i]-1)*drs[i]);
        maxrho=MAX(maxrho,static_cast<type0>(nrhos[i]-1)*drhos[i]);
        maxdr=MAX(maxdr,drs[i]);
        maxdrho=MAX(maxdrho,drhos[i]);
    }
    
    nr=static_cast<int>(maxr/maxdr+0.5);
    nrho=static_cast<int>(maxrho/maxdrho+0.5);
    dr=maxdr;
    drho=maxdrho;
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    allocate();
    
    type0 r,p,tmp0,tmp1;
    int k;
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        for(int i=0;i<nrho;i++)
        {
            r=static_cast<type0>(i)*drho;
            p=r/drhos[ityp];
            k=static_cast<int> (p);
            k=MIN(k,nrhos[ityp]-2);
            p-=k;
            p=MIN(p,1.0);
            
            F_arr[ityp][i][0]=interpolate(tmp_F[ityp],nrs[ityp],p,k);
        }
    }
    
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        for(int i=0;i<nr;i++)
        {
            r=static_cast<type0>(i)*dr;
            p=r/drs[ityp];
            k=static_cast<int> (p);
            k=MIN(k,nrs[ityp]-2);
            p-=k;
            p=MIN(p,1.0);
            rho_arr[type2rho[ityp][0]][i][0]=interpolate(tmp_rho[ityp],nrs[ityp],p,k);
        }
    }
    
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        for(int jtyp=0;jtyp<ityp+1;jtyp++)
        {
            for(int i=0;i<nr;i++)
            {
                
                r=static_cast<type0>(i)*dr;
                
                p=r/drs[ityp];
                k=static_cast<int> (p);
                k=MIN(k,nrs[ityp]-2);
                p-=k;
                p=MIN(p,1.0);
                tmp0=interpolate(tmp_zi[ityp],nrs[ityp],p,k);
                
                
                p=r/drs[jtyp];
                k=static_cast<int> (p);
                k=MIN(k,nrs[jtyp]-2);
                p-=k;
                p=MIN(p,1.0);
                tmp1=interpolate(tmp_zi[jtyp],nrs[jtyp],p,k);
                
                phi_r_arr[type2phi[ityp][jtyp]][i][0]=27.2*0.529*tmp0*tmp1;
            }
        }
    }
    
    
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        if(nrhos[ityp])
            delete [] tmp_F[ityp];
        
        if(nrs[ityp])
        {
            delete [] tmp_rho[ityp];
            delete [] tmp_zi[ityp];
        }
    }
    if(no_types)
    {
        delete [] tmp_F;
        delete [] tmp_rho;
        delete [] tmp_zi;
        delete [] drhos;
        delete [] drs;
        delete [] nrhos;
        delete [] nrs;
    }
    
}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void EAMFileReader::set_setfl()
{
    if(nfiles!=1)
        error->abort("for SetFL mode number of"
        " ff eam_dmd files should be one");

    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);

    ScriptReader::open_file(fp,files[0],"r");
    
    for(int i=0;i<4;i++)
        if(ScriptReader::read_line(fp,line)==-1)
            error->abort("%s file ended immaturely",files[0]);
    
    char** args=NULL;
    int args_cpcty=0;
    int nargs;
    
    nargs=ScriptReader::parse_line(line,args,args_cpcty);
    if(nargs<2 || nargs-1 < no_types)
        error->abort("invalid line in %s file: %s",files[0],line);
    
    int tot_no_types=nargs-1;
    
    int* type_ref;
    CREATE1D(type_ref,tot_no_types);
    for(int i=0;i<tot_no_types;i++)
        type_ref[i]=atom_types->find_type_exist(args[i+1]);
    
    int* tmp_type_ref;
    CREATE1D(tmp_type_ref,no_types);
    for(int i=0;i<no_types;i++)
        tmp_type_ref[i]=-1;
    for(int ityp=0;ityp<tot_no_types;ityp++)
        if(type_ref[ityp]!=-1)
            tmp_type_ref[type_ref[ityp]]=0;
    
    for(int i=0;i<no_types;i++)
        if(tmp_type_ref[i]==-1)
            error->abort("%s file does not contain "
            "parameters for element %s",files[0]
            ,atom_types->atom_names[i]);
    if(no_types)
        delete [] tmp_type_ref;
    
    
    if(ScriptReader::read_line(fp,line)==-1)
        error->abort("%s file ended immaturely",files[0]);
    nargs=ScriptReader::parse_line(line,args,args_cpcty);
    
    
    if(nargs!=5)
        error->abort("invalid line in %s file: %s",files[0],line);
    
    nrho=atoi(args[0]);
    nr=atoi(args[2]);
    drho=atof(args[1]);
    dr=atof(args[3]);
    
    if(nrho<5)
        error->abort("nrho in %s file should be larger than 5",files[0]);
    if(nr<5)
        error->abort("nr in %s file should be larger than 5",files[0]);
    if(drho<=0.0)
        error->abort("drho in %s file should be larger than 0.0",files[0]);
    if(dr<=0.0)
        error->abort("dr in %s file should be larger than 0.0",files[0]);
    
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    
    int ipos=0;
    int tot=tot_no_types*(nrho+nr)+tot_no_types*(tot_no_types+1)*nr/2;
    type0* tmp;
    CREATE1D(tmp,tot);
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(ScriptReader::read_line(fp,line)==-1)
            error->abort("%s file ended immaturely",files[0]);
        int lim=(ityp+1)*(nrho+nr);
        while (ipos<lim)
        {
            if(ScriptReader::read_line(fp,line)==-1)
                error->abort("%s file ended immaturely",files[0]);
            
            nargs=ScriptReader::parse_line(line,args,args_cpcty);
            for(int i=0;i<nargs;i++)
            {
                tmp[ipos]=atof(args[i]);
                ipos++;
            }
        }
        if(ipos!=lim)
            error->abort("unknown line in file %s:",files[0],line);
    }
    
    while (ipos<tot)
    {
        if(ScriptReader::read_line(fp,line)==-1)
            error->abort("%s file ended immaturely",files[0]);
        
        nargs=ScriptReader::parse_line(line,args,args_cpcty);
        for(int i=0;i<nargs;i++)
        {
            tmp[ipos]=atof(args[i]);
            ipos++;
        }
    }
    
    if(ipos!=tot)
        error->abort("unknown line in file %s:",files[0],line);
    
    delete [] line;
    if(args_cpcty)
        delete [] args;
    
    if(atoms->my_p==0)
        fclose(fp);
    
    
    allocate();
    
    
    int component;
    
    ipos=0;
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(type_ref[ityp]==-1)
        {
            ipos+=nrho+nr;
        }
        else
        {
            component=type_ref[ityp];
            
            for(int i=0;i<nrho;i++)
                F_arr[component][i][0]=tmp[ipos+i];
            
            ipos+=nrho;
            
            component=type2rho[type_ref[ityp]][type_ref[0]];
            for(int i=0;i<nr;i++)
                rho_arr[component][i][0]=tmp[ipos+i];
            
            ipos+=nr;
        }
    }
    
    
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(type_ref[ityp]==-1)
        {
            ipos+=nr*(ityp+1);
        }
        else
        {
            for(int jtyp=0;jtyp<ityp+1;jtyp++)
            {
                if(type_ref[jtyp]==-1)
                {
                    ipos+=nr;
                }
                else
                {
                    component=type2phi[type_ref[ityp]][type_ref[jtyp]];
                    
                    for(int i=0;i<nr;i++)
                        phi_r_arr[component][i][0]=tmp[ipos+i];
                    
                    ipos+=nr;
                }
            }
        }
    }
    
    delete [] tmp;
    
    if(tot_no_types)
        delete [] type_ref;
    
}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void EAMFileReader::set_fs()
{
    if(nfiles!=1)
        error->abort("one file is needed for eam");
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    
    ScriptReader::open_file(fp,files[0],"r");
    
    for(int i=0;i<4;i++)
        if(ScriptReader::read_line(fp,line)==-1)
            error->abort("%s file ended immaturely",files[0]);
    
    char** args=NULL;
    int args_cpcty=0;
    int nargs;
    
    nargs=ScriptReader::parse_line(line,args,args_cpcty);
    if(nargs<2 || nargs-1<no_types)
        error->abort("invalid line in %s file: %s",files[0],line);
    
    int tot_no_types=nargs-1;
    
    
    
    int* type_ref;
    CREATE1D(type_ref,tot_no_types);
    for(int i=0;i<tot_no_types;i++)
        type_ref[i]=atom_types->find_type_exist(args[i+1]);
    
    int* tmp_type_ref;
    CREATE1D(tmp_type_ref,no_types);
    for(int i=0;i<no_types;i++)
        tmp_type_ref[i]=-1;
    for(int ityp=0;ityp<tot_no_types;ityp++)
        if(type_ref[ityp]!=-1)
            tmp_type_ref[type_ref[ityp]]=0;
    
    for(int i=0;i<no_types;i++)
        if(tmp_type_ref[i]==-1)
            error->abort("%s file does not contain "
            "parameters for element %s",files[0]
            ,atom_types->atom_names[i]);
    
    if(no_types)
        delete [] tmp_type_ref;
    
    
    if(ScriptReader::read_line(fp,line)==-1)
        error->abort("%s file ended immaturely",files[0]);
    nargs=ScriptReader::parse_line(line,args,args_cpcty);
    
    if(nargs!=5)
        error->abort("invalid line in %s file: %s",files[0],line);
    
    nrho=atoi(args[0]);
    nr=atoi(args[2]);
    drho=atof(args[1]);
    dr=atof(args[3]);
    
    
    if(nrho<5)
        error->abort("nrho in %s file should be larger than 5",files[0]);
    if(nr<5)
        error->abort("nr in %s file should be larger than 5",files[0]);
    if(drho<=0.0)
        error->abort("drho in %s file should be larger than 0.0",files[0]);
    if(dr<=0.0)
        error->abort("dr in %s file should be larger than 0.0",files[0]);
    
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    
    int ipos;
    int tot=tot_no_types*nrho+tot_no_types*tot_no_types*nr
    +tot_no_types*(tot_no_types+1)*nr/2;
    int tot_loc=nrho+tot_no_types*nr;
    int ipos_loc;
    
    type0* tmp;
    CREATE1D(tmp,tot);
    
    
    ipos=0;
    
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(ScriptReader::read_line(fp,line)==-1)
            error->abort("%s file ended immaturely",files[0]);
        ipos_loc=0;
        while (ipos_loc<tot_loc)
        {
            if(ScriptReader::read_line(fp,line)==-1)
                error->abort("%s file ended immaturely",files[0]);
            
            nargs=ScriptReader::parse_line(line,args,args_cpcty);
            
            if(ipos_loc+nargs>tot_loc)
                error->abort("%s file ended immaturely",files[0]);
            
            for(int i=0;i<nargs;i++)
            {
                tmp[ipos+ipos_loc]=atof(args[i]);
                ipos_loc++;
            }
        }
        ipos+=ipos_loc;
    }
    
    tot_loc=tot_no_types*(tot_no_types+1)*nr/2;
    ipos_loc=0;
    while (ipos_loc<tot_loc)
    {
        if(ScriptReader::read_line(fp,line)==-1)
            error->abort("%s file ended immaturely",files[0]);
        
        nargs=ScriptReader::parse_line(line,args,args_cpcty);
        
        if(ipos_loc+nargs>tot_loc)
            error->abort("%s file ended immaturely",files[0]);
        
        for(int i=0;i<nargs;i++)
        {
            tmp[ipos+ipos_loc]=atof(args[i]);
            ipos_loc++;
        }
    }
    
    delete [] line;
    if(args_cpcty)
        delete [] args;
    if(atoms->my_p==0)
        fclose(fp);
    
    allocate();
    
    
    int component;
    
    ipos=0;
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(type_ref[ityp]==-1)
        {
            ipos+=nr*tot_no_types+nrho;
        }
        else
        {
            component=type_ref[ityp];
            
            for(int i=0;i<nrho;i++)
                F_arr[component][i][0]=tmp[ipos+i];
            
            ipos+=nrho;
            
            for(int jtyp=0;jtyp<tot_no_types;jtyp++)
            {
                if(type_ref[jtyp]==-1)
                {
                    ipos+=nr;
                }
                else
                {
                    component=type2rho[type_ref[ityp]][type_ref[jtyp]];
                    
                    for(int i=0;i<nr;i++)
                        rho_arr[component][i][0]=tmp[ipos+i];
                    ipos+=nr;
                }
                
            }
        }
    }
    
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(type_ref[ityp]==-1)
        {
            ipos+=nr*(ityp+1);
        }
        else
        {
            for(int jtyp=0;jtyp<ityp+1;jtyp++)
            {
                if(type_ref[jtyp]==-1)
                {
                    ipos+=nr;
                }
                else
                {
                    component=type2phi[type_ref[ityp]][type_ref[jtyp]];
                    
                    for(int i=0;i<nr;i++)
                        phi_r_arr[component][i][0]=tmp[ipos+i];
                    
                    ipos+=nr;
                }
            }
        }
    }
    
    delete [] tmp;
    if(tot_no_types)
        delete [] type_ref;
    
}
/*--------------------------------------------
 clean up the arrays
 --------------------------------------------*/
void EAMFileReader::deallocate()
{
    if(allocated==0)
        return;
    
    DEL_2D(type2phi);
    DEL_2D(type2rho);
    DEL_2D(type2rho_pair_ij);
    DEL_2D(type2rho_pair_ji);
    DEL_2D(type2phi_pair_ij);
    DEL_2D(type2phi_pair_ji);
    delete [] cut;
    delete [] cut_sq;
    
    type2phi=NULL;
    type2rho=NULL;
    type2rho_pair_ij=NULL;
    type2rho_pair_ji=NULL;
    type2phi_pair_ij=NULL;
    type2phi_pair_ji=NULL;
    cut=NULL;
    cut_sq=NULL;
    
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
        {
            for(int i=0;i<nr;i++)
                delete [] phi_r_arr[ityp][i];
            delete [] phi_r_arr[ityp];
        }
        delete [] phi_r_arr;
        
        for(int ityp=0;ityp<no_types;ityp++)
        {
            for(int i=0;i<nr;i++)
                delete [] rho_arr[ityp][i];
            delete [] rho_arr[ityp];
        }
        delete [] rho_arr;
        
        
        for(int ityp=0;ityp<no_types;ityp++)
        {
            for(int i=0;i<nrho;i++)
                delete [] F_arr[ityp][i];
            delete [] F_arr[ityp];
        }
        delete [] F_arr;
    }
    else if(eam_mode==FINNIS_FL)
    {
        for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
        {
            for(int i=0;i<nr;i++)
                delete [] phi_r_arr[ityp][i];
            delete [] phi_r_arr[ityp];
        }
        delete [] phi_r_arr;
        
        for(int ityp=0;ityp<no_types*no_types;ityp++)
        {
            for(int i=0;i<nr;i++)
                delete [] rho_arr[ityp][i];
            delete [] rho_arr[ityp];
        }
        delete [] rho_arr;
        
        for(int ityp=0;ityp<no_types;ityp++)
        {
            for(int i=0;i<nrho;i++)
                delete [] F_arr[ityp][i];
            delete [] F_arr[ityp];
        }
        delete [] F_arr;
    }
    
    allocated=0;
}
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void EAMFileReader::allocate()
{
    CREATE1D(cut,no_types*(no_types+1)/2);
    CREATE1D(cut_sq,no_types*(no_types+1)/2);
    CREATE_2D(type2phi,no_types,no_types);
    CREATE_2D(type2rho,no_types,no_types);
    CREATE_2D(type2rho_pair_ij,no_types,no_types);
    CREATE_2D(type2rho_pair_ji,no_types,no_types);
    CREATE_2D(type2phi_pair_ij,no_types,no_types);
    CREATE_2D(type2phi_pair_ji,no_types,no_types);
    
    for(int ityp=0;ityp<no_types;ityp++)
        for(int jtyp=0;jtyp<no_types;jtyp++)
        {
            type2phi[ityp][jtyp]=COMP(ityp,jtyp);
            type2rho[ityp][jtyp]=jtyp*no_types+ityp;
        }

    stride=0;
    
    for(int ityp=0;ityp<no_types;ityp++)
        for(int jtyp=0;jtyp<no_types;jtyp++)
        {
            if(ityp==jtyp)
            {
                type2rho_pair_ij[ityp][jtyp]=type2rho_pair_ji[jtyp][ityp]=stride;
                stride++;
                type2phi_pair_ij[ityp][jtyp]=type2phi_pair_ji[jtyp][ityp]=stride;
                stride++;
            }
            else
            {
                
                type2rho_pair_ij[ityp][jtyp]=stride;
                stride++;
                type2rho_pair_ji[jtyp][ityp]=stride;
                stride++;
                type2phi_pair_ij[ityp][jtyp]=type2phi_pair_ji[jtyp][ityp]=stride;
                stride++;
                
            }
        }
    
    
    int rho_phi_ncomp=0;
    int F_ncomp=0;
    
    if(mode==MD_mode)
    {
        rho_phi_ncomp=7;
        F_ncomp=7;
    }
    else if(mode==DMD_mode)
    {
        rho_phi_ncomp=4;
        F_ncomp=5;
    }
    
    
    CREATE1D(phi_r_arr,no_types*(no_types+1)/2);
    for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
        CREATE1D(phi_r_arr[ityp],nr);
    
    for(int ityp=0;ityp<no_types*(no_types+1)/2;ityp++)
        for(int i=0;i<nr;i++)
            CREATE1D(phi_r_arr[ityp][i],rho_phi_ncomp);
    
    CREATE1D(rho_arr,no_types*no_types);
    rho_r_arr=rho_arr;
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        for(int ityp=0;ityp<no_types;ityp++)
            CREATE1D(rho_arr[ityp],nr);
        for(int ityp=0;ityp<no_types;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(rho_arr[type2rho[ityp][0]][i],rho_phi_ncomp);
        
        
        for(int ityp=0;ityp<no_types;ityp++)
            for(int jtyp=0;jtyp<no_types;jtyp++)
                rho_arr[type2rho[ityp][jtyp]]
                =rho_arr[type2rho[ityp][0]];
    }
    else if(eam_mode==FINNIS_FL)
    {
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            CREATE1D(rho_arr[ityp],nr);
        for(int ityp=0;ityp<no_types*no_types;ityp++)
            for(int i=0;i<nr;i++)
                CREATE1D(rho_arr[ityp][i],rho_phi_ncomp);
    }
    
    CREATE1D(F_arr,no_types);
    for(int ityp=0;ityp<no_types;ityp++)
        CREATE1D(F_arr[ityp],nrho);
    
    for(int ityp=0;ityp<no_types;ityp++)
        for(int i=0;i<nrho;i++)
            CREATE1D(F_arr[ityp][i],F_ncomp);
    

    allocated=1;
}
/*--------------------------------------------
 set the arrays
 --------------------------------------------*/
void EAMFileReader::set_arrays()
{
    if(mode==MD_mode)
    {
        for(int itype=0;itype<no_types;itype++)
            interpolate_7(nrho,drho,F_arr[itype]);
        
        for(int itype=0;itype<no_types;itype++)
            for(int jtype=0;jtype<itype+1;jtype++)
                interpolate_7(nr,dr,phi_r_arr[type2phi[itype][jtype]]);
        
        if(eam_mode==FUNC_FL || eam_mode==SET_FL)
        {
            for(int itype=0;itype<no_types;itype++)
                interpolate_7(nr,dr,rho_arr[type2rho[itype][0]]);
        }
        else if(eam_mode==FINNIS_FL)
        {
            for(int itype=0;itype<no_types;itype++)
                for(int jtype=0;jtype<no_types;jtype++)
                    interpolate_7(nr,dr,rho_arr[type2rho[itype][jtype]]);
        }
        
    }
    else if(mode==DMD_mode)
    {
        for(int itype=0;itype<no_types;itype++)
            interpolate_5(nrho,drho,F_arr[itype]);
        
        for(int itype=0;itype<no_types;itype++)
            for(int jtype=0;jtype<itype+1;jtype++)
                interpolate_4(nr,dr,phi_r_arr[type2phi[itype][jtype]]);
        
        if(eam_mode==FUNC_FL || eam_mode==SET_FL)
        {
            for(int itype=0;itype<no_types;itype++)
                for(int k=0;k<nr;k++)
                    rho_arr[type2rho[itype][0]][k][0]*=static_cast<type0>(k)*dr;
            
            for(int itype=0;itype<no_types;itype++)
                interpolate_4(nr,dr,rho_arr[type2rho[itype][0]]);
        }
        else if(eam_mode==FINNIS_FL)
        {
         
            for(int itype=0;itype<no_types;itype++)
                for(int jtype=0;jtype<no_types;jtype++)
                    for(int k=0;k<nr;k++)
                        rho_arr[type2rho[itype][jtype]][k][0]*=static_cast<type0>(k)*dr;
          
            for(int itype=0;itype<no_types;itype++)
                for(int jtype=0;jtype<no_types;jtype++)
                    interpolate_4(nr,dr,rho_arr[type2rho[itype][jtype]]);
        }
    }
    
    rc=(static_cast<type0>(nr)-1.0)*dr;
    rho_max=(static_cast<type0>(nrho)-1.0)*drho;
    type0 r_cut;
    for(int itype=0;itype<no_types;itype++)
        for(int jtype=itype;jtype<no_types;jtype++)
        {
            type0 ** phi_arr_=phi_r_arr[type2phi[itype][jtype]];
            type0 ** rho_arr_0=rho_arr[type2rho[itype][jtype]];
            type0 ** rho_arr_1=rho_arr[type2rho[jtype][itype]];
            int nr_=nr-1;
            while(phi_arr_[nr_][0]==0.0 &&
                  rho_arr_0[nr_][0]==0.0 &&
                  rho_arr_1[nr_][0]==0.0)
            {
                nr_--;
            }
            r_cut=static_cast<type0>((nr_))*dr;
            cut[COMP(itype,jtype)]=r_cut;
            cut_sq[COMP(itype,jtype)]=r_cut*r_cut;
        }

}
/*--------------------------------------------
 interpolate for a single array
 --------------------------------------------*/
type0 EAMFileReader::interpolate(type0* arr
,int n,type0 p,int k)
{
    
    type0 coef0,coef1,coef2,coef3,tmp;
    
    coef0=arr[k];
    
    if(k==0)
    {
        coef1=arr[1]-arr[0];
        tmp=0.5*(arr[2]-arr[0]);
        coef2=(-arr[0]+2.0*arr[1]-arr[2])/2.0;
        coef3=(arr[0]-2.0*arr[1]+arr[2])/2.0;
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }
    else if(k==1)
    {
        coef1=0.5*(arr[2]-arr[0]);
        tmp=((arr[0]-arr[4])+8.0*(arr[3]-arr[1]))/12.0;
        coef2=(11.0*arr[0]-28.0*arr[1]+24.0*arr[2]-8.0*arr[3]+arr[4])/12.0;
        coef3=(-5.0*arr[0]+16.0*arr[1]-18.0*arr[2]+8.0*arr[3]-arr[4])/12.0;
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }
    else if(k==n-2)
    {
        coef1=0.5*(arr[n-1]-arr[n-3]);
        tmp=arr[n-1]-arr[n-2];
        coef2=arr[n-3]-3.0*arr[n-2]+2.0*arr[n-1];
        coef3=0.5*(-arr[n-3]+3.0*arr[n-2]-2.0*arr[n-1]);
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }
    else if(k==n-1)
    {
        coef1=arr[n-1]-arr[n-2];
        return coef0+coef1*p;
    }
    else
    {
        coef1=((arr[k-2]-arr[k+2])+
               8.0*(arr[k+1]-arr[k-1]))/12.0;
        tmp=((arr[k-1]-arr[k+3])+
             8.0*(arr[k+2]-arr[k]))/12.0;
        
        coef2=(-2.0*arr[k-2]+15.0*arr[k-1]-28.0*arr[k]+21.0*arr[k+1]-6.0*arr[k+2])/12.0;
        coef3=(arr[k-2]-7.0*arr[k-1]+16.0*arr[k]-17.0*arr[k+1]+7.0*arr[k+2])/12.0;
        return ((coef3*p+coef2)*p+coef1)*p+coef0;
    }
}
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void EAMFileReader::interpolate_4(int n
,type0 delta,type0** spline)
{
    spline[0][1]=spline[1][0]-spline[0][0];
    spline[1][1]=0.5*(spline[2][0]-spline[0][0]);
    spline[n-2][1]=0.5*(spline[n-1][0]-spline[n-3][0]);
    spline[n-1][1]=spline[n-1][0]-spline[n-2][0];
    
    for(int i=2;i<n-2;i++)
        spline[i][1]=((spline[i-2][0]-spline[i+2][0])+
        8.0*(spline[i+1][0]-spline[i-1][0]))/12.0;
    
    for(int i=0;i<n-1;i++)
    {
        spline[i][2]=3.0*(spline[i+1][0]-spline[i][0])-
        2.0*spline[i][1]-spline[i+1][1];
        spline[i][3]=spline[i][1]+spline[i+1][1]-
        2.0*(spline[i+1][0]-spline[i][0]);
    }
    
    spline[n-1][2]=0.0;
    spline[n-1][3]=0.0;

}
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void EAMFileReader::interpolate_5(int n
,type0 delta,type0** spline)
{
    spline[0][1]=spline[1][0]-spline[0][0];
    spline[1][1]=0.5*(spline[2][0]-spline[0][0]);
    spline[n-2][1]=0.5*(spline[n-1][0]-spline[n-3][0]);
    spline[n-1][1]=spline[n-1][0]-spline[n-2][0];
    
    for(int i=2;i<n-2;i++)
        spline[i][1]=((spline[i-2][0]-spline[i+2][0])+
        8.0*(spline[i+1][0]-spline[i-1][0]))/12.0;
    
    for(int i=0;i<n-1;i++)
    {
        spline[i][2]=3.0*(spline[i+1][0]-spline[i][0])-
        2.0*spline[i][1]-spline[i+1][1];
        spline[i][3]=spline[i][1]+spline[i+1][1]-
        2.0*(spline[i+1][0]-spline[i][0]);
    }
    
    
    spline[n-1][2]=0.0;
    spline[n-1][3]=0.0;
    
    
    
    for(int i=0;i<n-1;i++)
    {
        spline[i][4]=(spline[i+1][2]-spline[i][2])/6.0-0.5*spline[i][3];
    }
    spline[n-1][4]=0.0;

}
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void EAMFileReader::interpolate_7(int n,type0 delta
,type0** spline)
{
    spline[0][1]=spline[1][0]-spline[0][0];
    spline[1][1]=0.5*(spline[2][0]-spline[0][0]);
    spline[n-2][1]=0.5*(spline[n-1][0]-spline[n-3][0]);
    spline[n-1][1]=spline[n-1][0]-spline[n-2][0];
    
    for(int i=2;i<n-2;i++)
        spline[i][1]=((spline[i-2][0]-spline[i+2][0])+
            8.0*(spline[i+1][0]-spline[i-1][0]))/12.0;
    
    for(int i=0;i<n-1;i++)
    {
        spline[i][2]=3.0*(spline[i+1][0]-spline[i][0])-
        2.0*spline[i][1]-spline[i+1][1];
        spline[i][3]=spline[i][1]+spline[i+1][1]-
        2.0*(spline[i+1][0]-spline[i][0]);
    }
    
    spline[n-1][2]=0.0;
    spline[n-1][3]=0.0;
    for(int i=0;i<n;i++)
    {
        spline[i][4]=spline[i][1]/delta;
        spline[i][5]=2.0*spline[i][2]/delta;
        spline[i][6]=3.0*spline[i][3]/delta;
    }
        
}
