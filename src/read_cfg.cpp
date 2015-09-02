/*--------------------------------------------
 Created by Sina on 07/23/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <stdlib.h>
#include "read_cfg.h"
#include "atom_types.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Read_cfg::Read_cfg(MAPP* mapp,int narg,char** args)
:Read(mapp)
{
    
    if(narg!=3)
        error->abort("read cfg should only have 1 argument");
    file_name=args[2];
    if(atoms->dimension!=3)
        error->abort("read cfg works only with 3 dimensional boxes");
    
    curr_id=0;
    
    
    // setup the defaults
    basic_length=1.0;
    R=1.0;
    entry_count=6;
    ext_cfg=0;
    vel_chk=1;

    CREATE2D(H_x,3,3);
    CREATE2D(H_x_d,3,3);
    CREATE2D(H0,3,3);
    CREATE2D(eta,3,3);
    CREATE2D(eta_sq,3,3);
    CREATE2D(trns,3,3);
    
    M3ZERO(H_x);
    M3ZERO(H_x_d);
    M3ZERO(H0);
    M3ZERO(eta);
    M3ZERO(trns);
    
    trns[0][0]=trns[1][1]=trns[2][2]=1.0;
    H0[0][0]=H0[1][1]=H0[2][2]=1.0;
    // end of defaults
    
    
    CREATE1D(line,MAXCHAR);
    if(atoms->my_p_no==0)
    {
        cfgfile=fopen(file_name,"r");
        if(cfgfile==NULL)
            error->abort("file %s not found",file_name);
    }
    
    
    fpos_t pos;
    header_cmplt=0;
    while (!header_cmplt)
    {
        if(atoms->my_p_no==0)
        {
            fgetpos(cfgfile,&pos);
            fgets(line,MAXCHAR,cfgfile);
        }
        MPI_Bcast(line,MAXCHAR,MPI_CHAR,0,world);
        read_header();
    }
    if(atoms->my_p_no==0)
        fsetpos(cfgfile,&pos);
    
    set_box();
    
    atoms->auto_grid_proc();
 
    if(mapp->mode==MD_mode)
    {
        x_n=atoms->add<type0>(1, 3,"x");
        id_n=atoms->add<int>(0, 1,"id");
        type_n=atoms->add<int>(1, 1,"type");
        /*
        x_n=atoms->find("x");
        id_n=atoms->find("id");
        type_n=atoms->find("type");
        */
        if(vel_chk)
        {
            x_d_n=atoms->add<type0>(0, 3,"x_d");
            vec_list=new VecLst(mapp,4,x_n,x_d_n,type_n,id_n);
            
            ch_id=atoms->vectors[x_n]->byte_size;
            ch_type=ch_id+atoms->vectors[id_n]->byte_size;
            ch_x_d=ch_type+atoms->vectors[type_n]->byte_size;
            tmp_buff_size=6;
            CREATE1D(tmp_buff,tmp_buff_size);
        }
        else
        {
            
            vec_list=new VecLst(mapp,3,x_n,type_n,id_n);
            
            ch_id=atoms->vectors[x_n]->byte_size;
            ch_type=ch_id+atoms->vectors[id_n]->byte_size;
            tmp_buff_size=3;
            CREATE1D(tmp_buff,tmp_buff_size);
        }

    }
    else if(mapp->mode==DMD_mode)
    {
        if((entry_count-3)%2!=0)
            error->abort("entry_count-3 in %s file should be divisible by 2",file_name);
        dmd_no_types=(entry_count-3)/2;
        
        
        x_n=atoms->add<type0>(1,3+dmd_no_types,"x");
        id_n=atoms->add<int>(0, 1,"id");
        c_n=atoms->add<type0>(1,dmd_no_types,"c");
        vec_list=new VecLst(mapp,3,x_n,id_n,c_n);
        
        ch_id=atoms->vectors[x_n]->byte_size;
        ch_c=ch_id+atoms->vectors[id_n]->byte_size;
        
        tmp_buff_size=3+2*dmd_no_types;
        CREATE1D(tmp_buff,tmp_buff_size);
        
        if(ext_cfg==0)
            error->abort("based on entry_count in %s file only %d types of atoms should be present",file_name,dmd_no_types);
    }
    
    CREATE1D(ch_buff,vec_list->byte_size);
    
    atom_cmplt=0;
    while (!atom_cmplt)
    {
        if(atoms->my_p_no==0)
            fgets(line,MAXCHAR,cfgfile);
        MPI_Bcast(line,MAXCHAR,MPI_CHAR,0,world);
        if(mapp->mode==MD_mode)
            read_atom_md();
        else if(mapp->mode==DMD_mode)
            read_atom_dmd();
    }
    
    if(mapp->mode==DMD_mode)
    {
        if(atom_types->no_types!=dmd_no_types)
            error->abort("%s file should conatin properties of %d atom types",file_name,dmd_no_types);
    }
    
    
    int loc_no=atoms->natms;
    int tot_no=0;
    MPI_Allreduce(&loc_no,&tot_no,1,MPI_INT,MPI_SUM,world);
    
    if (tot_no!=atoms->tot_natms)
        error->abort("the number of atoms dont match"
                     " recheck your cfg file: %d %d",
                     atoms->tot_natms,tot_no);
    
    if(atoms->my_p_no==0)
        fclose(cfgfile);
    
    atoms->s2x(atoms->natms);
    
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Read_cfg::~Read_cfg()
{
    delete vec_list;
    delete [] line;
    
    
    for(int i=0;i<3;i++)
    {
        delete [] eta_sq[i];
        delete [] eta[i];
        delete [] H0[i];
        delete [] H_x[i];
        delete [] H_x_d[i];
        delete [] trns[i];
    }
    
    delete [] eta_sq;
    delete [] eta;
    delete [] H0;
    delete [] H_x;
    delete [] H_x_d;
    delete [] trns;
    
    delete [] ch_buff;
    
    if(tmp_buff_size)
        delete [] tmp_buff;
}
/*--------------------------------------------
 reads the header of the cfg file
 --------------------------------------------*/
void Read_cfg::read_header()
{
    char* command;
    int narg = mapp->hash_remover(line,command);
    type0 tmp;
    int icmp,jcmp,tmpno;
    char* strtmp1;
    char* strtmp2;
    char* strtmp3;
    
    CREATE1D(strtmp1,MAXCHAR);
    CREATE1D(strtmp2,MAXCHAR);
    CREATE1D(strtmp3,MAXCHAR);
    
    
    if(strcmp(command,".NO_VELOCITY.")==0)
    {
        vel_chk=0;
        ext_cfg=1;
        entry_count-=3;
    }
    else if(sscanf(command,"Transform(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
    {
        icmp--;
        jcmp--;
        if (icmp>2 || icmp<0)
            error->abort("wrong component in %s file for Transform(%d,%d)",file_name,icmp+1,jcmp+1);
        if(jcmp>2 || jcmp<0)
            error->abort("wrong component in %s file for Transform(%d,%d)",file_name,icmp+1,jcmp+1);
        trns[icmp][jcmp]=tmp;
    }
    else if(sscanf(command,"eta(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
    {
        icmp--;
        jcmp--;
        if (icmp>2 || icmp<0)
            error->abort("wrong component in %s file for eta(%d,%d)",file_name,icmp+1,jcmp+1);
        if(jcmp>2 || jcmp<0)
            error->abort("wrong component in %s file for eta(%d,%d)",file_name,icmp+1,jcmp+1);
        eta[icmp][jcmp]=tmp;
    }
    else if(sscanf(command,"entry_count = %d",&tmpno)==1)
    {
        entry_count=tmpno;
        ext_cfg=1;
        int mincomp=3+(3*vel_chk);
        if (entry_count < mincomp)
            error->abort("entry_count in %s should at least be equal to %d",file_name,mincomp);
    }
    else if(sscanf(command,"H0(%d,%d) = %lf A",&icmp,&jcmp,&tmp)==3)
    {
        icmp--;
        jcmp--;
        if (icmp>2 || icmp<0)
            error->abort("wrong component in %s file for H(%d,%d)",file_name,icmp+1,jcmp+1);
        if(jcmp>2 || jcmp<0)
            error->abort("wrong component in %s file for H(%d,%d)",file_name,icmp+1,jcmp+1);
        H0[icmp][jcmp]=tmp;
    }
    else if(sscanf(command,"R = %lf %s",&tmp,strtmp1)==2)
    {
        R=tmp;
    }
    else if(sscanf(command,"Number of particles = %d",&tmpno)==1)
        atoms->tot_natms=tmpno;
    else if(sscanf(command,"auxiliary[%d] = %s [%s]",&icmp,strtmp1,strtmp2)==3)
    {
        int mincomp=3+(3*vel_chk);
        if(icmp+mincomp+1>entry_count)
            error->abort("wrong component in %s file for auxiliary[%d], %d+%d+1 > entry_count",file_name,icmp,mincomp,icmp);
    }
    else if(sscanf(command,"A = %lf Angstrom (basic length-scale)",&tmp)==1)
    {
        if(tmp<=0.0)
            error->abort("A in %s file should be greater than 0.0",file_name);
        basic_length=tmp;
    }
    else
    {
        if (narg==8&&ext_cfg==0)
            header_cmplt=1;
        else if (narg==1&&ext_cfg)
            header_cmplt=1;
        else if (narg==0){}
        else
            error->warning("invalid line in %s file: %s",file_name,command);
    }
    
    if(narg)
        delete [] command;
    delete [] strtmp3;
    delete [] strtmp2;
    delete [] strtmp1;    
    
}
/*--------------------------------------------
 calculates H from H0, Transform, and eta;
 ** make sure H is zeroed before;
 for now we disregard eta;
 remeber to fix it later;
 --------------------------------------------*/
void Read_cfg::set_box()
{
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            for (int k=0;k<3;k++)
                H_x[i][j]+=H0[i][k]*trns[k][j];
    
    int chk=1;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            if(eta[i][j]!=0.0)
                chk=0;
    if(chk==0)
    {
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                eta[i][j]*=2.0;
        for(int i=0;i<3;i++)
            eta[i][i]++;
        
        M3sqroot(eta,eta_sq);
        
        M3EQV(H_x,H0);
        M3ZERO(H_x);
        for (int i=0;i<3;i++)
            for (int j=0;j<3;j++)
                for (int k=0;k<3;k++)
                    H_x[i][j]+=H0[i][k]*eta_sq[k][j];
    }
    
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
        {
            H_x[i][j]*=basic_length;
            atoms->H[i][j]=H_x[i][j];
        }
    M3EQV(H_x,H_x_d);
    
    if (M3DET(H_x)==0.0)
        error->abort("determinant of H in %s file is 0.0",file_name);
    
    XMath* xmath= new XMath(mapp);
    xmath->square2lo_tri(H_x,atoms->H);
    delete xmath;
    
    M3EQV(atoms->H,H_x);
    M3INV_TRI_LOWER(atoms->H,atoms->B);
    
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
void Read_cfg::read_atom_md()
{
    char** arg;
    type0 mass;
    int narg=mapp->parse_line(line,arg);
    
    if(atoms->my_p_no==0)
    {
        if(feof(cfgfile))
            atom_cmplt=1;
    }
    MPI_Bcast(&atom_cmplt,1,MPI_INT,0,world);
    if(atom_cmplt)
    {
        if(narg)
        {
            for(int i=0;i<narg;i++)
                delete [] arg[i];
            delete [] arg;
        }
        return;
    }
    
    
    if (narg!=8 && ext_cfg==0)
        error->abort("invalid line in %s file: %s",file_name,line);
    if (ext_cfg && !(narg==1 || narg==entry_count))
        error->abort("invalid line in %s file: %s",file_name,line);
    
    if(ext_cfg)
    {
        if(narg==1)
        {
            mass=static_cast<type0>(atof(arg[0]));
            if(narg)
            {
                for(int i=0;i<narg;i++)
                    delete [] arg[i];
                delete [] arg;
            }
            if(atoms->my_p_no==0)
                fgets(line,MAXCHAR,cfgfile);
            MPI_Bcast(line,MAXCHAR,MPI_CHAR,0,world);
            
            narg=mapp->parse_line(line,arg);
            if(narg!=1)
                error->abort("invalid line in %s file: %s",file_name,line);
            if(mass<=0.0)
                error->abort("mass of %s %s file (%lf) should be greater than 0.0",arg[0],file_name,line,mass);
            last_type=
            atom_types->add_type(mass,arg[0]);
        }
        else if(narg==entry_count)
        {
            if (vel_chk)
            {
                for (int i=0;i<6;i++)
                    tmp_buff[i]=0.0;
                
                for (int i=0;i<3;i++)
                    tmp_buff[i]=static_cast<type0>(atof(arg[i]));
                for (int i=3;i<6;i++)
                    tmp_buff[i]=static_cast<type0>(atof(arg[i]))*R;
                add_atom_read_x(last_type);
            }
            else
            {
                for (int i=0;i<3;i++)
                    tmp_buff[i]=0.0;
                
                for (int i=0;i<3;i++)
                    tmp_buff[i]=static_cast<type0>(atof(arg[i]));
                add_atom_read_x(last_type);
            }
            
        }
        else
            error->abort("invalid line in %s file: %s",file_name,line);
    }
    else
    {
        if(narg==8)
        {
            mass=static_cast<type0>(atof(arg[0]));
            last_type=
            atom_types->add_type(mass,arg[1]);
            if(mass<=0.0)
                error->abort("mass of %s %s file (%lf) should be greater than 0.0",arg[1],file_name,line,mass);
 
            for (int i=0;i<3;i++)
                tmp_buff[i]=static_cast<type0>(atof(arg[i+2]));
            for (int i=3;i<6;i++)
                tmp_buff[i]=static_cast<type0>(atof(arg[i+2]))*R;
            add_atom_read_x(last_type);

        }
        else
            error->abort("invalid line in %s file: %s",file_name,line);
    }
    
    
    for(int i=0;i<narg;i++)
        delete [] arg[i];
    if(narg)
        delete [] arg;
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
void Read_cfg::read_atom_dmd()
{
    char** arg;
    type0 mass;
    int narg=mapp->parse_line(line,arg);
    
    if(atoms->my_p_no==0)
    {
        if(feof(cfgfile))
            atom_cmplt=1;
    }
    MPI_Bcast(&atom_cmplt,1,MPI_INT,0,world);
    if(atom_cmplt)
    {
        if(narg)
        {
            for(int i=0;i<narg;i++)
                delete [] arg[i];
            delete [] arg;
        }
        return;
    }
    
    

    if (ext_cfg && !(narg==1 || narg==entry_count))
        error->abort("invalid line in %s file: %s",file_name,line);
    

    if(narg==1)
    {
        mass=static_cast<type0>(atof(arg[0]));
        if(narg)
        {
            for(int i=0;i<narg;i++)
                delete [] arg[i];
            delete [] arg;
        }
        if(atoms->my_p_no==0)
            fgets(line,MAXCHAR,cfgfile);
        MPI_Bcast(line,MAXCHAR,MPI_CHAR,0,world);
        
        narg=mapp->parse_line(line,arg);
        if(narg!=1)
            error->abort("invalid line in %s file: %s",file_name,line);
        if(mass<=0.0)
            error->abort("mass of %s %s file (%lf) should be greater than 0.0",arg[0],file_name,line,mass);
        last_type=atom_types->add_type(mass,arg[0]);
        if(last_type>=dmd_no_types)
            error->abort("based on entry_count in %s file only %d types of atoms should be present",file_name,dmd_no_types);
    }
    else if(narg==entry_count)
    {
        for (int i=0;i<entry_count;i++)
            tmp_buff[i]=static_cast<type0>(atof(arg[i]));
        for(int i=3;i<3+dmd_no_types;i++)
            if(tmp_buff[i]<=0.0)
                error->abort("values of alpha vector in %s file should be greater than 0.0",file_name);
        for(int i=3+dmd_no_types;i<3+2*dmd_no_types;i++)
            if(tmp_buff[i]<-1.0 || tmp_buff[i]>1.0)
                error->abort("values of c vector in %s file should be between 0.0 & 1.0",file_name);
        
        add_atom_read_x();
        
    }
    else
        error->abort("invalid line in %s file: %s",file_name,line);

    
    for(int i=0;i<narg;i++)
        delete [] arg[i];
    if(narg)
        delete [] arg;
}
/*--------------------------------------------
 addatom_read_x
 --------------------------------------------*/
void Read_cfg::add_atom_read_x(int t)
{
    
    for(int i=0;i<3;i++)
    {
        while(tmp_buff[i]>=1)
            tmp_buff[i]--;
        while(tmp_buff[i]<0)
            tmp_buff[i]++;
    }
    
    if(vel_chk)
    {
        type0* x_d;
        CREATE1D(x_d,3);
        
        V3ZERO(x_d);
        for (int i=0;i<3;i++)
            for (int j=0;j<3;j++)
                x_d[j]+=tmp_buff[i+3]*H_x_d[i][j]*R;
        for (int i=0;i<3;i++)
            tmp_buff[i+3]=x_d[i];
        
        delete [] x_d;
    }
    
    
    for(int i=0;i<3;i++)
        if(!(atoms->s_lo[i]<=tmp_buff[i]
             && tmp_buff[i]<atoms->s_hi[i]))
        {
            curr_id++;
            return;
        }
    
    memcpy(ch_buff,tmp_buff,3*sizeof(type0));
    memcpy(&ch_buff[ch_type],&t,sizeof(int));
    memcpy(&ch_buff[ch_id],&curr_id,sizeof(int));
    
    if(vel_chk)
        memcpy(&ch_buff[ch_x_d],&tmp_buff[3],3*sizeof(type0));
    
    atoms->unpack_read(ch_buff,1,vec_list);
    curr_id++;
}
/*--------------------------------------------
 addatom_read_x
 --------------------------------------------*/
void Read_cfg::add_atom_read_x()
{
    for(int i=0;i<3;i++)
    {
        while(tmp_buff[i]>=1.0)
            tmp_buff[i]--;
        while(tmp_buff[i]<0.0)
            tmp_buff[i]++;
    }
    
    for(int i=0;i<3;i++)
        if(!(atoms->s_lo[i]<=tmp_buff[i]
             && tmp_buff[i]<atoms->s_hi[i]))
        {
            curr_id++;
            return;
        }
    
    memcpy(ch_buff,tmp_buff,(3+dmd_no_types)*sizeof(type0));
    memcpy(&ch_buff[ch_id],&curr_id,sizeof(int));
    memcpy(&ch_buff[ch_c],&tmp_buff[3+dmd_no_types]
    ,dmd_no_types*sizeof(type0));
    
    atoms->unpack_read(ch_buff,1,vec_list);
    curr_id++;
}
/*--------------------------------------------
 calculates square root of 3x3 matrix
 ref: L. P. Franca
 An Algorithm to Compute The Square Root of
 a 3x3 Positive Definite Matrix
 Computers Math. Applic. Vol. 18, No. 5,
 pp. 459-466, 1989
 --------------------------------------------*/
void Read_cfg::M3sqroot(type0** A,type0** Asq)
{
    type0 IA=0;
    for(int i=0;i<3;i++)
        IA+=A[i][i];
    type0 IIA=0;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            IIA-=A[i][j]*A[j][i];
    
    IIA+=IA*IA;
    IIA*=0.5;
    type0 IIIA=M3DET(A);
    type0 k=IA*IA-3*IIA;
    if(k<0.0)
        error->abort("eta in %s should be positive definite",file_name);
    if(k<TOLERANCE)
    {
        if(IA<=0.0)
            error->abort("eta in %s should be positive definite",file_name);
        M3ZERO(Asq);
        for(int i=0;i<3;i++)
            Asq[i][i]=sqrt(IA/3.0);
        return;
    }
    
    type0 l=IA*(IA*IA -4.5*IIA)+13.5*IIIA;
    
    type0 temp=l/(k*sqrt(k));
    if(temp>1.0||temp<-1.0)
        error->abort("eta in %s should be positive definite",file_name);
    
    type0 phi=acos(temp);
    type0 lambda=sqrt((1.0/3.0)*(IA+2*sqrt(k)*cos(phi/3.0)));
    
    type0 IIIAsq=sqrt(IIIA);
    type0 y=-lambda*lambda+IA+2*IIIAsq/lambda;
    if(y<0.0)
        error->abort("eta in %s should be positive definite",file_name);
    type0 IAsq=lambda+sqrt(y);
    type0 IIAsq=0.5*(IAsq*IAsq-IA);
    
    type0 coef0=IAsq*IIAsq-IIIAsq;
    if(coef0==0)
        error->abort("eta in %s should be positive definite",file_name);
    coef0=1.0/coef0;
    
    M3ZERO(Asq);
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
                Asq[i][j]-=coef0*A[i][k]*A[k][j];
    
    type0 coef1=coef0*(IAsq*IAsq-IIAsq);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            Asq[i][j]+=coef1*A[i][j];
    
    type0 coef2=coef0*IAsq*IIIAsq;
    for(int i=0;i<3;i++)
        Asq[i][i]+=coef2;
    
}
