/*--------------------------------------------
 Created by Sina on 07/23/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <stdlib.h>
#include "read_cfg.h"
#include "atom_types.h"
#include "xmath.h"
#include "error.h"
#include "memory.h"

using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Read_cfg::Read_cfg(MAPP* mapp,int nargs,char** args)
:Read(mapp)
{
    
    if(nargs!=3)
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

    CREATE_2D(H_x,3,3);
    CREATE_2D(H_x_d,3,3);
    CREATE_2D(H0,3,3);
    CREATE_2D(eta,3,3);
    CREATE_2D(trns,3,3);
    
    M3ZERO(H_x);
    M3ZERO(H_x_d);
    M3ZERO(H0);
    M3ZERO(eta);
    M3ZERO(trns);
    
    trns[0][0]=trns[1][1]=trns[2][2]=1.0;
    H0[0][0]=H0[1][1]=H0[2][2]=1.0;
    // end of defaults
    
    
    FILE* fp=NULL;
    mapp->open_file(fp,file_name,"r");
    
    read_header(fp);
    set_box();
    atoms->auto_grid();
    set_vec_lst();
    read_atom(fp);
    
    if(atoms->my_p==0)
        fclose(fp);
    set_vecs();
    
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Read_cfg::~Read_cfg()
{
    
    DEL_2D(eta);
    DEL_2D(H0);
    DEL_2D(H_x);
    DEL_2D(H_x_d);
    DEL_2D(trns);
    
    if(ch_buff_sz)
        delete [] ch_buff;
    
    if(tmp_buff_size)
        delete [] tmp_buff;
    
    if(nvecs)
        delete [] vec_lst;
}
/*--------------------------------------------
 reads the header of the cfg file
 --------------------------------------------*/
void Read_cfg::read_header(FILE* fp)
{
    int nargs=0;
    char* line;
    CREATE1D(line,MAXCHAR);
    
    int icmp,jcmp,tmpno;
    int header_cmplt=0;
    type0 tmp;
    fpos_t pos;
    if(atoms->my_p==0)
        fgetpos(fp,&pos);
    
    while(!header_cmplt && mapp->read_line(fp,line)!=-1)
    {
        
        nargs=mapp->hash_remover(line);
        if((nargs==8 && ext_cfg==0) || (nargs==1 && ext_cfg))
        {
            header_cmplt=1;
            if(atoms->my_p==0)
                fsetpos(fp,&pos);
            continue;
        }
        if(atoms->my_p==0)
            fgetpos(fp,&pos);
        
        if(nargs==0)
            continue;
        
        
        if(sscanf(line,"Number of particles = %d",&tmpno)==1)
        {
            if(tmpno<=0)
                error->abort("Number of particles in %s file should be greater than 0",file_name);
            
            atoms->tot_natms=tmpno;
        }
        else if(sscanf(line,"A = %lf Angstrom (basic length-scale)",&tmp)==1)
        {
            if(tmp<=0.0)
                error->abort("A in %s file should be greater than 0.0",file_name);
            if(tmp!=1.0)
            {
                error->warning("A (basic length-scale) is not equal to 1, mapp will set it to 1");
                tmp=1.0;
            }

            basic_length=tmp;
        }
        else if(strcmp(line,".NO_VELOCITY.")==0)
        {
            vel_chk=0;
            ext_cfg=1;
            entry_count-=3;
        }
        else if(sscanf(line,"entry_count = %d",&tmpno)==1)
        {
            entry_count=tmpno;
            ext_cfg=1;
            int mincomp=3+(3*vel_chk);
            if(entry_count < mincomp)
                error->abort("entry_count in %s should at least be equal to %d",file_name,mincomp);
        }
        else if(sscanf(line,"R = %lf %*s",&tmp)==1)
        {
            R=tmp;
        }
        else if(sscanf(line,"auxiliary[%d] = %*s [%*s]",&icmp)==1)
        {
            int mincomp=3+(3*vel_chk);
            if(icmp+mincomp+1>entry_count)
                error->abort("wrong component in %s file for auxiliary[%d], %d+%d+1 > entry_count",file_name,icmp,mincomp,icmp);
        }
        else if(sscanf(line,"Transform(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if(icmp>2 || icmp<0)
                error->abort("wrong component in %s file for Transform(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp>2 || jcmp<0)
                error->abort("wrong component in %s file for Transform(%d,%d)",file_name,icmp+1,jcmp+1);
            trns[icmp][jcmp]=tmp;
        }
        else if(sscanf(line,"eta(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if(icmp>2 || icmp<0)
                error->abort("wrong component in %s file for eta(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp>2 || jcmp<0)
                error->abort("wrong component in %s file for eta(%d,%d)",file_name,icmp+1,jcmp+1);
            eta[icmp][jcmp]=tmp;
        }
        else if(sscanf(line,"H0(%d,%d) = %lf A",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if(icmp>2 || icmp<0)
                error->abort("wrong component in %s file for H(%d,%d)",file_name,icmp+1,jcmp+1);
            if(jcmp>2 || jcmp<0)
                error->abort("wrong component in %s file for H(%d,%d)",file_name,icmp+1,jcmp+1);
            H0[icmp][jcmp]=tmp;
        }
        else
        {
            error->abort("invalid line in %s file: %s|",file_name,line);
        }
    }
    delete [] line;
    
    if(!header_cmplt)
        error->abort("file %s ended unexpectedly",file_name);
}
/*--------------------------------------------
 calculates H from H0, Transform, and eta;
 ** make sure H is zeroed before;
 for now we disregard eta;
 remeber to fix it later;
 --------------------------------------------*/
void Read_cfg::set_box()
{
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
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

        type0** eta_sq;
        
        CREATE_2D(eta_sq,3,3);
        XMath* xmath=new XMath();
        if(xmath->M3sqroot(eta,eta_sq)==0)
            error->abort("eta in %s should be positive definite",file_name);
        delete xmath;
        
        M3EQV(H_x,H0);
        M3ZERO(H_x);
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                for(int k=0;k<3;k++)
                    H_x[i][j]+=H0[i][k]*eta_sq[k][j];
        
        DEL_2D(eta_sq);
        
    }
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        {
            H_x[i][j]*=basic_length;
            atoms->H[i][j]=H_x[i][j];
        }
    M3EQV(H_x,H_x_d);
    
    if(M3DET(H_x)==0.0)
        error->abort("determinant of H in %s file is 0.0",file_name);
    
    XMath* xmath= new XMath();
    xmath->square2lo_tri(H_x,atoms->H);
    delete xmath;
    
    M3EQV(atoms->H,H_x);
    M3INV_LT(atoms->H,atoms->B);
    
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
void Read_cfg::read_atom(FILE* fp)
{
    char** args=NULL;
    int nargs;
    int args_cpcty=0;
    type0 mass=0.0;
    int mass_flag=1;
    md_type atm_type=0;
    int type_defined=0;
    
    char* line;
    CREATE1D(line,MAXCHAR);
    
    if(mapp->mode==MD_mode)
    {
        while(mapp->read_line(fp,line)!=-1)
        {
            nargs=mapp->parse_line(line,args,args_cpcty);
            
            if(nargs==0)
                continue;
            
            if(ext_cfg)
            {
                if(nargs!=1 && nargs!=entry_count)
                    error->abort("invalid line in %s file: %s",file_name,line);
                if(nargs==1)
                {
                    if(mass_flag)
                    {
                        mass=static_cast<type0>(atof(args[0]));
                        mass_flag=0;
                    }
                    else
                    {
                        atm_type=atom_types->add_type(mass,args[0]);
                        if(mass<=0.0)
                            error->abort("mass of %s %s file (%lf) should be greater than 0.0",args[0],file_name,line,mass);
                        mass_flag=1;
                        type_defined=1;
                    }
                }
                else
                {
                    if(type_defined==0)
                        error->abort("line %s in file %s comes before any element was defined",line,file_name);
                    
                    for(int i=0;i<3*(1+vel_chk);i++)
                        tmp_buff[i]=static_cast<type0>(atof(args[i]));
                    for(int i=3;i<3*(1+vel_chk);i++)
                        tmp_buff[i]*=R;
                    add_atom_read_x(atm_type);
                }
            }
            else
            {
                if(nargs!=8)
                    error->abort("invalid line in %s file: %s",file_name,line);
                mass=static_cast<type0>(atof(args[0]));
                atm_type=atom_types->add_type(mass,args[1]);
                
                for(int i=0;i<3*(1+vel_chk);i++)
                    tmp_buff[i]=static_cast<type0>(atof(args[i+2]));
                for(int i=3;i<3*(1+vel_chk);i++)
                    tmp_buff[i]*=R;
                add_atom_read_x(atm_type);
            }
        }
    }
    else if(mapp->mode==DMD_mode)
    {
        while(mapp->read_line(fp,line)!=-1)
        {
            nargs=mapp->parse_line(line,args,args_cpcty);
            
            if(nargs==0)
                continue;

            if(nargs!=1 && nargs!=entry_count)
                error->abort("invalid line in %s file: %s",file_name,line);
            if(nargs==1)
            {
                if(mass_flag)
                {
                    mass=static_cast<type0>(atof(args[0]));
                    mass_flag=0;
                }
                else
                {
                    atm_type=atom_types->add_type(mass,args[0]);
                    if(mass<=0.0)
                        error->abort("mass of %s %s file (%lf) should be greater than 0.0",args[0],file_name,line,mass);
                    mass_flag=1;
                }
            }
            else
            {
                for(int i=0;i<entry_count;i++)
                    tmp_buff[i]=static_cast<type0>(atof(args[i]));
                for(int i=3;i<3+dmd_no_types;i++)
                    if(tmp_buff[i]<0.0)
                        error->abort("values of alpha vector in %s file should be greater than 0.0",file_name);
                for(int i=3+dmd_no_types;i<3+2*dmd_no_types;i++)
                    if(tmp_buff[i]<-1.0 || tmp_buff[i]>1.0)
                        error->abort("values of c vector in %s file should be between 0.0 & 1.0",file_name);
                add_atom_read_x();
            }
        }
    }
    
    if(MAXCHAR)
        delete [] line;
    if(args_cpcty)
        delete [] args;
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
        type0 x_d[3];
        
        V3ZERO(x_d);
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                x_d[j]+=tmp_buff[i+3]*H_x_d[i][j]*R;
        for(int i=0;i<3;i++)
            tmp_buff[i+3]=x_d[i];
    }
    
    
    for(int i=0;i<3;i++)
        if(!(atoms->s_lo[i]<=tmp_buff[i]
             && tmp_buff[i]<atoms->s_hi[i]))
        {
            curr_id++;
            return;
        }
    
    memcpy(ch_buff,tmp_buff,3*sizeof(type0));
    if(mapp->mode==MD_mode)
    {
        md_type tt=static_cast<md_type>(t);
        memcpy(&ch_buff[ch_type],&tt,sizeof(md_type));
    }
    else if(mapp->mode==DMD_mode)
    {
        dmd_type tt=static_cast<md_type>(t);
        memcpy(&ch_buff[ch_type],&tt,sizeof(dmd_type));
    }
    
    memcpy(&ch_buff[ch_id],&curr_id,sizeof(int));
    
    if(vel_chk)
        memcpy(&ch_buff[ch_x_d],&tmp_buff[3],3*sizeof(type0));
    
    atoms->insert(ch_buff,vec_lst,nvecs,1);
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
    
    atoms->insert(ch_buff,vec_lst,nvecs,1);
    curr_id++;
}
/*--------------------------------------------
 addatom_read_x
 --------------------------------------------*/
void Read_cfg::set_vec_lst()
{
    if(mapp->mode==MD_mode)
    {
        mapp->x=new Vec<type0>(atoms,3);
        mapp->id=new Vec<int>(atoms,1);
        mapp->type=new Vec<md_type>(atoms,1);


        if(vel_chk)
        {
            if(mapp->x_d!=NULL)
                delete mapp->x_d;
            mapp->x_d=new Vec<type0>(atoms,3);
            
            ch_id=mapp->x->byte_sz;
            ch_type=ch_id+mapp->id->byte_sz;
            ch_x_d=ch_type+mapp->type->byte_sz;
            tmp_buff_size=6;
            CREATE1D(tmp_buff,tmp_buff_size);
            
            nvecs=4;
            CREATE1D(vec_lst,nvecs);
            vec_lst[0]=mapp->x;
            vec_lst[1]=mapp->id;
            vec_lst[2]=mapp->type;
            vec_lst[3]=mapp->x_d;
        }
        else
        {
            ch_id=mapp->x->byte_sz;
            ch_type=ch_id+mapp->id->byte_sz;
            tmp_buff_size=3;
            CREATE1D(tmp_buff,tmp_buff_size);

            nvecs=3;
            CREATE1D(vec_lst,nvecs);
            vec_lst[0]=mapp->x;
            vec_lst[1]=mapp->id;
            vec_lst[2]=mapp->type;
            
        }

    }
    else if(mapp->mode==DMD_mode)
    {
        if((entry_count-3)%2!=0)
            error->abort("entry_count-3 in %s file should be divisible by 2",file_name);
        dmd_no_types=(entry_count-3)/2;

        mapp->x=new Vec<type0>(atoms,3+dmd_no_types);
        mapp->id=new Vec<int>(atoms,1);

        if(mapp->c!=NULL)
            delete mapp->c;
        mapp->c=new Vec<type0>(atoms,dmd_no_types);
        
        ch_id=mapp->x->byte_sz;
        ch_c=ch_id+mapp->id->byte_sz;
        
        tmp_buff_size=3+2*dmd_no_types;
        CREATE1D(tmp_buff,tmp_buff_size);

        nvecs=3;
        CREATE1D(vec_lst,nvecs);
        vec_lst[0]=mapp->x;
        vec_lst[1]=mapp->id;
        vec_lst[2]=mapp->c;
        
        if(ext_cfg==0)
            error->abort("read cfg supports only extended format for DMD mode");
    }
    
    ch_buff_sz=0;
    for(int ivec=0;ivec<nvecs;ivec++)
        ch_buff_sz+=vec_lst[ivec]->byte_sz;
        
    CREATE1D(ch_buff,ch_buff_sz);
}
/*--------------------------------------------
 addatom_read_x
 --------------------------------------------*/
void Read_cfg::set_vecs()
{
        if(mapp->mode==DMD_mode)
    {
        if(atom_types->no_types!=dmd_no_types)
            error->abort("%s file should conatin properties of %d atom types",file_name,dmd_no_types);
    }
    
    int loc_no=atoms->natms;
    int tot_no=0;
    MPI_Allreduce(&loc_no,&tot_no,1,MPI_INT,MPI_SUM,world);
    
    if(tot_no!=atoms->tot_natms)
        error->abort("the number of atoms dont match"
        " recheck your cfg file: %d %d",
        atoms->tot_natms,tot_no);
    
    atoms->s2x_lcl();

    if(mapp->mode==DMD_mode)
    {
        type0* c=mapp->c->begin();
        
        int cdim_lcl=0,cdim;
        
        for(int iatm=0;iatm<atoms->natms;iatm++)
        {
            int d=0;
            for(int j=0;j<dmd_no_types;j++)
                if(c[iatm*dmd_no_types+j]!=-1.0)
                    d++;
            cdim_lcl=MAX(d,cdim_lcl);
        }
        MPI_Allreduce(&cdim_lcl,&cdim,1,MPI_INT,MPI_MAX,world);
        
        if(mapp->ctype!=NULL)
            delete mapp->ctype;
        mapp->ctype=new Vec<dmd_type>(atoms,dmd_no_types);
        
        dmd_type* ctype=mapp->ctype->begin();
        type0* x=mapp->x->begin();
        
        dmd_type icurs;
        
        for(int iatm=0;iatm<atoms->natms;iatm++)
        {
            for(int j=0;j<dmd_no_types;j++)
                ctype[iatm*dmd_no_types+j]=j;

            icurs=0;
            for(int j=0;j<dmd_no_types;j++)
                if(c[iatm*dmd_no_types+j]>=0.0)
                {
                    swap(c[iatm*dmd_no_types+j],c[iatm*dmd_no_types+icurs]);
                    swap(x[iatm*(dmd_no_types+3)+3+j],x[iatm*(dmd_no_types+3)+3+icurs]);
                    swap(ctype[iatm*dmd_no_types+j],ctype[iatm*dmd_no_types+icurs]);
                    icurs++;
                }
        }
        
        mapp->ctype->change_dimension(cdim);
        mapp->x->change_dimension(0.0,dmd_no_types,dmd_no_types-cdim);
        mapp->c->change_dimension(-1.0,dmd_no_types,dmd_no_types-cdim);
        
        mapp->c->assign_print_format("%e ");
    }
}

