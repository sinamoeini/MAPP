#include <stdlib.h>
#include "atom_types.h"
#include "write_cfg.h"
#define TRESH 0.00005
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Write_cfg::Write_cfg(MAPP* mapp,int narg
                     ,char** arg) : Write(mapp)
{
    usr_nabled=0;
    sorting=0;
    if(narg<4)
        error->abort("wrong write cfg command");
    
    write_step_tally=atoi(arg[2]);
    int lngth=static_cast<int>(strlen(arg[3]))+1;
    CREATE1D(file_name,lngth);
    for(int i=0;i<lngth;i++)
        file_name[i]=arg[3][i];
    
    
    
    
    if(mapp->mode==DMD_mode)
    {
        no_vecs=2;
        CREATE1D(vec_list,no_vecs);
        vec_list[0]=atoms->find("x");
        vec_list[1]=atoms->find("c");
        c_n=atoms->find("c");
        CREATE1D(vec_name,no_vecs);
        
        lngth=static_cast<int>(strlen((char*)"x"))+1;
        CREATE1D(vec_name[0],lngth);
        memcpy(vec_name[0],(char*)"x",lngth*sizeof(char));
        
        lngth=static_cast<int>(strlen((char*)"c"))+1;
        CREATE1D(vec_name[1],lngth);
        memcpy(vec_name[1],(char*)"c",lngth*sizeof(char));
    }
    else
    {
        no_vecs=1;
        CREATE1D(vec_list,no_vecs);
        vec_list[0]=atoms->find("x");
        CREATE1D(vec_name,no_vecs);

        lngth=static_cast<int>(strlen((char*)"x"))+1;
        CREATE1D(vec_name[0],lngth);
        memcpy(vec_name[0],(char*)"x",lngth*sizeof(char));
        
    }
    
    
    id_n=-1;
    int iarg=4;
    int tmp;
    
    while(iarg<narg)
    {
        if(strcmp(arg[iarg],"sort")==0)
        {
            sorting=1;
            iarg++;
        }
        else if(strcmp(arg[iarg],"usr")==0)
        {
            usr_nabled=1;
            iarg++;
        }
        else
        {
            if(strcmp(arg[iarg],"id")==0 ||
               strcmp(arg[iarg],"type")==0 ||
               strcmp(arg[iarg],"c")==0 ||
               strcmp(arg[iarg],"x")==0)
                error->abort("%s cannot be included in auxilary",arg[iarg]);
            
            tmp=atoms->find_exist(arg[iarg]);

            lngth=static_cast<int>(strlen(arg[iarg]))+1;
            GROW(vec_name,no_vecs,no_vecs+1);
            CREATE1D(vec_name[no_vecs],lngth);
            memcpy(vec_name[no_vecs],arg[iarg],lngth*sizeof(char));
            
            GROW(vec_list,no_vecs,no_vecs+1);
            vec_list[no_vecs]=tmp;
            no_vecs++;
            iarg++;
        }
    }
    
    if(mapp->mode==MD_mode)
    {
        type_n=atoms->find("type");
        GROW(vec_list,no_vecs,no_vecs+1);
        vec_list[no_vecs]=type_n;
        
        lngth=static_cast<int>(strlen((char*)"type"))+1;
        GROW(vec_name,no_vecs,no_vecs+1);
        CREATE1D(vec_name[no_vecs],lngth);
        memcpy(vec_name[no_vecs],(char*)"type",lngth*sizeof(char));
        
        no_vecs++;
    }
    
    if(sorting)
    {
        id_n=atoms->find("id");
        GROW(vec_list,no_vecs,no_vecs+1);
        vec_list[no_vecs]=id_n;
        
        lngth=static_cast<int>(strlen((char*)"id"))+1;
        GROW(vec_name,no_vecs,no_vecs+1);
        CREATE1D(vec_name[no_vecs],lngth);
        memcpy(vec_name[no_vecs],(char*)"id",lngth*sizeof(char));
        no_vecs++;
    }
    
    x_dim=atoms->vectors[0]->dim;
    dim=atoms->dimension;
    if(usr_nabled && mapp->mode==DMD_mode)
        CREATE1D(clr_r,4);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Write_cfg::~Write_cfg()
{
    if(usr_nabled && mapp->mode==DMD_mode)
        delete [] clr_r;
    delete [] file_name;
    
    for(int i=0;i<no_vecs;i++)
        delete [] vec_name[i];
    if(no_vecs)
    {
        delete [] vec_list;
        delete [] vec_name;
    }
}
/*--------------------------------------------
 write file
 --------------------------------------------*/
void Write_cfg::write_file(int stp)
{
    if(mapp->mode==MD_mode)
        write_file_md(stp);
    else
        write_file_dmd(stp);
}
/*--------------------------------------------
 write file
 --------------------------------------------*/
void Write_cfg::write_file_md(int stp)
{
    
    for(int i=0;i<no_vecs;i++)
        atoms->vectors[vec_list[i]]->gather_dump();
    
    
    FILE* fp=NULL;
    FILE* fp_usr=NULL;
    
    if(atoms->my_p_no==0)
    {
        x2s(atoms->tot_natms);
        char* filename;
        CREATE1D(filename,MAXCHAR);
        sprintf (filename, "%s.%010d.cfg",file_name,stp);
        fp=fopen(filename,"w");
        delete [] filename;
        if(usr_nabled)
        {
            CREATE1D(filename,MAXCHAR);
            sprintf (filename, "%s.%010d.usr",file_name,stp);
            fp_usr=fopen(filename,"w");
            delete [] filename;
        }

        
        // write the header
        fprintf(fp,"Number of particles = %d\n",atoms->tot_natms);
        fprintf(fp,"A = 1 Angstrom (basic length-scale)\n");
        
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                fprintf(fp,"H0(%d,%d) = %lf A\n",i+1,j+1,atoms->H[i][j]);
        
        fprintf(fp,".NO_VELOCITY.\n");
        
        fprintf(fp,"entry_count = %d\n",tot_dim);
        
        
        // write the body
        if(sorting)
        {
            int icomp=0;
            
            for(int i=1;i<no_vecs-2;i++)
            {
                if(atoms->vectors[vec_list[i]]->dim>1)
                {
                    for(int idim=0;idim<atoms->vectors[vec_list[i]]->dim;idim++)
                        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n"
                                ,icomp++,atoms->vectors[vec_list[i]]->name,idim);
                }
                else
                {
                    fprintf(fp,"auxiliary[%d] = %s [reduced unit]\n"
                            ,icomp++,atoms->vectors[vec_list[i]]->name);
                }
                
            }
            
            int* id;
            int tot_natms=atoms->tot_natms;
            atoms->vectors[id_n]->ret_dump(id);
            int* sort;
            CREATE1D(sort,tot_natms);
            for(int i=0;i<tot_natms;i++)
                sort[id[i]]=i;
            
            int* type;
            atoms->vectors[type_n]->ret_dump(type);
            
            for(int itype=0;itype<atom_types->no_types;itype++)
            {
                fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                        ,atom_types->atom_names[itype]);
                for(int i=0;i<tot_natms;i++)
                {
                    if(type[sort[i]]==itype)
                    {
                        for(int j=0;j<no_vecs-2;j++)
                            atoms->vectors[vec_list[j]]->print_dump(fp,sort[i]);
                        fprintf(fp,"\n");
                        if(usr_nabled)
                        {
                            for(int j=0;j<4;j++)
                                fprintf(fp_usr,"%lf ",atom_types->clr_rad[itype][j]);
                            fprintf(fp_usr,"\n");
                        }
                        
                    }
                }
            }
            
            if(tot_natms)
                delete [] sort;
        }
        else
        {
            int icomp=0;
            
            for(int i=1;i<no_vecs-1;i++)
            {
                if(atoms->vectors[vec_list[i]]->dim>1)
                {
                    for(int idim=0;idim<atoms->vectors[vec_list[i]]->dim;idim++)
                        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n"
                                ,icomp++,atoms->vectors[vec_list[i]]->name,idim);
                }
                else
                {
                    fprintf(fp,"auxiliary[%d] = %s [reduced unit]\n"
                            ,icomp++,atoms->vectors[vec_list[i]]->name);
                }
                
            }
            int tot_natms=atoms->tot_natms;
            
            int* type;
            atoms->vectors[type_n]->ret_dump(type);
            
            for(int itype=0;itype<atom_types->no_types;itype++)
            {
                fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                        ,atom_types->atom_names[itype]);
                for(int i=0;i<tot_natms;i++)
                {
                    if(type[i]==itype)
                    {
                        for(int j=0;j<no_vecs-1;j++)
                            atoms->vectors[vec_list[j]]->print_dump(fp,i);
                        fprintf(fp,"\n");
                        if(usr_nabled)
                        {
                            for(int j=0;j<4;j++)
                                fprintf(fp_usr,"%lf ",atom_types->clr_rad[itype][j]);
                            fprintf(fp_usr,"\n");
                        }
                    }
                }
            }
        }
        
        fclose(fp);
        if(usr_nabled)
            fclose(fp_usr);
    }
    
    for(int i=0;i<no_vecs;i++)
        atoms->vectors[vec_list[i]]->del_dump();
    
}
/*--------------------------------------------
 write file
 --------------------------------------------*/
void Write_cfg::write_file_dmd(int stp)
{
    
    for(int i=0;i<no_vecs;i++)
        atoms->vectors[vec_list[i]]->gather_dump();
    
    
    FILE* fp=NULL;
    FILE* fp_usr=NULL;
    if(atoms->my_p_no==0)
    {
        x2s(atoms->tot_natms);
        char* filename;
        CREATE1D(filename,MAXCHAR);
        sprintf (filename, "%s.%010d.cfg",file_name,stp);
        fp=fopen(filename,"w");
        delete [] filename;
        
        if(usr_nabled)
        {
            CREATE1D(filename,MAXCHAR);
            sprintf (filename, "%s.%010d.usr",file_name,stp);
            fp_usr=fopen(filename,"w");
            delete [] filename;
        }
        
        
        // write the header
        fprintf(fp,"Number of particles = %d\n",atoms->tot_natms);
        fprintf(fp,"A = 1 Angstrom (basic length-scale)\n");
        
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                fprintf(fp,"H0(%d,%d) = %lf A\n",i+1,j+1,atoms->H[i][j]);
        
        fprintf(fp,".NO_VELOCITY.\n");
        
        fprintf(fp,"entry_count = %d\n",tot_dim);
        
        
        // write the body
        if(sorting)
        {
            int icomp=0;
            int xx_dim=x_dim-3;
            if(xx_dim==1)
            {
                fprintf(fp,"auxiliary[%d] = alpha [reduced unit]\n",icomp++);
            }
            else if(xx_dim>1)
            {
                for(int idim=0;idim<xx_dim;idim++)
                    fprintf(fp,"auxiliary[%d] = alpha_%d [reduced unit]\n",icomp++,idim);
            }
            
            
            for(int i=1;i<no_vecs-1;i++)
            {
                if(atoms->vectors[vec_list[i]]->dim>1)
                {
                    for(int idim=0;idim<atoms->vectors[vec_list[i]]->dim;idim++)
                        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n"
                                ,icomp++,atoms->vectors[vec_list[i]]->name,idim);
                }
                else
                {
                    fprintf(fp,"auxiliary[%d] = %s [reduced unit]\n"
                            ,icomp++,atoms->vectors[vec_list[i]]->name);
                }
                
            }
            
            type0* c;
            atoms->vectors[c_n]->ret_dump(c);
            int* id;
            int tot_natms=atoms->tot_natms;
            atoms->vectors[id_n]->ret_dump(id);
            int* sort;
            CREATE1D(sort,tot_natms);
            for(int i=0;i<tot_natms;i++)
                sort[id[i]]=i;
            
            for(int itype=0;itype<atom_types->no_types;itype++)
                fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                        ,atom_types->atom_names[itype]);
            
            
            
            for(int i=0;i<tot_natms;i++)
            {
                for(int j=0;j<no_vecs-1;j++)
                    atoms->vectors[vec_list[j]]->print_dump(fp,sort[i]);
                fprintf(fp,"\n");
                if(usr_nabled)
                {
                    find_type(&c[sort[i]*atom_types->no_types],TRESH,clr_r);
                    for(int j=0;j<3;j++)
                        fprintf(fp_usr,"%f ",clr_r[j]);
                    fprintf(fp_usr,"%f\n",clr_r[3]);
                }
            }
            
            if(tot_natms)
                delete [] sort;
        }
        else
        {
            int icomp=0;
            int xx_dim=x_dim-3;
            if(xx_dim==1)
            {
                fprintf(fp,"auxiliary[%d] = alpha [reduced unit]\n",icomp++);
            }
            else if(xx_dim>1)
            {
                for(int idim=0;idim<xx_dim;idim++)
                    fprintf(fp,"auxiliary[%d] = alpha_%d [reduced unit]\n",icomp++,idim);
            }
            
            for(int i=1;i<no_vecs;i++)
            {
                if(atoms->vectors[vec_list[i]]->dim>1)
                {
                    for(int idim=0;idim<atoms->vectors[vec_list[i]]->dim;idim++)
                        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n"
                                ,icomp++,atoms->vectors[vec_list[i]]->name,idim);
                }
                else
                {
                    fprintf(fp,"auxiliary[%d] = %s [reduced unit]\n"
                            ,icomp++,atoms->vectors[vec_list[i]]->name);
                }
            }
            int tot_natms=atoms->tot_natms;
            
            for(int itype=0;itype<atom_types->no_types;itype++)
                fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                        ,atom_types->atom_names[itype]);
            type0* c;
            atoms->vectors[c_n]->ret_dump(c);
            
            for(int i=0;i<tot_natms;i++)
            {
                for(int j=0;j<no_vecs;j++)
                    atoms->vectors[vec_list[j]]->print_dump(fp,i);
                fprintf(fp,"\n");
                if(usr_nabled)
                {
                    find_type(&c[i*atom_types->no_types],TRESH,clr_r);
                    for(int j=0;j<3;j++)
                        fprintf(fp_usr,"%f ",clr_r[j]);
                    fprintf(fp_usr,"%f\n",clr_r[3]);
                }
            }
            

        }
        
        fclose(fp);
        if(usr_nabled)
            fclose(fp_usr);
    }
    
    for(int i=0;i<no_vecs;i++)
        atoms->vectors[vec_list[i]]->del_dump();
    
}
/*--------------------------------------------
 transform x 2 s
 --------------------------------------------*/
void Write_cfg::x2s(int no)
{
    
    type0* x;
    atoms->vectors[0]->ret_dump(x);
    type0** B=atoms->B;
    
    int icomp=0;
    for(int i=0;i<no;i++)
    {
        for(int j=0;j<dim;j++)
        {
            x[icomp+j]=x[icomp+j]*B[j][j];
            for(int k=j+1;k<dim;k++)
                x[icomp+j]+=B[k][j]*x[icomp+k];
            
        }
        
        for(int j=0;j<dim;j++)
        {
            while(x[icomp+j]<0.0)
                x[icomp+j]++;
            while(x[icomp+j]>=1.0)
                x[icomp+j]--;
            
        }
        icomp+=x_dim;
    }
    
}
/*--------------------------------------------
 find a type of atom based on tolerance
 and other shit it is for dmd mode
 --------------------------------------------*/
int Write_cfg::find_type(type0* c
,type0 tol,type0* clr_r)
{
    type0 c_m_eff,rat,c_v=1.0,c_max=0.0;
    int ityp=-1;
    
    for(int i=0;i<atom_types->no_types;i++)
    {
        if(c[i]>=0.0)
        {

            c_v-=c[i];
            if(ityp==-1)
            {
                ityp=i;
                c_max=c[i];
            }
            else
            {
                if(c[i]>c_max)
                {
                    ityp=i;
                    c_max=c[i];
                }
            }
        }
    }

    c_m_eff=1.0-0.5*(1.0+c_v-c_max);
    
    if(c_m_eff<tol)
    {
        rat=(tol-c_m_eff)/tol;
        
        rat=c_m_eff/tol;
        for(int i=0;i<3;i++)
            clr_r[i]=atom_types->clr_rad[ityp][i]*rat+1.0-rat;
        
        clr_r[3]=atom_types->clr_rad[ityp][3];
        return -1;
    }
    else
    {
        /*
        for(int i=0;i<3;i++)
            clr_r[i]=1.0;

        clr_r[3]=atom_types->clr_rad[ityp][3];
         */
        for(int i=0;i<4;i++)
            clr_r[i]=atom_types->clr_rad[ityp][i];
        return ityp;
    }
    
}
/*--------------------------------------------
 init()
 --------------------------------------------*/
void Write_cfg::init_indv()
{
    for(int i=0;i<no_vecs;i++)
        vec_list[i]=atoms->find(vec_name[i]);
    
    //chek wether two vectors are the same
    
    for(int i=0;i<no_vecs;i++)
        for(int j=i+1;j<no_vecs;j++)
            if(vec_list[i]==vec_list[j])
                error->abort("duplicate "
                             "atomic vectos in write "
                             "cfg cannot be the same");
    
    tot_dim=0;
    for(int i=0;i<no_vecs;i++)
        tot_dim+=atoms->vectors[vec_list[i]]->dim;
    
    if(mapp->mode==MD_mode)
        tot_dim-=atoms->vectors[type_n]->dim;
    if(sorting)
        tot_dim-=atoms->vectors[id_n]->dim;

}

