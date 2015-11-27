#include <stdlib.h>
#include "atom_types.h"
#include "write_cfg.h"
#include "error.h"
#include "memory.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Write_cfg::Write_cfg(MAPP* mapp,int nargs
,char** args) : Write(mapp)
{
    thresh=1.0;
    strech=1.0;
    usr_nabled=0;
    sorting=0;
    ndump_vecs=0;
    ndump_vecs_names=0;
    dim=3;
    
    if(mapp->mode==MD_mode)
    {
        GROW(dump_vecs,ndump_vecs,ndump_vecs+1);
        dump_vecs[ndump_vecs]=mapp->type;
        ndump_vecs++;
        GROW(dump_vecs,ndump_vecs,ndump_vecs+1);
        dump_vecs[ndump_vecs]=mapp->x;
        ndump_vecs++;
    }
    else if(mapp->mode==DMD_mode)
    {
        GROW(dump_vecs,ndump_vecs,ndump_vecs+1);
        dump_vecs[ndump_vecs]=mapp->x;
        ndump_vecs++;
        GROW(dump_vecs,ndump_vecs,ndump_vecs+1);
        dump_vecs[ndump_vecs]=mapp->c;
        ndump_vecs++;
    }
    

    
    write_step_tally=atoi(args[2]);
    int lngth=static_cast<int>(strlen(args[3]))+1;
    CREATE1D(file_name,lngth);
    memcpy(file_name,args[3],lngth*sizeof(char));

    int iarg=4;
    while(iarg<nargs)
    {
        if(strcmp(args[iarg],"sort")==0)
        {
            sorting=1;
            iarg++;
        }
        else if(strcmp(args[iarg],"usr")==0)
        {
            usr_nabled=1;
            iarg++;
        }
        else if(strcmp(args[iarg],"thresh")==0)
        {
            iarg++;
            if(iarg==nargs)
                error->abort("thresh in write cfg should be followed by it's value");
            thresh=atof(args[iarg]);
            if(thresh<0.0 || thresh>1.0)
                error->abort("thresh in write cfg should be between 0.0 & 1.0");
            iarg++;
        }
        else if(strcmp(args[iarg],"stretch")==0)
        {
            iarg++;
            if(iarg==nargs)
                error->abort("stretch in write cfg should be followed by it's value");
            strech=atof(args[iarg]);
            if(strech<0.0)
                error->abort("stretch in write cfg should be greater than 0.0");
            iarg++;
        }
        else
        {
            add_to_dump_vec_nams(args[iarg]);
            iarg++;
        }
    }
    
    if(sorting)
    {
        vec** dump_vecs_;
        CREATE1D(dump_vecs_,ndump_vecs+1);
        for(int ivec=0;ivec<ndump_vecs;ivec++)
            dump_vecs_[ivec+1]=dump_vecs[ivec];
        dump_vecs_[0]=mapp->id;
        if(ndump_vecs)
            delete [] dump_vecs;
        dump_vecs=dump_vecs_;
        ndump_vecs++;

    }


    if(usr_nabled || mapp->mode==DMD_mode)
        CREATE1D(clr_r,4);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Write_cfg::~Write_cfg()
{
    if(ndump_vecs_names)
    {
        for(int ivec=0;ivec<ndump_vecs_names;ivec++)
            delete [] dump_vecs_names[ivec];
        delete [] dump_vecs_names;
    }
    
    if(usr_nabled || mapp->mode==DMD_mode)
        delete [] clr_r;
    
    delete [] file_name;
    
    if(ndump_vecs)
        delete [] dump_vecs;
}
/*--------------------------------------------
 write file
 --------------------------------------------*/
void Write_cfg::write_file(int stp)
{
    if(mapp->mode==MD_mode)
        write_file_md(stp);
    else if(mapp->mode==DMD_mode)
        write_file_dmd(stp);
}
/*--------------------------------------------
 write file
 --------------------------------------------*/
void Write_cfg::write_file_md(int stp)
{
    for(int ivec=0;ivec<ndump_vecs;ivec++)
        dump_vecs[ivec]->gather_dump();
        
    FILE* fp=NULL;
    FILE* fp_usr=NULL;
    
    if(atoms->my_p==0)
    {
        x2s(atoms->tot_natms);
        char* filename;
        CREATE1D(filename,MAXCHAR);
        sprintf (filename, "%s.%010d.cfg",file_name,stp);
        fp=fopen(filename,"w");
        if(fp==NULL)
            error->abort("Cannot create file %s for write cfg",filename);
        delete [] filename;
        if(usr_nabled)
        {
            CREATE1D(filename,MAXCHAR);
            sprintf (filename, "%s.%010d.usr",file_name,stp);
            fp_usr=fopen(filename,"w");
            if(fp_usr==NULL)
                error->abort("Cannot create file %s for write cfg",filename);
            delete [] filename;
        }

        
        // write the header
        fprintf(fp,"Number of particles = %d\n",atoms->tot_natms);
        fprintf(fp,"A = %lf Angstrom (basic length-scale)\n",strech);
        
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                fprintf(fp,"H0(%d,%d) = %lf A\n",i+1,j+1,atoms->H[i][j]);
        
        fprintf(fp,".NO_VELOCITY.\n");
        
        fprintf(fp,"entry_count = %d\n",tot_dim);
        
        
        // write the body
        if(sorting)
        {
            int icomp=0;
            
            for(int ivec=3;ivec<ndump_vecs;ivec++)
            {
                if(dump_vecs[ivec]->dim>1)
                {
                    for(int idim=0;idim<dump_vecs[ivec]->dim;idim++)
                        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n"
                            ,icomp++,dump_vecs[ivec]->name,idim);
                }
                else
                {
                    fprintf(fp,"auxiliary[%d] = %s [reduced unit]\n"
                            ,icomp++,dump_vecs[ivec]->name);
                }
                
            }
            
            int tot_natms=atoms->tot_natms;
            int* id=mapp->id->begin_dump();
            int* sort;
            CREATE1D(sort,tot_natms);
            for(int i=0;i<tot_natms;i++)
                sort[id[i]]=i;
            
            md_type* type=mapp->type->begin_dump();
            md_type itype=-1;
            int iatm;
            for(int i=0;i<tot_natms;i++)
            {
                iatm=sort[i];
                if(type[iatm]!=itype)
                {
                    itype=type[iatm];
                    fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                    ,atom_types->atom_names[itype]);
                }
                
                for(int ivec=2;ivec<ndump_vecs;ivec++)
                    dump_vecs[ivec]->print_dump(fp,iatm);
                fprintf(fp,"\n");
                if(usr_nabled)
                {
                    for(int j=0;j<4;j++)
                        fprintf(fp_usr,"%lf ",atom_types->clr_rad[itype][j]);
                    fprintf(fp_usr,"\n");
                }
            }
            
            
            if(tot_natms)
                delete [] sort;
        }
        else
        {
            int icomp=0;
            
            for(int ivec=2;ivec<ndump_vecs;ivec++)
            {
                if(dump_vecs[ivec]->dim>1)
                {
                    for(int idim=0;idim<dump_vecs[ivec]->dim;idim++)
                        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n"
                                ,icomp++,dump_vecs[ivec]->name,idim);
                }
                else
                {
                    fprintf(fp,"auxiliary[%d] = %s [reduced unit]\n"
                            ,icomp++,dump_vecs[ivec]->name);
                }
                
            }
            int tot_natms=atoms->tot_natms;
            
            md_type* type=mapp->type->begin_dump();
            md_type itype=-1;
            for(int iatm=0;iatm<tot_natms;iatm++)
            {
                if(type[iatm]!=itype)
                {
                    itype=type[iatm];
                    fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                            ,atom_types->atom_names[itype]);
                }
                
                for(int ivec=1;ivec<ndump_vecs;ivec++)
                    dump_vecs[ivec]->print_dump(fp,iatm);
                fprintf(fp,"\n");
                if(usr_nabled)
                {
                    for(int j=0;j<4;j++)
                        fprintf(fp_usr,"%lf ",atom_types->clr_rad[itype][j]);
                    fprintf(fp_usr,"\n");
                }
            }
            
        }
        
        fclose(fp);
        if(usr_nabled)
            fclose(fp_usr);
    }
    
    for(int ivec=0;ivec<ndump_vecs;ivec++)
        dump_vecs[ivec]->del_dump();
    
}
/*--------------------------------------------
 write file
 --------------------------------------------*/
void Write_cfg::write_file_dmd(int stp)
{
    
    for(int ivec=0;ivec<ndump_vecs;ivec++)
        dump_vecs[ivec]->gather_dump(mapp->ctype);

    
    
    FILE* fp=NULL;
    FILE* fp_usr=NULL;
    if(atoms->my_p==0)
    {
        x2s(atoms->tot_natms);

        char* filename;
        CREATE1D(filename,MAXCHAR);
        sprintf (filename, "%s.%010d.cfg",file_name,stp);
        fp=fopen(filename,"w");
        if(fp==NULL)
            error->abort("Cannot create file %s for write cfg",filename);
        delete [] filename;
        
        if(usr_nabled)
        {
            CREATE1D(filename,MAXCHAR);
            sprintf (filename, "%s.%010d.usr",file_name,stp);
            fp_usr=fopen(filename,"w");
            if(fp_usr==NULL)
                error->abort("Cannot create file %s for write cfg",filename);
            delete [] filename;
        }
        
        
        // write the header
        fprintf(fp,"Number of particles = %d\n",atoms->tot_natms);
        fprintf(fp,"A = %lf Angstrom (basic length-scale)\n",strech);
        
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                fprintf(fp,"H0(%d,%d) = %lf A\n",i+1,j+1,atoms->H[i][j]);
        
        fprintf(fp,".NO_VELOCITY.\n");
        
        fprintf(fp,"entry_count = %d\n",tot_dim);
        
        
        int icomp=0;
        int xx_dim=mapp->x->orig_dim-3;
        if(xx_dim==1)
        {
            fprintf(fp,"auxiliary[%d] = alpha [reduced unit]\n",icomp++);
            fprintf(fp,"auxiliary[%d] = c [reduced unit]\n",icomp++);
        }
        else if(xx_dim>1)
        {
            for(int idim=0;idim<xx_dim;idim++)
                fprintf(fp,"auxiliary[%d] = alpha_%d [reduced unit]\n",icomp++,idim);
            for(int idim=0;idim<xx_dim;idim++)
                fprintf(fp,"auxiliary[%d] = c_%d [reduced unit]\n",icomp++,idim);
        }
        // write the body

        if(sorting)
        {
            for(int ivec=3;ivec<ndump_vecs;ivec++)
            {
                if(dump_vecs[ivec]->orig_dim>1)
                {
                    for(int idim=0;idim<dump_vecs[ivec]->orig_dim;idim++)
                        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n"
                            ,icomp++,dump_vecs[ivec]->name,idim);
                }
                else
                {
                    fprintf(fp,"auxiliary[%d] = %s [reduced unit]\n"
                            ,icomp++,dump_vecs[ivec]->name);
                }
                
            }
            
            for(int itype=0;itype<atom_types->no_types;itype++)
                fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                        ,atom_types->atom_names[itype]);
            
            int tot_natms=atoms->tot_natms;
            type0* c=mapp->c->begin_dump();
            int c_dim=mapp->c->orig_dim;
            int* id=mapp->id->begin_dump();
            int* sort;
            
            CREATE1D(sort,tot_natms);
            for(int i=0;i<tot_natms;i++)
                sort[id[i]]=i;
            
            int itype=atom_types->no_types-1;
            int jtype;
            int iatm;
            for(int i=0;i<tot_natms;i++)
            {
                iatm=sort[i];
                jtype=find_type(&c[iatm*c_dim],clr_r);
                if(jtype!=itype)
                {
                    itype=jtype;
                    fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                    ,atom_types->atom_names[itype]);
                }
                
                for(int ivec=1;ivec<ndump_vecs;ivec++)
                    dump_vecs[ivec]->print_dump(fp,iatm);
                fprintf(fp,"\n");
                if(usr_nabled)
                {
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
            for(int ivec=2;ivec<ndump_vecs;ivec++)
            {
                if(dump_vecs[ivec]->orig_dim>1)
                {
                    for(int idim=0;idim<dump_vecs[ivec]->orig_dim;idim++)
                        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n"
                            ,icomp++,dump_vecs[ivec]->name,idim);
                }
                else
                {
                    fprintf(fp,"auxiliary[%d] = %s [reduced unit]\n"
                            ,icomp++,dump_vecs[ivec]->name);
                }
                
            }
            
            for(int itype=0;itype<atom_types->no_types;itype++)
                fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                        ,atom_types->atom_names[itype]);
            
            int tot_natms=atoms->tot_natms;
            type0* c=mapp->c->begin_dump();
            int c_dim=mapp->c->orig_dim;
            
            
            int itype=atom_types->no_types-1;
            int jtype;
            for(int iatm=0;iatm<tot_natms;iatm++)
            {
                jtype=find_type(&c[iatm*c_dim],clr_r);
                if(jtype!=itype)
                {
                    itype=jtype;
                    fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                            ,atom_types->atom_names[itype]);
                }
                
                for(int ivec=0;ivec<ndump_vecs;ivec++)
                    dump_vecs[ivec]->print_dump(fp,iatm);
                fprintf(fp,"\n");
                if(usr_nabled)
                {
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
    
    for(int ivec=0;ivec<ndump_vecs;ivec++)
        dump_vecs[ivec]->del_dump();
    
}
/*--------------------------------------------
 transform x 2 s
 --------------------------------------------*/
void Write_cfg::x2s(int no)
{
    
    type0* x=mapp->x->begin_dump();
    int x_dim=mapp->x->orig_dim;
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
,type0* clr_r)
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
    
    if(c_m_eff<thresh)
    {
        rat=(thresh-c_m_eff)/thresh;
        
        rat=c_m_eff/thresh;
        for(int i=0;i<3;i++)
            clr_r[i]=atom_types->clr_rad[ityp][i]*rat+1.0-rat;
        
        clr_r[3]=atom_types->clr_rad[ityp][3];
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
    }
    return ityp;
    
}
/*--------------------------------------------
 init()
 --------------------------------------------*/
void Write_cfg::init()
{
    for(int ivec=0;ivec<ndump_vecs_names;ivec++)
    {
        vec* v=atoms->find_vec(dump_vecs_names[ivec]);
        GROW(dump_vecs,ndump_vecs,ndump_vecs+1);
        dump_vecs[ndump_vecs]=v;
        ndump_vecs++;
    }
    
    tot_dim=0;
    for(int ivec=0;ivec<ndump_vecs;ivec++)
        tot_dim+=dump_vecs[ivec]->orig_dim;

    if(mapp->mode==MD_mode)
        tot_dim-=mapp->type->dim;
    
    if(sorting)
        tot_dim-=mapp->id->dim;
    
    Write::init();

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Write_cfg::fin()
{
    Write::fin();
    
    int ndump_vecs_=ndump_vecs-ndump_vecs_names;
    vec** dump_vecs_;
    CREATE1D(dump_vecs_,ndump_vecs_);
    memcpy(dump_vecs_,dump_vecs,ndump_vecs_*sizeof(vec*));
    if(ndump_vecs)
        delete [] dump_vecs;
    
    dump_vecs=dump_vecs_;
    ndump_vecs=ndump_vecs_;
}
/*--------------------------------------------
 init()
 --------------------------------------------*/
void Write_cfg::add_to_dump_vec_nams(const char* name)
{
    GROW(dump_vecs_names,ndump_vecs_names,ndump_vecs_names+1);
    int length=static_cast<int>(strlen(name)+1);
    CREATE1D(dump_vecs_names[ndump_vecs_names],length);
    memcpy(dump_vecs_names[ndump_vecs_names],name,length*sizeof(char));
    ndump_vecs_names++;
    
    
}


