#include "thermo_dynamics.h"
#include "error.h"
#include "memory.h"
#include "atoms.h"
using namespace MAPP_NS;
/*--------------------------------------------
 init the quantity
 --------------------------------------------*/
void ThermoQuantity::init(const char* name)
{
    hdr_name_lngth=1+static_cast<int>(strlen(name));
    value=0.0;
    hdr_name= new char[hdr_name_lngth];
    for(int i=0;i<hdr_name_lngth;i++)
        hdr_name[i]=name[i];
    mod_hdr_name_lngth=0;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ThermoQuantity::ThermoQuantity()
{
    hdr_name_lngth=0;
    mod_hdr_name_lngth=0;
}
/*--------------------------------------------
 copy constructor
 --------------------------------------------*/
ThermoQuantity::ThermoQuantity(ThermoQuantity& other)
{
    hdr_name_lngth=other.hdr_name_lngth;
    hdr_name= new char[hdr_name_lngth];
    memcpy(hdr_name,other.hdr_name,hdr_name_lngth*sizeof(char));
    mod_hdr_name_lngth=other.mod_hdr_name_lngth;
    mod_hdr_name=new char[mod_hdr_name_lngth];
    memcpy(mod_hdr_name,other.mod_hdr_name,mod_hdr_name_lngth*sizeof(char));
    value=other.value;
}
/*--------------------------------------------
 move constructor
 --------------------------------------------*/
ThermoQuantity::ThermoQuantity(ThermoQuantity&& other)
{
    hdr_name_lngth=other.hdr_name_lngth;
    other.hdr_name_lngth=0;
    hdr_name=other.hdr_name;
    other.hdr_name=NULL;
    mod_hdr_name_lngth=other.mod_hdr_name_lngth;
    other.mod_hdr_name_lngth=0;
    mod_hdr_name=other.mod_hdr_name;
    other.mod_hdr_name=NULL;
    value=other.value;
}
/*--------------------------------------------
 copy operator
 --------------------------------------------*/
ThermoQuantity& ThermoQuantity::operator=(ThermoQuantity& other)
{
    this->hdr_name_lngth=other.hdr_name_lngth;
    this->hdr_name= new char[hdr_name_lngth];
    memcpy(this->hdr_name,other.hdr_name,this->hdr_name_lngth*sizeof(char));
    this->mod_hdr_name_lngth=other.mod_hdr_name_lngth;
    this->mod_hdr_name=new char[this->mod_hdr_name_lngth];
    memcpy(this->mod_hdr_name,other.mod_hdr_name,this->mod_hdr_name_lngth*sizeof(char));
    this->value=other.value;
    return *this;
}
/*--------------------------------------------
 move operator
 --------------------------------------------*/
ThermoQuantity& ThermoQuantity::operator=(ThermoQuantity&& other)
{
    this->hdr_name_lngth=other.hdr_name_lngth;
    other.hdr_name_lngth=0;
    this->hdr_name=other.hdr_name;
    other.hdr_name=NULL;
    this->mod_hdr_name_lngth=other.mod_hdr_name_lngth;
    other.mod_hdr_name_lngth=0;
    this->mod_hdr_name=other.mod_hdr_name;
    other.mod_hdr_name=NULL;
    this->value=other.value;
    return *this;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ThermoQuantity::~ThermoQuantity()
{
    if(hdr_name_lngth)
        delete [] hdr_name;
    if(mod_hdr_name_lngth)
        delete [] mod_hdr_name;
}
/*--------------------------------------------
 modify for header
 --------------------------------------------*/
void ThermoQuantity::mod(int lngth)
{
    if(mod_hdr_name_lngth)
        delete [] mod_hdr_name;
    mod_hdr_name_lngth=lngth+1;
    
    mod_hdr_name=new char[mod_hdr_name_lngth];
    
    int no=mod_hdr_name_lngth-hdr_name_lngth;
    int lft_n;
    lft_n=no/2;
    memset(mod_hdr_name,' ',mod_hdr_name_lngth-1);
    mod_hdr_name[mod_hdr_name_lngth-1]='\0';
    for(int i=0;i<hdr_name_lngth-1;i++)
    {
        mod_hdr_name[i+lft_n]=hdr_name[i];
    }
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ThermoDynamics::ThermoDynamics(int nargs,char** args)
{
    step_name_lngth=15;
    mod_lngth=precision+9;    
    
    if(__dim__!=3)
        Error::abort("the thermodynamics "
        "works only with box dimension 3");

    no_quantities=nargs;
    quantities=new ThermoQuantity[no_quantities];
    
    
    for(int i=0;i<no_quantities;i++)
        quantities[i].init(args[i]);
    
    for(int i=0;i<no_quantities;i++)
    {
        quantities[i].mod(mod_lngth);
    }


    CREATE1D(step_name,step_name_lngth+1);
    int no=step_name_lngth-strlen("Step")+1;
    int tmp_lngth=static_cast<int>(strlen("Step"));
    int left_no=no/2;
    memset(step_name,' ',step_name_lngth);
    for(int i=0;i<tmp_lngth;i++)
        step_name[i+left_no]="Step"[i];
    step_name[step_name_lngth]='\0';
    
    
    
    int tmp0=7
    +static_cast<int>(log(static_cast<type0>(mod_lngth-1))/log(10))
    +static_cast<int>(log(static_cast<type0>(precision))/log(10));

    int tmp1=7
    +static_cast<int>(log(static_cast<type0>(step_name_lngth-2))/log(10));
    
    CREATE1D(qform,tmp0);
    CREATE1D(sform,tmp1);
    sprintf(qform,"%%%d.%de ",mod_lngth-1,precision);
    sprintf(sform," %%0%dd ",step_name_lngth-2);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ThermoDynamics::~ThermoDynamics()
{

    delete [] sform;
    delete [] qform;
    delete [] step_name;
    delete [] quantities;
}
/*--------------------------------------------
 add quantites
 --------------------------------------------*/
int ThermoDynamics::add(const char* q)
{
    ThermoQuantity* quantities_=new ThermoQuantity[no_quantities+1];
    for(int i=0;i<no_quantities;i++)
        quantities_[i]=std::move(quantities[i]);
    if(no_quantities)
        delete [] quantities;
    quantities=quantities_;
    quantities[no_quantities].init(q);
    quantities[no_quantities].mod(mod_lngth);
    no_quantities++;
    return no_quantities-1;
}
/*--------------------------------------------
 print output
 --------------------------------------------*/
void ThermoDynamics::hdr_print()
{
    if(atoms->my_p==0)
    {
        fprintf(output,"\n");
        fprintf(output," ");
        for(int i=0;i<step_name_lngth;i++)
            fprintf(output,"-");
        
        for(int i=0;i<no_quantities;i++)
        {
            fprintf(output,"-");
            for(int j=0;j<mod_lngth;j++)
                fprintf(output,"-");
        }
        fprintf(output,"\n");
        
        fprintf(output,"|");
        fprintf(output,"%s",step_name);
        for(int i=0;i<no_quantities;i++)
        {
            fprintf(output,"|");
            fprintf(output,"%s",quantities[i].mod_hdr_name);
        }
        fprintf(output,"|\n");
        
        
        
        fprintf(output,"|");
        for(int i=0;i<step_name_lngth;i++)
            fprintf(output,"-");
        
        for(int i=0;i<no_quantities;i++)
        {
            fprintf(output,"|");
            for(int j=0;j<mod_lngth;j++)
                fprintf(output,"-");
        }
        fprintf(output,"|\n");
    }
}
/*--------------------------------------------
 print output
 --------------------------------------------*/
void ThermoDynamics::val_print()
{

    print_step=step_no+step_tally;
    if(atoms->my_p==0)
    {
        fprintf(output,"|");
        
        fprintf(output,sform,step_no);
        
        for(int i=0;i<no_quantities;i++)
        {
            fprintf(output,"|");
            fprintf(output,qform,quantities[i].value);
        }
        fprintf(output,"|\n");
    }
}
/*--------------------------------------------
 print output
 --------------------------------------------*/
void ThermoDynamics::thermo_print()
{
    if(print_step!=step_no)
        return;
    print_step=step_no+step_tally;
    val_print();
}
/*--------------------------------------------
 print output
 --------------------------------------------*/
void ThermoDynamics::tail_print()
{
    if(atoms->my_p==0)
    {
        fprintf(output," ");
        for(int i=0;i<step_name_lngth;i++)
            fprintf(output,"-");
        
        for(int i=0;i<no_quantities;i++)
        {
            fprintf(output,"-");
            for(int j=0;j<mod_lngth;j++)
                fprintf(output,"-");
        }

        fprintf(output,"\n");
    }
    
}
/*--------------------------------------------
 initiated before a run
 --------------------------------------------*/
void ThermoDynamics::init()
{
    init_step=step_no;
    hdr_print();
    val_print();
    print_step=step_no+step_tally;
}
/*--------------------------------------------
 finish after a run
 --------------------------------------------*/
void ThermoDynamics::fin()
{
    if(step_no!=print_step-step_tally)
        val_print();
    tail_print();
}
/*--------------------------------------------
 initiated before a run
 --------------------------------------------*/
int ThermoDynamics::test_prev_step()
{
    if(step_no+1==print_step)
        return 1;
    else
        return 0;
}
/*--------------------------------------------
 update stress values
 --------------------------------------------*/
void ThermoDynamics::update(int qstrt,int qlngth
,type0* values)
{
    for(int i=0;i<qlngth;i++)
        quantities[i+qstrt].value=values[i];

}
/*--------------------------------------------
 update PE
 --------------------------------------------*/
void ThermoDynamics::update(int qindx,type0 value)
{
    quantities[qindx].value=value;
}

