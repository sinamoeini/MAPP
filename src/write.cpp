#include "write.h"
#include "timer.h"
#include "MAPP.h"
#include "global.h"
#include "comm.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Write::Write():
atom_types(mapp->atom_types),
world(comm->world)
{
    last_write_step=-1;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Write::~Write()
{
    
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void Write::init()
{
    if(last_write_step!=step_no)
    {
        mapp->timer->start(WRITE_TIME_mode);
        write_file(step_no);
        mapp->timer->stop(WRITE_TIME_mode);
        last_write_step=step_no;
    }
    write_step=step_no+write_step_tally;
}
/*--------------------------------------------
 write the file
 --------------------------------------------*/
void Write::write()
{
    if(write_step!=step_no)
        return;
    mapp->timer->start(WRITE_TIME_mode);
    write_file(step_no);
    mapp->timer->stop(WRITE_TIME_mode);
    last_write_step=step_no;
    write_step=step_no+write_step_tally;
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void Write::fin()
{

    if(last_write_step!=step_no)
    {
        mapp->timer->start(WRITE_TIME_mode);
        write_file(step_no);
        mapp->timer->stop(WRITE_TIME_mode);
    }
    last_write_step=step_no;
}
