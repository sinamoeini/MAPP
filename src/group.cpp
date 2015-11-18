#include "group.h"
#include "atoms.h"
#include "xmath.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Group::Group(MAPP* mapp):InitPtrs(mapp)
{
    group_list_size=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Group::~Group()
{

}
/*--------------------------------------------
 read line and broadcast
 --------------------------------------------*/
void Group::update_ids()
{

}

