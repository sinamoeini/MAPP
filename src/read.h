#ifndef __MAPP__read__
#define __MAPP__read__
#include <mpi.h>
namespace MAPP_NS
{
    class Read
    {
    private:
    protected:
        class AtomTypes*& atom_types;
        MPI_Comm& world;
    public:
        Read();
        virtual ~Read();

    };
}
#endif

