#ifndef __MAPP__write__
#define __MAPP__write__
#include <mpi.h>
namespace MAPP_NS
{
    class Write
    {
    private:
    protected:
        int write_step_tally;
        int write_step;
        int last_write_step;
        class AtomTypes*& atom_types;
        MPI_Comm& world;
    public:
        Write();
        virtual ~Write();
        virtual void write_file(int)=0;
        virtual void init();
        void write();
        virtual void fin();
    };
}
#endif
