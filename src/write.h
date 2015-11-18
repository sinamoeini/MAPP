#ifndef __MAPP__write__
#define __MAPP__write__

#include "init.h"

namespace MAPP_NS {
    class Write: protected InitPtrs
    {
    private:
    protected:
        int write_step_tally;
        int write_step;
        int last_write_step;
    public:
        Write(MAPP *);
        virtual ~Write();
        virtual void write_file(int)=0;
        virtual void init();
        void write();
        virtual void fin();
    };
    
}
#endif
