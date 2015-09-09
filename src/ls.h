#ifndef __MAPP__ls__
#define __MAPP__ls__
#include "init.h"
#include "atoms.h"
#include "xmath.h"
namespace MAPP_NS
{
    enum
    {
        LS_S,
        LS_F_DOWNHILL,
        LS_F_GRAD0,
        LS_MIN_ALPHA,
        
        MIN_S_TOLERANCE,
        MIN_F_MAX_ITER,
        
        B_S,
        B_F_MAX_ALPHA,
        B_F_DOWNHILL
    };
       
    
    class LineSearch : protected InitPtrs
    {
    private:
    protected:
        
        type0 epsilon;
        type0 epsilon_3_4;
        
        type0 max_dx;
        type0 golden;
        
        type0 energy(type0);
        type0 energy(type0,type0&);

        
        int bracket(type0,type0,type0&,type0&,type0&,type0&,type0&,type0&);
        
        void init_manip(type0&,type0&,type0&);
        void reset();
        
        XMath* xmath;
    public:
        type0** H_prev;
        type0** B_prev;
        type0** h_H;
        type0** f_H;
        type0* nrgy_strss;
        
        VecLst* vecs_comm;
        
        // for x mode
        int x_prev_n;
        int h_n;
        int f_n;
        int dim;
        int x_dim;
        int chng_box;
        int affine;
        
        type0 prev_val;
        
        LineSearch(MAPP*);
        virtual ~LineSearch();
        virtual int line_min(type0&,type0&,int)=0;
        
        //void test(type0,type0,type0);
    };

}


#endif 
