#ifndef __MAPP__line_search__
#define __MAPP__line_search__

#include <iostream>
#include "atoms.h"
namespace MAPP_NS {
    enum {LS_S,LS_F_DOWNHILL,LS_F_GRAD0,LS_F_ALPHAMIN,MIN_F_MAX_ITER,MIN_F_TOLERANCE};
    
    class LineSearch : protected InitPtrs
    {
    private:
    protected:
        int x_n;
        int x_prev_n;
        int f_n;
        int dim;
        int x_dim;
        int dof_n;
        
        type0 inner_f_h();
        type0 inner_f_h_s();
        void normalize_h();
        void normalize_h_s();
        type0 energy(type0);
        VecLst* vecs_comm;
        type0 d_max;
        type0 s_max;
        type0** N;
        type0** M;
        
        int** H_dof;
        void ccc(type0**,type0**,type0**,type0);
    public:
        type0** H_prev;
        type0** B_prev;
        type0** h_H;
        type0** f_H;
        int chng_box;
        ThermoDynamics* thermo;
        LineSearch(MAPP *,VecLst*);
        virtual ~LineSearch()=0;
        virtual int line_min(type0&,type0&)=0;
        int h_n;
    };
    
    class LineSearch_BackTrack : public LineSearch
    {
    private:
    protected:
        type0 c,rho,alpha_max,alpha_min;
    public:
        LineSearch_BackTrack(MAPP *,VecLst*);
        ~LineSearch_BackTrack();
        int line_min(type0&,type0&);
        int line_min_s(type0&,type0&);
    };

}

#endif /* defined(__MAPP__line_search__) */
