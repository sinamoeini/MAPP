#ifdef DMD_Style
    DMDStyle(DMD_bdf,bdf)
#else
#ifndef __MAPP__dmd_bdf__
#define __MAPP__dmd_bdf__
#include "dmd.h"
namespace MAPP_NS
{
    class DMD_bdf:public DMDImplicit
    {
    private:
        type0* t;
        type0* l;
        
        void update_A_bar();

        type0** A;
        type0** A_bar;
        
        type0* dy;
        type0* z;
        
        
        type0 lo_err_fac[2];
        type0 hi_err_fac[2];
        
        bool lo_avail,hi_avail;
        type0 eta[3];
        

        
        void eta_calc();
        
        Vec<type0>** vecs_1;

        void update_z();
        
        
        class DMD_bdf_sub
        {
        private:
        protected:
            type0*& t;
            type0& dt;
            
            type0& beta;
            type0& beta_inv;
            
            type0& err_fac;
            type0(& lo_err_fac)[2];
            type0(& hi_err_fac)[2];
            
            int& q;
            int& dq;
            int& const_q;
            int& q_max;
            
            type0*& l;
            type0**& A_bar;
            
            
            
            void l_calc(int);

            
        public:
            DMD_bdf_sub(DMD_bdf* bdf):
            t(bdf->t),
            dt(bdf->dt),
            
            beta(bdf->beta),
            beta_inv(bdf->beta_inv),
            
            err_fac(bdf->err_fac),
            lo_err_fac(bdf->lo_err_fac),
            hi_err_fac(bdf->hi_err_fac),
            
            q(bdf->q),
            dq(bdf->dq),
            const_q(bdf->const_q),
            q_max(bdf->q_max),
            
            l(bdf->l),
            A_bar(bdf->A_bar)
            {}
            virtual ~DMD_bdf_sub(){}
            
            virtual void beta_calc()=0;
            virtual void prep_A_l()=0;
            virtual void err_fac_calc()=0;
        };
        
        class DMD_bdf_f_flc:public DMD_bdf_sub
        {
        private:
        protected:
        public:
            DMD_bdf_f_flc(DMD_bdf* bdf):DMD_bdf_sub(bdf)
            {}
            ~DMD_bdf_f_flc(){}
            void beta_calc();
            void prep_A_l();
            void err_fac_calc();
        };
        
        class DMD_bdf_f_vc:public DMD_bdf_sub
        {
        private:
        protected:
        public:
            DMD_bdf_f_vc(DMD_bdf* bdf):DMD_bdf_sub(bdf)
            {}
            ~DMD_bdf_f_vc(){}
            void beta_calc();
            void prep_A_l();
            void err_fac_calc();
        };
        
        class DMD_bdf_y_flc:public DMD_bdf_sub
        {
        private:
        protected:
        public:
            DMD_bdf_y_flc(DMD_bdf* bdf):DMD_bdf_sub(bdf)
            {}
            ~DMD_bdf_y_flc(){}
            void beta_calc();
            void prep_A_l();
            void err_fac_calc();
        };
        
        class DMD_bdf_y_vc:public DMD_bdf_sub
        {
        private:
        protected:
        public:
            DMD_bdf_y_vc(DMD_bdf* bdf):DMD_bdf_sub(bdf)
            {}
            ~DMD_bdf_y_vc(){}
            void beta_calc();
            void prep_A_l();
            void err_fac_calc();
        };
        
        DMD_bdf_sub* bdf_eng;
    protected:
        // stuff for book keeping
        
        void allocate();
        void deallocate();
        
        //functions to be accessed by run()
        void restart();
        void update_for_next();
        bool interpolate();
        void interpolate_fail();
        void ord_dt(type0&);
        void err_calc();
        void start();
        void err_fac_calc();
    public:
        DMD_bdf(MAPP*,int,char**);
        ~DMD_bdf();
    };
}

#endif 
#endif