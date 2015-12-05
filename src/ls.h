#ifndef __MAPP__ls__
#define __MAPP__ls__
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
    /*
    const char *ls_msg[]=
    {
        "success",
        "not downhill direction",
        "gradient is 0",
        "reached minimum alpha",

        "tolerance reached",
        "maximum iteration",
        
        "bracketing success",
        "bracketing reached maximum alpha",
        "bracketing not downhill",        
    };*/
    
    template<class Func>
    class LineSearch
    {
    private:
    protected:
        
        type0 epsilon;
        type0 epsilon_3_4;
        
        type0 max_dx;
        type0 golden;
        type0 prev_val,h_norm;
        int bracket(type0,type0,type0&,type0&,type0&,type0&,type0&,type0&);
        Func* func;
    public:

        LineSearch(MAPP*);
        virtual ~LineSearch();
        virtual int line_min(type0&,type0&,int)=0;
        virtual void init(Func*);
        void test(type0,type0,type0);
    };

}
#include <limits>
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
template<class Func>
LineSearch<Func>::LineSearch(MAPP* mapp)
{
    max_dx=1.0;
    prev_val=0.0;
    golden=0.5+0.5*sqrt(5.0);
    epsilon=std::numeric_limits<type0>::epsilon();
    epsilon_3_4=pow(epsilon,0.75);
    
    prev_val=-1.0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<class Func>
LineSearch<Func>::~LineSearch()
{
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
template<class Func>
void LineSearch<Func>::init(Func* func_)
{
    func=func_;
}
/*--------------------------------------------
 bracketing routine
 --------------------------------------------*/
template<class Func>
int LineSearch<Func>::bracket(type0 dfa,type0 max_a,type0& a,
type0& b,type0& c,type0& fa,type0& fb,type0& fc)
{
    type0 u,fu,r,q,ulim;
   
    q=1.0e-14/(max_a*fabs(dfa));
    if(0.01<q && q<1.0)
    {
        b=q*max_a;
    }
    else if(q<=1.0)
    {
        b=max_a;
    }
    else
        b=0.01*max_a;
    
    fb=func->F(b);
    
    if(fb>=fa)
    {
        c=b;
        fc=fb;
        int iter=20;
        type0 r=(fa-fc)/(c*dfa);
        b=c*0.5/(1.0+r);
        while(fb>=fa && iter && b>epsilon)
        {
            fb=func->F(b);
            if(fb<fa)
                continue;
            c=b;
            fc=fb;
            iter--;
            r=(fa-fc)/(c*dfa);
            b=c*0.5/(1.0+r);
        }
        
        if(fb>=fa)
            return B_F_DOWNHILL;
        
        return B_S;
    }

    fc=fb;
    while (fb>=fc)
    {
        
        c=b+golden*(b-a);
        if(c>=max_a)
        {
            c=max_a;
            fc=func->F(c);
            return B_S;
        }
        
        fc=func->F(c);
        if(fc>fb)
            continue;
        
        ulim=MIN(b+(golden+2.0)*(b-a),max_a);
        ulim=max_a;
        
        r=(b-a)*(fb-fc);
        q=(b-c)*(fb-fa);
        
        u=b-((b-c)*q-(b-a)*r)/(2.0*(q-r));
        if(b<u && u<c)
        {
            fu=func->F(u);
            if(fu<fc)
            {
                a=b;
                b=u;
                fa=fb;
                fb=fu;
                return B_S;
            }
            else if (fu>fb)
            {
                c=u;
                fc=fu;
                return B_S;
            }
            
            a=b;
            b=c;
            
            fa=fb;
            fb=fc;
        }
        else if(u>c)
        {
            
            u=MIN(u,ulim);
            fu=func->F(u);
            a=b;
            b=c;
            c=u;
            
            fa=fb;
            fb=fc;
            fc=fu;
            
            if(fu>fc)
            {
                return B_S;
            }
        }
        else
        {
            a=b;
            b=c;
            
            fa=fb;
            fb=fc;
        }
    }
    
    return B_S;
}
/*--------------------------------------------
 reset to initial position
 --------------------------------------------*/
template<class Func>
void LineSearch<Func>::test(type0 fa,type0 dfa,type0 max_a)
{
    int no=100;
    type0 frac=1.0e-2*max_a;
    type0 fu,u=0.0;
    fa=func->F(0);
    
    printf("dfa %e\n",dfa);
    printf("u fu f_x u*dfa\n");

    for(int i=0;i<no;i++)
    {
        fu=func->F(u);
        printf("%22.20lf %22.20lf %22.20lf \n",u,fu-fa,u*dfa);
        u+=frac;
    }
    
}



#endif 
