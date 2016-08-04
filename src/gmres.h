#ifndef __MAPP__gmres__
#define __MAPP__gmres__

#include "init.h"
/*--------------------------------------------
 
 --------------------------------------------*/

namespace MAPP_NS
{
    
    template<typename T0,class C0>
    class GMRES_m
    {
    private:
        const int m;
        const int n;
        C0*& kernel;
        T0** Q;
        T0** H;
        T0* b_hat;
        
        T0* c;
        T0* s;
        
        T0* y;
        
        T0 tol;
    protected:
    public:
        

        GMRES_m(int n_,int m_,T0 tol_,C0*& kernel_):
        m(m_),
        n(n_),
        kernel(kernel_),
        tol(tol_)
        {
            
            Q=new T0*[m+1];
            *Q=new T0[(m+1)*n];
            for(int i=1;i<m+1;i++)
                Q[i]=Q[i-1]+n;

            
            H=new T0*[m];
            *H=new T0[m*(m+1)/2];
            for(int i=1;i<m+1;i++)
                H[i]=H[i-1]+i;
            

            
            b_hat=new T0[m+1];
            c=new T0[m];
            s=new T0[m];
            y=new T0[m];
            
        }
        
        ~GMRES_m()
        {
            delete [] *Q;
            delete [] Q;
            delete [] *H;
            delete [] H;
            delete [] b_hat;
            delete [] c;
            delete [] s;
            delete [] y;
        }
        
        
        void solve(T0* b,T0* x)
        {
            T0 last_h=0.0,tmp_0,b_norm,err,err0;
            b_norm=norm(b);
            int max_iter=7;
            
            b_hat[0]=b_norm;
            for(int i=1;i<m+1;i++)
                b_hat[i]=0.0;
            
            for(int i=0;i<n;i++)
                Q[0][i]=b[i]/b_hat[0];
            
            
            err=err0=1.0;
            
            while (err>tol && max_iter)
            {
                for(int i=0;i<m && err>tol;i++)
                {
                    kernel->prod(Q[i],Q[i+1]);
                    
                    for(int j=0;j<i+1;j++)
                    {
                        H[i][j]=inner(Q[j],Q[i+1]);
                        add(-H[i][j],Q[j],Q[i+1]);
                    }
                    
                    last_h=norm(Q[i+1]);
                    mul(1.0/last_h,Q[i+1]);
                    
                    for(int j=0;j<i;j++)
                    {
                        tmp_0=c[j]*H[i][j]-s[j]*H[i][j+1];
                        H[i][j+1]=s[j]*H[i][j]+c[j]*H[i][j+1];
                        H[i][j]=tmp_0;
                    }
                    
                    tmp_0=sqrt(H[i][i]*H[i][i]+last_h*last_h);
                    
                    c[i]=H[i][i]/tmp_0;
                    s[i]=-last_h/tmp_0;
                    
                    H[i][i]=c[i]*H[i][i]-s[i]*last_h;
                    b_hat[i+1]=s[i]*b_hat[i];
                    b_hat[i]*=c[i];
                    
                    err=fabs(b_hat[i+1])/b_norm;
                    err0=err;
                    if(err<tol)
                    {
                        solve_y(i+1);
                        add_eq(y,Q,i+1,x);
                    }
                    
                }
                
                if(err<tol)
                    continue;
                
                solve_y(m);
                add_eq(y,Q,m,x);

                
                max_iter--;
                
                if (max_iter==0)
                    continue;
                

                tmp_0=fabs(b_hat[m]);
                
                if(b_hat[m]>0.0)
                    b_hat[m]=1.0;
                else
                    b_hat[m]=-1.0;

                for(int i=m-1;i>-1;i--)
                {
                    b_hat[i]=s[i]*b_hat[i+1];
                    b_hat[i+1]*=c[i];
                }
                
                add_(b_hat,Q,m+1,Q[0]);

                b_hat[0]=tmp_0;
                for(int i=1;i<m+1;i++)
                    b_hat[i]=0.0;
                
                
            }

        }
        
        
        
        T0 inner(T0* v0,T0* v1)
        {
            T0 ans=0.0;
            for(int i=0;i<n;i++)
                ans+=v0[i]*v1[i];
            return ans;
        }
        
        void mul(T0 s,T0* v)
        {
            for(int i=0;i<n;i++)
                v[i]*=s;
        }
        
        void add(T0 s,T0* v0,T0* v1)
        {
            for(int i=0;i<n;i++)
                v1[i]+=s*v0[i];
        }
        
        void add_(T0* coefs,T0** vecs,int nvecs,T0* v)
        {
            T0 tmp;
            for(int i=0;i<n;i++)
            {
                tmp=0.0;
                for(int ivec=0;ivec<nvecs;ivec++)
                    tmp+=coefs[ivec]*vecs[ivec][i];
                v[i]=tmp;
            }
        }
        
        void add_eq(T0* coefs,T0** vecs,int nvecs,T0* v)
        {
            for(int i=0;i<n;i++)
                for(int ivec=0;ivec<nvecs;ivec++)
                    v[i]+=coefs[ivec]*vecs[ivec][i];
        }
        
        T0 norm(T0* v)
        {
            return sqrt(inner(v,v));
        }
        
        void solve_y(int nvec)
        {
            for(int i=nvec-1;i>-1;i--)
            {
                y[i]=b_hat[i];
                for(int j=i+1;j<nvec;j++)
                    y[i]-=H[j][i]*y[j];
                y[i]/=H[i][i];
            }
        }
        
    };
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    template<typename T0,class C0>
    class GMRES: protected InitPtrs
    {
    private:
        const int m;
        const int dim;
        int n;
        C0& kernel;
        Vec<T0>** vecs;
        T0** Q;
        T0** H;
        T0* b_hat;
        
        T0* cos;
        T0* sin;
        
        T0* y;
        T0* ans_lcl;
        
        
        T0 calc(int ivec)
        {
            for(int j=0;j<ivec+1;j++)
                ans_lcl[j]=0.0;
            
            for(int i=0;i<n;i++)
                for(int j=0;j<ivec+1;j++)
                    ans_lcl[j]+=Q[ivec+1][i]*Q[j][i];
            
            MPI_Allreduce(ans_lcl,H[ivec],ivec+1,MPI_TYPE0,MPI_SUM,world);
            
            T0 ans_lcl_=0.0,ans;
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<ivec+1;j++)
                    Q[ivec+1][i]-=H[ivec][j]*Q[j][i];
                
                ans_lcl_+=Q[ivec+1][i]*Q[ivec+1][i];
            }
            
            MPI_Allreduce(&ans_lcl_,&ans,1,MPI_TYPE0,MPI_SUM,world);
            ans=1.0/sqrt(ans);
            
            for(int i=0;i<n;i++)
                Q[ivec+1][i]*=ans;
            
            return 1.0/ans;
        }
        
        
        T0 solve_y(int nvecs,type0* x)
        {
            for(int i=nvecs-1;i>-1;i--)
            {
                y[i]=b_hat[i];
                for(int j=i+1;j<nvecs;j++)
                    y[i]-=H[j][i]*y[j];
                y[i]/=H[i][i];
            }
            
            
            for(int i=0;i<n;i++)
                for(int ivec=0;ivec<nvecs;ivec++)
                    x[i]+=y[ivec]*Q[ivec][i];
            
            T0 norm=0.0;
            for(int ivec=0;ivec<nvecs;ivec++)
                norm+=y[ivec]*y[ivec];
            
            return sqrt(norm);
        }
        
        
        void restart()
        {
            T0 tmp_0=fabs(b_hat[m]);
            
            
            if(b_hat[m]>0.0)
                b_hat[m]=1.0;
            else
                b_hat[m]=-1.0;
            
            for(int i=m-1;i>-1;i--)
            {
                b_hat[i]=sin[i]*b_hat[i+1];
                b_hat[i+1]*=cos[i];
            }
            
            for(int i=0;i<n;i++)
            {
                Q[0][i]*=b_hat[0];
                for(int ivec=1;ivec<m+1;ivec++)
                    Q[0][i]+=b_hat[ivec]*Q[ivec][i];
            }
            
            b_hat[0]=tmp_0;
            for(int i=1;i<m+1;i++)
                b_hat[i]=0.0;
        }
        
        void start(T0* b)
        {
            T0 ans_lcl=0.0,norm,inv_norm;
            for(int i=0;i<n;i++)
                ans_lcl+=b[i]*b[i];
            
            MPI_Allreduce(&ans_lcl,&norm,1,MPI_TYPE0,MPI_SUM,world);
            norm=sqrt(norm);
            inv_norm=1.0/norm;
            
            for(int i=0;i<n;i++)
                Q[0][i]=inv_norm*b[i];
            
            b_hat[0]=norm;
        }
        
    protected:
    public:
        
        
        GMRES(MAPP* mapp,int m_,int dim_,C0& kernel_):
        InitPtrs(mapp),
        m(m_),
        dim(dim_),
        kernel(kernel_)
        {
            
            Q=new T0*[m+1];
            vecs=new Vec<T0>*[m+1];
            for(int ivec=0;ivec<m+1;ivec++)
                vecs[ivec]=new Vec<T0>(atoms,dim);
            
            H=new T0*[m+1];
            *H=new T0[m*(m+1)/2+1];
            for(int i=1;i<m+1;i++)
                H[i]=H[i-1]+i;
            
            refresh();
            
            b_hat=new T0[m+1];
            cos=new T0[m];
            sin=new T0[m];
            y=new T0[m];
            ans_lcl=new T0[m+1];
            
        }
        
        ~GMRES()
        {
            for(int ivec=0;ivec<m+1;ivec++)
                delete vecs[ivec];
            delete [] vecs;
            delete [] Q;
            
            delete [] *H;
            delete [] H;
            delete [] b_hat;
            delete [] cos;
            delete [] sin;
            delete [] y;
        }
        
        void refresh()
        {
            for(int ivec=0;ivec<m+1;ivec++)
                Q[ivec]=vecs[ivec]->begin();
            n=atoms->natms*dim;
        }
        
        
        T0 solve(T0 tol,T0* b,T0* x,int& iter,T0& norm)
        {
            T0 last_h=0.0,tmp_0,res_norm,b_norm;
            int max_iter=1;
            iter=0;
            start(b);
            
            b_norm=b_hat[0];
            
            res_norm=b_norm;
            
            while (res_norm>tol && max_iter)
            {
                for(int i=0;i<m;i++)
                {
                    atoms->update(vecs[i]);
                    kernel(vecs[i],vecs[i+1]);
                    iter++;
                    
                    last_h=calc(i);
                    
                    for(int j=0;j<i;j++)
                    {
                        tmp_0=cos[j]*H[i][j]-sin[j]*H[i][j+1];
                        H[i][j+1]=sin[j]*H[i][j]+cos[j]*H[i][j+1];
                        H[i][j]=tmp_0;
                    }
                    
                    tmp_0=sqrt(H[i][i]*H[i][i]+last_h*last_h);
                    
                    cos[i]=H[i][i]/tmp_0;
                    sin[i]=-last_h/tmp_0;
                    
                    H[i][i]=cos[i]*H[i][i]-sin[i]*last_h;
                    b_hat[i+1]=sin[i]*b_hat[i];
                    b_hat[i]*=cos[i];
                    
                    res_norm=fabs(b_hat[i+1]);

                    
                    if(res_norm<tol)
                    {
                        norm=solve_y(i+1,x);
                        return fabs(b_hat[i+1]);
                    }
                }
                
                
                
                norm=solve_y(m,x);
                
                max_iter--;
                
                if(max_iter==0)
                    return fabs(b_hat[m]);
                restart();
            }
            return b_hat[m];
        }
        
        
        
        
        
        
        
        
    };
    
    
    
    
    
}


#endif
