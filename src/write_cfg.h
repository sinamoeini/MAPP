#ifdef Write_Style
WriteStyle(Write_cfg,cfg)
#else
#ifndef __MAPP__write_cfg__
#define __MAPP__write_cfg__

#include "write.h"
namespace MAPP_NS
{
    class Write_cfg: public Write
    {
    private:
    protected:
        char* file_name;
        int* vec_list;
        char** vec_name;
        int lst_compltd;
        int no_vecs;
        int id_n;
        int type_n;
        int sorting;
        int tot_dim;
        int x_dim;
        int dim;
        int c_n;
        void cmplt_lst();
        void write_file_dmd(int);
        void write_file_md(int);
        void x2s(int);
    public:
        Write_cfg(MAPP *,int,char**);
        ~Write_cfg();
        void write_file(int);
    };
    
}
#endif
#endif
