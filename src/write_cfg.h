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
        vec** dump_vecs;
        int ndump_vecs;
        char** dump_vecs_names;
        int ndump_vecs_names;
        
        char* file_name;
        int sorting;
        int usr_nabled;
        int tot_dim;
        int dim;
        int c_n;
        void write_file_dmd(int);
        void write_file_md(int);
        void x2s(int);
        int find_type(type0*,type0*);
        type0* clr_r;
        type0* tmp_arr;
        type0 thresh;
        type0 strech;
        void add_to_dump_vec_nams(const char*);
        void open_write(int,FILE*&,FILE*&);
    public:
        Write_cfg(MAPP *,int,char**);
        ~Write_cfg();
        void write_file(int);
        void init();
        void fin();
    };
    
}
#endif
#endif
