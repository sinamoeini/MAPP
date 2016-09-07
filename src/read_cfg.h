#ifdef Read_Style
    ReadStyle(Read_cfg,cfg)
#else
#ifndef __MAPP__read_cfg__
#define __MAPP__read_cfg__
#include "read.h"
#include "type_def.h"
#include <stdio.h>
namespace MAPP_NS
{
    class Read_cfg : public Read
    {
    private:
        
        type0** H0;
        type0** eta;
        type0** trns;
        type0** H_x;
        type0** H_x_d;
                
        type0 basic_length;
        type0 R;
        int entry_count;
        int ext_cfg;
        int vel_chk;
        
        void read_header(FILE*);
        void read_atom(FILE*);
        void set_box();
        void add_atom_read_x(int);
        void add_atom_read_x();
        
        int ch_x_d,ch_type,ch_id,ch_c;
        int dmd_no_types;
        
        int curr_id;
        byte* ch_buff;
        int ch_buff_sz;
        class vec** vec_lst;
        int nvecs;
        type0* tmp_buff;
        int tmp_buff_size;
        char* file_name;
        
        void set_vec_lst();
        void set_vecs();
        
    protected:
    public:
        Read_cfg(int,char**);
        ~Read_cfg();
    };
}
#endif
#endif

