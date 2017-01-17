//
//  comm.hpp
//  MAPP
//
//  Created by Sina on 11/21/16.
//  Copyright © 2016 MIT. All rights reserved.
//

#ifndef __MAPP__comm__
#define __MAPP__comm__

#include <mpi.h>
#include "global.h"
/*----------------------------------------------------------------------------------------------------------------
 _____   _____       ___  ___       ___  ___   _   _   __   _   _   _____       ___   _____   _   _____   __   _  
/  ___| /  _  \     /   |/   |     /   |/   | | | | | |  \ | | | | /  ___|     /   | |_   _| | | /  _  \ |  \ | | 
| |     | | | |    / /|   /| |    / /|   /| | | | | | |   \| | | | | |        / /| |   | |   | | | | | | |   \| | 
| |     | | | |   / / |__/ | |   / / |__/ | | | | | | | |\   | | | | |       / / | |   | |   | | | | | | | |\   | 
| |___  | |_| |  / /       | |  / /       | | | |_| | | | \  | | | | |___   / /  | |   | |   | | | |_| | | | \  | 
\_____| \_____/ /_/        |_| /_/        |_| \_____/ |_|  \_| |_| \_____| /_/   |_|   |_|   |_| \_____/ |_|  \_|
 
 ----------------------------------------------------------------------------------------------------------------*/
namespace MAPP_NS
{
    class Communication
    {
    private:
        static int get_rank(MPI_Comm& comm)
        {
            int i;
            MPI_Comm_rank(comm,&i);
            return i;
        }
        static int get_size(MPI_Comm& comm)
        {
            int i;
            MPI_Comm_size(comm,&i);
            return i;
        }
        int tot_n;
        /*
         * total number of nodes
         * set in: constructor
         */
        
        int my_n;
        /*
         * my node id
         * set in: constructor
         */
        
        int* p_per_n;
        /*
         * number of processors per each node
         * 1st input: node number
         * allocated in: constructor
         * deallocated in: destructor
         * set in: constructor
         */
        
        int** n_p_grid;
        /*
         * the processor ids in each node
         * 1st input: node number
         * 2nd input: ith processor in that node
         * 0 <= i < p_per_n[node]
         * n_p_grid[i][j] id of jth proc in ith node
         * allocated in: constructor
         * deallocated in: destructor
         * set in: constructor
         */
        
        
    protected:
    public:
        MPI_Comm& world;
        
        const int tot_p;
        /*
         * total number of procs. in this simulation
         * set in: constructor
         */
        const int my_p;
        /*
         * the id of my proc.
         * set in: constructor
         */
        
        int tot_p_grid[__dim__];
        /*
         * tot_p_grid[idm]: number of processors in idim dimension
         * allocated in: constructor
         * deallocated in: destructor
         * set in: auto_grid() / man_grid ()
         */
        
        int my_loc[__dim__];
        /*
         * the location of my proc. in proc. grid (domain)
         * my_loc[idim] my location in dimension idim
         * allocated in: constructor
         * deallocated in: destructor
         * set in: auto_grid() / man_grid ()
         */
        
        int neigh_p[__dim__][2];
        /*
         * the number of my neighbor processors
         * 1st input: dimension
         * 2nd input: direction; 0 for behind me, 1 for infront of me
         * neigh_p[idim][idir]: id of my neighbor proc. in my idim and idir direction
         * allocated in: constructor
         * deallocated in: destructor
         * set in: auto_grid() / man_grid ()
         */
        
        type0 s_lo[__dim__];
        type0 s_hi[__dim__];
        
        
        Communication(MPI_Comm&);
        ~Communication();
        
        void auto_grid(type0(&)[__dim__][__dim__]);
        /*
         * does an automatic griding of procs.
         * gets H as an input
         */
        
        void __man_grid__(int(&)[__dim__]);
        /*
         * does a manual griding of procs.
         * gets the grid as input
         */
        
        //move this to comm
        bool grid_established;
        
        type0 skin;
        unsigned long xchng_id;
        void man_grid(int (&)[__dim__]);
        void auto_grid();
        void reset();
    };
}
#endif
