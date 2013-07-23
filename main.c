// MD Simulator
// By Jacob Wagner

#include "mpi.h"
#include "headers.h"

//global variables
Control System;
Prototypes* Prototype;
Particles* Particle;
NeighborList* Neighbor;

#include "prototypes.h"
#include "rngs.h"
#include "rngs.c"
#include "cpu_time.c"
#include "closing.c"
#include "err.c"
#include "forces.c"
#include "initialize.c"
#include "integrate.c"
#include "output.c"

int main(int argc, char **argv)
{
	//double ctime;
	//declare variables
    int moreTimeStep;
	int numprocs, rank, namelen;

    moreTimeStep = 1;
    System.counter = 0;
    
    //set-up MPI stuff
    MPI_Comm universe = MPI_COMM_WORLD;
	// Initialize
	MPI_Init(&argc, &argv);
 	MPI_Comm_size(universe, &numprocs);
  	MPI_Comm_rank(universe, &rank);
	double start_time = MPI_Wtime();
	
    initialize_system();//sets stepLimit
    System.pressure /= 1.6605;
    wrap_PBC();

	double init_time = MPI_Wtime();
    //ctime = cpu_time();
    printf("initialization was finished\n");

    while(moreTimeStep)
    {
        if((System.counter % 10) == 0)
        {
            build_neighbor_list();
            //printf("build neighbor list\n");
            printf("counter %d\n", System.counter);

        }

        interval_output();

        //printf("calc forces\n"); /**/
        calc_forces();
        //printf("integrate\n"); /**/
        time_step();

        System.counter ++;
        if(System.counter>=System.n_steps) moreTimeStep = 0;

    }

    
    double out_time = MPI_Wtime();
    VMD_output();

	//ctime = cpu_time() - ctime;

	//printf("CPU time was %16.6f secs.\n",ctime);
    printf("\nin closing\n");
	FreeMemory();
    printf("end program.\n");
	double end_time = MPI_Wtime();
	double total_time = end_time - start_time;
	printf("total run time was %lf seconds with %d cores\n", total_time, numprocs);
	printf("INIT: %lf\n", init_time - start_time);
	printf("FORCE,INTEGRATE,NEIGHBOR: %lf\n", out_time - init_time);
	printf("OUT: %lf\n", end_time - out_time);
	//MPI_finalize

   MPI_Finalize(); 
	return 0;

}
