#include <stdlib.h>   /* malloc */
#include <stdio.h>    /* prints */
#include <mpi.h>      /* timings */
#include <string.h>

/* Set time and space extents */
#define TIME_INTERVAL 2.0
#define SPACE_INTERVAL 100.0
#define HEAT_SOURCE 1.0

/*---------------------------------------------
This code (heat.c) solves the discretized heat equation,
Ut-Uxx=HEAT_SOURCE, on the intervals

0<x<SPACE_INTERVAL and 0<t<TIME_INTERVAL

with zero boundary codinitions and initial condition.
It uses an explicit Euler scheme, and checks that the
CFL stability condition is not violated.

The code uses MPI functions but only coded to run on
one process.

The number of spatial unknowns (N) and number of time
steps (TIMESTEPS) can be set from the command line by
using -n and -t each followed by the appropriate number.
Additionally, to print the computed solution after a
certain time step use -debug follwed by the time number
of the time step you want to print at.  For example:

If you compile in ~/ASSN4, naming the executable heat

mpicc -o heat heat.c

Then to run with 20 spatial unknowns for 100 time steps
the last line of your subimssion script (.sh) should be:

 mpiexec -n 1 ~/ASSN4/heat  -n 20 -t 100

To do the same run but print the solution at the end:

 mpiexec -n 1 ~/ASSN4/heat  -n 20 -t 100 -debug 100

To do the first time step, stop, and print:

 mpiexec -n 1 ~/ASSN4/heat  -n 20 -t 100 -debug 1

-----------------------------------------------*/
int main(int argc, char* argv[]) {

    /* declarations */
    double h, hSquared, k;
    double *uNew, *uOld;
    double startTime, endTime;
    int ierr,numprocs, myid;
    int i,timeCounter;

    /* default values for problem size, number of time steps, and printouts */
    int N=4;
    int TIMESTEPS=10;
    int DEBUG=0;
    int arg_index=1;

    /* declarations for optional debug printouts */
    FILE       *fp;
    char        filename[100];

    /* intialize MPI */
    ierr = MPI_Init(&argc, &argv);
    if (ierr != MPI_SUCCESS) {
        printf("Error in MPI_Init = %i\n",ierr);
        return -1;
    }
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if (ierr != MPI_SUCCESS) {
        printf("Error in MPI_Comm_size = %i\n",ierr);
        return -1;
    }
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (ierr != MPI_SUCCESS) {
        printf("Error in MPI_Comm_rank = %i\n",ierr);
        return -1;
    }

    /* Read optional command line arguments to overide defaults and allow prints*/
    while (arg_index < argc)
    {
        /* Optionally redefine local problem size from command line */
        if ( strcmp(argv[arg_index], "-n") == 0 )
        {
            arg_index++;
            N = atoi(argv[arg_index++]);
        }
            /* Optionally redefine number of timesteps from command line */
        else if ( strcmp(argv[arg_index], "-t") == 0 )
        {
            arg_index++;
            TIMESTEPS = atoi(argv[arg_index++]);
        }
            /* Read optional debug print flag from command line */
        else if ( strcmp(argv[arg_index], "-debug") == 0 )
        {
            arg_index++;
            DEBUG = atoi(argv[arg_index++]);
        }
    }

    /* error checking: this code is setup to run properly only on one processor */
    if (numprocs > 1) {
        printf("Error serial code only");
        return(-1);
    }

    /* error checking: this code only allows prints for small problems (N<101) */
    if ((N > 101) & (DEBUG!=0)) {
        printf("Error because debugging prints only available for N < 101");
        return(-1);
    }

    /* time step size */
    k=TIME_INTERVAL/TIMESTEPS;
    /* grid step size */
    h=SPACE_INTERVAL/(N+1);
    hSquared=h*h;

    /* error checking: Courant, Friedrichs, Lewy conditon means we expect stable
       results only if 2k < h^2.  Comment out return to allow possible debugging with
       a single time step, but numerical results are almost certainly inaccurate.*/
    if (2*k > hSquared) {
        printf("Error because CFL condition violated, decrease N or increase TIMESTEPS");
        return(-1);
    }

    /* attempt to allocate memory for new and old temperature and zero out if successful*/
    if ( (uNew=(double *)malloc((N+2)*sizeof(double))) == NULL ) {
        printf("Failed memory allocation for uNew");
        return(-1);
    }
    if ( (uOld=(double *)malloc((N+2)*sizeof(double))) == NULL ) {
        printf("Failed memory allocation for uOld");
        return(-1);
    }
    for (i=0; i < N+2; i++)
    {
        uNew[i] = 0.0;
        uOld[i] = 0.0;
    }

    /* start timer for code execution time */
    startTime=MPI_Wtime();

    /* time stepping loop */
    for (timeCounter=0; timeCounter < TIMESTEPS; timeCounter++)
    {
        /* for each gridpoint, compute new temperature */
        for (i=1; i < N+1; i++)
        {
            uNew[i]=uOld[i]+k*HEAT_SOURCE+k*(uOld[i+1]-2*uOld[i]+uOld[i-1])/hSquared;
        }
        /* replace old temps to get ready for next timestep */
        for (i=1; i < N+1; i++)
        {
            uOld[i]=uNew[i];
        }
        if (DEBUG == timeCounter+1) break;
    }

    endTime=MPI_Wtime();

    /* optional printout of solution, first check if file is available */
    if (DEBUG != 0)
    {
        sprintf(filename,"ZZZOutputSolution.%05d",myid);
        if ((fp = fopen(filename, "w")) == NULL)
        {
            printf("Error: can't open output file %s\n", filename);
            DEBUG=0;
        }
    }
    if (DEBUG != 0)
    {
        for (i=1; i < N+1; i++)
        {
            fprintf(fp,"u(%i)=%f\n",i,uNew[i]);
        }
        fclose(fp);
    }

    if (myid==0)
    {
        printf("N=%d, P=%d, and Tsteps= %d\n",N,numprocs,TIMESTEPS);
        printf("runtime is=%.16f\n",endTime-startTime);
    }


    free(uNew);
    free(uOld);

    MPI_Finalize();
}


