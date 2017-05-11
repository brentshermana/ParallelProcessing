/* Jim E. Jones
   FIT Mathematics
   Math 4280
   12 March 2016 */


/* Inclusions */
#include <omp.h>      /* OpenMP */
#include <stdlib.h>   /* malloc */
#include <stdio.h>    /* prints */

/* Set time and space extents */
#define ROD_LENGTH 1.0
#define HEAT_SOURCE 1.0

/*---------------------------------------------

This code (jacobi.c) solves the discretized steady
state heat equation,
-Uxx=HEAT_SOURCE, on the intervals

0<x<ROD_LENGTH

with zero boundary conditions at left and right endpoints
It uses finite differences and a Jacobi method to find
an approximate solution to the equations begining with
a zero initial guess.

The code uses OpenMP functions but only coded to run on
one process.

The number of spatial unknowns (N) and number of jacobi
steps (NUMSTEPS) can be set from the command line by
using -n and -s each followed by the appropriate number.

To set number of threads use -p followed by number.


gcc -fopenmp -o jacobi jacobi.c

Then to run with 20 spatial unknowns for 100 jacobi steps
and 2 threads, type:

./jacobi  -n 20 -s 100 -p 2

-----------------------------------------------*/
int main(int argc, char* argv[]) {

  /* declarations */
  double h, hSquared, k;
  double *uNew, *uOld;
  int threadCount=1, myid=0;
  int i,stepCounter;

  /* default values for problem size, number of time steps, and printouts */
  int N=4;
  int NUMSTEPS=10;
  int arg_index=1;

  /* timers */
  double startTime, endTime;

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
      else if ( strcmp(argv[arg_index], "-s") == 0 )
      {
         arg_index++;
         NUMSTEPS = atoi(argv[arg_index++]);
      }
      /* Optionally redefine number of threads */
      else if ( strcmp(argv[arg_index], "-p") == 0 )
      {
         arg_index++;
         threadCount = atoi(argv[arg_index++]);
      }
  }

  /* set number of threads */
  omp_set_num_threads(threadCount);
  /* each thread print "hello' */
#pragma omp parallel private(threadCount, myid)
     {
       /* Obtain thread number */
       myid = omp_get_thread_num();
       printf("Hello World from thread = %d\n", myid);

       /* Only master thread does this */
       if (myid == 0)
       {
         threadCount = omp_get_num_threads();
         printf("Number of threads = %d\n", threadCount);
       }
    }  /* All threads join master thread and disband */


  /* grid step size */
  h=ROD_LENGTH/(N+1);
  hSquared=h*h;

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
  startTime = omp_get_wtime();


  /* time stepping loop */
  for (stepCounter=0; stepCounter < NUMSTEPS; stepCounter++)
    {
       /* for each gridpoint, compute new temperature */
#pragma omp parallel for
       for (i=1; i < N+1; i++)
       {
          uNew[i]=0.5*(hSquared*HEAT_SOURCE+uOld[i+1]+uOld[i-1]);
       }
       /* replace old temps to get ready for next timestep */
#pragma omp parallel for
       for (i=1; i < N+1; i++)
       {
          uOld[i]=uNew[i];
       }
    }

  /* stop timer */
  endTime = omp_get_wtime();


  printf("N=%d, P=%d, and steps= %d\n",N,threadCount,NUMSTEPS);

  /* output solution time */
  printf("     runtime = %.16e\n",endTime-startTime);

  free(uNew);
  free(uOld);
}


