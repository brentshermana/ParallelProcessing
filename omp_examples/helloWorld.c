/* Jim E. Jones
   FIT Mathematics
   Math 4280
   10 March 2015 */

/* Inclusions */
#include <omp.h>      /* OpenMP */
#include <stdlib.h>   /* malloc */
#include <stdio.h>    /* prints */

int main (int argc, char *argv[]) 
{
  int threadCount=1;
  int threadID;
  int arg_index=1;


  /* Read optional command line arguments to overide defaults and allow prints*/
  while (arg_index < argc)
  {
    /* Optionally redefine number of threads from command line */
    if ( strcmp(argv[arg_index], "-threads") == 0 )
    {
       arg_index++;
       threadCount = atoi(argv[arg_index++]);
    }
  }

  /* set number of threads */
  omp_set_num_threads(threadCount);


  /* Begin parallel block, each thread prints hello */

  #pragma omp parallel private(threadID)
  {
    /* Obtain thread number */
    threadID = omp_get_thread_num();
    printf("Hello World from thread = %d\n", threadID);
  
    /* Only master thread does this print */
    #pragma omp master
      printf("Number of threads = %d\n", threadCount);
      threadCount = omp_get_num_threads();
    #pragma end master`

  }  /* End of parallel block*/
}
