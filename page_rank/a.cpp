#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#define DAMPING_FACTOR 0.15
#define MATVEC_ITERATIONS 1000

//two arguments
//1: number of threads
//2: number of web pages

int main (int argc, char* argv[]) {
   //copied from jacobi.c
   int argIndex = 1;
   int id = 0;
   int numThreads = atoi(argv[argIndex++]);
   omp_set_num_threads(numThreads);

   #pragma omp parallel private(numThreads, id)
   {
       /* Obtain thread number */
       id = omp_get_thread_num();

       /* Only master thread does this */
       if (id == 0)
       {
         numThreads = omp_get_num_threads();
         printf("Number of threads = %d\n", numThreads);
       }
   }  /* All threads join master thread and disband */

   //-----------------------------------------------------
   const int pages = atoi(argv[argIndex++]);

   //won't bother parallelizing initialization... it's not the bottleneck
   double ** matrix = new double*[pages];
   for (int i = 0; i < pages; i++) {
      matrix[i] = new double[pages];
   }


   //each page i is connected to i+1 and i-1, unless on either end
   for (int page = 0; page < pages; page++) { //from
      int before = page - 1;
      int after = page + 1;

      double edgeWeight = 0.5;
      if (before == -1 || after == pages) { //for first and last, edge weight are 1
         edgeWeight = 1.0;
      }
      for (int adj = 0; adj < pages; adj++) { // to
         matrix[adj][page] = 0;

         if (adj == before) {
            matrix[adj][page] += edgeWeight * (1.0 - DAMPING_FACTOR);
         }
         if (adj == after) {
            matrix[adj][page] += edgeWeight * (1.0 - DAMPING_FACTOR);
         }
         matrix[adj][page] += DAMPING_FACTOR * (1.0 / pages);
      }
   }

   printf("created the matrix");

   //TEST: PRINT THE MATRIX:
   /*
   for (int row = 0; row < pages; row++) {
      for (int col = 0; col < pages; col++) {
         printf("%.3f ", matrix[row][col]);
      }
      printf("\n");
   }
   */


   double rank[pages];
   //inital estimate assumes all pages are equal
   for (int i = 0; i < pages; i++) {
      rank[i] = 1.0 / pages;
   }

   double startTime;
   if (id == 0)
      startTime = omp_get_wtime();

   double newRank[pages];
   for (int i = 0; i < MATVEC_ITERATIONS; i++) {

      //multiply matrix by vector:
      #pragma omp parallel for //each thread owns a matrix row, so we don't have to worry about reduction, etc.
      for (int row = 0; row < pages; row++) {
         newRank[row] = 0.0;
         for (int col = 0; col < pages; col++) {
            newRank[row] += rank[col] * matrix[row][col];
         }
      }
      //copy over:
      #pragma omp parallel for
      for (int j = 0; j < pages; j++) {
         rank[j] = newRank[j];
      }
   }

   double elapsed;
   if (id == 0)
      elapsed = omp_get_wtime() - startTime;

   int minIndex = 0;
   int maxIndex = 0;
   double sum = 0;
   //print out the result:
   for (int i = 0; i < pages; i++) {
      if (rank[i] < rank[minIndex]) {
         minIndex = i;
      }
      if (rank[i] > rank[maxIndex]) {
         maxIndex = i;
      }
      sum += rank[i];
      printf("The rank for page %d is: %.5f\n", i, rank[i]);
   }
   printf("===============================\n");
   printf("Sum: %.8f\n", sum);
   printf("Duration: %.7f\n", elapsed);
   printf("Highest ranked page: %d, with rank %.4f\n", maxIndex, rank[maxIndex]);
   printf("Lowest ranked page: %d, with rank %.4f\n", minIndex, rank[minIndex]);
}
