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
   int pages = atoi(argv[argIndex++]);

   //won't bother parallelizing initialization... it's not the bottleneck

   //for Sparse Row Storage
   int nnz = (pages-1) * 2;
   int rows = pages;
   int cols = pages;
   double A[nnz]; //assumption: pages > 1
   int Ai[rows+1];
   Ai[rows] = nnz; //last index of Ai is the total number of nonzeroes
   int Aj[nnz];
   //each page i is connected to i+1 and i-1, unless on either end
   int aCounter = 0;
   for (int to = 0; to < pages; to++) { //ROW
      Ai[to] = aCounter; //mark the beginning of the row
      for (int from = 0; from < pages; from++) { //COL
         int before = from - 1;
         int after = from + 1;

         if (to == before || to == after) {
            double edgeWeight = 0.5;
            if (before == -1 || after == pages) {
               edgeWeight = 1.0;
            }

            A[aCounter] = edgeWeight * (1.0 - DAMPING_FACTOR);
            Aj[aCounter] = from;
            aCounter++;
         }

      }
   }





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
      #pragma omp parallel for
      for (int row = 0; row < pages; row++) {
         newRank[row] = 0.0;

         int rowCounter = Ai[row];
         int rowEnd = Ai[row+1];
         for (int col = 0; col < pages; col++) {
            double matrixValue = 1.0 / pages * DAMPING_FACTOR;
            if (rowCounter < rowEnd && Aj[rowCounter] == col) { //nonzero field
               matrixValue += A[rowCounter++];
            }
            newRank[row] += rank[col] * matrixValue;
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



   //print out the result:
   int minIndex = 0;
   int maxIndex = 0;
   double sum = 0;

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
