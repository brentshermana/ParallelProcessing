#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define ITERATIONS 1000
#define DAMPING_FACTOR 0.15

int max(int a, int b) {
   if (a > b) return a;
   return b;
}
int min (int a, int b) {
   if (a < b) return a;
   return b;
}
bool inBounds(int x, int lo, int hi) {
   return x <= hi && x >= lo;
}

//Adapted from Stack Overflow (lost the link)
double **allocContiguousDouble(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);
    return array;
}

int main(int argc, char* argv[]) {
   /* Initialize MPI and get number of processes and my number or rank*/
   MPI_Init(&argc, &argv);
   int procs, id;
   MPI_Comm_size(MPI_COMM_WORLD, &procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   int numPages = atoi(argv[1]);

   if (numPages < procs) {
      if (id == 0)
         printf("number of procs > matrix rows, run again with fewer procs or more rows");
      return 0;
   }

   //does each process have access to the command line args? better not to assume
   MPI_Bcast(&numPages, 1, MPI_INT, 0, MPI_COMM_WORLD);


   //ratio of pages:procs. guaranteed to be an integer
   int ratio = numPages/procs;

   int startRow = id*ratio;
   int numRows = ratio;

   //first index: the row number. second, the index of the matrix
   double ** matrixRows= allocContiguousDouble(numRows, numPages);
   double vectorParts[numRows]; //the portion of the vector this thread is responsible for

   for (int i = 0; i < numRows; i++) {
      vectorParts[i] = 1.0 / numPages;
   }

   //dampening matrix:
   for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numPages; j++) {
         matrixRows[i][j] = (1.0/numPages) * DAMPING_FACTOR;
      }
   }
   //links between pages
   for (int to = startRow; to < startRow + numRows; to++) { //row = to index
      int from1 = to+1;
      if (inBounds(from1, 0, numPages-1)) { //from1 is valid page index
         double edgeWeight = 0.5;
         if (from1 == 0 || from1 == numPages-1) { //from1 only has one edge
            edgeWeight = 1.0;
         }
         matrixRows[to - startRow][from1] += edgeWeight * (1.0 - DAMPING_FACTOR);
      }

      int from2 = to-1;
      if (inBounds(from2, 0, numPages-1)) { //from1 is valid page index
         double edgeWeight = 0.5;
         if (from2 == 0 || from2 == numPages-1) { //from1 only has one edge
            edgeWeight = 1.0;
         }
         matrixRows[to-startRow][from2] += edgeWeight * (1.0-DAMPING_FACTOR);
      }
   }

   double startTime;
   if (id == 0)
      startTime = MPI_Wtime();

   double fullVector[numPages];
   //now run the calculations:
   for (int iter = 0; iter < ITERATIONS; iter++) {
      MPI_Allgather(vectorParts, numRows, MPI_DOUBLE, fullVector, numRows, MPI_DOUBLE, MPI_COMM_WORLD);

      for (int row = 0; row < numRows; row++) { //iterate through assigned rows
         vectorParts[row] = 0;
         for (int col = 0; col < numPages; col++) {
            vectorParts[row] += matrixRows[row][col] * fullVector[col];
         }
      }
   }

   double elapsed;
   if (id == 0)
      elapsed = MPI_Wtime() - startTime;


   MPI_Gather(vectorParts, numRows, MPI_DOUBLE, fullVector, numRows, MPI_DOUBLE, 0, MPI_COMM_WORLD);


   int minIndex = 0;
   int maxIndex = 0;
   double sum = 0;
   if (id == 0) {
      for (int i = 0; i < numPages; i++) {
         if (fullVector[i] < fullVector[minIndex]) {
            minIndex = i;
         }
         if (fullVector[i] > fullVector[maxIndex]) {
            maxIndex = i;
         }
         sum += fullVector[i];
         printf("Page %3d has rank %.3f\n", i, fullVector[i]);
      }
   }
   printf("===============================\n");
   printf("Sum: %.8f\n", sum);
   printf("Duration: %.7f\n", elapsed);
   printf("Highest ranked page: %d, with rank %.4f\n", maxIndex, fullVector[maxIndex]);
   printf("Lowest ranked page: %d, with rank %.4f\n", minIndex, fullVector[minIndex]);


   MPI_Finalize();
}
