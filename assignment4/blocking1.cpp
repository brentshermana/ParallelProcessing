#include <stdlib.h>   /* malloc */
#include <stdio.h>    /* prints */
#include <mpi.h>      /* timings */
#include <cstring>
/* Set time and space extents */
#define TIME_INTERVAL 2.0
#define SPACE_INTERVAL 100.0
#define HEAT_SOURCE 1.0

using namespace std;

int nps;
int pn;
int N;
int myN;

/*---------------------------------------------This code (heat.c) solves the discretized heat equation, Ut-Uxx=HEAT_SOURCE, on the intervals
0<x<SPACE_INTERVAL and 0<t<TIME_INTERVAL
with zero boundary codinitions and initial condition. It uses an explicit Euler scheme, and checks that the CFL stability condition is not violated.
The code uses MPI functions but only coded to run on one process.
The number of spatial unknowns (N) and number of time steps (TIMESTEPS) can be set from the command line by using -n and -t each followed by the appropriate number. Additionally, to print the computed solution after a certain time step use -debug follwed by the time number of the time step you want to print at.  For example:
If you compile in ~/ASSN4, naming the executable heat mpicc -o heat heat.c Then to run with 20 spatial unknowns for 100 time steps
the last line of your subimssion script (.sh) should be:  mpiexec -n 1 ~/ASSN4/heat  -n 20 -t 100 To do the same run but print the solution at the end:  mpiexec -n 1 ~/ASSN4/heat  -n 20 -t 100 -debug 100 To do the first time step, stop, and print:  mpiexec -n 1 ~/ASSN4/heat  -n 20 -t 100 -debug 1 -----------------------------------------------*/
int main(int argc, char *argv[]) {


    /* default values for problem size, number of time steps, and printouts */

    N = 4;
    int TIMESTEPS = 10;
    bool DEBUG = false;
    int arg_index = 1;
    /* declarations for optional debug printouts */ //  FILE *fp;
    // char filename[100];


    /* intialize MPI */
    int ierr = MPI_Init(&argc, &argv);
    if (ierr != MPI_SUCCESS) {
        printf("Error in MPI_Init = %i\n", ierr);
        return -1;
    }
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nps);
    if (ierr != MPI_SUCCESS) {
        printf("Error in MPI_Comm_size = %i\n", ierr);
        return -1;
    }
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &pn);
    if (ierr != MPI_SUCCESS) {
        printf("Error in MPI_Comm_rank = %i\n", ierr);
        return -1;
    }

    /* Read optional command line arguments to override defaults and allow prints*/
    while (arg_index < argc) {       /* Optionally redefine local problem size from command line */
        if (strcmp(argv[arg_index], "-n") == 0) {
            arg_index++;
            N = atoi(argv[arg_index++]);
        }       /* Optionally redefine number of timesteps from command line */
        else if (strcmp(argv[arg_index], "-t") == 0) {
            arg_index++;
            TIMESTEPS = atoi(argv[arg_index++]);
        }       /* Read optional debug print flag from command line */
        else if (strcmp(argv[arg_index], "-debug") == 0) {
            DEBUG = true;
            arg_index++;
        }
        else if (strcmp(argv[arg_index], "-times") == 0) {
            DEBUG = true;
            arg_index++;
        }
        else {
            arg_index++;
        }
    }

    double k = TIME_INTERVAL / TIMESTEPS;   /* grid step size */
    myN = N/nps;
    double h = SPACE_INTERVAL / (N + 1);
    double hSquared = h * h;
    /*
     * error checking: Courant, Friedrichs, Lewy conditon means we expect stable
     * results only if 2k < h^2.  Comment out return to allow possible debugging with
     * a single time step, but numerical results are almost certainly inaccurate.
     */
    if (2 * k > hSquared) {
        printf("Error because CFL condition violated, decrease N or increase TIMESTEPS");
        return (-1);
    }

    double* newV = new double[myN+2];
    double* oldV = new double[myN+2];

    for (int i = 0; i < myN + 2; i++) {
        newV[i] = 0.0;
        oldV[i] = 0.0;
    }
    bool leftP = pn > 0;
    bool rightP = pn < nps-1;
    /* start timer for code execution time */
    /* time stepping loop */
    double startTime = MPI_Wtime();

    for (int timeCounter = 0; timeCounter < TIMESTEPS; timeCounter++) {
        /* for each gridpoint, compute new temperature */
        for (int i = 1; i < myN + 1; i++) {
            newV[i] = oldV[i] + k * HEAT_SOURCE + k * (oldV[i + 1] - 2 * oldV[i] + oldV[i - 1]) / hSquared;
        }
        /* replace old temps to get ready for next timestep */
        for (int i = 1; i < myN + 1; i++) {
            oldV[i] = newV[i];
        }

        if (pn % 2 == 1) {
            //receive
            if (leftP) {
                MPI_Recv(&oldV[0], 1, MPI_DOUBLE, pn-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rightP) {
                MPI_Recv(&oldV[myN+1], 1, MPI_DOUBLE, pn+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            //send
            if (rightP) {
                MPI_Send(&newV[myN], 1, MPI_DOUBLE, pn+1, 0, MPI_COMM_WORLD);
            }
            if (leftP) {
                MPI_Send(&(newV[1]), 1, MPI_DOUBLE, pn-1, 0, MPI_COMM_WORLD);
            }
        }
        else {
            //send
            if (rightP) {
                MPI_Send(&newV[myN], 1, MPI_DOUBLE, pn+1, 0, MPI_COMM_WORLD);
            }
            if (leftP) {
                MPI_Send(&(newV[1]), 1, MPI_DOUBLE, pn-1, 0, MPI_COMM_WORLD);
            }
            //receive
            if (rightP) {
                MPI_Recv(&oldV[myN+1], 1, MPI_DOUBLE, pn+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (leftP) {
                MPI_Recv(&oldV[0], 1, MPI_DOUBLE, pn-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

    double endTime = MPI_Wtime();
    if (pn == 0) {
        printf("N=%d, P=%d, and Tsteps= %d\n", N, nps, TIMESTEPS);
        printf("runtime is=%.16f\n", endTime - startTime);
    }
    if (DEBUG) {
        double* all = new double[N];
        MPI_Gather(&oldV[1], myN, MPI_DOUBLE, all, myN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (pn == 0) {
            for (int i = 0; i < N; i++) {
                printf("%.6f ", all[i]);
            }
            printf("\n");
        }

    }

    delete[] newV;
    delete[] oldV;
    MPI_Finalize();
}