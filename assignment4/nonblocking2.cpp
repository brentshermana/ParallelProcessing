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

double k;
double hSquared;

double mathStuff(double left, double center, double right) {
    return center + k * HEAT_SOURCE + k * (right - 2 * center + left) / hSquared;
}
int ltog(int l) {
    return l + (myN*pn) - 1;
}
int gtol(int g) {
    return g - (myN*pn) + 1;
}
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
        else {
            arg_index++;
        }
    }

    k = TIME_INTERVAL / TIMESTEPS;   /* grid step size */
    myN = N/nps;
    double h = SPACE_INTERVAL / (N + 1);
    hSquared = h * h;
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
        //first, compute gridpoints which must be sent:
        newV[1] = mathStuff(oldV[0], oldV[1], oldV[2]);
        newV[myN] = mathStuff(oldV[myN-1], oldV[myN], oldV[myN+1]);

        //make IO calls as soon as possible
        MPI_Request lSend = MPI_REQUEST_NULL;
        MPI_Request rRec = MPI_REQUEST_NULL;
        MPI_Request rSend = MPI_REQUEST_NULL;
        MPI_Request lRec = MPI_REQUEST_NULL;
        if (leftP) {
            MPI_Isend(&(newV[1]), 1, MPI_DOUBLE, pn-1, 0, MPI_COMM_WORLD, &lSend);
            MPI_Irecv(&newV[0], 1, MPI_DOUBLE, pn-1, 0, MPI_COMM_WORLD, &lRec);
        }
        if (rightP) {
            MPI_Irecv(&newV[myN+1], 1, MPI_DOUBLE, pn+1, 0, MPI_COMM_WORLD, &rRec);
            MPI_Isend(&newV[myN], 1, MPI_DOUBLE, pn+1, 0, MPI_COMM_WORLD, &rSend);
        }

        //below - do as much work as possible between communications:
        //compute other points
        for (int i = 2; i < myN; i++) {
            newV[i] = mathStuff(oldV[i-1], oldV[i], oldV[i+1]);
        }
        //copy over items besides those which are asynchronously received
        for (int i = 1; i < myN + 1; i++) {
            oldV[i] = newV[i];
        }

        //wait for transactions to finish
        if (leftP) {
            MPI_Wait(&lSend, MPI_STATUS_IGNORE);
            MPI_Wait(&lRec, MPI_STATUS_IGNORE);
        }
        if (rightP) {
            MPI_Wait(&rSend, MPI_STATUS_IGNORE);
            MPI_Wait(&rRec, MPI_STATUS_IGNORE);
        }

        //finally, copy over values from IO:
        oldV[0] = newV[0];
        oldV[myN+1] = newV[myN+1];
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