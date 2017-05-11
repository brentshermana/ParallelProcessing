all: omp mpi

omp: a b

mpi: c

a: a.cpp
	gcc -fopenmp -o a a.cpp -lstdc++

b: b.cpp
	gcc -fopenmp -o b b.cpp -lstdc++

c: c.cpp
	mpicc -o c c.cpp

clean: rm *.o a b c
