#CC=g++
CFLAGS=-Wall -fopenmp
#-std=c++11


all: fs prs

fs: nodemap.cpp fringesearch.cpp
	g++ $(CFLAGS)  nodemap.cpp fringesearch.cpp fs_main.cpp -o fs

prs: nodemap.cpp fringesearch.cpp
	g++ $(CFLAGS)  nodemap.cpp fringesearch.cpp prs_main.cpp -o prs

clean:
	rm fs prs
