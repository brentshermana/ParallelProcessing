#include "nodemap.h"

#ifndef FRINGESEARCH_H
#define FRINGESEARCH_H

#include <list>
#include <cstdlib>

using namespace std;

struct fs {
	int iterations;

	int threshold;

	MetaMap * mmap;

	int increment;

	coord start;

	coord * goals;

	int goalsFound;

	int numGoals;

	list<Node*> * now; // = new list<Node*>();
	list<Node*> * later; // = new list<Node*>();

	list<Node*> ** paths; //as many paths as there are goals
};

int manhattan(int, int, int, int);

int manhattan(coord a, coord b);

fs* buildFS(MetaMap* mmap, int increment, coord start, coord* goals, int numGoals);

int fsearch(fs* fs, int maxIterations);

list<coord>* hlsearch(MetaMap * mmap, coord start, coord goal);

//list<Node*>* getPath(fs* fs, Node* end);
list<Node*>* getPath(Node* end);

#endif
