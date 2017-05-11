#include <cstdlib>
#include <cstddef>
#include <stdio.h>
#include <stdint.h>

#include <iostream>

#include "fringesearch.h"

using namespace std;

int manhattan(int a, int b, int ai, int bi) {
	return abs(a-ai) + abs(b-bi);
}

int manhattan(coord a, coord b) {
	return manhattan(a.x, a.y, b.x, b.y);
}

fs* buildFS(MetaMap* mmap, int increment, coord start, coord* goals, int numGoals) {
	fs* search = (fs*)malloc(sizeof(fs));
	search->mmap = mmap;
	search->increment = increment;
	search->start = start;
	search->goals = goals;
	search->numGoals = numGoals;
	search->goalsFound = 0;
	search->iterations = 0;

	int minH = INT_MAX;
	for (int g = 0; g < numGoals; g++) {
		int temp = manhattan(start, goals[g]);
		if (minH > temp) minH = temp;
	}
	search->threshold = minH;

	search->now = new list<Node*>();
	search->later = new list<Node*>();
	search->paths = new list<Node*>*[numGoals];
	for (int i = 0; i < numGoals; i++) {
		//all paths start null, will be set once a path is found
		search->paths[i] = NULL;
	}

	Node* origin = getNode(search->mmap->real, search->start.x, search->start.y);
	origin->owner = origin;
	origin->cost = 0;
	search->now->push_front(origin);

	return search;
}

list<Node*>* getPath(Node* end) {
	list<Node*>* path = new list<Node*>();
	Node* temp = end;
	while (temp != NULL) {
		path->push_front(temp);
		temp = temp->parent;
	}
	return path;
}

int fsearch(fs* fs, int maxIterations) {
	//cout << "fsearch call" << endl << flush;
	int lastListSwap = -1;
	for (int i = 0; i < maxIterations; i++) {
		//cout << "it " << i << endl << flush;
		//if we've exanded through the current threshold
		if (fs->now->empty()) {
			//cout << "FS Reached End" << endl << flush;
			if (lastListSwap == i) {             // the current fs instance is out of nodes
				//cout << "Ret -1" << endl << flush;
				return -1;
			}
			else {
				//cout << "Swap" << endl << flush;
				lastListSwap = i;

				fs->threshold += fs->increment;

				if (fs->later == NULL) {
					cout << "WARN - LATER LIST IS NULL" << endl << flush;
				}

				delete fs->now;
				fs->now = fs->later;
				fs->later = new list<Node*>();

				if (fs->now == NULL) {
					cout << "WARN - NOW LIST IS NULL" << endl << flush;
				}
			}
			i--;
		}
		else {
			Node* n = fs->now->front();

			int f = INT_MAX;      //will be the smallest f value of candidate goals
			for (int g = 0; g < fs->numGoals; g++) {
				//if we already have a path to a given goal, don't use
				// it as a candidate for the best heuristic value
				if (fs->paths[g] != NULL) continue;

				int fTemp = n->cost + manhattan(n->coordinate, fs->goals[g]);
				if (fTemp < f) f = fTemp;
			}

			if (f > fs->threshold) {
				fs->later->push_back(n);
			}
			else {
				//cout << "expand node: " << n->coordinate.x << " " << n->coordinate.y << endl << flush;
				//expand children
				int nx = n->coordinate.x;
				int ny = n->coordinate.y;
				int x[] = {nx+1, nx-1, nx, nx  };
				int y[] = {ny,   ny, ny+1, ny-1};
				for (int i = 0; i < 4; i++) {
					//if coordinate pair is valid
					if (x[i] >= 0 && x[i] < fs->mmap->real->cols && y[i] >= 0 && y[i] < fs->mmap->real->cols) {
						Node* child = getNode(fs->mmap->real, x[i], y[i]);



						if (!isBlocked(child)) {

							omp_lock_t * childLock = lockFor(fs->mmap, child->coordinate);

							omp_set_lock(childLock);
							//if the child isn't owned by another process
							if (child->owner == NULL || child->owner->coordinate == fs->start) {
								if (child->cost > n->cost+1) { //and the path we've found to it is best so far
									child->owner = n->owner; //make sure we own it
									omp_unset_lock(childLock); //setting the owner ensures that this instance owns the node
									child->parent = n;
									child->cost = n->cost+1;                                 //set its cost appropriately
									fs->now->push_back(child);
									//cout << "push child: " << child->coordinate.x << " " << child->coordinate.y << endl << flush;
								}
								else {
									omp_unset_lock(childLock);
								}
							}
							else {
								omp_unset_lock(childLock); //no more need for synchronization
								//we've found a path to another process
								for (int g = 0; g < fs->numGoals; g++) {
									if (fs->goals[g] == child->owner->coordinate && fs->paths[g] == NULL) {
										fs->goalsFound++;
										list<Node*>* pathToN = getPath(n);
										pathToN->push_back(child);
										fs->paths[g] = pathToN;
									}
								}
							}

							//TODO: UNLOCK this node
						}
						else {
							//printf("BLOCKED\n");
						}
					}
				}
			}
			fs->now->pop_front();
		}
	}
	return 0; //zero indicates the pathfinding instance isn't out of nodes
}

list<coord>* hlsearch(MetaMap * mmap, coord start, coord goal) {
	cout << "Enter HL Search" << endl << flush;

	coord hlStart = bigToLittle(mmap, start);
	coord hlGoal = bigToLittle(mmap, goal);

	MetaMap fake_mmap = {
		mmap->meta, //real
		mmap->meta, //meta
		mmap->locks,
		1
	};
//fs* buildFS(MetaMap* mmap, int increment, coord start, coord* goals, int numGoals)
	coord goalClaimerGoals[] = {hlStart};
	fs * goalClaimer = buildFS(
		&fake_mmap,
		1, //increment (minimized, to ensure a good path)
		hlGoal, //start for fs
		&(goalClaimerGoals[0]), //goal for fs
		1 //numGoals
	);
	coord goals[] = {hlGoal};
	fs * search = buildFS(
		&fake_mmap,
		1,
		hlStart,
		&(goals[0]),
		1
	);

	cout << "HLSearch: Completed setup. About to run search" << endl << flush;

	int arbitraryIteration = 10;
	while (search->goalsFound == 0) {
		if (fsearch(search, arbitraryIteration) == -1) {
			if (search->goalsFound == 1) {
				break;
			}
			else {
				return NULL; //null retval means no high level path
			}
		}
	}
	//found a high level path if this point is reached

	//FIRST restore mmap's meta component to its original form
	Map * meta = mmap->meta;
	int metaNodes = meta->rows * meta->cols;
	for (int i = 0; i < metaNodes; i++) {
		meta->nodes[i].parent = NULL;
		meta->nodes[i].owner = NULL;
		meta->nodes[i].cost = INT_MAX;
	}

	//finally, construct the path from fs:
	list<coord> * ret = new list<coord>();
	list<Node*>::iterator it = search->paths[0]->begin();
	while (it != search->paths[0]->end()) {
		ret->push_back((*it)->coordinate);
		it++;
	}

	return ret;
}
