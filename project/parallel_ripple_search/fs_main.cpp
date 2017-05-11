#include "nodemap.h"
#include "fringesearch.h"

#include <math.h>

#include <stdio.h>
#include <iostream>

using namespace std;

int main() {
	MapParams params = {
		16, //sidelen
		.0, //obsratio
		.0 / log2(16) //change
	};
	MetaMap* mmap = buildMap(
	    params,
	    1, //seed
	    4 //max sidelength for hl
	    );

	coord goals[] = {{15,15}}; //one goal
	fs* search = buildFS(
	    mmap,
	    3, //increment
	    {0,0}, //start
	    goals, //the goals
	    1 //numGoals
	    );

	//necessary for another instance to claim the goal (otherwise it won't be recognized as a goal)
	coord otherGoal[] = {{0,0}};
	fs* other = buildFS(
	    mmap,
	    3,
	    {15,15},
	    otherGoal,
	    1
	    );


	while (search->goalsFound < 1) {
		printf("iterate\n");
		if (fsearch(search, 10) == -1) break;
	}
	//printf("Found Path:\n");
	if (search->goalsFound == 1) {

		cout << "found path!" << endl << flush;
		list<Node*>::iterator it = search->paths[0]->begin();
		while (it != search->paths[0]->end()) {
			Node* n = *it;
			//printf("%d, %d\n", n->coordinate.x, n->coordinate.y);
			cout << "node " << n->coordinate.x << " " << n->coordinate.y << endl << flush;
			it++;
		}

	}
	else {
		cout << "No Path Found" << endl <<flush;
	}
	printf("Done.\n");

}
