#include <iostream>
#include <list>
#include <stdlib.h>
#include <cstddef>
#include "nodemap.h"

using namespace std;

int manhattan(int a, int b, int ai, int bi) {
	return abs(a-ai) + abs(b-bi);
}


int main() {

	Map * map = buildMap();

	//randomly select a start point:
	int startX;
	int startY;
	do {
		startX = randRange(0, gridSize);
		startY = randRange(0, gridSize);
	} while (blocked(startX, startY, grid));
	cout << startX << ", " << startY << endl;

	//randomly select an end point:
	int endX;
	int endY;
	do {
		endX = randRange(0, gridSize);
		endY = randRange(0, gridSize);
	} while (blocked(startX, startY, grid) && endX != startX);
	cout << endX << ", " << endY << endl;

	//now, fringe search:
	int inc = 3;
	int threshold = manhattan(startX, startY, endX, endY);
	list<Node*> * now = new list<Node*>();
	list<Node*> * later = new list<Node*>();
	Node * goalNode;
	bool done = false;

	Node * start = new Node;
	start->x = startX;
	start->y = startY;
	start->parent = NULL;
	start->cost = 0;

	now->push_front( start );
	closed[startX][startY] = true;
	while (now->size() > 0 && !done) {
		threshold += inc;
		for (list<Node*>::iterator it = now->begin(); it != now->end() && !done; ) {
			Node * n = *it;
			int f = n->cost + manhattan(n->x, n->y, endX, endY);
			if (f > threshold) {
				//defer for later
				later->push_back(n);
			}
			else if (n->x == endX && n->y == endY) {
				goalNode = n;
				done = true;
			}
			else {
				//expand children
				int x[] = {n->x+1, n->x-1, n->x, n->x};
				int y[] = {n->y,   n->y, n->y+1, n->y-1};
				for (int i = 0; i < 4; i++) {
					if (validIndex(x[i], gridSize) && validIndex(y[i], gridSize) && !closed[x[i]][y[i]]) {
						Node* child = new Node;
						child->x = x[i];
						child->y = y[i];
						child->cost = n->cost+1;
						child->parent = n;

						cout << "Child of " << child->parent->x << ", " << child->parent->y << " is " << x[i] << ", " << y[i] << endl;

						closed[x[i]][y[i]] = true;
						now->push_back(child);
					}
				}
			}
			it = now->erase(it);
		}
		delete now;
		now = later;
		later = new list<Node*>();
	}

	//print results
	if (!done) {
		cout << "No Path" << endl;
	}
	else {
		cout << "Path Found" << endl;
		list<Node> path = list<Node>();
		Node * temp = goalNode;
		while (true) {
			path.push_front(*temp);
			if (temp->parent == NULL) break;
			else {
				Node * parent = temp->parent;
				cout << "Parent [" << parent->x << "," << parent->y << "]" << endl;
				cout << "Temp   [" << temp->x << "," << temp->y << "]" << endl;
				if (parent->x == temp->x && parent->y == temp->y) return 0;
				temp = parent;
			}
		}
		while (path.size() > 0) {
			Node next = path.front();
			path.pop_front();
			cout << "[" << next.x << "," << next.y << "]" << " ";
		}
		cout << endl;
	}

	return 0;
}
