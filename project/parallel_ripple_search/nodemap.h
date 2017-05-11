#include <limits.h>
#include <cstddef>
#include <omp.h>

#ifndef NODEMAP_H
#define NODEMAP_H

struct coord {
	int x;
	int y;
};
bool operator==(const coord& l, const coord& r);

struct MapParams {
	int sidelength;
	double obsratio;
	double change;
};

struct Node {
	Node(coord _coordinate) {
		coordinate = _coordinate;
		//default values:
		blocked = 0;
		cost = INT_MAX;
		parent = NULL;
		owner = NULL;
	}
	coord coordinate;
	int blocked;
	int cost;
	Node* parent;
	Node* owner;
};

struct Map {
	Node* nodes;
	int rows;
	int cols;
	double percchange;
};

struct MetaMap {
	Map* real;
	Map* meta;
	omp_lock_t* locks;
	int factor; // >=1
};

struct Bounds {
	int row;
	int col;
	int rowlength;
	int collength;
	double perc;
};

omp_lock_t * lockFor(MetaMap* mmap, coord globalCoord);

int littleToBigI(MetaMap* mmap, coord little);
coord littleToBig(MetaMap* mmap, coord little);
int bigToLittleI(MetaMap* mmap, coord big);
coord bigToLittle(MetaMap* mmap, coord big);

bool validX(Map& map, int x);
bool validY(Map& map, int y);

int indexOf(Map& map, int x, int y);

Node * getNode(Map& map, int x, int y);
Node * getNode(Map* map, int x, int y);

bool isBlocked(Map& map, int x, int y);
bool isBlocked(Node* node);

Map* initializeMap(MapParams&);
struct Bounds* initializeBounds(MapParams&);

void obsFiller(struct Map& map, struct Bounds& bounds);
void printMap(Map& map);
void saveFile(Map& map, char* filename);

MetaMap* buildMap(MapParams params, int seed, int maxHighLevelSideLen);
Map* highLevelMap(Map& map, int newSideLen, double occupancyThreshold);
Map* occupancy(Map&);

#endif
