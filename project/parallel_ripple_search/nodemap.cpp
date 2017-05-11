#include "stdlib.h"
#include "stdio.h"
#include <stdio.h>
#include "string.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>

#include "nodemap.h"

using namespace std;

bool operator==(const coord& l, const coord& r) {
	return l.x == r.x && l.y ==r.y;
}

double genRand() {
	int percent = rand() % 100;
	return percent == 0 ? 1.0 : ((double) percent) / 100.0;
}
double genRandSign() {
	double val = (rand() % 100) >= 50 ? 1.0 : -1.0;
	return val;
}

Node * getNode(Map& map, int x, int y) {
	return &(map.nodes[indexOf(map, x, y)]);
}
Node * getNode(Map* map, int x, int y) {
	return getNode(*map, x, y);
}

bool isBlocked(Map& map, int x, int y) {
	return getNode(map, x, y)->blocked;
}
bool isBlocked(Node* node) {
	return node->blocked;
}


int indexOf(Map& map, int x, int y) {
	return x * map.cols + y;
}

bool validX(Map& map, int x) {
	return x >=0 && x < map.cols;
}
bool validY(Map& map, int y) {
	return y >=0 && y < map.rows;
}

int littleToBigI(MetaMap* mmap, coord little) {
	coord bigCoord = littleToBig(mmap, little);
	return indexOf(*(mmap->real), bigCoord.x, bigCoord.y);
}

coord littleToBig(MetaMap* mmap, coord little) {
	int x = mmap->factor * little.x;
	int y = mmap->factor * little.y;
	coord ret = {x, y};
	return ret;
}

coord bigToLittle(MetaMap* mmap, coord big) {
	int x = big.x / mmap->factor;
	int y = big.y / mmap->factor;
	coord ret = {x, y};
	return ret;
}
int bigToLittleI(MetaMap* mmap, coord big) {
	coord littleCoord = bigToLittle(mmap, big);
	return indexOf(*(mmap->meta), littleCoord.x, littleCoord.y);
}

omp_lock_t * lockFor(MetaMap* mmap, coord globalCoord) {
	int lockIndex = bigToLittleI(mmap, globalCoord);
	return &(mmap->locks[lockIndex]);
}

/* belongs in main
   void getParams(MapParams& params, int argc, char** argv) {
    int arg = 0;

    params.sidelength = 20;
    params.obsratio = 0.3;

    while (arg < argc) {
        if (strcmp(argv[arg], "--size") == 0) {
            arg++;
            params.sidelength = atoi(argv[arg]);
        } else if (strcmp(argv[arg], "--obsratio") == 0) {
            params.obsratio = atof(argv[++arg]);
        } else if (strcmp(argv[arg], "--change") == 0) {
            params.change = atof(argv[++arg]);
        }
        arg++;
    }

    params.change = 0.4 / log2(params.sidelength);
   }
 */

MetaMap* buildMap(MapParams params, int seed, int maxHighLevelSideLen) {
	srand(seed);

	//make the real map:
	struct Map* map = initializeMap(params);
	struct Bounds* bounds = initializeBounds(params);
	obsFiller(*map, *bounds);

	int sum = 0;
	for (int i = 0; i < params.sidelength; i++) {
		for (int j = 0; j < params.sidelength; j++) {
			sum += isBlocked(*map, i, j);
		}
	}

	double average = ((double)sum) / (params.sidelength*params.sidelength);

	//IMPORTANT: if we don't make the threshold higher than the average, then the concentration
	// in the HL map will be about 0.5, which is rather high (resulting in fewer paths)
	double cutoff = average * 1.2;//max(0.03, average*1.1);

	//make the high level map:
	if (maxHighLevelSideLen > params.sidelength) maxHighLevelSideLen = params.sidelength;
	Map* highLevel = highLevelMap(*map, maxHighLevelSideLen, cutoff);

	//create the locks:
	int numLocks = maxHighLevelSideLen * maxHighLevelSideLen;
	omp_lock_t * locks = new omp_lock_t[numLocks];
	for (int i = 0; i < numLocks; i++) {
		omp_init_lock(&(locks[i]));
	}

	MetaMap* mmap = (MetaMap*)malloc(sizeof(MetaMap));
	mmap->factor = map->cols / highLevel->cols;
	mmap->meta = highLevel;
	mmap->real = map;
	mmap->locks = locks;

	return mmap;
}

Map* initializeMap(MapParams& params) {
	int sidelength = params.sidelength;
	double change = params.change;
	struct Map* map = (Map*) malloc(sizeof(struct Map));
	map->rows = sidelength;
	map->cols = sidelength;
	map->nodes = (Node*) malloc(sidelength*sidelength * sizeof(Node));
	for (int x = 0; x < sidelength; x++) {
		for (int y = 0; y < sidelength; y++) {
			int index = indexOf(*map, x, y);
			coord nc = {x,y};
			map->nodes[index] = Node(nc);
		}
	}
	map->percchange = change;
	return map;
}

//uses pixelation to reduce the map from its current dimensions to a new, smaller dimension
Map* highLevelMap(Map& map, int newSideLen, double occupancyThreshold) {
	int highLevelSquares = newSideLen * newSideLen;
	if (highLevelSquares%4  != 0)
		printf("Warn: highLevelSquares is not a multiple of four");
	Map* hl = (Map*) malloc(sizeof(struct Map));
	hl->nodes = (Node*) malloc(highLevelSquares * sizeof(Node));
	hl->rows = newSideLen;
	hl->cols = hl->rows;
	for (int x = 0; x < hl->cols; x++) {
		for (int y = 0; y < hl->rows; y++) {
			coord nc = {x,y};
			hl->nodes[indexOf(*hl, x, y)] = Node(nc);
		}
	}
	hl->percchange = 0.0; //TODO: BRENT DOESN'T KNOW WHAT THIS MEANS... so it's zero

	int step = map.cols/newSideLen;
	//for each high level coordinate i,j
	for (int i = 0; i < hl->rows; i++) {
		for (int j = 0; j < hl->cols; j++) {
			int sum = 0;
			int xlo = i*step;
			int xhi = xlo + step;
			int ylo = j*step;
			int yhi = ylo + step;
			for (int x = xlo; x < xhi; x++) {
				for (int y = ylo; y < yhi; y++) {
					//TODO: assumption on the following line, that map values are all ones (occupied) or zeroes (empty)
					sum+=getNode(map, x, y)->blocked;
				}
			}
			double fraction = ((double)sum)/(step*step);
			int value = 0;
			if (fraction > occupancyThreshold) {
				value = 1;
			}

			getNode(*hl, i, j)->blocked = value;
			//setMapValue(*hl, i, j, value);
		}
	}
	return hl;
}

struct Bounds* initializeBounds(MapParams& params) {
	int sidelength = params.sidelength;
	double obsratio = params.obsratio;
	struct Bounds* bounds = (Bounds*) malloc(sizeof(struct Bounds));
	bounds->row = 0;
	bounds->col = 0;
	bounds->rowlength = sidelength;
	bounds->collength = sidelength;
	bounds->perc = obsratio;

	return bounds;
}

void obsFiller(struct Map& map, struct Bounds& bounds) {
	//printf("\nRow: %d, Col %d, RowL %d, ColL %d\n", bounds.row, bounds.col, bounds.rowlength, bounds.collength);


	if (bounds.rowlength <= 1 && bounds.collength <= 1) { //base case
		if (bounds.row < map.rows && bounds.col < map.cols && bounds.row >=0 && bounds.col >= 0) {
			//printf("Drawing row: %d, col %d. Percent %f\n", bounds.row, bounds.col, bounds.perc);
			int val = genRand() < bounds.perc ? 1 : 0;
			getNode(map, bounds.row, bounds.col)->blocked = val;
		}
		else {
			return; //bounds are... out of bounds
		}
	} else {
		int rowlength = bounds.rowlength / 2;
		int collength = bounds.collength / 2;

		int rowstart;
		int colstart;
		int rowuse;
		int coluse;

//        printf("Row: %d, Col: %d, Original: %d, Halved: %d\n", bounds.row, bounds.col, bounds.sidelength, newsidelength);

		//struct Bounds* upperLeft = (Bounds*) malloc(sizeof(struct Bounds));

		struct Bounds upperLeft = {
			bounds.row,
			bounds.col,
			rowlength,
			collength,
			bounds.perc + map.percchange * genRand() * genRandSign()
		};

/*
		upperLeft->rowlength = rowlength;
		upperLeft->collength = collength;
		upperLeft->row = bounds.row;
		upperLeft->col = bounds.col;
		upperLeft->perc = bounds.perc + map.percchange * genRand() * genRandSign();
		*/
		obsFiller(map, upperLeft);
		//free(upperLeft);

		if (bounds.collength > 1) {

			colstart = bounds.col + collength;
			if (bounds.collength % 2 != 0) {
//                if (bounds.col % 2 == 0 && collength > 1) {
//                    coluse = collength + 2;
//                } else

				if (bounds.collength > 2) {
					coluse = collength + 1;
				} else {
					coluse = collength;
				}
			} else {
				coluse = collength;
			}

			rowstart = bounds.row;
			rowuse = rowlength;

/*
			struct Bounds* upperRight = (Bounds*) malloc(sizeof(struct Bounds));

			upperRight->rowlength = rowuse;
			upperRight->collength = coluse;
			upperRight->row = rowstart;
			upperRight->col = colstart;
*/
			struct Bounds upperRight = {
				rowstart,
				colstart,
				rowuse,
				coluse,
				bounds.perc+map.percchange*genRand()*genRandSign()
			};

			//upperRight->perc = bounds.perc + map.percchange * genRand() * genRandSign();
			obsFiller(map, upperRight);
			//free(upperRight);

		}

//        if (bounds.rowlength == 24 && bounds.collength == 24) {
//            printf("Ready to start real deal\n");
//        }
		if (bounds.rowlength > 1 && bounds.collength > 1) {

//            if (bounds.rowlength == 24 && bounds.collength == 24) {
//                printf("Ready to start real deal\n");
//            }
			colstart = bounds.col + collength;
			rowstart = bounds.row + rowlength;

			if (bounds.collength % 2 != 0) {
//                if (bounds.col % 2 == 0 && collength > 1) {
//                    coluse = collength + 2;
//                } else

				if (bounds.collength > 2) {
					//printf("Bounds incremented %d\n", collength + 1);
					coluse = collength + 1;
				} else {
					coluse = collength;
				}
			} else {
				coluse = collength;
			}

			if (bounds.rowlength % 2 != 0) {
//                if (bounds.row % 2 == 0 && rowlength > 1) {
//                    rowuse = rowlength + 2;
//                } else

				if (bounds.rowlength > 2) {
					//printf("Bounds incremented %d\n", rowlength + 1);
					rowuse = rowlength + 1;
				} else {
					rowuse = rowlength;
				}
			} else {
				rowuse = rowlength;
			}

/*
			struct Bounds* lowerRight = (Bounds*) malloc(sizeof(struct Bounds));
			lowerRight->rowlength = rowuse;
			lowerRight->collength = coluse;
			lowerRight->row = rowstart;
			lowerRight->col = colstart;
			lowerRight->perc = bounds.perc + map.percchange * genRand() * genRandSign();
			*/

			struct Bounds lowerRight = {
				rowstart,
				colstart,
				rowuse,
				coluse,
				bounds.perc + map.percchange * genRand() * genRandSign()
			};
			obsFiller(map, lowerRight);
			//free(lowerRight);

		}

		if (bounds.rowlength > 1) {
			colstart = bounds.col;
			coluse = collength;

			rowstart = bounds.row + rowlength;
			if (rowlength % 2 != 0) {
//                if (bounds.row % 2 == 0 && rowlength > 1) {
//                    rowuse = rowlength + 2;
//                } else

				if (bounds.rowlength > 2) {
					rowuse = rowlength + 1;
				} else {
					rowuse = rowlength;
				}
			} else {
				rowuse = rowlength;
			}
/*
			struct Bounds* lowerLeft = (Bounds*) malloc(sizeof(struct Bounds));
			lowerLeft->rowlength = rowuse;
			lowerLeft->collength = coluse;
			lowerLeft->perc = bounds.perc + map.percchange * genRand() * genRandSign();
			lowerLeft->row = rowstart;
			lowerLeft->col = colstart;
			*/


			struct Bounds lowerLeft = {
				rowstart,
				colstart,
				rowuse,
				coluse,
				bounds.perc + map.percchange * genRand() * genRandSign()
			};
			obsFiller(map, lowerLeft);
			//free(lowerLeft);
		}
	}
}

void printMap(Map& map) {
	for (int i = 0; i < map.rows; i++) {
		for (int j = 0; j < map.rows; j++) {
			printf("%d ", isBlocked(map, i, j));
		}
		printf("\n");
	}
}

void saveFile(Map& map, char* filename) {

	FILE* fp = fopen(filename, "w+");

	for (int i = 0; i < map.rows; i++) {
		for (int j = 0; j < map.cols; j++) {
			fprintf(fp, "%d ", isBlocked(map, i, j));
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	printf("Closing fp\n");
}

/*
   Map* occupancy(Map& course) {
    Map* occ = (Map*) malloc(sizeof(Map));
    int* map = (int*) malloc(course.rows * course.cols * sizeof(int));
    occ->map = map;
    occ->rows = course.rows;
    occ->cols = course.cols;


    for (int i = 0; i < course.rows; i++) {
        for (int j = 0; j < course.cols; j++) {
            Node* n = getNode(occ, i, j);
            n->blocked = 0;
            //setMapValue(*occ, i, j, 0);
        }
    }

    return occ;
 */
