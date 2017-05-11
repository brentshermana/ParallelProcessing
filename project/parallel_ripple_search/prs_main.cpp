#include "nodemap.h"
#include "fringesearch.h"
#include <stdbool.h>

#include <omp.h>

#include <math.h>

#include <stdio.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    //parse CLI options:

    //first: The side length of the map (side length):
    int mapSideLen = atoi(argv[1]);
    //second: The starting obstacle ratio
    double obsRatio = atof(argv[2]);
    //third: The side length of the high level map
    int hlSideLen = atoi(argv[3]);
    //fourth: The seed
    int seed = atoi(argv[4]);
    //fifth: The number of threads
    int threads = atoi(argv[5]);
    if (threads < 3) {
        cout << "Must assign at least three threads -- the manager, and the two essential cores" << endl << flush;
        return 0;
    }
    int maxThreads = omp_get_max_threads();
    if (maxThreads < threads) {
        cout << "omp only allows " << maxThreads << " Threads" << endl << flush;
        return 0;
    }
    omp_set_num_threads(threads);

    cout << "args:" << endl <<
        "sideLen " << mapSideLen << endl <<
        "ratio " << obsRatio << endl <<
        "hl side len " << hlSideLen << endl <<
        "seed " << seed << endl <<
        "threads " << threads << endl << flush;

    cout << "got args successfully" << endl << flush;

	//singly threaded setup:
    double change = obsRatio / log2(1.0 * mapSideLen);
	MapParams params = {
		mapSideLen,
		obsRatio,
		change     //change
	};
	MetaMap* mmap = buildMap(
	    params,
	    seed,
	    hlSideLen
	);

    cout << "Constructed map" << endl <<flush;

    //start in bottom left corner, end in top right
    coord start = {1,1};
    coord goal = {mapSideLen - 2, mapSideLen - 2};

    coord hlStart = bigToLittle(mmap, start);
    coord hlGoal = bigToLittle(mmap, goal);

    //manually unblock start and goal:
    getNode(mmap->real, start.x, start.y)->blocked = 0;
    getNode(mmap->real, goal.x, goal.y)->blocked = 0;
    getNode(mmap->meta, hlGoal.x, hlGoal.y)->blocked = 0;
    getNode(mmap->meta, hlStart.x, hlStart.y)->blocked = 0;

    cout << "About to HL Search" << endl << flush;

    //Run pathfinding on higher level graph
    list<coord>* hlPath = hlsearch(mmap, start, goal);

    if (hlPath == NULL) {

        printMap(*mmap->real);
        cout << "\n" << "\n";
        printMap(*mmap->meta);

        cout << "High Level search didn't yield a path. Try a different seed" << endl << flush;
        return 0;
    }

    cout << "HL Search Complete" << endl << flush;

    list<coord>::iterator it = hlPath->begin();
    while (it != hlPath->end()) {
        coord current = *it;
        //cout << current.x << " " << current.y << endl << flush;
        it++;
    }

    cout << "Assigning Cores" << endl << flush;

    int cores = threads-1; //slave cores, 0->(threads-2)

    coord coreStartPoints[cores];

    int step = hlPath->size() / cores; // deliberate rounding down

    //assign the start points of each core
    coreStartPoints[cores-1] = littleToBig(mmap, hlPath->back()); //goal core
    for (int c = 0; c < cores-1; c++) { //start core and all others
        coreStartPoints[c] = littleToBig(mmap, hlPath->front());
        for (int s = 0; s < step; s++) {
            hlPath->pop_front();
        }
    }

    for (int c = 0; c < cores; c++) {
        coord crd = coreStartPoints[c];
        getNode(mmap->real, crd.x, crd.y)->blocked = 0; //manually ensure no core is assigned to a blocked grid
        cout << "core " << c << " assigned " << crd.x << " " << crd.y << endl;
    }

    fs** searchInstances = new fs*[cores];

    cout << "Initializing fs instances" << endl << flush;
    for (int c = 0; c < cores; c++) {
        coord* goals;
        int numGoals;
        if (c==0) {
            goals = new coord[1];
            numGoals = 1;
            goals[0] = coreStartPoints[c+1];
        }
        else if (c==cores-1) {
            numGoals = 1;
            goals = new coord[1];
            goals[0] = coreStartPoints[c-1];
        }
        else {
            numGoals = 2;
            goals = new coord[2];
            goals[0] = coreStartPoints[c-1];
            goals[1] = coreStartPoints[c+1];
        }

        searchInstances[c] = buildFS(
            mmap,
            1,
            coreStartPoints[c],
            goals,
            numGoals
        );
    }

    int* masterHalt = new int[cores]; //for master to tell slaves to stop
    int* slaveAck = new int[cores];  //for slaves to signal that they've stopped
    int* totallyDone = new int[cores]; //indicates that a slave has found all neighbors
    int* outOfNodes = new int[cores];

    int master = threads-1; //master core
    for (int i = 0; i < cores; i++) {
        masterHalt[i] = 0;
        slaveAck[i] = 0;
        totallyDone[i] = 0;
        outOfNodes[i] = 0;
    }

    double startTime = omp_get_wtime();
    #pragma omp parallel // this is where the magic happens
    {
        int id = omp_get_thread_num();

        if (id == master) { //master coordinates other cores
            int totallyDoneCount = 0;
            while (totallyDoneCount < cores) {
                //wait for rendevous
                /*
                int rendevous = 0;
                while (!rendevous) {
                    int count = 0;
                    for (int c = 0; c < cores; c++) {
                        int ack, tdone;
                        #pragma omp atomic read
                        ack = slaveAck[c];
                        #pragma omp atomic read
                        tdone = totallyDone[c];

                        if (tdone || ack) count++;
                    }
                    rendevous = (count == cores);
                    if (rendevous) cout << rendevous << endl;
                    cout << "Count: " << count << " Cores " << cores << endl;
                }
                */
                //act on each core
                for (int c = 0; c < cores; c++) {
                    int done;
                    #pragma omp atomic read
                    done = totallyDone[c];
                    if (done) continue;

                    int ack;
                    #pragma omp atomic read
                    ack = slaveAck[c];
                    if (ack) { //slave is waiting to be checked up on
                        fs* cInst = searchInstances[c];
                        int badStatus;
                        #pragma omp atomic read
                        badStatus = outOfNodes[c];
                        if (cInst->goalsFound == cInst->numGoals) {
                            cout << "Master: slave " << c << " successfully finished" << endl << flush;
                            //mark slave as finished
                            totallyDoneCount++;
                            #pragma omp atomic write
                            totallyDone[c] = 1;
                        }
                        else if (badStatus) {
                            cout << "Master: slave " << c << " FAILED" << endl << flush;
                            totallyDoneCount++;
                            #pragma omp atomic write
                            totallyDone[c] = 1;
                        }
                        #pragma omp atomic write
                        masterHalt[c] = 0;
                        #pragma omp atomic write
                        slaveAck[c] = 0;
                    }
                }
            }
        }
        else {
            //cout << "I am slave " << id << endl;
            fs* mySearch = searchInstances[id];
            int done = 0;
            while (!done) {
                //check for master interrupts
                int haltByMaster;
                #pragma omp atomic read
                haltByMaster = masterHalt[id];
                if (haltByMaster) {
                    #pragma omp atomic write
                    slaveAck[id] = 1; //signal that you're waiting
                    continue;
                }

                //check if waiting on master to be attended to
                int waitingOnMaster;
                #pragma omp atomic read
                waitingOnMaster = slaveAck[id];
                if (waitingOnMaster) {
                    continue;
                }

                //cout << "core " << id << " searching" << endl;
                int retStatus = fsearch(mySearch, 2000);
                #pragma omp atomic write
                outOfNodes[id] = retStatus;
                #pragma omp atomic write
                slaveAck[id] = 1;

                //check for validation from master that the core is done
                #pragma omp atomic read
                done = totallyDone[id];
            }

        }
    }

    double totalTime = omp_get_wtime() - startTime;

    cout << "End Parallel Section" << endl << flush;

    cout << "Constructing Master Path" << endl << flush;
    list<coord> * masterList = new list<coord>();

    cout << "Iterating from start node" << endl << flush;
    list<Node*>::iterator nit = searchInstances[0]->paths[0]->begin();
    while (nit != searchInstances[0]->paths[0]->end()) {
        masterList->push_back((*nit)->coordinate);
        nit++;
    }

    for (int i = 1; i < cores-1; i++) {
        cout << "Intermediary Node " << i << endl << flush;
        //get all coordinates leading to i
        Node* bridge = getNode(*mmap->real, masterList->back().x, masterList->back().y);
        list<Node*> * toI = getPath(bridge);
        toI->pop_back();
        while (!toI->empty()) {
            masterList->push_back(toI->back()->coordinate);
            toI->pop_back();
        }

        //get all coordinates owned by i which lead to i+1
        list<Node*> * toNext = searchInstances[i]->paths[1];
        nit = toNext->begin();
        while (nit != toNext->end()) {
            masterList->push_back((*nit)->coordinate);
            nit++;
        }
    }

    cout << "Goal Node" << endl << flush;
    Node* bridge = getNode(*mmap->real, masterList->back().x, masterList->back().y);
    list<Node*> * toGoal = getPath(bridge);
    toGoal->pop_back();
    while (toGoal->size() > 0) {
        coord c = toGoal->back()->coordinate;
        masterList->push_back(c);
        toGoal->pop_back();
    }

    cout << "Master path produced:" << endl;

    list<coord>::iterator mit = masterList->begin();
    coord prev = masterList->front();
    while (mit != masterList->end()) {
        coord c = *mit;
        if (manhattan(prev, c) > 1) {
            cout << "Master Path Integrity Check Fail" << endl;
            return 0;
        }
        prev = c;
        //cout << c.x << " " << c.y << endl;
        mit++;
    }

    cout << "Cores: " << threads << endl;
    cout << "Time: " << totalTime << endl;

}
