#include "enumPerfMatchings.h"
#include "../matchmaker/matchmaker.h"
#include "time.h"
#include <omp.h>


/************************************************
 *  Helper functions
 ***********************************************/

// Return a psuedo-random integer in {lo, ..., hi} (so both lo and hi inclusive)
int randIntBetween(int lo, int hi) {
    return lo + rand() % (hi - lo + 1);
    // rand() Returns a pseudo-random integer between 0 and RAND_MAX.
}



/************************************************
 *  Component parity functions
 ***********************************************/

// Based on nauty
// Returns the number of components of the graph g of
// order numVertices
int numComponents(graph *g, int maxm, int numVertices) {
    /* Find number of components of undirected graph */
    int* myqueue = (int*)malloc(numVertices*sizeof(int));
    if (myqueue == NULL) {
        fprintf(stderr,"Error allocating memory for myqueue.\n");
        exit(-1);
    }
    bool* myvisited = (bool*)calloc(numVertices, sizeof(bool));
    if (myvisited == NULL) {
        fprintf(stderr,"Error allocating memory for myvisited.\n");
        exit(-1);
    }

    int nc = 0;

    for (int v = 0; v < numVertices; ++v) {
        if (myvisited[v]) continue;
        ++nc;
        myqueue[0] = v;

        myvisited[v] = true;

        int head = 0;
        int tail = 1;
        while (head < tail) {
            int w = myqueue[head++];
            FOREACH(nbr, GRAPHROW(g, w, maxm), maxm) {
                if (myvisited[nbr] == false) {
                    myvisited[nbr] = true;
                    myqueue[tail++] = nbr;
                }
            }
        }
    }

    free(myqueue);
    free(myvisited);

    return nc;
}

// Returns the parity of the number of components when deleting matching from the graph g
// of order numVertices. The matching contains the vertices numbered according to the
// MatchMaker library, which can be renumbered with inverseRenumPart1 and inverseRenumPart2.
int getParity(graph *g, int numVertices, int maxm, int *matching, int *inverseRenumPart1, int *inverseRenumPart2) {
    int parity;

    #pragma omp critical(useGraph)
    {
        for(int j=0; j<numVertices/2; j++) {
            REMOVEONEEDGE(g, inverseRenumPart2[j], inverseRenumPart1[matching[j]], maxm);
        }
        parity = numComponents(g, maxm, numVertices) % 2;
        for(int j=0; j<numVertices/2; j++) {
            ADDONEEDGE(g, inverseRenumPart2[j], inverseRenumPart1[matching[j]], maxm);
        }
    }

    return parity;
}

// Returns true if there is 1 component when deleting matching from the graph g of order
// numVertices; otherwise returns false. The matching contains the vertices numbered according to the
// MatchMaker library, which can be renumbered with inverseRenumPart1 and inverseRenumPart2.
bool is2FactorHamiltonian(graph *g, int numVertices, int maxm, int *matching, int *inverseRenumPart1,
                          int *inverseRenumPart2) {
    bool still2FactorHamiltonian = true;

    #pragma omp critical(useGraph)
    {
        for(int j=0; j<numVertices/2; j++) {
            REMOVEONEEDGE(g, inverseRenumPart2[j], inverseRenumPart1[matching[j]], maxm);
        }
        if (numComponents(g, maxm, numVertices) != 1) {
            still2FactorHamiltonian = false;
        }
        for(int j=0; j<numVertices/2; j++) {
            ADDONEEDGE(g, inverseRenumPart2[j], inverseRenumPart1[matching[j]], maxm);
        }
    }
    return still2FactorHamiltonian;
}



/************************************************
 *  Printing functions
 ***********************************************/

// Prints the perfect matching perfMatching, numbered numMatching of a graph of order numVertices.
void printPerfMatching(unsigned long long* numMatching, int* perfMatching, int numVertices) {
    printf("%llu: ", *numMatching);
    for(int j=0; j<numVertices/2; j++) {
        printf("(%d,%d) ", perfMatching[2*j], perfMatching[2*j+1]);
    }
    printf("\n");
}


// Prints the 2-factors when deleting the perfect matching perfMatching, numbered numMatching
// of the graph g of order numVertices.
void print2Factor(graph *g, unsigned long long* numMatching, int* perfMatching, int numVertices, int maxm) {
    printf("%llu: ", *numMatching);

    for(int j=0; j<numVertices/2; j++) {
        REMOVEONEEDGE(g, perfMatching[2*j], perfMatching[2*j+1], maxm);
    }

    int* myqueue = (int*)malloc(numVertices*sizeof(int));
    if (myqueue == NULL) {
        fprintf(stderr,"Error allocating memory for myqueue.\n");
        exit(-1);
    }
    bool* myvisited = (bool*)calloc(numVertices, sizeof(bool));
    if (myvisited == NULL) {
        fprintf(stderr,"Error allocating memory for myvisited.\n");
        exit(-1);
    }

    for (int v = 0; v < numVertices; ++v) {
        if (myvisited[v]) continue;

        printf("( ");
        myqueue[0] = v;
        myvisited[v] = true;
        int head = 0;
        int tail = 1;
        while (head < tail) {
            int w = myqueue[head++];
            printf("%d ", w);

            FOREACH(nbr, GRAPHROW(g, w, maxm), maxm) {
                if (myvisited[nbr] == false) {
                    myvisited[nbr] = true;
                    myqueue[tail++] = nbr;
                }
            }
        }
        printf(") ");
    }

    free(myqueue);
    free(myvisited);

    printf("\n");
}


// Prints the sizes of the 2-factors when deleting the perfect matching perfMatching, numbered numMatching
// of the graph g of order numVertices.
void print2FactorSize(graph *g, unsigned long long* numMatching, int* perfMatching, int numVertices, int maxm) {
    printf("%llu: ", *numMatching);

    for(int j=0; j<numVertices/2; j++) {
        REMOVEONEEDGE(g, perfMatching[2*j], perfMatching[2*j+1], maxm);
    }

    int* myqueue = (int*)malloc(numVertices*sizeof(int));
    if (myqueue == NULL) {
        fprintf(stderr,"Error allocating memory for myqueue.\n");
        exit(-1);
    }
    bool* myvisited = (bool*)calloc(numVertices, sizeof(bool));
    if (myvisited == NULL) {
        fprintf(stderr,"Error allocating memory for myvisited.\n");
        exit(-1);
    }

    printf("( ");

    for (int v = 0; v < numVertices; ++v) {
        if (myvisited[v]) continue;

        int ctrVertices = 0;
        myqueue[0] = v;
        myvisited[v] = true;
        int head = 0;
        int tail = 1;
        while (head < tail) {
            int w = myqueue[head++];
            ctrVertices++;

            FOREACH(nbr, GRAPHROW(g, w, maxm), maxm) {
                if (myvisited[nbr] == false) {
                    myvisited[nbr] = true;
                    myqueue[tail++] = nbr;
                }
            }
        }
        printf("%d ", ctrVertices);

    }

    printf(")\n");

    free(myqueue);
    free(myvisited);
}


/************************************************
 *  Perfect matching enumeration functions
 ***********************************************/

bool diffParityFoundRec;

// Recursively generate perfect matchings.
void enumMatchingsRec(graph *g, int numVertices, unsigned long long * matching_count, int *parity,
                 int maxm, bool* done, bool hamiltonian, bool* visitedVertices,
                 int* perfMatching, int numEdgesMatching, bool* inPartition1,
                 bool printPerfMatchings, bool print2Factors, bool print2FactorSizes) {

    if (diffParityFoundRec) {
        return;
    }

    bool quit = false;
    #pragma omp critical(end)
    {
        if (*done) { quit = true; }
    }
    if (quit) return;

    int currVertex = 0;
    for (; currVertex < numVertices; currVertex++) {
        if (!visitedVertices[currVertex] && inPartition1[currVertex]) break;
    }

    if (currVertex == numVertices) {
        int sameParity = 1;
        (*matching_count)++;

        if (printPerfMatchings) {
            printPerfMatching(matching_count, perfMatching, numVertices);
        } else if (print2Factors) {
            print2Factor(g, matching_count, perfMatching, numVertices, maxm);
        } else if (print2FactorSizes) {
            print2FactorSize(g, matching_count, perfMatching, numVertices, maxm);
        }

        #pragma omp critical(useGraph)
        {
            for(int j=0; j<numVertices/2; j++) {
                REMOVEONEEDGE(g, perfMatching[2*j], perfMatching[2*j+1], maxm);
            }
            if (hamiltonian) {
                if (numComponents(g, maxm, numVertices) != 1) {
                    sameParity = 0;
                }
            } else {
                int newParity = numComponents(g, maxm, numVertices) % 2;
                if (*parity == -1) {
                    // *parity == -1 if this is the first matching found and the parity
                    // has not been determined yet. So, we set *parity = newParity.
                    *parity = newParity;
                } else if (newParity != *parity) {
                    sameParity = 0;
                }
            }
            for(int j=0; j<numVertices/2; j++) {
                ADDONEEDGE(g, perfMatching[2*j], perfMatching[2*j+1], maxm);
            }
        }
        if (sameParity == 0) {
            diffParityFoundRec = true;
        }
        return;
    }

    visitedVertices[currVertex] = true;
    set* nbrs = malloc(maxm*sizeof(set));
    #pragma omp critical(useGraph)
    {
        memcpy(nbrs,GRAPHROW(g, currVertex, maxm),sizeof(set)*maxm);
    }
    FOREACH(neigh, nbrs, maxm) {
        if (!visitedVertices[neigh]) {
            perfMatching[numEdgesMatching*2] = currVertex;
            perfMatching[numEdgesMatching*2+1] = neigh;
            visitedVertices[neigh] = true;
            enumMatchingsRec(g, numVertices, matching_count, parity, maxm,
                                        done, hamiltonian, visitedVertices, perfMatching, numEdgesMatching+1,
                                        inPartition1, printPerfMatchings, print2Factors, print2FactorSizes);
            visitedVertices[neigh] = false;
        }
    }
    visitedVertices[currVertex] = false;
    free(nbrs);
}

// Function to call to enumerate perfect matchings of graph g, till one is found leading to a number of 2-factors
// of different parity or not equal to 1 (if hamiltonian == true).
// Returns 1 if all enumerated perfect metchings lead to the same parity or one 2-factor (if hamiltonian == true).
// Returns 0 if not all enumerated perfect metchings lead to the same parity or one 2-factor (if hamiltonian == true)
// or another thread finished execution.
int enumMatchings(graph *g, int * col_ids, int * col_ptrs, int numVertices, unsigned long long * matching_count,
                int *inverseRenumPart1, int *inverseRenumPart2, int maxm, bool genAll,
                bool* done, bool hamiltonian, bool* inPartition1, bool printPerfMatchings, bool print2Factors,
                bool print2FactorSizes) {
    if (genAll) {
        diffParityFoundRec = false;
        bool* visitedVertices = (bool*)calloc(numVertices, sizeof(bool));
        if (visitedVertices == NULL) {
            fprintf(stderr,"Error allocating memory for visitedVertices.\n");
            exit(-1);
        }
        int* perfMatching = (int*)malloc(2*numVertices * sizeof(int));
        if (perfMatching == NULL) {
            fprintf(stderr,"Error allocating memory for perfMatching.\n");
            exit(-1);
        }

        *matching_count = 0;
        int parity = -1;

        enumMatchingsRec(g, numVertices, matching_count, &parity, maxm, done,
                                             hamiltonian, visitedVertices, perfMatching, 0, inPartition1,
                                             printPerfMatchings, print2Factors, print2FactorSizes);

        free(visitedVertices);
        free(perfMatching);
        return diffParityFoundRec ? 0 : 1;
    }


    int* match = (int*)malloc(numVertices/2 * sizeof(int));
    if (match == NULL) {
        fprintf(stderr,"Error allocating memory for match.\n");
        exit(-1);
    }
    int* row_match = (int*)malloc(numVertices/2 * sizeof(int));
    if (row_match == NULL) {
        fprintf(stderr,"Error allocating memory for row_match.\n");
        exit(-1);
    }

    srand(time(NULL));

    int parity = -1;
    bool quit = false;
    while (!quit) {
        for(int i=0; i<numVertices/2; i++)
            row_match[i] = -1;

        // matching id = 1 -> DFS
        // cheap id = 5 -> randomized Karp-Sipser
        matching(col_ptrs, col_ids, match, row_match, numVertices/2, numVertices/2,
         1, 3, 1.0);
        (*matching_count)++;

        if (hamiltonian) {
            if (!is2FactorHamiltonian(g, numVertices, maxm, row_match, inverseRenumPart1, inverseRenumPart2)) {
                free(match);
                free(row_match);
                return 0;
            }
        } else {
            int newParity = getParity(g, numVertices, maxm, row_match, inverseRenumPart1, inverseRenumPart2);
            if (parity == -1)
                // *parity == -1 if this is the first matching found and the parity
                // has not been determined yet. So, we set *parity = newParity.
                parity = newParity;
            else if (newParity != parity) {
                free(match);
                free(row_match);
                return 0;
            }
        }

        #pragma omp critical(end)
        {
            if (*done) { quit = true; }
        }
    }
    free(match);
    free(row_match);
    return 0;
}
