#include <stdio.h>
#include <getopt.h>
#include <omp.h>

#include "nauty/gtools.h"
#include "nauty/gutils.h"
#include "nauty/nauconnect.h"
#include "perfMatchings/enumPerfMatchings.h"


/************************************************
 *  --- 2-factor parity checker ---
 ***********************************************/

#define USAGE \
"\nUsage: ./2FactorParityChecker [-l# -n# -N#] [-H -s] [-p | -f | -F] [-h]\n\n"

#define HELPTEXT \
"Checks if cubic biparite graphs are pseudo 2-factor isomorphic or 2-factor hamiltonian.\n\
It reads these cubic graphs in graph6 or sparse6 format from stdin. Input graphs\n\
needs to be cubic, but need not be bipartite; non-bipartite graphs are filtered out.\n\
Output is sent to stdout; error messages are sent to stderr. This \n\
implementation uses parallelisation using openmp, so make sure it is \n\
installed.\n\
\n\
Underneath are the optional arguments.\n\
\n\
    -c, --checkConnectedness\n\
        first check whether the graphs are essentially 4-edge-\n\
        connected. Significantly increases the runtime, so only\n\
        enable this when needed.\n\
    -h, --help\n\
        print this help message\n\
    -H, --hamiltonian\n\
        check if the graphs are 2-factor hamiltonian. When -H\n\
        is not set, it checks if the graph is pseudo 2-factor\n\
        isomorphic.\n\
    -l#, --minLineNumber=#\n\
        only starts testing from the #'th graph received\n\
        from stdin\n\
    -n#, --minOrder=#\n\
        only check the graphs in the input that have at least\n\
        # vertices\n\
    -N#, --maxOrder=#\n\
        only check the graphs in the input that have at most\n\
        # vertices. In case you check graphs with more than\n\
        10 000 vertices, then you need to set this argument!\n\
    -s, --stopFound\n\
        terminate the program as soon as a pseudo 2-factor\n\
        isomorphic / 2-factor hamiltonian graph is found\n\
    -p, --printPerfMatchings\n\
        prints the perfect matchings of the graph\n\
    -f, --print2Factors\n\
        prints the 2-factors of the graph\n\
    -F, --print2FactorSizes\n\
        prints the sizes of the 2-factors of the graph\n"


/************************************************
 *  Global variables
 ***********************************************/

bool *inPartition1;
int maxm;
bool hamiltonian = false;
bool printPerfMatchings = false;
bool print2Factors = false;
bool print2FactorSizes = false;

#define NUMIDS 4


// Assumes graph is 3-edge-connected.
// Returns true if the given graph g is essentially 4-edge-connected; otherwise returns false
// and prints a nontrivial 3-edge cut.
bool isEssentially4EdgeConnected(graph *g, int numVertices) {
    for(int v11 = 0; v11 < numVertices; v11++) {
        FOREACH(v12, GRAPHROW(g, v11, maxm),maxm) {
            if (v12 > v11) {
                REMOVEONEEDGE(g, v11, v12, maxm);
                for(int v21 = 0; v21 < numVertices; v21++) {
                    FOREACH(v22, GRAPHROW(g, v21, maxm),maxm) {
                        if (v22 > v21) {
                            // Check if {v11,v12} and {v21,v22} are not the same edge
                            if (v11 == v21 && v12 == v22)
                                continue;
                            REMOVEONEEDGE(g, v21, v22, maxm);
                            for(int v31 = 0; v31 < numVertices; v31++) {
                                // Check if v31 is not incident to both {v11,v12} and {v21,v22},
                                // because otherwise this would be a trivial cut
                                // are not all connected to the same vertex.
                                if (v31 == v11 && v31 == v21 ||
                                    v31 == v11 && v31 == v22 ||
                                    v31 == v12 && v31 == v21 ||
                                    v31 == v12 && v31 == v22)
                                    continue;
                                FOREACH(v32, GRAPHROW(g, v31, maxm),maxm) {
                                    if (v32 > v31) {
                                        // Check if v33 is not incident to both {v11,v12} and {v21,v22},
                                        // because otherwise this would be a trivial cut
                                        // are not all connected to the same vertex.
                                        if (v32 == v11 && v32 == v21 ||
                                            v32 == v11 && v32 == v22 ||
                                            v32 == v12 && v32 == v21 ||
                                            v32 == v12 && v32 == v22)
                                            continue;
                                        REMOVEONEEDGE(g, v31, v32, maxm);
                                        if (!isthisconnected(g, maxm, numVertices, 1, false)) {
                                            printf("{%d,%d} {%d,%d} {%d,%d} is a non-trivial edge-cut\n",
                                                   v11, v12, v21, v22, v31, v32);
                                            return false;
                                        }
                                        ADDONEEDGE(g, v31, v32, maxm);
                                    }
                                }
                            }
                        }
                        ADDONEEDGE(g, v21, v22, maxm);
                    }
                }
                ADDONEEDGE(g, v11, v12, maxm);
            }
        }
    }
    return true;
}


// Assumes the graph g of order numVertices is bipartite. Sets inPartition1[v] = true if
// vertex v is in the same partite set as vertex 0; otherwise sets inPartition1[v] = false.
void setPartition1(graph *g, int numVertices) {
    // Always set the partition of g as the one with vertex 0
    inPartition1[0] = true;

    bool* vis = (bool*)calloc(numVertices, sizeof(bool));
    if (vis == NULL) {
        fprintf(stderr,"Error allocating memory for vis.\n");
        exit(-1);
    }
    int* q = (int*)malloc(numVertices * sizeof(int));
    if (q == NULL) {
        fprintf(stderr,"Error allocating memory for q.\n");
        exit(-1);
    }
    int front = 0, back = 0;
    vis[0] = true;
    q[back++] = 0;
    while (front < back) {
        int currV = q[front++];
        FOREACH(nbr, GRAPHROW(g, currV, maxm),maxm) {
            if (!vis[nbr]) {
                inPartition1[nbr] = inPartition1[currV] ? false : true;
                vis[nbr] = true;
                q[back++] = nbr;
            }
        }
    }
    free(vis);
    free(q);
}

// Returns true if cubic bipartite graph g of order numVertices
// is pseudo 2-factor isomorphic or 2-factor hamiltonian (if
// global variable hamiltionian == true); otherwise returns false.
bool haveSameParity(graph *g, int numVertices) {
    int* col_ptrs = (int*)malloc(numVertices * sizeof(int));
    if (col_ptrs == NULL) {
        fprintf(stderr,"Error allocating memory for col_ptrs.\n");
        exit(-1);
    }
    int r = 3; // degree = 3
    for (int i = 0; i < numVertices; i++) {
        col_ptrs[i] = i*r;
    }

    // Do renumbering
    int* renumbering = (int*)malloc(numVertices * sizeof(int));
    if (renumbering == NULL) {
        fprintf(stderr,"Error allocating memory for renumbering.\n");
        exit(-1);
    }
    int* inverseRenumPart1 = (int*)malloc(numVertices/2 * sizeof(int));
    if (inverseRenumPart1 == NULL) {
        fprintf(stderr,"Error allocating memory for inverseRenumPart1.\n");
        exit(-1);
    }
    int* inverseRenumPart2 = (int*)malloc(numVertices/2 * sizeof(int));
    if (inverseRenumPart2 == NULL) {
        fprintf(stderr,"Error allocating memory for inverseRenumPart2.\n");
        exit(-1);
    }
    int ctrPart1 = 0;
    int ctrPart2 = 0;
    for (int i = 0; i < numVertices; i++) {
        if (inPartition1[i]) {
            renumbering[i] = ctrPart1;
            inverseRenumPart1[ctrPart1] = i;
            ctrPart1++;
        } else {
            renumbering[i] = ctrPart2;
            inverseRenumPart2[ctrPart2] = i;
            ctrPart2++;
        }
    }

    // Make col_ids
    int* col_ids = (int*)malloc((numVertices*r/2) * sizeof(int));
    if (col_ids == NULL) {
        fprintf(stderr,"Error allocating memory for col_ids.\n");
        exit(-1);
    }

    for (int v = 0; v < numVertices; v++) {
        if (inPartition1[v]) {
            int ctrNbrs = 0;
            FOREACH(nbr, GRAPHROW(g, v, maxm),maxm) {
                col_ids[renumbering[v]*r+ctrNbrs] = renumbering[nbr];
                ctrNbrs++;
            }
        }
    }


    bool done = false;
    int allSameParity;
    unsigned long long ctrPerfMatchings = 0;

    if (printPerfMatchings || print2Factors || print2FactorSizes) {
        // In case perfect matchings, the 2-factors or 2-factor sizes need to
        // be printed, only oen thread is usd for enumeration of the
        // perfect matchings.
        allSameParity = enumMatchings(g, col_ids, col_ptrs, numVertices, &ctrPerfMatchings,
                                       inverseRenumPart1, inverseRenumPart2, maxm,
                                       true, &done,
                                       hamiltonian, inPartition1, printPerfMatchings, print2Factors,
                                       print2FactorSizes);
    } else {
        unsigned long long matching_count[NUMIDS] = { 0 };

        int allSameParityArr[NUMIDS] = { -2 };
        int idxResultParity = -1;

        #pragma omp parallel shared(done)
        {
            #pragma omp for
            for (int i = 0; i < NUMIDS; i++) {
                if (done) continue;

                // When i == 0, we will try to generate all perfect matchings.
                // When i != 0 we will generate perfect matchings until a different parity is found
                // and ohter threads are stopped as well.

                // Perform the actual computation
                allSameParityArr[i] = enumMatchings(g, col_ids, col_ptrs, numVertices, &matching_count[i],
                                           inverseRenumPart1, inverseRenumPart2, maxm,
                                           i == 0, &done,
                                           hamiltonian, inPartition1,
                                           false, false, false);

                #pragma omp critical(end)
                {
                    if (!done) {
                        idxResultParity = i;
                        done = true; // Set the shared flag
                    }
                }
            }
        }
        allSameParity = allSameParityArr[idxResultParity];
        ctrPerfMatchings = matching_count[idxResultParity];
    }

    free(col_ids);
    free(col_ptrs);
    free(renumbering);
    free(inverseRenumPart1);
    free(inverseRenumPart2);

    if (allSameParity == 0) {
        printf("Number of matchings until different parity found = %llu\n", ctrPerfMatchings);
    }
    if (allSameParity == 1) {
        if (print2Factors || print2FactorSizes) {
            printf("Total number of 2-factors = %llu\n", ctrPerfMatchings);
        } else {
            printf("Total number of matchings = %llu\n", ctrPerfMatchings);
        }
    }

    return allSameParity == 1 ? true : false;
}



// main function
int main(int argc, char ** argv) {

    int minOrder = 1;
    int minLineNumber = 1;
    int maxn = 10000;

    bool stopFoundPseudo2Factor = false;
    bool checkConnectedness = false;

    // Parsing arguments based on: faultCost (https://github.com/JarneRenders/faultCost)

    int opt;
    while (true) {
        int option_index = 0;
        static struct option long_options[] =
        {
            {"minOrder", required_argument, NULL, 'n'},
            {"maxOrder", required_argument, NULL, 'N'},
            {"minLineNumber", required_argument, NULL, 'l'},
            {"hamiltonian", no_argument, NULL, 'H'},
            {"stopFound", no_argument, NULL, 's'},
            {"checkConnectedness", no_argument, NULL, 'c'},
            {"help", no_argument, NULL, 'h'},
            {"printPerfMatchings", no_argument, NULL, 'p'},
            {"print2Factors", no_argument, NULL, 'f'},
            {"print2FactorSizes", no_argument, NULL, 'F'}
        };
        opt = getopt_long(argc, argv, "cfFhHl:n:N:ps", long_options, &option_index);
        if (opt == -1) break;
        switch(opt) {
            case 'l':
                minLineNumber = (int) strtol(optarg, NULL, 10);
                break;
            case 'n':
                minOrder = (int) strtol(optarg, NULL, 10);
                break;
            case 'N':
                maxn = (int) strtol(optarg, NULL, 10);
                break;
            case 'c':
                checkConnectedness = true;
                break;
            case 'H':
                hamiltonian = true;
                break;
            case 's':
                stopFoundPseudo2Factor = true;
                break;
            case 'p':
                printPerfMatchings = true;
                break;
            case 'f':
                print2Factors = true;
                break;
            case 'F':
                print2FactorSizes = true;
                break;
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case '?':
                fprintf(stderr,"Error: Unknown option: %c\n", optopt);
                fprintf(stderr, "%s", USAGE);
                fprintf(stderr,"Use %s --help for more detailed instructions.\n", argv[0]);
                return 1;
            }
    }

    if (printPerfMatchings && print2Factors ||
        printPerfMatchings && print2FactorSizes ||
        print2Factors && print2FactorSizes) {
        fprintf(stderr, "Only allowed to print either the perfect matchings or the 2-factors or the 2-factor's sizes.\n");
        return 1;
    }

    int numVertices;

    //  Start looping over lines of stdin.
    char *graphString = NULL;
    size_t size;
    int ctrAllGraphs = 0;
    int ctrTested = 0;
    int ctrAllSameParity = 0;

    maxm = SETWORDSNEEDED(maxn);

    inPartition1 = (bool*)malloc(maxn * sizeof(bool));
    if (inPartition1 == NULL) {
        fprintf(stderr,"Memory allocation inPartition1 failed\n");
        return -1;
    }

    DYNALLSTAT(graph,g,g_sz);

    while(getline(&graphString, &size, stdin) != -1) {
        ctrAllGraphs++;
        if (ctrAllGraphs < minLineNumber) continue;

        numVertices = graphsize(graphString);

        if (numVertices < minOrder) continue;
        if (numVertices > maxn) continue;

        ctrTested++;

        maxm = SETWORDSNEEDED(numVertices);
        nauty_check(WORDSIZE,maxm,numVertices,NAUTYVERSIONID);

        printf("n=%d\n", numVertices);
        printf("line num=%d\n", ctrAllGraphs);

        DYNALLOC2(graph,g,g_sz,maxm,numVertices,"malloc");

        stringtograph(graphString,g,maxm);

        int theGirth = girth(g, maxm, numVertices);
        if (!hamiltonian && theGirth < 6) continue;

        bool isBipartite = isbipartite(g, maxm, numVertices);
        if (!isBipartite) {
            printf("not bipartite\n");
            continue;
        }

        if (checkConnectedness) {
            bool is3EdgeConnected = isthisedgeconnected(g, maxm, numVertices, 3);
            if (!is3EdgeConnected) {
                printf("not 3-edge connected\n");
                continue;
            }
            bool essentially4EdgeConnected = isEssentially4EdgeConnected(g, numVertices);
            if (!essentially4EdgeConnected) {
                printf("not essentially 4-edge connected\n");
                continue;
            }
        }

        setPartition1(g, numVertices);

        bool allSameParity = haveSameParity(g, numVertices);

        if (allSameParity) {
            ctrAllSameParity++;
            printf("girth = %d\n", theGirth);
            if (hamiltonian) {
                printf("All 2-factors are hamiltonian.\n\n");
            } else {
                printf("All 2-factors have the same parity.\n\n");
            }
            if (stopFoundPseudo2Factor) break;
        }
    }

    printf("Tested %d graph(s)\n", ctrTested);
    if (hamiltonian) {
        printf("Found %d graph(s) where all 2-factors are hamiltonian cycles\n", ctrAllSameParity);
    } else {
        printf("Found %d graph(s) where all 2-factors have the same parity\n", ctrAllSameParity);
    }

    free(graphString);
    free(inPartition1);

    return 0;
}






