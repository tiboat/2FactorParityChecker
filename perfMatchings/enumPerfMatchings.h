
#ifndef ENUMPERFMATCHINGS_H
#define ENUMPERFMATCHINGS_H

#include <stdbool.h>
#include "../nauty/nauty.h"

#define REMOVEONEEDGE(g,i,j,macromaxm)\
DELELEMENT(GRAPHROW(g,i,macromaxm),j); DELELEMENT(GRAPHROW(g,j,macromaxm),i)

//  Macro's for nauty representation
#define FOREACH(element,nautySet,macromaxm)\
for(int element = nextelement((nautySet),macromaxm,-1); (element) >= 0;\
(element) = nextelement((nautySet),macromaxm,(element)))

// Function to call to enumerate perfect matchings of graph g, till one is found leading to a number of 2-factors
// of different parity or not equal to 1 (if hamiltonian == true).
// Returns 1 if all enumerated perfect metchings lead to the same parity or one 2-factor (if hamiltonian == true).
// Returns 0 if not all enumerated perfect metchings lead to the same parity or one 2-factor (if hamiltonian == true)
// or another thread finished execution.
int enumMatchings(graph *g, int * col_ids, int * col_ptrs, int n, unsigned long long * matching_count,
                int *inverseRenumPart1, int *inverseRenumPart2, int maxm, bool genAll,
                bool *done, bool hamiltonian, bool* inPartition1, bool printPerfMatchings, bool print2Factors,
                bool print2FactorSizes);

#endif //ENUMPERFMATCHINGS_H
