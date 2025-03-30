# 2-factor parity checker

This repository contains code and data related to the manuscript "TODO". This code allows to check cubic bipartite graphs on being pseudo 2-factor isomorphic or 2-factor hamiltonian.

## Data

In the directory `data/grayGraph` one can find data related to the Gray graph, which is pseudo 2-factor isomorphic.

- `adjacency_list.txt` contains the adjacency list of the Gray graph.
- `gray.g6` contains the Gray graph in graph6 format. See [here](https://users.cecs.anu.edu.au/~bdm/data/formats.txt) for more information on the graph6 format.
- `2-factors.txt` lists the 2-factors of the Gray graph.
- `2-factor_sizes.txt` lists the different sizes of the 2-factors of the Gray graph.
- `perfect_matchings.txt` lists the perfect matchings of the Gray graph.

## Code

### Compilation
Compile with the make command.

```
make
```

### Execution
This helptext can be found by executing `./2FactorParityChecker -h`.

Usage: `./2FactorParityChecker [-l# -n# -N#] [-H -s] [-p | -f | -F] [-h]`


Checks if cubic biparite graphs are pseudo 2-factor isomorphic or 2-factor hamiltonian.
It reads these cubic graphs in graph6 or sparse6 format from stdin. Input graphs
needs to be cubic, but need not be bipartite; non-bipartite graphs are filtered out.
Output is sent to stdout; error messages are sent to stderr. This
implementation uses parallelisation using openmp, so make sure it is
installed.

Underneath are the optional arguments.

    -c, --checkConnectedness
        first check whether the graphs are essentially 4-edge-
        connected. Significantly increases the runtime, so only
        enable this when needed.
    -h, --help
        print this help message
    -H, --hamiltonian
        check if the graphs are 2-factor hamiltonian. When -H
        is not set, it checks if the graph is pseudo 2-factor
        isomorphic.
    -l#, --minLineNumber=#
        only starts testing from the #'th graph received
        from stdin
    -n#, --minOrder=#
        only check the graphs in the input that have at least
        # vertices
    -N#, --maxOrder=#
        only check the graphs in the input that have at most
        # vertices. In case you check graphs with more than
        10 000 vertices, then you need to set this argument!
    -s, --stopFound
        terminate the program as soon as a pseudo 2-factor
        isomorphic / 2-factor hamiltonian graph is found
    -p, --printPerfMatchings
        prints the perfect matchings of the graph
    -f, --print2Factors
        prints the 2-factors of the graph
    -F, --print2FactorSizes
        prints the sizes of the 2-factors of the graph




## Resources

Parts of the code are taken from or based on the following resources.

- [MatchMaker](https://gitlab.inria.fr/bora-ucar/matchmaker): a library used for finding perfect matchings in bipartite graphs using various algorithms. The needed files are included in the directory `matchmaker`.
- [nauty](https://pallini.di.uniroma1.it/): functionality regarding the manipulation of graphs and computing graph properties. The needed files of nauty version 2.8.9 are included in the directory `nauty`.
- [faultCost](https://github.com/JarneRenders/faultCost/): used the same style of flags given to the executable.


## Author

Tibo Van den Eede

E-mail: tibo [dot] vandeneede [at] kuleuven [dot] be