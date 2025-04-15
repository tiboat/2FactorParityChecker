#!/usr/bin/env python


"""
Lift construction script with two options:
  1. Lift the Gray graph using (Z_{3},+).
  2. Lift the theta graph (two vertices and 3 parallel edges) using (Z_{15},+).
  
Usage (example):
    python genAllLiftsGrayOrTheta.py --graph-type theta --min-girth 6
    python genAllLiftsGrayOrTheta.py --graph-type gray  --min-girth 10

If --graph-type is "gray", the program uses the Gray graph.
If --graph-type is "theta", the program uses the theta graph.

Author: Tibo Van den Eede
Code based on output of Chat-GPT
"""

import networkx as nx
import argparse
import copy
import sys

# --- Helper function: check for short cycles using minimum cycle basis ---
def has_short_cycle(G, threshold):
    """
    Returns True if the undirected graph G has any cycle of length less than threshold.
    """
    return nx.girth(G) < threshold

# --- Recursive backtracking function ---
def assign_non_tree_arcs(non_tree_arcs, current_index, current_lift, group, threshold):
    """
    Recursively assign a group element from 'group' (list of integers) to each non-tree arc.
    For each non-tree arc (u,v) and assignment c, we add, for every g in group, an edge:
         ((u, g), (v, (g + c) mod n)).
    Prune if the current lifted graph already contains a cycle with length less than threshold.
    
    Returns a list of valid lifted graphs.
    """
    if current_index == len(non_tree_arcs):
        # All non-tree arcs have been assigned.
        if not has_short_cycle(current_lift, threshold):
            print(nx.to_graph6_bytes(current_lift, header=False).decode("utf-8"), end="")
        return

    # Get the current arc to assign; arc is a tuple (u, v)
    arc = non_tree_arcs[current_index]
    for c in group:
        # Make a deep copy to extend this branch
        lift_copy = copy.deepcopy(current_lift)
        # For every g in group, add an edge for the current arc with assignment c:
        # edge: ((u, g), (v, (g + c) mod len(group)))

        parallel_edge = False
        for g in group:
            if not lift_copy.has_edge((arc[0], g), (arc[1], (g + c) % len(group))):
                lift_copy.add_edge((arc[0], g), (arc[1], (g + c) % len(group)))
            else:
                # If parallel edges appear -> prune
                parallel_edge = True
                break
        
        if parallel_edge or has_short_cycle(lift_copy, threshold):
            continue

        # Recurse on the next non-tree arc
        assign_non_tree_arcs(non_tree_arcs, current_index + 1, lift_copy, group, threshold)
       

# --- Build the lifted graph base for spanning tree arcs ---
def build_lift_base(G_nodes, spanning_tree_arcs, group):
    """
    Create a lifted graph with vertex set = { (v, g) : v in nodes and g in group }.
    For each spanning tree arc (u, v), add for each g in group an edge from (u, g) to (v, g).
    Returns an undirected graph.
    """
    lift = nx.Graph()
    # Add all vertices: (v, g) for each v in V and g in Z_n.
    for v in G_nodes:
        for g in group:
            lift.add_node((v, g))
    # For each spanning tree arc (u, v), add edges ((u, g), (v, g))
    for u, v in spanning_tree_arcs:
        for g in group:
            lift.add_edge((u, g), (v, g))
    return lift

# --- Function to build the starting graph based on the chosen option ---
def build_starting_graph(graph_type):
    """
    Depending on graph_type ("gray" or "theta"), returns:
      - the base graph G (undirected),
      - a directed version DG with an arbitrary orientation, and
      - the group to use (list of integers).
      
    For "gray": we attempt to read graph6 data (here using a placeholder).
    For "theta": we build the theta graph, a multigraph on two vertices with 3 parallel edges;
                 we then assign all edges the direction (0,1).
    """
    if graph_type == "gray":
        gray_graph_bytes = b"us??????????B??o@C?GG@A?GO?B??CG?AG??S??@G??B????_???G???A????G???A?????O???@????A????A????O??????G????@??????_@????_A????B?????B??????CG?????a??????Q?????AG??????@_?????@_???????g??????B??O?????AOO?????@_A??????A`??????@OC?????c??_????@_??"
        try:
            G = nx.from_graph6_bytes(gray_graph_bytes)
        except Exception as e:
            print("Error reading Gray graph from graph6 data:", e)
            sys.exit(1)
        # For undirected graphs, assign a direction arbitrarily (here: from smaller to larger).
        DG = nx.DiGraph()
        for u, v in G.edges():
            if u < v:
                DG.add_edge(u, v)
            else:
                DG.add_edge(v, u)
        group = list(range(3))  # Z_3
    elif graph_type == "theta":
        # Build the theta graph.
        # Use a MultiDiGraph to allow parallel edges.
        DG = nx.MultiDiGraph()
        DG.add_node(0)
        DG.add_node(1)
        # Add 3 parallel edges, all directed from 0 to 1.
        for _ in range(3):
            DG.add_edge(0, 1)
        # For spanning tree purposes, we create a simple undirected version where
        # parallel edges are replaced by a single edge.
        G_simple = nx.Graph()
        G_simple.add_edge(0, 1)
        # The spanning tree of a 2-vertex graph is the single edge (0,1).
        group = list(range(15))  # Z_{15}
        # In DG, all edges are (0,1) – we will later treat one as spanning-tree edge and the others as non-tree.
        G = G_simple  # For consistency, use the simple version for node set.
    else:
        print("Unknown graph type:", graph_type)
        sys.exit(1)
    return G, DG, group

# --- Main function ---
def main():
    parser = argparse.ArgumentParser(
        description="Lift Construction for Two Graph Options: 'gray' (with Z_{3}) or 'theta' (with Z_{15})."
    )
    parser.add_argument("--graph-type", choices=["gray", "theta"], default="gray",
                        help="Choose the base graph: 'gray' for the Gray graph, 'theta' for the theta graph.")
    parser.add_argument("--min-girth", type=int, default=10,
                        help="Minimum allowed cycle length in the lifted graph (girth threshold).")
    args = parser.parse_args()

    # Build the starting graph based on the chosen flag.
    G, DG, group = build_starting_graph(args.graph_type)

    # Compute a spanning tree.
    # For graphs with parallel edges we work on the underlying simple graph.
    if args.graph_type == "theta":
        # For the theta graph, the spanning tree is the single edge (0,1).
        T = nx.DiGraph()
        T.add_edge(0, 1)
    else:
        # For the Gray graph, compute a spanning tree on the undirected G.
        T_undirected = nx.minimum_spanning_tree(G)
        T = nx.DiGraph()
        for u, v in T_undirected.edges():
            # Orient to agree with DG if possible.
            if DG.has_edge(u, v):
                T.add_edge(u, v)
            else:
                T.add_edge(v, u)

    # Determine non-tree arcs: edges in DG (ignoring parallel edges in theta case)
    non_tree_arcs = []
    if args.graph_type == "theta":
        non_tree_arcs.append((0, 1))
        non_tree_arcs.append((0, 1))
    else:
        # Gray graph case.
        for u, v in DG.edges():
            if not T.has_edge(u, v):
                non_tree_arcs.append((u, v))
    

    # Build the initial lifted graph using the spanning tree arcs.
    lift_base = build_lift_base(G.nodes(), list(T.edges()), group)

    # Process non-tree arcs recursively.
    assign_non_tree_arcs(non_tree_arcs, 0, lift_base, group, args.min_girth)


if __name__ == '__main__':
    main()
