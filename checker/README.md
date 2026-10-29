# Checker

`check_supplement.py` checks every item of the supplemental data, using only the files of the folder `supplemental_data`. It first compares every file with its SHA-256 checksum in `SHA256SUMS.txt`. Then, for each item, it rebuilds every graph from its record and checks every witness by its kind, together with the structure of the item.

The seven kinds of witness (Appendix A of the paper) are checked in `witnesses.py`. An embedding is a rotation system whose facial walks are traced, its genus following from Euler's formula. K5 branch sets are five disjoint connected vertex sets with an edge between every two. An isomorphism is a bijection carrying edges exactly onto edges. A minor model is a connected vertex set for each vertex of the target graph, with an edge between the sets of adjacent vertices. A cut is a set of at most two vertices whose deletion disconnects the graph. A Kuratowski subgraph is a set of edges forming a subdivision of K5 or K3,3. A certificate that a graph has no K5 minor writes it as a subgraph of a gluing of planar graphs and copies of V8 along cliques of at most three vertices.

Each module checks one group of items, and its opening comment lists what it checks.

| Module | Items |
|---|---|
| `items_obstructions.py` | O1–O4, TO1, TO2, N1 |
| `items_witnesses.py` | W1–W21 |
| `items_rotations.py` | R1–R11 |
| `items_bridges.py` | E1, E2, E3, E4, E7, E8, E12 |
| `items_coverings.py` | E5, E13 |
| `items_enumerations_a.py` | E6, E9, E10, E11, E14–E19 |
| `items_enumerations_b.py` | E20–E34 |

Where a verdict names only the edge whose deletion or contraction gives a nontoroidal proper minor, the nontoroidality of that minor comes from the genus computation in the folder `genus`, and the checker counts these minors. The DRAT proofs of E4, E5 and E13 are checked with drat-trim. The checker finds drat-trim through the environment variable `DRAT_TRIM` or on the PATH, and otherwise compiles it from `tools/drat-trim.c` with the C compiler `cc` the first time it is needed. Without drat-trim it reports these proofs as not checked.

## Running it

The checker needs Python 3 and the package networkx. From the top folder of the repository:

```
cd checker
python3 check_supplement.py ../supplemental_data
```

It prints one line for each item and ends with `all checks passed`, or with a list of failures. The whole run takes about half an hour on one core of a recent computer. R8, R4, E5 and TO2 take the longest. One item, or a few, can be checked alone:

```
python3 check_supplement.py ../supplemental_data --items E4,TO2
```

If nauty's program `labelg` is on the PATH, or its path is given by the environment variable `LABELG`, the checker also recomputes the canonical forms recorded in O2, O4 and N1.
