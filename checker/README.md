# Checker

`check_supplement.py` checks every item of the supplemental data, using only the files of the folder `supplemental_data`. It first compares every file with its SHA-256 checksum in `SHA256SUMS.txt`. Then, for each item, it rebuilds every graph from its record and checks every witness by its kind, together with the structure of the item.

The seven kinds of witness (Appendix A of the paper) are checked in `witnesses.py`. An embedding is a rotation system whose facial walks are traced, its genus following from Euler's formula. K5 branch sets are five disjoint connected vertex sets with an edge between every two. An isomorphism is a bijection carrying edges exactly onto edges. A minor model is a connected vertex set for each vertex of the target graph, with an edge between the sets of adjacent vertices. A cut is a set of at most two vertices whose deletion disconnects the graph. A Kuratowski subgraph is a set of edges forming a subdivision of K5 or K3,3. A certificate that a graph has no K5 minor writes it as a subgraph of a gluing of planar graphs and copies of V8 along cliques of at most three vertices.

Each module checks one group of items, and its opening comment lists what it checks.

| Module | Items |
|---|---|
| `items_obstructions.py` | O1–O4, TO1, TO2, original N1 format |
| `items_new_obstructions.py` | expanded N1 with 120 graphs |
| `items_witnesses.py` | W1–W16 |
| `items_rotations.py` | R1–R11 |
| `items_rotation_R12.py` | R12 |
| `items_bridges.py` | E1, E2, E3, E6, E7, E11 |
| `items_coverings.py` | E4, E12 |
| `items_enumerations_a.py` | E5, E8, E9, E10, E13–E18 |
| `items_enumerations_b.py` | E19–E29 |
| `items_three_root_completions.py` | E30 |

Where a verdict names only the edge whose deletion or contraction gives a nontoroidal proper minor, the nontoroidality of that minor comes from the genus computation in the folder `genus`, and the checker counts these minors. The DRAT proofs of E3, E4 and E12 are checked with drat-trim. The checker finds drat-trim through the environment variable `DRAT_TRIM` or on the PATH, and otherwise compiles it from `tools/drat-trim.c` with the C compiler `cc` the first time it is needed. Without drat-trim it reports these proofs as not checked.

## Running it

The checker needs Python 3 and the package networkx. From the top folder of the repository:

```
cd checker
python3 check_supplement.py ../supplemental_data
```

It prints one line for each item and ends with `all checks passed`, or with a list of failures. The whole run takes about half an hour on one core of a recent computer. R8, R4, E4 and TO2 take the longest. One item, or a few, can be checked alone:

```
python3 check_supplement.py ../supplemental_data --items E3,TO2
```

If nauty's program `labelg` is on the PATH, or its path is given by the environment variable `LABELG`, the checker also recomputes the canonical forms recorded in O2, O4 and N1.

The expanded N1 check validates all edge-deletion embeddings, all contraction embeddings for its 22 minor obstructions, 119 K5 models and T74’s certificate in O2. It also checks the recorded genus-two computations for the 120 graphs and the nontoroidal contraction for each of the other 98 graphs, including input hashes, exact-search settings and completed output. These checks establish the records’ consistency without repeating the exhaustive searches. The generator’s `--recompute` option repeats those searches. Absence from Myrvold’s full 250,815-graph topological collection is checked by direct isomorphism tests after invariant filtering and, when available, by nauty canonical forms.
