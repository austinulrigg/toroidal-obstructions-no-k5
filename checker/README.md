# Checking the supplemental data

From the root of the repository, run the checker with Python 3 and NetworkX. The second command checks only the named items.

```sh
python3 checker/check_supplement.py supplemental_data
python3 checker/check_supplement.py supplemental_data --items E3,TO2
```

The checker first compares the files in `supplemental_data` with `SHA256SUMS.txt`, even when checking a single item. For each selected item it then checks the graphs, witnesses and counts specified in the corresponding module below. It prints a result for each item and lists any failures. A full run takes about half an hour on one core.

| Module | Items |
| --- | --- |
| `items_obstructions.py` | O1–O4, TO1, TO2 and the original N1 format |
| `items_new_obstructions.py` | N1 (120 graphs) |
| `items_witnesses.py` | W1–W16 |
| `items_rotations.py`, `items_rotation_R12.py` | R1–R11 and R12, respectively |
| `items_bridges.py` | E1, E2, E3, E6, E7, E11 |
| `items_coverings.py` | E4, E12 |
| `items_enumerations_a.py` | E5, E8–E10, E13–E18 |
| `items_enumerations_b.py` | E19–E29 |
| `items_three_root_completions.py` | E30 |

The common witness checks are in [witnesses.py](witnesses.py). For an embedding, the checker traces the facial walks of the given rotation system and computes its genus by Euler’s formula. For a minor, it checks that the branch sets are nonempty, connected and pairwise disjoint, with an edge between the sets for each edge of the target graph. The [data guide](../supplemental_data/README.md) describes the other witnesses.

Some enumeration records give only the edge whose deletion or contraction yields a nontoroidal minor. The checker verifies the edge operation. The genus of the minor is determined by the separate [genus computations](../genus/README.md).

For E3, E4 and E12, the checker verifies the DRAT proof with `drat-trim`. It uses the executable named by `DRAT_TRIM`, searches for `drat-trim` on the search path, or compiles the [bundled source](../tools/drat-trim.c) with `cc`. The source is distributed under the [drat-trim MIT license](https://github.com/marijnheule/drat-trim). If no executable is available, the item reports that its DRAT proof was skipped. The checker also recomputes the canonical forms in O2, O4 and N1 when nauty’s `labelg` is on the search path or its path is given in `LABELG`.

For N1, the checker verifies all edge-deletion embeddings, all edge-contraction embeddings for the 22 minor obstructions, the 119 K5 minor models and the certificate in O2 for T74. It checks the inputs, settings and outputs of the recorded genus-two computations for the 120 graphs and one contraction of each of the other 98 graphs. The [N1 generator](../programs/README.md#the-additional-obstructions-in-n1) can repeat these searches. To compare N1 with all 250,815 entries in the Myrvold–Woodcock topological collection, the checker tests isomorphism for graphs with matching invariants.
