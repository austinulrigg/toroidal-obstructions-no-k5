# Supplemental data

The paper cites each item by its label. The [index](../PAPER-REFERENCES.md) gives the corresponding theorem, lemma or section.

| Items | Contents |
| --- | --- |
| [O1–O4](no_K5_obstructions/) | The 74 minor obstructions, their certificates and the comparisons with the Myrvold–Woodcock collection |
| [W1–W16](witnesses/) | Individual witnesses used in the proofs |
| [R1–R12](rotations/) | Rotation systems and facial walks |
| [E1–E30](enumerations/) | Candidate graphs, verdicts and counts |
| [TO1 and TO2](no_K5_topological_obstructions/) | The 381 topological obstructions and the vertex splits used to find them |
| [N1](new_obstructions/) | The 120 obstructions in the author’s catalogue absent from the Myrvold–Woodcock topological collection |

The records are JSON. Larger files have the extension `.json.gz`. E3, E4 and E12 also contain DIMACS formulas and DRAT proofs. [SHA256SUMS.txt](SHA256SUMS.txt) gives the checksum of each file in this folder except the manifest itself.

## Graphs and vertex names

A graph lists its vertices in `vertices` and its edges as pairs in `edges`. O1 numbers the vertices from 0 in graph6 order. The proof records use the paper’s labels wherever possible.

In Section 4, the twins are $p,q$ on $\{s_1,s_2,s_3\}$ and $r,z$ on $\{t_1,t_2,t_3\}$. In Sections 5 and 6, $p,q$ have neighborhood $\{a,b,c\}$, and $r,z$ have neighborhood $\{a,b,d\}$ or $\{a,d,e\}$, respectively. The labels for additional twins start at `p3` and `r3`. Interior vertices of a bridge are usually `u0`, `u1`, …, or `i0`, `i1`, … in some Section 4 records. A label such as `v_abd` denotes a vertex with neighborhood $\{a,b,d\}$. Each enumeration gives its conventions in `input.json`.

## Witnesses

| Type | What the record gives |
| --- | --- |
| `embedding` | The neighbors of each vertex in cyclic order and the genus. After the directed edge $uv$, a facial walk takes $vw$, where $w$ follows $u$ in the order at $v$. Euler’s formula gives the genus of each component. |
| `minor model` | A nonempty connected branch set for each vertex of the named minor. The sets are disjoint, with an edge between the sets of adjacent vertices. |
| `K5 branch sets` | Five disjoint nonempty connected sets with an edge between every two. |
| `isomorphism` | The target graph and a vertex bijection preserving adjacency in both directions. |
| `cut` | At most two vertices whose deletion disconnects the graph. |
| `Kuratowski subgraph` | The edges of a subdivision of $K_5$ or $K_{3,3}$. |

For a graph with no $K_5$ minor, the certificate gives a decomposition along cliques of order at most three. A `clique-sum` record lists the separator and the vertices of each part. Each part is induced on its listed vertices, with the separator completed to a clique. The leaves give plane rotation systems (`planar`) or isomorphisms to the Wagner graph (`V8`). Wagner’s theorem then gives the result for the original graph.

## Enumerations

Four files occur in every enumeration. `input.json` describes the construction and its parameters. `candidates.json` lists the graphs or bridge forms, `verdicts.json` gives the result for each candidate by ID, and `summary.json` gives the counts. E29 also records its base graphs in `bases.json.gz`.

In E1–E29, a toroidal verdict gives an embedding of genus zero or one. A cut, a $K_5$ minor or an isomorphism to a listed graph has the corresponding witness above. E6 and E7 also identify bridges from the paper and record the outcomes of their planarity tests. When a graph violates a hypothesis, the record specifies which one, for example two distinct twin triples meeting in two vertices.

For a proper nontoroidal minor, the witness names the edge deleted or contracted. Some records give an allowed reduction of the bridge. Depending on the case, nontoroidality follows from a minor model, a failed embedding test or an [exact genus computation](../genus/README.md). A failed test identifies a Kuratowski subdivision or a bridge that fits in none of the available faces.

In E3, E4 and E12 the candidates are the graphs used to cover the allowed combinations. `cover.cnf` asserts that some allowed combination is contained in no recorded toroidal graph and contains no recorded nontoroidal graph. The DRAT proof establishes that this formula is unsatisfiable. E4 stores its proof as `cover.drat.gz`.

E30 gives the orientable genus for each of the eight choices of attachment edges on each of the three graphs in Lemma 3.10. Every record includes a rotation system, and each nonplanar case also includes a Kuratowski subdivision.

## Rotation tables

Each rotation system is a choice of one cyclic order at every vertex. The table lists these choices in vertex order, starting at index 0, with the last vertex varying fastest. For rotation system $i$, `rows[i]` is the index of its entry in `distinct rows`.

An entry gives the genus and the attachment sequences of the facial walks. Each sequence records every appearance of an attachment in walk order, starting at its least cyclic shift. R1–R3 also describe the region along a specified path in a torus embedding. The `claims` field records the checks made from the table.

[R8](rotations/R8/) contains 21,220,992 rotation systems for twelve exteriors. For the largest exterior, the 17,694,720 row indices are stored consecutively in `rows_1.json.gz` and `rows_2.json.gz`.

## The additional obstructions in N1

N1 compares the author’s November 11, 2026 catalogue with the 250,815 graphs in the [Myrvold–Woodcock topological collection](https://webhome.cs.uvic.ca/~wendym/torus/all_obs.txt.gz). The source collection is included as `N1_myrvold_topological.txt.gz`. Of the 120 graphs absent from that collection, 22 are minor obstructions. N1.1–N1.13 retain their original labels. N1.14 is $\mathfrak T_{74}$. Each of the other 119 graphs has a recorded $K_5$ minor.

Each record gives the catalogue identifier, graph6 string, edges, connectivity, a genus-two embedding and torus embeddings for every edge deletion. For the 22 minor obstructions, N1 also gives a torus embedding for every edge contraction. For each of the other 98 graphs, it gives an edge whose contraction has genus two. A contraction merges the second endpoint into the first and keeps the remaining labels.

`N1_sources.json.gz` preserves the original catalogue records, the release manifest and the genus-search inputs, settings, outputs and checksums. The rotations give the upper bounds on genus. The lower bounds come from completed PAGE or multi_genus computations. The [checker](../checker/README.md) verifies the records, matches the vertex labels after contraction and compares the graphs with the source collection by isomorphism. The [N1 generator](../programs/README.md) can also rerun the genus computations.
