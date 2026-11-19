# Supplemental data for "Toroidal Obstructions with No K5 Minor"

This folder holds every computation the paper cites. An item cited as \supp{X} in the text is stored under the name X in the subfolder given by its letters:

- `no_K5_obstructions`: O1–O4, the 74 obstructions and their certificates.
- `witnesses`: W1–W16, single certificates cited in proofs.
- `rotations`: R1–R12, rotation tables.
- `enumerations`: E1–E30, one folder per enumeration.
- `no_K5_topological_obstructions`: TO1 and TO2, for the topological obstructions.
- `new_obstructions`: N1, 120 topological obstructions absent from Myrvold’s collection, including 22 minor obstructions.

The data files are JSON, except the Myrvold source collection `new_obstructions/N1_myrvold_topological.txt.gz` and the propositional formulas and their proofs in E3, E4 and E12 (`cover.cnf` in DIMACS format and `cover.drat` in DRAT format). A file larger than 25 MB is stored gzip-compressed, with the extension `.json.gz` (or `.drat.gz`); these are the three files of `rotations/R8`, the candidates and verdicts of E4, E8 and E20, the candidates of E25, `E29/bases.json.gz` and `E4/cover.drat.gz`. `new_obstructions/N1_sources.json.gz` also compresses the original website evidence. Any gzip tool opens them, and Python reads the compressed JSON files with `json.load(gzip.open(path, 'rt'))`.

## Graphs and vertex names

A graph is written as `{"vertices": [...], "edges": [[u, v], ...]}`, using the vertex names of the text. Where the text names no vertex, the file says what its names mean: each enumeration in the `vertex names` field of its `input.json`, each other item in its `what` or `note` field. The general conventions are these. In Section 4 the twins are p, q on {s1, s2, s3} and r, z on {t1, t2, t3}, and where the proof writes a = s1, b = t1 and so on, the file does too. In Sections 5 and 6 the twins are p, q on {a, b, c} and r, z on {a, b, d} or {a, d, e}; further vertices with these neighbourhoods are p3, p4, ... and r3, r4, .... Interior vertices of a bridge are u0, u1, ... (i0, i1, ... in some Section 4 items). A vertex that the text describes only by its neighbourhood is named by it, for example v_abd for a vertex with neighbourhood {a, b, d}. The graphs T1–T74 in O1 have vertices 0, ..., n − 1, numbered as in their graph6 strings (the format of nauty).

## The seven kinds of witness

Each certificate is a JSON object whose `type` says which kind it is.

- An embedding, `{"type": "embedding", "genus": g, "rotation": {v: [neighbours of v in cyclic order]}}`. Its facial walks are traced by the rule: after the dart u→v comes v→w, where w follows u in the cyclic order at v. Euler's formula then gives the genus.
- A K5 minor, `{"type": "K5 branch sets", "branch sets": [five vertex lists]}`: the sets are disjoint and connected, with an edge between every two.
- An isomorphism, `{"type": "isomorphism", "onto": name, "map": {v: image}}`: a bijection onto a graph of O1 or another named graph, carrying edges exactly onto edges.
- A minor model, `{"type": "minor model", "of": "Ti", "branch sets": {vertex of Ti: vertex list}}`: one connected set for each vertex of Ti, with an edge between two sets whenever the corresponding vertices are adjacent.
- A cut, `{"type": "cut", "vertices": [...]}`: at most two vertices whose deletion disconnects the graph.
- A Kuratowski subgraph, `{"type": "Kuratowski subgraph", "edges": [...]}`: the edges of a subdivision of K5 or K3,3, so the graph containing them is nonplanar.
- A certificate that a graph has no K5 minor: a tree of `{"type": "clique-sum", "separator": [...], "parts": [{"vertices": [...], "certificate": ...}]}` whose leaves are `{"type": "planar", "rotation": ...}` (a plane rotation system) or `{"type": "V8", "isomorphism": ...}`. Each part is the subgraph induced by its vertices with the separator made a clique, and separators have at most three vertices. The graph is a subgraph of a gluing of planar graphs and copies of V8 along cliques of order at most three, so it has no K5 minor by Wagner's theorem.

## Enumerations

Each folder E1–E30 has four files. `input.json` states the case: the passage of the paper, the construction and its parameters, the vertex names and what each verdict means in this case. `candidates.json` lists the candidates, each with an `id`. `verdicts.json` gives one verdict and witness for each candidate. `summary.json` gives the counts, including the numbers quoted in the text. In E3, E4 and E12 the candidates are recorded graphs, and `cover.cnf` with `cover.drat` (gzip-compressed for E4) is a propositional formula asserting that some combination of the case is not covered by them, with a DRAT proof of unsatisfiability that drat-trim accepts.

E1–E29 use the following verdicts.

- `toroidal`: an embedding of genus at most one.
- `K5 minor`: branch sets, in the graph or in a named minor of it.
- `cut of at most two vertices`: the cut.
- `isomorphic to Ti`, or to a bridge given in the text (E6, E7): an isomorphism. In E6 and E7 the record that the bridge passes every test of the case is kept beside it.
- `excluded by a hypothesis of the case`: the fact used, for example two twin triples meeting twice, with its own certificate where it needs one.
- `proper nontoroidal minor`: the edge deleted or contracted (`{"type": "proper minor", "operation": "delete" or "contract", "edge": [u, v]}`; in E1, E6 and E7, `{"type": "allowed reduction", "operation": "delete edge" or "contract edge", "arguments": [u, v]}` acting on the bridge). Under `why the minor is nontoroidal`, or in the named failing test, the file records why that minor is nontoroidal: a test of the case that fails, with a Kuratowski subgraph for each failure or the bridge that fits no face, or a minor model of some Ti. Where nothing is recorded, the minor contains a graph whose genus is in the genus computations of Appendix A.

E30 gives all eight choices of edges among the attachments of each of T(z), S(x;y) and Q(x;y,z). Each of the 24 verdicts gives the computed orientable genus and a witness for its planarity decision.

## Rotation tables

Each table lists every rotation system of its graph, with its genus and the attachments met by each facial walk, in order. For each graph, `cyclic orders` lists, for every vertex, all cyclic orders of its neighbours, each starting with its alphabetically first neighbour. Rotation system number i (counting from 0) takes, at each vertex, the cyclic order whose position is the corresponding digit of i in the mixed radix given by the numbers of cyclic orders, in the order of `vertices`, with the last vertex varying fastest. Since many rotation systems give the same row, the distinct rows are listed once under `distinct rows`, and `rows` gives, for every rotation system in turn, the position of its row in that list. A row is the genus followed by the attachment sequences of the facial walks: each sequence lists the attachments the walk passes through, one entry per passage, starting from its lexicographically least cyclic shift, and the faces are sorted. In R1–R3 a row also gives the region along a path (the faces on the two sides of the path glued along it). Each table ends with `claims`: every statement the text makes about the table, the data computed from the rows, and the result of the check. R8, the largest table (21,220,992 rotation systems of twelve exteriors), is a folder: `R8.json.gz` is the table, and the rows of its largest exterior (17,694,720 rotation systems) are stored in `rows_1.json.gz` and `rows_2.json.gz`, in order, 9,000,000 rows to a part, so that no file is larger than 25 MB.

## The additional obstructions of N1

N1 lists the 120 topological obstructions in the website release of November 11, 2026 that are absent from all 250,815 graphs in Myrvold’s topological collection. Exactly 22 are minor obstructions. N1.1–N1.13 retain their original names, graph labels and certificates. N1.14 is T74. Each of the other 119 graphs has explicit K5 branch sets.

For each graph, N1 records its website identifier, graph6 string, edges, connectivity, genus-two embedding and a torus embedding of every edge deletion. For each minor obstruction it also records a torus embedding of every edge contraction. Each of the other 98 graphs has an edge contraction with a recorded completed genus-two computation. Contraction records in N1 merge the second endpoint into the first and retain the other vertex labels. The original website evidence relabels contracted graphs consecutively, and the checker checks this correspondence.

`N1_sources.json.gz` contains the original website evidence JSON text, its SHA-256 and Git-blob provenance, the original thirteen records, K5 branch sets, the published release manifest and the comparison results. The completed PAGE or multi_genus searches include their exact input checksums, program checksums, search settings and output. A genus-two rotation proves an upper bound. The lower bound comes from the recorded completed exact search. The checker checks the identity and consistency of these records, and `python3 programs/new_obstructions/make_N1.py /tmp/N1-rebuilt --recompute` reruns their required genus computations with multi_genus.

`N1_myrvold_topological.txt.gz` is the full source collection downloaded from https://webhome.cs.uvic.ca/~wendym/torus/all_obs.txt.gz. Its SHA-256 is `3772744c0af7aa4677b9fa8dc19195dfc5e04f81837411f4b4bc7eb921e83e4c`. The checker compares every N1 graph against this collection by invariants and direct isomorphism tests, and also by canonical forms when nauty’s `labelg` is available.

## Checksums

`SHA256SUMS.txt` gives the SHA-256 checksum of every other file in the folder.
