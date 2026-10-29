# Supplemental data for "Toroidal Obstructions with No K5 Minor"

This folder holds every computation the paper cites. An item cited as \supp{X} in the text is stored under the name X in the subfolder given by its letters:

- `no_K5_obstructions`: O1–O4, the 74 obstructions and their certificates.
- `witnesses`: W1–W21, single certificates cited in proofs.
- `rotations`: R1–R12, rotation tables.
- `enumerations`: E1–E34, one folder per enumeration.
- `no_K5_topological_obstructions`: TO1 and TO2, for the topological obstructions.
- `new_obstructions`: N1, thirteen toroidal minor obstructions that are not in the collection of Myrvold and Woodcock.

Every file is JSON, except the propositional formulas and their proofs in E4, E5 and E13 (`cover.cnf` in DIMACS format and `cover.drat` in DRAT format). A file larger than 25 MB is stored gzip-compressed, with the extension `.json.gz` (or `.drat.gz`); these are the three files of `rotations/R8`, the candidates and verdicts of E5, E9 and E22, the candidates of E30, `E34/bases.json.gz` and `E5/cover.drat.gz`. Any gzip tool opens them, and Python reads them with `json.load(gzip.open(path, 'rt'))`.

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

Each folder E1–E34 has four files. `input.json` states the case: the passage of the paper, the construction and its parameters, the vertex names and what each verdict means in this case. `candidates.json` lists the candidates, each with an `id`. `verdicts.json` gives one verdict and witness for each candidate. `summary.json` gives the counts, including the numbers quoted in the text. In E4, E5 and E13 the candidates are recorded graphs, and `cover.cnf` with `cover.drat` (gzip-compressed for E5) is a propositional formula asserting that some combination of the case is not covered by them, with a DRAT proof of unsatisfiability that drat-trim accepts.

The verdicts are these.

- `toroidal`: an embedding of genus at most one.
- `K5 minor`: branch sets, in the graph or in a named minor of it.
- `cut of at most two vertices`: the cut.
- `isomorphic to Ti`, or to a bridge given in the text (E7, E8): an isomorphism. In E7 and E8 the record that the bridge passes every test of the case is kept beside it.
- `excluded by a hypothesis of the case`: the fact used, for example two twin triples meeting twice, with its own certificate where it needs one.
- `proper nontoroidal minor`: the edge deleted or contracted (`{"type": "proper minor", "operation": "delete" or "contract", "edge": [u, v]}`; in E1, E7 and E8, `{"type": "allowed reduction", "operation": "delete edge" or "contract edge", "arguments": [u, v]}` acting on the bridge). Under `why the minor is nontoroidal`, or in the named failing test, the file records why that minor is nontoroidal: a test of the case that fails, with a Kuratowski subgraph for each failure or the bridge that fits no face, or a minor model of some Ti. Where nothing is recorded, the minor contains a graph whose genus is in the genus computations of Appendix A.

## Rotation tables

Each table lists every rotation system of its graph, with its genus and the attachments met by each facial walk, in order. For each graph, `cyclic orders` lists, for every vertex, all cyclic orders of its neighbours, each starting with its alphabetically first neighbour. Rotation system number i (counting from 0) takes, at each vertex, the cyclic order whose position is the corresponding digit of i in the mixed radix given by the numbers of cyclic orders, in the order of `vertices`, with the last vertex varying fastest. Since many rotation systems give the same row, the distinct rows are listed once under `distinct rows`, and `rows` gives, for every rotation system in turn, the position of its row in that list. A row is the genus followed by the attachment sequences of the facial walks: each sequence lists the attachments the walk passes through, one entry per passage, starting from its lexicographically least cyclic shift, and the faces are sorted. In R1–R3 a row also gives the region along a path (the faces on the two sides of the path glued along it). Each table ends with `claims`: every statement the text makes about the table, the data computed from the rows, and the result of the check. R8, the largest table (21,220,992 rotation systems of twelve exteriors), is a folder: `R8.json.gz` is the table, and the rows of its largest exterior (17,694,720 rotation systems) are stored in `rows_1.json.gz` and `rows_2.json.gz`, in order, 9,000,000 rows to a part, so that no file is larger than 25 MB.

## The thirteen obstructions of N1

N1 gives thirteen toroidal minor obstructions, N1.1, ..., N1.13, that are not in the collection of Myrvold and Woodcock. They were found in searches during this research. Each has a K5 minor, so none of them is among T1–T74. For each graph, N1 gives its graph6 string and edges, an embedding of genus two, a torus embedding of every edge deletion and of every edge contraction, K5 branch sets, its connectivity and its canonical graph6 string (nauty's `labelg`). Since every edge deletion and edge contraction is toroidal, every proper minor is toroidal. That the graph itself is not toroidal is the genus computation recorded under `genus` (genus two by multi_genus, as bundled with graph-genus 0.1.3), which `programs/new_obstructions/make_N1.py` repeats. Each graph is isomorphic to no entry of the collection in O4: N1 lists the entries with the same invariants (orders, sizes, degrees and triangles), and there are none.

## Checksums

`SHA256SUMS.txt` gives the SHA-256 checksum of every other file in the folder.
