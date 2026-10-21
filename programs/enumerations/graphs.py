"""Graph routines shared by the enumeration programs.

Graphs are networkx.Graph objects. Embeddings are rotation systems, stored
as a dict sending each vertex to the cyclic list of its neighbors, together
with the face walks that the rotation system determines.

Every genus computation goes through PAGE (the graph_genus package).
"""

import gzip
import os
import random
import subprocess
import tempfile
from itertools import combinations
from pathlib import Path

import networkx as nx
import graph_genus


# ----------------------------------------------------------------------------
# Reading and building graphs

def from_graph6(s):
    """The graph with the given graph6 string, on vertices 0, ..., n-1."""
    return nx.from_graph6_bytes(s.encode())


def read_myrvold(path):
    """The graphs of a Myrvold-Woodcock database file.

    Each line holds the order n and then the n(n-1)/2 entries above the
    diagonal of the adjacency matrix, read row by row: (0,1), (0,2), ...,
    (0,n-1), (1,2), and so on.
    """
    out = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            n, bits = line.split()
            n = int(n)
            pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
            assert len(pairs) == len(bits)
            G = nx.empty_graph(n)
            G.add_edges_from(p for p, b in zip(pairs, bits) if b == "1")
            out.append(G)
    return out


def delete_edge(G, e):
    """G - e."""
    H = G.copy()
    H.remove_edge(*e)
    return H


def contract_edge(G, e):
    """G/e as a simple graph: the endpoints become the first endpoint,
    loops are deleted and parallel edges are merged."""
    u, v = e
    return nx.contracted_nodes(G, u, v, self_loops=False, copy=True)


def k33_disjoint_union():
    """Two disjoint copies of K_{3,3}."""
    return nx.disjoint_union(nx.complete_bipartite_graph(3, 3),
                             nx.complete_bipartite_graph(3, 3))


def k33_sharing_a_vertex():
    """Two copies of K_{3,3} with exactly one vertex in common."""
    A = nx.complete_bipartite_graph(3, 3)
    B = nx.relabel_nodes(nx.complete_bipartite_graph(3, 3),
                         {i: i + 5 for i in range(6)})  # vertex 5 is shared
    return nx.compose(A, B)


def v8():
    """The Wagner graph V_8: an 8-cycle with its four long diagonals."""
    G = nx.cycle_graph(8)
    G.add_edges_from((i, i + 4) for i in range(4))
    return G


def isomorphism(G, H):
    """A vertex map from G onto H that sends edges to edges, or None."""
    if (G.number_of_nodes() != H.number_of_nodes()
            or G.number_of_edges() != H.number_of_edges()):
        return None
    matcher = nx.algorithms.isomorphism.GraphMatcher(G, H)
    if matcher.is_isomorphic():
        return {(int(k) if isinstance(k, int) else k): (int(v) if isinstance(v, int) else v)
                for k, v in matcher.mapping.items()}
    return None


# ----------------------------------------------------------------------------
# Embeddings

def faces_of_rotation(rotation):
    """The face walks of a rotation system.

    A face is traced dart by dart: arriving at v along the edge uv, leave
    along the edge vw, where w follows u in the rotation at v. Each face is
    returned as the list of vertices met in order.
    """
    successor = {}
    for v, nbrs in rotation.items():
        k = len(nbrs)
        for i, u in enumerate(nbrs):
            successor[(v, u)] = nbrs[(i + 1) % k]
    unused = {(u, v) for v, nbrs in rotation.items() for u in nbrs}
    faces = []
    while unused:
        start = min(unused)
        walk, dart = [], start
        while True:
            unused.discard(dart)
            u, v = dart
            walk.append(u)
            dart = (v, successor[(v, u)])
            if dart == start:
                break
        faces.append(walk)
    return faces


def genus_of_rotation(G, rotation):
    """The genus of the surface on which a rotation system of the connected
    graph G embeds it, by Euler's formula V - E + F = 2 - 2g."""
    if G.number_of_edges() == 0:
        return 0
    f = len(faces_of_rotation(rotation))
    twice = 2 - G.number_of_nodes() + G.number_of_edges() - f
    assert twice % 2 == 0
    return twice // 2


def _page_with_time_limit(adjacency, seconds):
    """PAGE called the way graph_genus.embed calls it (package version 0.1.3),
    but stopped after the given number of seconds, raising
    subprocess.TimeoutExpired."""
    from graph_genus import _core
    max_degree = max(len(x) for x in adjacency)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "graph.adj"
        _core._write_page_adjacency_file(adjacency, path, max_degree)
        env = os.environ.copy()
        env.update({"STDOUT": "1", "PBN": "1", "S": "0",
                    "DEG": str(max_degree), "ADJ": str(path)})
        proc = subprocess.run([_core._tool_path("page")], stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, env=env, timeout=seconds)
    if proc.returncode != 0:
        raise RuntimeError("PAGE exited with status %d" % proc.returncode)
    return _core._parse_page_solution(proc.stdout.decode(), proc.stderr.decode(), adjacency)


def _page_block(B, seconds=None):
    """PAGE on one block B (two-connected, at least three vertices)."""
    verts = sorted(B.nodes())
    index = {v: i for i, v in enumerate(verts)}
    adjacency = [[index[u] for u in sorted(B.neighbors(v))] for v in verts]
    if seconds is None:
        g, rot = graph_genus.embed(adjacency)
    else:
        g, rot = _page_with_time_limit(adjacency, seconds)
    rotation = {verts[i]: [verts[j] for j in rot[i]] for i in range(len(verts))}
    for v in verts:
        assert sorted(rotation[v]) == sorted(B.neighbors(v)), "bad rotation"
    return g, rotation


def page_embedding(G, seconds=None):
    """The genus of G from PAGE, with an embedding of that genus.

    PAGE is run on each block of G with at least three vertices. The genus of
    a graph is the sum of the genera of its blocks (Battle, Harary, Kodama and
    Youngs), and an embedding of G of that genus is obtained by joining the
    block rotations at each cut vertex, one block's neighbors after another's.
    The face count of the joined rotation system is checked against the sum.

    The embedding is returned as one rotation system and face list for each
    connected component of G. With a time limit in seconds, a PAGE call that
    runs longer raises subprocess.TimeoutExpired.
    """
    total, parts = 0, []
    for comp in sorted(nx.connected_components(G), key=min):
        H = G.subgraph(comp)
        verts = sorted(H.nodes())
        rotation = {v: [] for v in verts}
        g_comp = 0
        for block in sorted(nx.biconnected_components(H), key=min):
            B = H.subgraph(block)
            if B.number_of_nodes() == 2:  # a bridge
                u, v = sorted(block)
                rotation[u].append(v)
                rotation[v].append(u)
                continue
            g, rot = _page_block(B, seconds)
            g_comp += g
            for v in block:
                rotation[v].extend(rot[v])
        if genus_of_rotation(H, rotation) != g_comp:
            raise RuntimeError("the joined rotation system does not have PAGE's genus")
        total += g_comp
        parts.append({"vertices": verts, "genus": g_comp,
                      "rotation": {str(v): rotation[v] for v in verts},
                      "faces": faces_of_rotation(rotation)})
    return total, parts


def page_genus(G):
    return page_embedding(G)[0]


# ----------------------------------------------------------------------------
# Connectivity and twins

def connectivity(G):
    """Vertex connectivity; 0 for a disconnected graph."""
    if not nx.is_connected(G):
        return 0
    return nx.node_connectivity(G)


def twin_triples(G):
    """The twin triples of G, each with its twins.

    Twins are distinct vertices of degree three with the same neighborhood.
    A twin triple is an independent triple that is the neighborhood of at
    least two twins.
    """
    by_nbhd = {}
    for v in G.nodes():
        if G.degree(v) == 3:
            by_nbhd.setdefault(tuple(sorted(G.neighbors(v))), []).append(v)
    triples = {}
    for nbhd, vs in by_nbhd.items():
        independent = not any(G.has_edge(x, y) for x, y in combinations(nbhd, 2))
        if len(vs) >= 2 and independent:
            triples[nbhd] = sorted(vs)
    return dict(sorted(triples.items()))


def case_of(G, triples):
    """The case of the proof of Theorem 1.1 that G falls under, computed from
    its connectivity and twin triples."""
    k = connectivity(G)
    if k <= 1:
        return "connectivity at most one"
    if k == 2:
        return "connectivity two"
    if isomorphism(G, nx.complete_bipartite_graph(3, 7)) is not None:
        return "K3,7"
    ts = [set(t) for t in triples]
    if len(ts) < 2:
        return "three-connected with fewer than two twin triples"
    shared = [len(a & b) for a, b in combinations(ts, 2)]
    if 0 in shared:
        return "Theorem 4.1"
    if 2 in shared:
        return "Theorem 5.1"
    if all(s == 1 for s in shared):
        return "Theorem 6.1"
    return "unclassified"


def row_of_theorem_61(G, triples):
    """The row of Theorem 6.1 in which G is listed: the first of abcde, bcde,
    acde (up to symmetry) for which some two twin pairs on distinct twin
    triples abc, ade give a component of G - V(Psi_3) with that
    neighborhood, and "none of these" otherwise. Also returns a witness."""
    found = {}
    for A, B in combinations(list(triples), 2):
        (a,) = set(A) & set(B)
        others = (set(A) | set(B)) - {a}
        for pa in combinations(triples[A], 2):
            for pb in combinations(triples[B], 2):
                psi = set(A) | set(B) | set(pa) | set(pb)
                rest = G.subgraph(set(G.nodes()) - psi)
                for comp in nx.connected_components(rest):
                    nbhd = {u for v in comp for u in G.neighbors(v)} - set(comp)
                    if nbhd == set(A) | set(B):
                        row = "abcde"
                    elif nbhd == others:
                        row = "bcde"
                    elif a in nbhd and len(nbhd & others) == 3 and len(nbhd) == 4:
                        row = "acde"
                    else:
                        continue
                    found.setdefault(row, {
                        "twin triples": [list(A), list(B)],
                        "twin pairs": [list(pa), list(pb)],
                        "component": sorted(comp),
                        "neighborhood": sorted(nbhd)})
    for row in ("abcde", "bcde", "acde"):
        if row in found:
            return row, found[row]
    return "none of these", None


# ----------------------------------------------------------------------------
# K5 minors

def is_k5_model(G, sets):
    """True when the five vertex sets are disjoint, each induces a connected
    subgraph of G, and an edge of G joins every two of them."""
    if len(sets) != 5:
        return False
    allv = [v for s in sets for v in s]
    if len(allv) != len(set(allv)):
        return False
    if not all(s and nx.is_connected(G.subgraph(s)) for s in sets):
        return False
    return all(any(G.has_edge(x, y) for x in a for y in b)
               for a, b in combinations(sets, 2))


def find_k5_minor(G, trials=300, deletion_trials=0, seed=0):
    """Branch sets of a K5 minor of the connected graph G, or None if the
    search finds none. A None answer proves nothing.

    Each trial contracts the edges of G in a random order, skipping edges
    inside one set, until five sets remain, and tests whether every two of
    them are joined by an edge. Each deletion trial first deletes a random
    set of vertices and then does the same in what remains, if connected.
    """
    rng = random.Random(seed)
    sets = _contract_search(G, trials, rng)
    if sets:
        return sets
    nodes = list(G.nodes())
    for _ in range(deletion_trials):
        keep = set(nodes) - set(rng.sample(nodes, rng.randint(1, len(nodes) - 5)))
        H = G.subgraph(keep)
        if nx.is_connected(H):
            sets = _contract_search(H, 1, rng)
            if sets:
                return sets
    return None


def _contract_search(G, trials, rng):
    edges = list(G.edges())
    nodes = list(G.nodes())
    for _ in range(trials):
        parent = {v: v for v in nodes}

        def find(v):
            while parent[v] != v:
                parent[v] = parent[parent[v]]
                v = parent[v]
            return v
        rng.shuffle(edges)
        k = len(nodes)
        for u, v in edges:
            if k == 5:
                break
            a, b = find(u), find(v)
            if a != b:
                parent[a] = b
                k -= 1
        if k != 5:
            continue
        groups = {}
        for v in nodes:
            groups.setdefault(find(v), []).append(v)
        sets = sorted((sorted(s) for s in groups.values()), key=min)
        if is_k5_model(G, sets):
            return sets
    return None


# ----------------------------------------------------------------------------
# No K5 minor: Wagner decomposition

def _pieces_are_minors(G, S, comps):
    """True when every piece G[C + S] + K[S] is a minor of G.

    This holds when |S| <= 1, when |S| = 2 and some other component is
    adjacent to both vertices of S (its path between them gives the edge),
    and when |S| = 3 and at least two other components are adjacent to all
    of S (contracting one onto a vertex of S gives two edges of the
    triangle, and the other gives a path for the third).
    """
    if len(S) <= 1:
        return True
    full = [c for c in comps
            if set(S) <= {u for v in c for u in G.neighbors(v)}]
    need = 1 if len(S) == 2 else 2
    for c in comps:
        if sum(1 for d in full if d is not c) < need:
            return False
    return True


_wagner_cache = {}
PLANARITY_SECONDS = 20


def wagner_decomposition(G):
    """A certificate that G has no K5 minor, or None if none is found.

    The certificate writes a supergraph of G as a clique-sum, along cliques of
    at most three vertices, of planar graphs (each with a plane rotation
    system from PAGE) and copies of V_8. A separator S splits G into the
    pieces G[C + S] plus all edges inside S, one for each component C of
    G - S. A clique-sum of graphs with no K5 minor has no K5 minor, and
    planar graphs and V_8 have none.

    Separators whose pieces are all minors of G are tried first. For those,
    a piece with a K5 minor would give G one, so the first such separator is
    kept and the search does not backtrack over it. Only a graph with no such
    separator is tested for planarity with PAGE or compared with V_8, and
    other separators are tried after that.
    """
    key = (frozenset(G.nodes()), frozenset(frozenset(e) for e in G.edges()))
    if key not in _wagner_cache:
        _wagner_cache[key] = _wagner(G)
    return _wagner_cache[key]


def _wagner(G):
    nodes = sorted(G.nodes())
    if G.number_of_nodes() <= 4:
        return {"type": "planar", "vertices": nodes,
                "embedding": page_embedding(G)[1]}
    if not nx.is_connected(G):
        pieces = [wagner_decomposition(G.subgraph(c).copy())
                  for c in sorted(nx.connected_components(G), key=min)]
        if any(p is None for p in pieces):
            return None
        return {"type": "sum", "separator": [], "pieces": pieces}
    good, other = [], []
    for k in (1, 2, 3):
        for S in combinations(nodes, k):
            rest = G.subgraph(set(nodes) - set(S))
            comps = sorted(nx.connected_components(rest), key=min)
            if len(comps) < 2:
                continue
            (good if _pieces_are_minors(G, S, comps) else other).append((S, comps))

    def split(S, comps):
        pieces = []
        for c in comps:
            P = G.subgraph(set(c) | set(S)).copy()
            P.add_edges_from(combinations(S, 2))
            cert = wagner_decomposition(P)
            if cert is None:
                return None
            pieces.append(cert)
        return {"type": "sum", "separator": list(S), "pieces": pieces}

    if good:
        # Every piece is a minor of G, so if one has no certificate, the
        # search gives up on G rather than trying other separators.
        return split(*good[0])
    m = isomorphism(G, v8())
    if m is not None:
        return {"type": "V8", "vertices": nodes, "map to V8": m}
    try:
        # A planar graph gets its genus-zero embedding from PAGE at once; a
        # PAGE call that runs past the limit only means no certificate here.
        g, emb = page_embedding(G, seconds=PLANARITY_SECONDS)
    except subprocess.TimeoutExpired:
        g = None
    if g == 0:
        return {"type": "planar", "vertices": nodes, "embedding": emb}
    for S, comps in other:
        cert = split(S, comps)
        if cert is not None:
            return cert
    return None


# ----------------------------------------------------------------------------
# All embeddings of a small graph, by plain enumeration of rotation systems

def from_edge_list(edges):
    G = nx.Graph()
    G.add_edges_from(tuple(e) for e in edges)
    return G


def number_of_rotation_systems(G):
    """The product over the vertices of (deg - 1)!."""
    from math import factorial
    out = 1
    for v in G.nodes():
        out *= factorial(max(G.degree(v) - 1, 0))
    return out


def all_rotation_systems(G):
    """Every rotation system of G, each exactly once: at each vertex the
    smallest neighbor comes first and the others follow in every order."""
    from itertools import permutations, product
    verts = sorted(G.nodes())
    choices = []
    for v in verts:
        nb = sorted(G.neighbors(v))
        choices.append([[nb[0]] + list(p) for p in permutations(nb[1:])] if nb else [[]])
    for combo in product(*choices):
        yield dict(zip(verts, combo))


def embeddings_up_to_genus(G, tally, max_genus=1):
    """Run through every rotation system of the connected graph G and yield
    (genus, rotation, faces) for those of genus at most max_genus. tally[0]
    is increased once for every rotation system tried, so that it can be
    compared with number_of_rotation_systems(G)."""
    n, m = G.number_of_nodes(), G.number_of_edges()
    for rot in all_rotation_systems(G):
        tally[0] += 1
        faces = faces_of_rotation(rot)
        g2 = 2 - n + m - len(faces)
        if g2 // 2 <= max_genus:
            yield g2 // 2, rot, faces


def canon_cyclic(seq):
    """A cyclic sequence up to rotation and reversal, as a string."""
    seq = list(seq)
    if not seq:
        return ""
    k = len(seq)
    cands = []
    for s in (seq, seq[::-1]):
        for i in range(k):
            cands.append(tuple(s[i:] + s[:i]))
    return " ".join(min(cands))


def restrict(walk, targets):
    return [v for v in walk if v in targets]


def glued_region(faces, u, v):
    """The region obtained by gluing the faces on the two sides of the edge
    uv along it. Returns a list of boundary walks: one walk for a disk, two
    when both sides are the same face (an annulus)."""
    darts_of = []
    for f in faces:
        darts_of.append([(f[i], f[(i + 1) % len(f)]) for i in range(len(f))])
    i1 = next(k for k, ds in enumerate(darts_of) if (u, v) in ds)
    i2 = next(k for k, ds in enumerate(darts_of) if (v, u) in ds)
    if i1 == i2:
        ds = darts_of[i1]
        a, b = ds.index((u, v)), ds.index((v, u))
        k = len(ds)
        seg1 = [ds[(a + 1 + t) % k] for t in range((b - a - 1) % k)]
        seg2 = [ds[(b + 1 + t) % k] for t in range((a - b - 1) % k)]
        return [[d[0] for d in seg1], [d[0] for d in seg2]], [i1]
    d1, d2 = darts_of[i1], darts_of[i2]
    a, b = d1.index((u, v)), d2.index((v, u))
    seg1 = [d1[(a + 1 + t) % len(d1)] for t in range(len(d1) - 1)]
    seg2 = [d2[(b + 1 + t) % len(d2)] for t in range(len(d2) - 1)]
    return [[d[0] for d in seg1] + [d[0] for d in seg2]], [i1, i2]


def automorphisms(G):
    matcher = nx.algorithms.isomorphism.GraphMatcher(G, G)
    return [dict(m) for m in matcher.isomorphisms_iter()]


def rotation_key(rot, verts):
    """A rotation system written canonically: each cyclic order started at
    its smallest neighbor."""
    out = []
    for v in verts:
        r = rot[v]
        if not r:
            out.append(())
            continue
        i = r.index(min(r))
        out.append(tuple(r[i:] + r[:i]))
    return tuple(out)


def embedding_class(rot, auts, verts):
    """The smallest key among the images of the rotation system under every
    automorphism, with and without reversing every cyclic order."""
    best = None
    for phi in auts:
        for mirror in (False, True):
            img = {}
            for v in verts:
                r = [phi[u] for u in rot[v]]
                img[phi[v]] = r[::-1] if mirror else r
            k = rotation_key(img, verts)
            if best is None or k < best:
                best = k
    return best


def k5_minor_by_partitions(G):
    """Branch sets of a K5 minor of the small connected graph G, or None if
    there is none, found by trying every partition of the vertices into five
    blocks. This decides the question: the unused vertices of a K5 model of
    a connected graph can always be added to an adjacent branch set, so a
    K5 minor exists exactly when some partition into five blocks has every
    block connected and every two blocks adjacent."""
    V = sorted(G.nodes())
    n = len(V)
    idx = {v: i for i, v in enumerate(V)}
    adj = [set(idx[u] for u in G[v]) for v in V]
    assign = [0] * n
    found = []

    def connected(block):
        block = set(block)
        start = next(iter(block))
        seen, stack = {start}, [start]
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y in block and y not in seen:
                    seen.add(y)
                    stack.append(y)
        return len(seen) == len(block)

    def rec(i, k):
        if found or n - i < 5 - k:
            return
        if i == n:
            if k != 5:
                return
            blocks = [[j for j in range(n) if assign[j] == b] for b in range(5)]
            if all(connected(b) for b in blocks) and all(
                    any(y in adj[x] for x in blocks[s] for y in blocks[t])
                    for s, t in combinations(range(5), 2)):
                found.append(sorted((sorted(V[j] for j in b) for b in blocks), key=min))
            return
        for b in range(min(k + 1, 5)):
            assign[i] = b
            rec(i + 1, max(k, b + 1))
    rec(0, 0)
    return found[0] if found else None
