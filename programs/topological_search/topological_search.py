"""The search in the proof of Theorem 8.1: the topological obstructions for the
torus with no K5 minor.

Start from T1, ..., T74 (item O1). Split every vertex of every graph found so
far in every way: the vertex v keeps its first neighbour and at least one
more, a new vertex gets at least two of the others, and v and the new vertex
are joined. Each split, up to isomorphism, gets one verdict:

  K5 minor       five branch sets (found by a SAT solver, which also shows
                 that the other graphs have no K5 minor);
  an edge deletion has a T_i minor
                 some G - e is nontoroidal; the program reduces G - e to a
                 minor that is a minor-order obstruction, which has no K5
                 minor and so is some T_i, and records its branch sets;
  kept           G has no K5 minor (a clique-sum of planar graphs and copies
                 of V8 along cliques of at most three vertices) and every
                 G - e is toroidal (a rotation system of genus at most one).
                 A kept graph is nontoroidal, since contracting its new edge
                 gives the graph it was split from, so it is a topological
                 obstruction; it is split in turn.

The search stops when a round keeps no new graph. TO1.json lists T1, ..., T74
and the kept graphs, TO2.json every split with its verdict, in the format of
the data folder. Whether a graph is toroidal is decided with multi_genus
(option l1, which stops once the genus is known to exceed one), and the torus
embeddings are found with PAGE, both from the package graph-genus; the
witnesses written are checked by the checker without any genus computation.

Usage:
    python topological_search.py ../../supplemental_data OUTPUT_FOLDER [--workers N]

The graphs T1, ..., T74 are read from no_K5_obstructions/O1.json of the data
folder. Then compare the result with the data folder:
    python compare.py ../../supplemental_data OUTPUT_FOLDER
"""
import argparse
import json
import os
import platform
import subprocess
import sys
import time
from collections import Counter
from itertools import combinations
from multiprocessing import Pool

import networkx as nx
import graph_genus


# ---------------------------------------------------------------- graphs as {vertex: set of neighbours}

def from_graph6(s):
    data = [ord(c) - 63 for c in s.strip()]
    n = data[0]
    bits = [(x >> k) & 1 for x in data[1:] for k in range(5, -1, -1)]
    adj = {v: set() for v in range(n)}
    k = 0
    for j in range(1, n):
        for i in range(j):
            if bits[k]:
                adj[i].add(j)
                adj[j].add(i)
            k += 1
    return adj


def to_graph6(adj):
    order = sorted(adj)
    n = len(order)
    bits = [1 if order[i] in adj[order[j]] else 0 for j in range(1, n) for i in range(j)]
    while len(bits) % 6:
        bits.append(0)
    return chr(n + 63) + ''.join(chr(int(''.join(map(str, bits[k:k + 6])), 2) + 63) for k in range(0, len(bits), 6))


def edges(adj):
    return sorted((u, v) for u in adj for v in adj[u] if u < v)


def n_edges(adj):
    return sum(len(s) for s in adj.values()) // 2


def copy(adj):
    return {v: set(s) for v, s in adj.items()}


def delete_edge(adj, e):
    u, v = e
    H = copy(adj)
    H[u].discard(v)
    H[v].discard(u)
    return H


def contract_edge(adj, e):
    """v merged into u; loops dropped, parallel edges merged."""
    u, v = e
    H = copy(adj)
    for w in H.pop(v):
        H[w].discard(v)
        if w != u:
            H[w].add(u)
            H[u].add(w)
    return H


def components(adj, vertices=None):
    vs = set(adj) if vertices is None else set(vertices)
    seen, out = set(), []
    for s in sorted(vs):
        if s in seen:
            continue
        comp, stack = {s}, [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y in vs and y not in seen:
                    seen.add(y)
                    comp.add(y)
                    stack.append(y)
        out.append(comp)
    return out


def nxg(adj):
    G = nx.Graph()
    G.add_nodes_from(adj)
    G.add_edges_from(edges(adj))
    return G


def v8():
    adj = {v: set() for v in range(8)}
    for i in range(8):
        for j in ((i + 1) % 8, (i + 4) % 8):
            adj[i].add(j)
            adj[j].add(i)
    return adj


class IsoIndex:
    """Graphs up to isomorphism: Weisfeiler-Lehman hash buckets, then VF2."""
    def __init__(self):
        self.buckets = {}

    def find(self, adj):
        G = nxg(adj)
        for key, H in self.buckets.get(nx.weisfeiler_lehman_graph_hash(G, iterations=3), []):
            if nx.is_isomorphic(G, H):
                return key
        return None

    def add(self, adj, key):
        G = nxg(adj)
        self.buckets.setdefault(nx.weisfeiler_lehman_graph_hash(G, iterations=3), []).append((key, G))


# ---------------------------------------------------------------- splits

def splits(g):
    """Every split of every vertex of g, in a fixed order, with the split data."""
    n = len(g)
    for v in sorted(g):
        nb = sorted(g[v])
        d = len(nb)
        for mask in range(1 << d):
            if not mask & 1 or not 2 <= bin(mask).count('1') <= d - 2:
                continue
            h = {x: set() for x in range(n + 1)}
            for x, y in edges(g):
                if v in (x, y):
                    w = y if x == v else x
                    t = v if mask & (1 << nb.index(w)) else n
                    h[t].add(w)
                    h[w].add(t)
                else:
                    h[x].add(y)
                    h[y].add(x)
            h[v].add(n)
            h[n].add(v)
            yield h, {'vertex': v, 'new vertex': n, 'neighbours moved to the new vertex': sorted(h[n] - {v})}


# ---------------------------------------------------------------- embeddings

def embedding(adj):
    """(genus, rotation) of a minimum-genus embedding. Vertices of degree at
    most one are removed (repeatedly) and put back afterwards, which does not
    change the genus; each remaining component is embedded by PAGE, and the
    genus is the sum over the components."""
    H = copy(adj)
    stripped = []
    queue = [v for v in sorted(H) if len(H[v]) <= 1]
    while queue:
        v = queue.pop()
        if v not in H or len(H[v]) > 1:
            continue
        nb = next(iter(H[v])) if H[v] else None
        stripped.append((v, nb))
        del H[v]
        if nb is not None:
            H[nb].discard(v)
            if len(H[nb]) <= 1:
                queue.append(nb)
    rot, genus = {}, 0
    for comp in components(H):
        order = sorted(comp)
        idx = {x: i for i, x in enumerate(order)}
        g, r = graph_genus.embed([[idx[y] for y in sorted(H[x])] for x in order])
        genus += g
        for i, nb in enumerate(r):
            rot[order[i]] = [order[j] for j in nb]
    for v, nb in reversed(stripped):
        rot[v] = [nb] if nb is not None else []
        if nb is not None:
            rot[nb] = rot.get(nb, []) + [v]
    return genus, rot


def multicode(adj):
    vs = sorted(adj)
    idx = {v: i + 1 for i, v in enumerate(vs)}
    out = [len(vs)]
    for v in vs[:-1]:
        out += sorted(idx[u] for u in adj[v] if idx[u] > idx[v])
        out.append(0)
    return bytes(out)


def multi_genus_program():
    p = os.path.join(os.path.dirname(graph_genus.__file__), 'bin',
                     'multi_genus' + ('.exe' if platform.system() == 'Windows' else ''))
    return os.environ.get('MULTI_GENUS_BIN', p)


def genus_above_one(graphs):
    """For each connected graph with at least two vertices, whether its genus
    is larger than one: multi_genus with the options l1 F L writes a byte 255
    when it starts each graph, followed by the graph itself when its genus is
    larger than one."""
    if not graphs:
        return []
    data = b''.join(multicode(G) for G in graphs)
    r = subprocess.run([multi_genus_program(), 'l1', 'F', 'L'], input=data, capture_output=True, check=True)
    chunks = r.stdout.split(b'\xff')
    assert chunks[0] == b'' and len(chunks) == len(graphs) + 1, 'unexpected output of multi_genus'
    return [bool(c) for c in chunks[1:]]


def nontoroidal(graphs):
    """For each graph, whether it is nontoroidal: two nonplanar components, or
    one component of genus larger than one."""
    out = [False] * len(graphs)
    jobs = []
    for i, G in enumerate(graphs):
        nonplanar = []
        for c in components(G):
            if len(c) >= 5:
                C = {v: G[v] & c for v in c}
                if not nx.check_planarity(nxg(C))[0]:
                    nonplanar.append(C)
        if len(nonplanar) >= 2:
            out[i] = True
        elif nonplanar:
            jobs.append((i, nonplanar[0]))
    for (i, _), f in zip(jobs, genus_above_one([C for _, C in jobs])):
        out[i] = f
    return out


def plane_rotation(adj):
    """A plane rotation system of a planar graph, or None."""
    ok, emb = nx.check_planarity(nxg(adj))
    if not ok:
        return None
    return {v: list(emb.neighbors_cw_order(v)) if adj[v] else [] for v in adj}


# ---------------------------------------------------------------- K5 minors

def is_k5_model(adj, sets):
    allv = [v for s in sets for v in s]
    if len(sets) != 5 or len(allv) != len(set(allv)) or not all(sets):
        return False
    if not all(len(components(adj, s)) == 1 for s in sets):
        return False
    return all(any(y in adj[x] for x in a for y in b) for a, b in combinations(sets, 2))


def k5_minor(adj):
    """Branch sets of a K5 minor, or None if the graph has none.

    In a connected graph, a K5 model extends to a partition of all the
    vertices into five connected, pairwise adjacent sets (an unused vertex
    joins a set it has a path to). For each component, a SAT solver looks for
    such a partition: every vertex in exactly one set, the sets numbered by
    their least vertex, an edge between every two sets. Connectivity is added
    as it is needed: when a set of a solution falls into pieces, then for each
    piece P, a vertex of P and a vertex of the set outside P may share the set
    only together with some neighbour of P. Every connected set satisfies these
    clauses, so an unsatisfiable formula shows that the component has no K5
    minor."""
    from pysat.solvers import Solver
    for comp in components(adj):
        if len(comp) < 5:
            continue
        V = sorted(comp)
        n = len(V)
        idx = {v: i for i, v in enumerate(V)}

        def x(i, c):
            return 1 + 5 * i + c
        top = 5 * n
        clauses = []
        for i in range(n):
            clauses.append([x(i, c) for c in range(5)])
            clauses += [[-x(i, a), -x(i, b)] for a, b in combinations(range(5), 2)]
            for c in range(1, 5):
                clauses.append([-x(i, c)] + [x(j, c - 1) for j in range(i)])
        for c in range(5):
            clauses.append([x(i, c) for i in range(n)])
        E = [(idx[u], idx[v]) for u, v in edges(adj) if u in comp]
        for a, b in combinations(range(5), 2):
            ws = []
            for p, q in E:
                for s_, t_ in ((p, q), (q, p)):
                    top += 1
                    ws.append(top)
                    clauses += [[-top, x(s_, a)], [-top, x(t_, b)]]
            clauses.append(ws)
        with Solver(name='g3', bootstrap_with=clauses) as S:
            while S.solve():
                m = set(l for l in S.get_model() if l > 0)
                bags = [[V[i] for i in range(n) if x(i, c) in m] for c in range(5)]
                cuts = []
                for c, bag in enumerate(bags):
                    pieces = components(adj, bag)
                    if len(pieces) <= 1:
                        continue
                    for piece in pieces:
                        first = min(piece)
                        other = min(set(bag) - piece)
                        boundary = {w for v in piece for w in adj[v] if w in comp} - piece
                        cuts.append([-x(idx[first], c), -x(idx[other], c)] + [x(idx[w], c) for w in sorted(boundary)])
                if not cuts:
                    assert is_k5_model(adj, bags)
                    return sorted(bags, key=min)
                for cl in cuts:
                    S.add_clause(cl)
    return None


def pieces_are_minors(adj, S, comps):
    """Each piece G[C + S] + K[S] is a minor of G: |S| <= 1, or |S| = 2 and
    another component meets both vertices of S, or |S| = 3 and two other
    components meet all of S."""
    if len(S) <= 1:
        return True
    full = [c for c in comps if set(S) <= {u for v in c for u in adj[v]}]
    need = 1 if len(S) == 2 else 2
    return all(sum(1 for d in full if d is not c) >= need for c in comps)


def no_k5_certificate(adj):
    """A clique-sum of planar graphs and copies of V8, along cliques of at most
    three vertices, containing adj; None if none is found. Separators whose
    pieces are minors of the graph are tried first."""
    rot = plane_rotation(adj)
    if rot is not None:
        return {'type': 'planar', 'rotation': {str(k): v for k, v in rot.items()}}
    if len(adj) == 8 and n_edges(adj) == 12:
        gm = nx.algorithms.isomorphism.GraphMatcher(nxg(adj), nxg(v8()))
        if gm.is_isomorphic():
            return {'type': 'V8', 'isomorphism': {str(k): v for k, v in gm.mapping.items()}}
    V = sorted(adj)
    good, other = [], []
    for k in range(0, 4):
        for S in combinations(V, k):
            comps = components(adj, [v for v in V if v not in S])
            if len(comps) < 2:
                continue
            (good if pieces_are_minors(adj, S, comps) else other).append((S, comps))

    def split(S, comps):
        parts = []
        for C in comps:
            vs = set(C) | set(S)
            P = {v: adj[v] & vs for v in vs}
            for x, y in combinations(S, 2):
                P[x].add(y)
                P[y].add(x)
            c = no_k5_certificate(P)
            if c is None:
                return None
            parts.append({'vertices': sorted(vs), 'certificate': c})
        return {'type': 'clique-sum', 'separator': list(S), 'parts': parts}

    if good:
        # its pieces are minors of the graph, so a piece without a certificate
        # means that no certificate is found for the graph either
        return split(*good[0])
    for S, comps in other:
        c = split(S, comps)
        if c is not None:
            return c
    return None


# ---------------------------------------------------------------- a minor T_i of a nontoroidal graph

CAT = None        # name -> graph, set in each worker
CAT_INDEX = None


def load_catalogue(root):
    d = json.load(open(os.path.join(root, 'no_K5_obstructions', 'O1.json')))
    return {r['name']: from_graph6(r['graph6']) for r in d['graphs']}


def init_worker(cat):
    global CAT, CAT_INDEX
    CAT = cat
    CAT_INDEX = IsoIndex()
    for name, g in cat.items():
        CAT_INDEX.add(g, name)


def minor_of_catalogue(H):
    """H is nontoroidal with no K5 minor. Delete or contract edges, one at a
    time, while the graph stays nontoroidal (the first operation that keeps it
    so, deletions before contractions, in the order of the edges); the result
    is a minor-order obstruction with no K5 minor, hence some T_i. Returns its
    name and the branch sets of the model in H (vertex numbers of T_i as in
    O1)."""
    G = {v: s for v, s in H.items() if s}
    blocks = {v: {v} for v in G}
    while True:
        es = edges(G)
        tries = [('delete', e, {v: s for v, s in delete_edge(G, e).items() if s}) for e in es] + \
                [('contract', e, contract_edge(G, e)) for e in es]
        flags = nontoroidal([t[2] for t in tries])
        hit = next((t for t, f in zip(tries, flags) if f), None)
        if hit is None:
            break
        op, (u, v), G2 = hit
        if op == 'delete':
            for x in set(G) - set(G2):
                del blocks[x]
        else:
            blocks[u] |= blocks.pop(v)
        G = G2
    name = CAT_INDEX.find(G)
    if name is None:
        raise RuntimeError('a minimal nontoroidal minor is not among T1, ..., T74')
    gm = nx.algorithms.isomorphism.GraphMatcher(nxg(G), nxg(CAT[name]))
    assert gm.is_isomorphic()
    model = {str(gm.mapping[v]): sorted(blocks[v]) for v in G}
    return name, dict(sorted(model.items(), key=lambda kv: int(kv[0])))


# ---------------------------------------------------------------- one verdict

def classify(h):
    sets = k5_minor(h)
    if sets is not None:
        return 'K5 minor', {'type': 'K5 branch sets', 'branch sets': sets}
    es = edges(h)
    flags = nontoroidal([delete_edge(h, e) for e in es])
    for e, f in zip(es, flags):
        if f:
            name, model = minor_of_catalogue(delete_edge(h, e))
            return 'an edge deletion has a T_i minor', {'type': 'minor model', 'edge deleted': list(e), 'of': name,
                                                        'branch sets': model}
    cert = no_k5_certificate(h)
    if cert is None:
        raise RuntimeError('no certificate found for a graph with no K5 minor: %s' % to_graph6(h))
    dels = []
    for e in es:
        g, rot = embedding(delete_edge(h, e))
        if g > 1:
            raise RuntimeError('PAGE and multi_genus disagree on %s' % to_graph6(h))
        dels.append({'edge': list(e), 'genus': g, 'rotation': {str(k): v for k, v in rot.items()}})
    return 'kept', {'no K5 certificate': cert, 'torus embeddings of every G - e': dels}


def classify_g6(s):
    return classify(from_graph6(s))


# ---------------------------------------------------------------- the search

WHAT1 = ('The 381 topological obstructions for the torus with no K5 minor (Theorem 8.1), as graph6 strings. The 74 '
         'that are T1,...,T74 are named.')
WHAT2 = ('The search in the proof of Theorem 8.1. Every split of every graph in TO1 is one of these 7,006 graphs up to '
         'isomorphism, and each of them is a split of a graph in TO1. Each graph is given as produced by its split: '
         'the vertex split, the new vertex, and the neighbours moved to it; contracting the edge between them gives '
         'the graph it was split from. Verdicts: kept (a certificate that it has no K5 minor, made of planar pieces '
         'and copies of V8 glued along cliques of order at most three, and a torus embedding of every edge deletion), '
         'K5 minor (branch sets), or an edge deletion containing some T_i as a minor (the edge and a minor model, '
         'vertex numbers of T_i as in O1). No genus computation is used.')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('data', help='the supplemental data folder (for T1, ..., T74 in O1)')
    ap.add_argument('out', help='folder to write TO1.json and TO2.json')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    cat = load_catalogue(a.data)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    index = IsoIndex()
    records = []
    listed = [(name, g) for name, g in cat.items()]          # T1..T74, then the kept graphs
    n_splits = 0
    frontier = list(listed)
    with Pool(a.workers, initializer=init_worker, initargs=(cat,)) as pool:
        rnd = 0
        while frontier:
            rnd += 1
            new = []
            for name, g in frontier:
                for h, sp in splits(g):
                    n_splits += 1
                    if index.find(h) is not None:
                        continue
                    k = len(records) + len(new)
                    index.add(h, k)
                    parent = name if name in cat else {'graph6': name}
                    new.append({'graph6': to_graph6(h), 'split from': parent, 'split': sp})
            verdicts = pool.map(classify_g6, [r['graph6'] for r in new], chunksize=4)
            frontier = []
            for r, (v, w) in zip(new, verdicts):
                r['verdict'], r['witness'] = v, w
                if v == 'kept':
                    frontier.append((r['graph6'], from_graph6(r['graph6'])))
            records += new
            listed += frontier
            print('round %d: %d new graphs, %d kept; %d graphs so far, %.0f s' % (
                rnd, len(new), len(frontier), len(records), time.time() - t0), flush=True)
    counts = Counter(r['verdict'] for r in records)
    orders = Counter(len(g) for _, g in listed)
    to1 = {'item': 'TO1', 'what': WHAT1, 'orders': {str(k): v for k, v in sorted(orders.items())},
           'graphs': [{'graph6': to_graph6(g), 'n': len(g), 'm': n_edges(g), **({'is': name} if name in cat else {})}
                      for name, g in sorted(listed, key=lambda x: (len(x[1]), n_edges(x[1])))]}
    to2 = {'item': 'TO2', 'what': WHAT2,
           'counts': {'graphs': len(records), 'splits of the %d graphs' % len(listed): n_splits,
                      'kept': counts['kept'], 'K5 minor': counts['K5 minor'],
                      'edge deletion with a T_i minor': counts['an edge deletion has a T_i minor']},
           'graphs': records}
    for name, d in (('TO1', to1), ('TO2', to2)):
        with open(os.path.join(a.out, name + '.json'), 'w') as f:
            json.dump(d, f, indent=1)
    print('%d topological obstructions; %d splits; %d graphs up to isomorphism: %s; %.0f s' % (
        len(listed), n_splits, len(records), dict(counts), time.time() - t0))


if __name__ == '__main__':
    main()
