"""The kinds of witness used in the supplemental data, and their checks.

A graph is a dictionary mapping each vertex to the set of its neighbours.
Every check raises CheckFailed (a subclass of AssertionError) with a message
when the witness is wrong, and returns normally when it is right.

The seven kinds (Appendix A of the paper and README.md of the data folder):
  embedding            a rotation system; faces traced, Euler's formula
  K5 branch sets       five disjoint connected sets, an edge between every two
  isomorphism          a bijection carrying edges exactly onto edges
  minor model          a connected set for each vertex of the target graph,
                       an edge between the sets of adjacent vertices
  cut                  at most two vertices whose deletion disconnects
  Kuratowski subgraph  the edges of a subdivision of K5 or K3,3
  no K5 certificate    a tree of clique-sums over at most three vertices whose
                       leaves are plane rotation systems or isomorphisms onto
                       V8 (Wagner's theorem)
"""
import gzip
import json
import os
from itertools import combinations


class CheckFailed(AssertionError):
    pass


def need(cond, *msg):
    if not cond:
        raise CheckFailed(' '.join(str(m) for m in msg) if msg else 'check failed')


# ---------------------------------------------------------------- files and graphs

def load(path):
    """Read a JSON file, or its .json.gz form."""
    if not os.path.exists(path) and os.path.exists(path + '.gz'):
        path += '.gz'
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as f:
        return json.load(f)


def find(path):
    """The path of a file stored either plain or gzip-compressed, or None."""
    for p in (path, path + '.gz'):
        if os.path.exists(p):
            return p
    return None


def graph(g):
    """{"vertices": [...], "edges": [[u, v], ...]} as an adjacency dictionary."""
    adj = {v: set() for v in g['vertices']}
    need(len(adj) == len(g['vertices']), 'repeated vertex')
    for u, v in g['edges']:
        need(u != v, 'loop at', u)
        need(u in adj and v in adj, 'edge with an unknown end', u, v)
        need(v not in adj[u], 'repeated edge', u, v)
        adj[u].add(v)
        adj[v].add(u)
    return adj


def copy(adj):
    return {v: set(s) for v, s in adj.items()}


def add_edges(adj, es):
    H = copy(adj)
    for u, v in es:
        need(u != v, 'loop at', u)
        H.setdefault(u, set()).add(v)
        H.setdefault(v, set()).add(u)
    return H


def edges(adj):
    return sorted((u, v) if str(u) <= str(v) else (v, u) for u in adj for v in adj[u] if str(u) < str(v))


def n_edges(adj):
    return sum(len(s) for s in adj.values()) // 2


def same_graph(a, b):
    return set(a) == set(b) and all(a[v] == b[v] for v in a)


def components(adj, vertices=None):
    vs = set(adj) if vertices is None else set(vertices)
    seen, comps = set(), []
    for s in vs:
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
        comps.append(comp)
    return comps


def delete_edge(adj, e):
    u, v = e
    need(u in adj and v in adj[u], 'the edge', u, v, 'is not in the graph')
    H = copy(adj)
    H[u].discard(v)
    H[v].discard(u)
    return H


def contract_edge(adj, e):
    """Contract uv: v is merged into u; loops are deleted and parallel edges merged."""
    u, v = e
    need(u in adj and v in adj[u], 'the edge', u, v, 'is not in the graph')
    H = copy(adj)
    for w in H.pop(v):
        H[w].discard(v)
        if w != u:
            H[w].add(u)
            H[u].add(w)
    return H


def delete_vertices(adj, vs):
    vs = set(vs)
    return {v: s - vs for v, s in adj.items() if v not in vs}


def from_graph6(s):
    s = s.strip()
    data = [ord(c) - 63 for c in s]
    need(data[0] < 63, 'graph6 with more than 62 vertices')
    n = data[0]
    bits = []
    for x in data[1:]:
        bits.extend((x >> k) & 1 for k in range(5, -1, -1))
    adj = {v: set() for v in range(n)}
    k = 0
    for j in range(1, n):
        for i in range(j):
            if bits[k]:
                adj[i].add(j)
                adj[j].add(i)
            k += 1
    return adj


def v8():
    """The Wagner graph V8 on 0, ..., 7: the 8-cycle with its four long diagonals."""
    adj = {v: set() for v in range(8)}
    for i in range(8):
        for j in ((i + 1) % 8, (i + 4) % 8):
            adj[i].add(j)
            adj[j].add(i)
    return adj


def k5():
    return {v: set(range(5)) - {v} for v in range(5)}


# ---------------------------------------------------------------- the seven kinds

def vertex_of(adj, x):
    """JSON object keys are strings; a graph whose vertices are integers (as
    T1, ..., T74 in O1) gets them back as integers."""
    if x in adj:
        return x
    if isinstance(x, str) and x.lstrip('-').isdigit() and int(x) in adj:
        return int(x)
    return x


def genus_of_rotation(adj, rot):
    """Genus of the embedding given by rot (vertex -> neighbours in cyclic
    order), summed over the components. After the dart u->v comes v->w, where
    w follows u in the cyclic order at v."""
    rot = {vertex_of(adj, k): [vertex_of(adj, x) for x in nb] for k, nb in rot.items()}
    need(set(rot) == set(adj), 'the rotation system does not list exactly the vertices of the graph')
    for v in adj:
        need(len(rot[v]) == len(set(rot[v])) and set(rot[v]) == adj[v], 'the rotation at', v,
             'is not a cyclic order of its neighbours')
    succ = {}
    for v, nb in rot.items():
        for i, u in enumerate(nb):
            succ[(v, u)] = nb[(i + 1) % len(nb)]
    total = 0
    for comp in components(adj):
        darts = {(u, v) for u in comp for v in adj[u]}
        m = len(darts) // 2
        if m == 0:
            continue
        faces = 0
        while darts:
            start = next(iter(darts))
            d = start
            faces += 1
            while True:
                darts.discard(d)
                d = (d[1], succ[(d[1], d[0])])
                if d == start:
                    break
        two_g = 2 - len(comp) + m - faces
        need(two_g >= 0 and two_g % 2 == 0, 'Euler characteristic of a component is not even')
        total += two_g // 2
    return total


def faces_of_rotation(adj, rot):
    """The facial walks of a rotation system, each a list of vertices."""
    genus_of_rotation(adj, rot)
    rot = {vertex_of(adj, k): [vertex_of(adj, x) for x in nb] for k, nb in rot.items()}
    succ = {}
    for v, nb in rot.items():
        for i, u in enumerate(nb):
            succ[(v, u)] = nb[(i + 1) % len(nb)]
    darts = {(u, v) for u in adj for v in adj[u]}
    out = []
    while darts:
        start = min(darts, key=lambda d: (str(d[0]), str(d[1])))
        d, walk = start, []
        while True:
            darts.discard(d)
            walk.append(d[0])
            d = (d[1], succ[(d[1], d[0])])
            if d == start:
                break
        out.append(walk)
    return out


def check_embedding(adj, w, genus=None, at_most=None):
    need(w.get('type') == 'embedding', 'not an embedding witness')
    g = genus_of_rotation(adj, w['rotation'])
    need(g == w['genus'], 'stated genus', w['genus'], 'but the faces give', g)
    if genus is not None:
        need(g == genus, 'genus', g, 'where', genus, 'is claimed')
    if at_most is not None:
        need(g <= at_most, 'genus', g, 'above', at_most)
    return g


def check_k5(adj, sets):
    sets = [list(s) for s in sets]
    allv = [v for s in sets for v in s]
    need(len(sets) == 5 and all(sets), 'K5 needs five nonempty branch sets')
    need(len(allv) == len(set(allv)), 'K5 branch sets meet')
    need(set(allv) <= set(adj), 'K5 branch set vertex not in the graph')
    for s in sets:
        need(len(components(adj, s)) == 1, 'K5 branch set not connected', s)
    for a, b in combinations(sets, 2):
        need(any(y in adj[x] for x in a for y in b), 'no edge between K5 branch sets', a, b)


def check_cut(adj, cut):
    need(1 <= len(cut) <= 2 and len(set(cut)) == len(cut) and set(cut) <= set(adj), 'not a set of one or two vertices')
    need(len(components(adj, set(adj) - set(cut))) >= 2, 'deleting', cut, 'does not disconnect the graph')


def check_iso(adj, target, m):
    """m maps the vertices of adj bijectively onto those of target, carrying
    edges exactly onto edges."""
    need(set(m) == set(adj), 'the map is not defined on exactly the vertices')
    need(len(set(m.values())) == len(m) and set(m.values()) == set(target), 'the map is not onto the target vertices')
    need(n_edges(adj) == n_edges(target), 'different numbers of edges')
    for u in adj:
        for v in adj[u]:
            need(m[v] in target[m[u]], 'the edge', u, v, 'is not carried onto an edge')


def check_minor_model(adj, target, blocks):
    """blocks: vertex of target -> list of vertices of adj."""
    need(set(blocks) == set(target), 'a branch set for each vertex of the target is needed')
    allv = [v for b in blocks.values() for v in b]
    need(len(allv) == len(set(allv)), 'branch sets meet')
    need(set(allv) <= set(adj), 'branch set vertex not in the graph')
    for x, b in blocks.items():
        need(b and len(components(adj, b)) == 1, 'branch set of', x, 'empty or not connected')
    for u in target:
        for v in target[u]:
            need(any(y in adj[x] for x in blocks[u] for y in blocks[v]), 'no edge between the branch sets of', u, v)


def check_kuratowski(adj, es):
    """The edges es lie in adj and form a subdivision of K5 or K3,3.
    Returns 'K5' or 'K3,3'."""
    sub = {}
    for u, v in es:
        need(u in adj and v in adj[u], 'Kuratowski edge', u, v, 'not in the graph')
        sub.setdefault(u, set()).add(v)
        sub.setdefault(v, set()).add(u)
    need(n_edges(sub) == len(es), 'repeated Kuratowski edge')
    need(len(components(sub)) == 1, 'Kuratowski subgraph not connected')
    branch = [v for v in sub if len(sub[v]) != 2]
    need(all(len(sub[v]) in (3, 4) for v in branch), 'a vertex of degree other than 2, 3, 4')
    links = {}
    for b in branch:
        for nxt in sub[b]:
            prev, cur = b, nxt
            while cur not in branch:
                (a,) = [x for x in sub[cur] if x != prev]
                prev, cur = cur, a
            k = frozenset((b, cur))
            need(len(k) == 2, 'a path returns to its start')
            links[k] = links.get(k, 0) + 1
    need(all(c == 2 for c in links.values()), 'two paths join the same branch vertices')
    if len(branch) == 5:
        need(all(len(sub[v]) == 4 for v in branch) and len(links) == 10, 'not a subdivision of K5')
        return 'K5'
    need(len(branch) == 6 and all(len(sub[v]) == 3 for v in branch) and len(links) == 9,
         'not a subdivision of K5 or K3,3')
    A = {branch[0]} | {v for v in branch if v != branch[0] and frozenset((branch[0], v)) not in links}
    B = set(branch) - A
    need(len(A) == 3 and all(frozenset((a, b)) in links for a in A for b in B), 'not a subdivision of K3,3')
    return 'K3,3'


def check_no_k5(adj, c):
    """A certificate that adj has no K5 minor: adj is a subgraph of a gluing of
    planar graphs and copies of V8 along cliques of order at most three."""
    t = c['type']
    if t == 'planar':
        need(genus_of_rotation(adj, c['rotation']) == 0, 'the rotation system of a planar part is not plane')
    elif t == 'V8':
        m = {vertex_of(adj, k): int(v) for k, v in c['isomorphism'].items()}
        check_iso(adj, v8(), m)
    elif t == 'clique-sum':
        S = [vertex_of(adj, x) for x in c['separator']]
        need(len(S) <= 3 and set(S) <= set(adj) and len(set(S)) == len(S), 'bad separator', S)
        parts = c['parts']
        need(len(parts) >= 2, 'a clique-sum needs two parts')
        vsets = [{vertex_of(adj, x) for x in p['vertices']} for p in parts]
        need(set().union(*vsets) == set(adj), 'the parts do not cover the graph')
        for a, b in combinations(vsets, 2):
            need(a & b == set(S), 'two parts meet outside the separator')
        owner = {}
        for i, vs in enumerate(vsets):
            need(vs - set(S), 'a part lies inside the separator')
            for v in vs - set(S):
                owner[v] = i
        for u in adj:
            for v in adj[u]:
                if u in owner and v in owner:
                    need(owner[u] == owner[v], 'an edge joins two parts', u, v)
        for p, vs in zip(parts, vsets):
            P = {v: adj[v] & vs for v in vs}
            for x, y in combinations(S, 2):
                P[x].add(y)
                P[y].add(x)
            check_no_k5(P, p['certificate'])
    else:
        raise CheckFailed('unknown kind of no-K5 certificate: %s' % t)


def as_vertex_of(T, x):
    """JSON object keys are strings; the vertices of T1, ..., T74 are integers."""
    if x in T:
        return x
    if isinstance(x, str) and x.isdigit() and int(x) in T:
        return int(x)
    return x


def check_witness(adj, w, targets):
    """Check a witness of one of the seven kinds against the graph adj.
    targets: name -> graph, for isomorphisms and minor models."""
    t = w.get('type')
    if t == 'embedding':
        return check_embedding(adj, w)
    if t == 'K5 branch sets':
        return check_k5(adj, w['branch sets'])
    if t == 'isomorphism':
        need(w['onto'] in targets, 'isomorphism onto an unknown graph', w['onto'])
        T = targets[w['onto']]
        return check_iso(adj, T, {k: as_vertex_of(T, x) for k, x in w['map'].items()})
    if t == 'minor model':
        need(w['of'] in targets, 'minor model of an unknown graph', w['of'])
        T = targets[w['of']]
        return check_minor_model(adj, T, {as_vertex_of(T, k): v for k, v in w['branch sets'].items()})
    if t == 'cut':
        return check_cut(adj, w['vertices'])
    if t == 'Kuratowski subgraph':
        return check_kuratowski(adj, w['edges'])
    if t == 'clique-sum' or t == 'planar' or t == 'V8':
        return check_no_k5(adj, w)
    raise CheckFailed('unknown kind of witness: %s' % t)


# ---------------------------------------------------------------- the obstructions

def catalogue(root):
    """T1, ..., T74 from item O1 (vertices 0, ..., n-1 in graph6 order)."""
    d = load(os.path.join(root, 'no_K5_obstructions', 'O1.json'))
    out = {}
    for r in d['graphs']:
        g = from_graph6(r['graph6'])
        need(len(g) == r['n'] and n_edges(g) == r['m'], 'O1: wrong size of', r['name'])
        if 'edges' in r:
            need(same_graph(g, graph({'vertices': list(range(r['n'])), 'edges': r['edges']})),
                 'O1: graph6 and edges of', r['name'], 'differ')
        out[r['name']] = g
    return out
