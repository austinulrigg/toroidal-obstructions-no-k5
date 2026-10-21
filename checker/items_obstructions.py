"""The obstructions and the topological search: O1-O4, TO1 and TO2.

O1  T1, ..., T74: graph6 and edge lists agree; the 74 graphs are pairwise
    nonisomorphic (checked directly, by isomorphism tests within classes of
    equal invariants).
O2  For each T_i: the genus-two embedding; the no-K5 certificate; a torus
    embedding of G - e and of G/e for every edge e; the connectivity and the
    twin triples, recomputed; the case of the proof of Theorem 1.1, matched
    with the connectivity and counted (20 of connectivity at most two, K3,7,
    and 15, 24 and 14 for Theorems 4.1, 5.1 and 6.1); the structure of T8 and
    T25; the canonical forms are distinct (and, if nauty's labelg is given by
    the environment variable LABELG, recomputed).
O3  The 68 graphs of connectivity two in the Myrvold-Woodcock collection:
    each equals its entry of the collection in O4; K5 branch sets for 50 and
    isomorphisms onto T_i for 18, these 18 being the T_i of connectivity two;
    pairwise nonisomorphic; and they are exactly the entries of the
    collection of connectivity two.
O4  The collection: the isomorphisms from T1, ..., T73 onto their entries,
    and that T74 is isomorphic to no entry.
TO1 The 381 graphs: their orders as stated, the 74 named ones isomorphic to
    T1, ..., T74, pairwise nonisomorphic.
TO2 The 7,006 graphs of the search: each is the split of its recorded parent
    as recorded (contracting the edge between the split vertex and the new
    vertex gives the parent exactly); each verdict's witness (a no-K5 certificate and a torus
    embedding of every edge deletion for a kept graph; K5 branch sets; or an
    edge whose deletion has a minor model of T_i); the counts; and that
    splitting every vertex of every graph of TO1 in every way gives exactly
    these graphs up to isomorphism, the kept ones being the graphs of TO1
    other than T1, ..., T74.
"""
import os
import shutil
import subprocess
from collections import Counter
from itertools import combinations

import networkx as nx

from witnesses import (catalogue, check_embedding, check_iso, check_k5, check_minor_model, check_no_k5, components,
                       contract_edge, delete_edge, delete_vertices, edges, from_graph6, load, n_edges, need,
                       same_graph, vertex_of)

D = 'no_K5_obstructions'


def nxg(adj):
    G = nx.Graph()
    G.add_nodes_from(adj)
    G.add_edges_from(edges(adj))
    return G


def invariant(adj):
    deg = {v: len(s) for v, s in adj.items()}
    tri = sum(1 for u in adj for v in adj[u] for w in adj[v] if w in adj[u]) // 6
    return (len(adj), n_edges(adj), tuple(sorted(deg.values())),
            tuple(sorted(tuple(sorted(deg[u] for u in adj[v])) for v in adj)), tri)


class IsoIndex:
    """Graphs up to isomorphism: classes by invariants, then VF2."""
    def __init__(self):
        self.buckets = {}

    def find(self, adj):
        G = nxg(adj)
        for k, H in self.buckets.get(invariant(adj), []):
            if nx.is_isomorphic(G, H):
                return k
        return None

    def add(self, adj, key):
        self.buckets.setdefault(invariant(adj), []).append((key, nxg(adj)))


def distinct(graphs, what):
    idx = IsoIndex()
    for k, adj in graphs:
        need(idx.find(adj) is None, what, k, 'is isomorphic to an earlier one')
        idx.add(adj, k)
    return idx


def connectivity(adj):
    """Vertex connectivity, up to 4 (4 means at least 4); 0 when disconnected."""
    V = list(adj)
    if len(components(adj)) > 1:
        return 0
    for k in (1, 2, 3):
        if len(V) <= k + 1:
            return k
        for S in combinations(V, k):
            if len(components(delete_vertices(adj, S))) > 1:
                return k
    return 4


def twin_triples(adj):
    groups = {}
    for v in adj:
        if len(adj[v]) == 3:
            groups.setdefault(frozenset(adj[v]), []).append(v)
    return {k: sorted(v) for k, v in groups.items() if len(v) >= 2}


def labelg():
    p = os.environ.get('LABELG') or shutil.which('labelg')
    return p if p and os.path.exists(p) else None


def canonical(g6s):
    r = subprocess.run([labelg(), '-q'], input='\n'.join(g6s) + '\n', capture_output=True, text=True, check=True)
    out = r.stdout.split()
    need(len(out) == len(g6s), 'labelg output')
    return out


def to_graph6(adj):
    order = sorted(adj)
    n = len(order)
    bits = [1 if order[i] in adj[order[j]] else 0 for j in range(1, n) for i in range(j)]
    while len(bits) % 6:
        bits.append(0)
    out = [n + 63] + [int(''.join(map(str, bits[k:k + 6])), 2) + 63 for k in range(0, len(bits), 6)]
    return ''.join(chr(c) for c in out)


def mw_graph(line):
    n, bits = line.split()
    n = int(n)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    need(len(bits) == len(pairs), 'a line of the collection has the wrong length')
    adj = {i: set() for i in range(n)}
    for (i, j), b in zip(pairs, bits):
        if b == '1':
            adj[i].add(j)
            adj[j].add(i)
    return adj


def ints(m):
    return {int(k): int(v) for k, v in m.items()}


def int_rotation(rot):
    return {int(k): [int(x) for x in v] for k, v in rot.items()}


def complete_bipartite(a, b):
    adj = {v: set() for v in range(a + b)}
    for x in range(a):
        for y in range(a, a + b):
            adj[x].add(y)
            adj[y].add(x)
    return adj


# ---------------------------------------------------------------- O1-O4

def check_O1(root):
    d = load(os.path.join(root, D, 'O1.json'))
    need([g['name'] for g in d['graphs']] == ['T%d' % i for i in range(1, 75)], 'O1 must list T1, ..., T74 in order')
    T = catalogue(root)          # also checks graph6 against the edge lists and the sizes
    distinct(sorted(T.items(), key=lambda kv: int(kv[0][1:])), 'O1:')
    return '74 graphs, graph6 and edges agree, pairwise nonisomorphic'


def check_O2(root):
    d = load(os.path.join(root, D, 'O2.json'))
    T = catalogue(root)
    need([g['name'] for g in d['graphs']] == ['T%d' % i for i in range(1, 75)], 'O2 must cover T1, ..., T74 in order')
    cases = Counter()
    n_emb = 0
    for g in d['graphs']:
        name = g['name']
        G = T[name]
        e = g['genus-two embedding']
        check_embedding(G, {'type': 'embedding', 'genus': e['genus'], 'rotation': int_rotation(e['rotation'])}, genus=2)
        es = {frozenset(x) for x in edges(G)}
        for key, op in (('G - e', delete_edge), ('G / e', contract_edge)):
            seen = set()
            for x in g[key]:
                u, v = int(x['edge'][0]), int(x['edge'][1])
                H = op(G, (u, v))
                check_embedding(H, {'type': 'embedding', 'genus': x['genus'], 'rotation': int_rotation(x['rotation'])},
                                at_most=1)
                seen.add(frozenset((u, v)))
                n_emb += 1
            need(seen == es and len(g[key]) == len(es), name, key, 'is not given exactly once for every edge')
        check_no_k5(G, g['no K5 certificate'])
        k = connectivity(G)
        stated = g['connectivity']
        want = {0: 'disconnected', 1: 'has a cut vertex'}.get(k, k)
        need(stated == want, name, 'connectivity stated', stated, 'but it is', want)
        tt = twin_triples(G)
        mine = sorted((sorted(t), tw) for t, tw in tt.items())
        theirs = sorted((sorted(x['triple']), sorted(x['twins'])) for x in g['twin triples'])
        need(mine == theirs, name, 'twin triples differ')
        case = g['case']
        cases[case] += 1
        if case.startswith('Theorem') or case == 'K3,7':
            need(k == 3, name, 'is in', case, 'but has connectivity', k)
        elif case == 'connectivity two':
            need(k == 2, name, 'connectivity two stated')
        else:
            need(case == 'connectivity at most one' and k <= 1, name, 'unexpected case', case)
        if case == 'K3,7':
            need(name == 'T1', 'only T1 is K3,7')
            A = sorted(v for v in G if len(G[v]) == 7)
            B = sorted(v for v in G if len(G[v]) == 3)
            check_iso(G, complete_bipartite(3, 7), dict(zip(A + B, range(10))))
        if 'structure' in g:
            m = ints(g['structure']['isomorphism'])
            target = {v: set() for v in range(len(G))}
            if name == 'T25':
                sides = [((0, 1, 2), (3, 4, 5)), ((6, 7, 8), (9, 10, 11))]
            else:
                need(name == 'T8', 'a structure is recorded only for T8 and T25')
                sides = [((0, 1, 2), (3, 4, 5)), ((5, 6, 7), (8, 9, 10))]
            for X, Y in sides:
                for x in X:
                    for y in Y:
                        target[x].add(y)
                        target[y].add(x)
            check_iso(G, target, m)
    need(dict(cases) == {'Theorem 4.1': 15, 'Theorem 5.1': 24, 'Theorem 6.1': 14, 'connectivity two': 18, 'K3,7': 1,
                         'connectivity at most one': 2}, 'case counts', dict(cases))
    can = [g['canonical graph6'] for g in d['graphs']]
    need(len(set(can)) == 74, 'O2: canonical forms not distinct')
    note = 'canonical forms distinct'
    if labelg():
        mine = canonical([to_graph6(T['T%d' % i]) for i in range(1, 75)])
        need(mine == can, 'O2: canonical forms differ from labelg')
        note = 'canonical forms recomputed with labelg'
    return '74 genus-two embeddings, %d torus embeddings of G-e and G/e, no-K5 certificates, connectivity, twin triples, ' \
           'cases %s; %s' % (n_emb, dict(cases), note)


def collection(root):
    d = load(os.path.join(root, D, 'O4.json'))
    lines = d['collection']['lines']
    need(len(lines) == d['collection']['entries'] == 17535, 'O4: the collection should have 17,535 entries')
    return d, [mw_graph(l) for l in lines]


def check_O3(root):
    d = load(os.path.join(root, D, 'O3.json'))
    T = catalogue(root)
    o2 = {g['name']: g for g in load(os.path.join(root, D, 'O2.json'))['graphs']}
    _, MW = collection(root)
    verdicts = Counter()
    iso_T = []
    for g in d['graphs']:
        adj = from_graph6(g['graph6'])
        e = g['entry']
        need(adj == MW[e - 1], 'O3: entry', e, 'is not the graph of that entry')
        need(len(adj) == g['n'] and n_edges(adj) == g['m'], 'O3: size of entry', e)
        w = g['witness']
        if g['verdict'] == 'K5 minor':
            check_k5(adj, w['branch sets'])
        else:
            tn = g['verdict'].split()[-1]
            need(w['type'] == 'isomorphism' and w['onto'] == tn, 'O3: entry', e)
            check_iso(adj, T[tn], ints(w['map']))
            iso_T.append(tn)
        verdicts[g['verdict'].split(' to ')[0]] += 1
    need(len(d['graphs']) == 68 and verdicts['K5 minor'] == 50 and len(iso_T) == 18, 'O3 counts', dict(verdicts))
    conn2 = sorted(n for n, g in o2.items() if g['case'] == 'connectivity two')
    need(sorted(iso_T) == conn2, 'O3: the 18 graphs without K5 minor are not the T_i of connectivity two')
    distinct([(g['entry'], from_graph6(g['graph6'])) for g in d['graphs']], 'O3: entry')
    mine = sorted(i + 1 for i, adj in enumerate(MW) if connectivity(adj) == 2)
    need(mine == sorted(g['entry'] for g in d['graphs']), 'O3: the entries of connectivity two differ')
    return '68 entries (all of connectivity two in the collection): 50 K5 branch sets, 18 isomorphisms onto the T_i ' \
           'of connectivity two; pairwise nonisomorphic'


def check_O4(root):
    d, MW = collection(root)
    T = catalogue(root)
    found = d['T graphs in the collection']
    need(sorted(found, key=lambda s: int(s[1:])) == ['T%d' % i for i in range(1, 74)] and
         d['T graphs not in the collection'] == ['T74'], 'O4: the T_i in the collection should be T1, ..., T73')
    for name, x in found.items():
        check_iso(T[name], MW[x['entry'] - 1], ints(x['isomorphism']))
    T74 = T['T74']
    inv = invariant(T74)
    same = [i for i, adj in enumerate(MW) if invariant(adj) == inv]
    for i in same:
        need(not nx.is_isomorphic(nxg(T74), nxg(MW[i])), 'O4: T74 is isomorphic to entry', i + 1)
    note = ''
    if labelg():
        can = canonical([to_graph6(adj) for adj in MW])
        need(can == d['collection']['canonical graph6 of each entry'], 'O4: canonical forms differ from labelg')
        need(canonical([to_graph6(T74)])[0] == d['canonical form of T74'] not in set(can), 'O4: canonical form of T74')
        note = '; canonical forms recomputed with labelg'
    return '17,535 entries; isomorphisms from T1, ..., T73 onto their entries; T74 isomorphic to none (%d entries ' \
           'share its invariants)%s' % (len(same), note)


# ---------------------------------------------------------------- TO1, TO2

def splits(g):
    """Every split of every vertex v of g: v keeps its first neighbour and at
    least one more, the new vertex n gets at least two of the others, and
    v and n are joined (both end with degree at least three)."""
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
            yield h


def check_TO1(root):
    d = load(os.path.join(root, 'no_K5_topological_obstructions', 'TO1.json'))
    T = catalogue(root)
    gs = d['graphs']
    need(len(gs) == 381, 'TO1 should list 381 graphs')
    orders = Counter()
    named = {}
    for g in gs:
        adj = from_graph6(g['graph6'])
        need(len(adj) == g['n'] and n_edges(adj) == g['m'], 'TO1: size of', g['graph6'])
        orders[len(adj)] += 1
        if 'is' in g:
            need(nx.is_isomorphic(nxg(adj), nxg(T[g['is']])), 'TO1:', g['graph6'], 'is not', g['is'])
            named[g['is']] = g['graph6']
    need(sorted(named, key=lambda s: int(s[1:])) == ['T%d' % i for i in range(1, 75)], 'TO1 must name T1, ..., T74 once each')
    need({int(k): v for k, v in d['orders'].items()} == dict(orders), 'TO1: orders')
    need(dict(sorted(orders.items())) == {10: 1, 11: 9, 12: 43, 13: 111, 14: 151, 15: 58, 16: 7, 17: 1},
         'TO1: the orders differ from the table of Theorem 8.1')
    distinct([(g['graph6'], from_graph6(g['graph6'])) for g in gs], 'TO1:')
    return '381 graphs, orders as in Theorem 8.1, T1, ..., T74 named, pairwise nonisomorphic'


def check_TO2(root):
    d = load(os.path.join(root, 'no_K5_topological_obstructions', 'TO2.json'))
    to1 = load(os.path.join(root, 'no_K5_topological_obstructions', 'TO1.json'))['graphs']
    T = catalogue(root)
    listed = {g['graph6'] for g in to1}
    recs = d['graphs']
    counts = Counter()
    n_rot = 0
    index = IsoIndex()
    for k, r in enumerate(recs):
        G = from_graph6(r['graph6'])
        a, b = r['split']['vertex'], r['split']['new vertex']
        need(b in G[a] and not (G[a] & G[b]), 'TO2', k, ': the split vertices are not adjacent with disjoint neighbourhoods')
        need(sorted(G[b] - {a}) == sorted(r['split']['neighbours moved to the new vertex']), 'TO2', k, ': moved neighbours')
        need(min(len(s) for s in G.values()) >= 3, 'TO2', k, ': a vertex of degree below three')
        p = r['split from']
        if isinstance(p, str):
            P = T[p]
        else:
            need(p['graph6'] in listed, 'TO2', k, ': its parent is not in TO1')
            P = from_graph6(p['graph6'])
        need(same_graph(contract_edge(G, (a, b)), P), 'TO2', k, ': contracting the split edge does not give the parent')
        w = r['witness']
        v = r['verdict']
        if v == 'kept':
            check_no_k5(G, w['no K5 certificate'])
            seen = set()
            for x in w['torus embeddings of every G - e']:
                u, z = x['edge']
                check_embedding(delete_edge(G, (u, z)), {'type': 'embedding', 'genus': x['genus'],
                                                          'rotation': x['rotation']}, at_most=1)
                seen.add(frozenset((u, z)))
                n_rot += 1
            need(seen == {frozenset(e) for e in edges(G)}, 'TO2', k, ': not every edge deletion is embedded')
        elif v == 'K5 minor':
            check_k5(G, w['branch sets'])
        else:
            need(v == 'an edge deletion has a T_i minor', 'TO2', k, ': unknown verdict', v)
            u, z = w['edge deleted']
            H = delete_edge(G, (u, z))
            Ti = T[w['of']]
            check_minor_model(H, Ti, {vertex_of(Ti, x): vs for x, vs in w['branch sets'].items()})
        counts[v] += 1
        need(index.find(G) is None, 'TO2', k, ': isomorphic to an earlier record')
        index.add(G, k)
    c = d['counts']
    need(len(recs) == c['graphs'] == 7006 and counts['kept'] == c['kept'] == 307 and
         counts['K5 minor'] == c['K5 minor'] == 4884 and
         counts['an edge deletion has a T_i minor'] == c['edge deletion with a T_i minor'] == 1815, 'TO2 counts', dict(counts))
    # every split of every graph of TO1 is a recorded graph, and every recorded graph is hit
    hit, n_splits = set(), 0
    for g in to1:
        for h in splits(from_graph6(g['graph6'])):
            n_splits += 1
            k = index.find(h)
            need(k is not None, 'TO2: a split of', g['graph6'], 'is not recorded')
            hit.add(k)
    need(hit == set(range(len(recs))), 'TO2: some recorded graph is not a split of a graph of TO1')
    need(n_splits == c['splits of the 381 graphs'], 'TO2: number of splits')
    kept = set()
    for g in to1:
        if 'is' in g:
            continue
        k = index.find(from_graph6(g['graph6']))
        need(k is not None and recs[k]['verdict'] == 'kept', 'TO2: a graph of TO1 is not a kept split')
        kept.add(k)
    need(kept == {k for k, r in enumerate(recs) if r['verdict'] == 'kept'}, 'TO2: kept graphs and TO1 differ')
    return '7,006 graphs, each the split of its parent; %d torus embeddings of edge deletions; %s; all %d splits of the ' \
           '381 graphs recorded, and the 307 kept graphs are the graphs of TO1 besides T1, ..., T74' % (
               n_rot, dict(counts), n_splits)


def check(root, item):
    return {'O1': check_O1, 'O2': check_O2, 'O3': check_O3, 'O4': check_O4, 'TO1': check_TO1,
            'TO2': check_TO2}[item](root)
