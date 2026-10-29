"""Make N1: thirteen toroidal minor obstructions that are not in the collection
of Myrvold and Woodcock.

The thirteen graphs were found in searches during this research. Each has a K5
minor, so none of them is among T1, ..., T74. For each graph G this program
records

  * a genus-two embedding of G, and the genus of G computed by multi_genus
    (option w writes a minimum genus embedding, checked here by tracing its
    faces), and with the option --page also by PAGE;
  * a torus embedding of G - e and of G/e for every edge e (for e = [u, v],
    v is merged into u, loops are deleted and parallel edges merged), each the
    minimum genus embedding written by multi_genus;
  * K5 branch sets;
  * its connectivity;
  * its canonical graph6 (nauty's labelg), and the entries of the collection
    (item O4) with the same invariants, none of which is isomorphic to G.

Since every edge deletion and every edge contraction of G is toroidal, every
proper minor of G is toroidal; with genus two, G is a toroidal minor
obstruction.

Usage, from the top folder of the repository:

    python programs/new_obstructions/make_N1.py OUTPUT_FOLDER [--page]

It writes OUTPUT_FOLDER/new_obstructions/N1.json. It needs graph-genus 0.1.3,
networkx and nauty's labelg (found through the environment variable LABELG, or
as labelg or nauty-labelg on the PATH).
"""
import json
import os
import random
import shutil
import subprocess
import sys
import time
from itertools import combinations
from pathlib import Path

import networkx as nx

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'genus'))
from common import decode_rotations, genus_of_rotation, multicode, tool, zero_indexed  # noqa: E402

# The thirteen graphs, in the order in which they were found.
GRAPHS = [
    'M??????UplPk{Wz_?',
    'M?????FAtsMgv?wK?',
    'M?????FAuSEg~?{K?',
    'L?CP@`Er@XRWoq',
    'L?CAGKtg_Qz@~A',
    'MD?@W@C?Od`MrOwS?',
    'L??acUoY?l@jv?',
    'L?_PIQOWKeL`B{',
    'L??Ge?rYFKBifc',
    'M?B?CcBGgod_OJ}_?',
    'LOC?HGINeFIkuW',
    'M???OJIhPsIgcoZ?_',
    'M???OKBS[wTC\\?x_?',
]


def from_graph6(s):
    G = nx.from_graph6_bytes(s.encode())
    return {v: set(G[v]) for v in G}


def edges(adj):
    return sorted((u, v) for u in adj for v in adj[u] if u < v)


def delete_edge(adj, e):
    u, v = e
    H = {x: set(s) for x, s in adj.items()}
    H[u].discard(v)
    H[v].discard(u)
    return H


def contract_edge(adj, e):
    """v is merged into u; loops are deleted and parallel edges merged."""
    u, v = e
    H = {x: set(s) for x, s in adj.items()}
    for w in H.pop(v):
        H[w].discard(v)
        if w != u:
            H[w].add(u)
            H[u].add(w)
    return H


def multi_genus_embeddings(graphs):
    """A minimum genus embedding of each connected graph, by multi_genus,
    with the vertices of each graph in sorted order."""
    data = b''.join(multicode(G) for G in graphs)
    p = subprocess.run([tool('multi_genus'), 'w', 'L'], input=data, capture_output=True, check=True)
    rots = decode_rotations(p.stdout)
    assert len(rots) == len(graphs)
    out = []
    for G, rot in zip(graphs, rots):
        vs = sorted(G)
        named = {vs[i]: [vs[u] for u in nb] for i, nb in enumerate(rot)}
        assert all(sorted(named[v]) == sorted(G[v]) for v in vs), 'embedding of another graph'
        g = genus_of_rotation(rot)
        out.append((g, {str(v): named[v] for v in vs}))
    return out


def page_genus(adj):
    import graph_genus as gg
    A = zero_indexed(adj)
    t = time.perf_counter()
    g, rot = gg.embed(A, algorithm='page')
    dt = time.perf_counter() - t
    assert genus_of_rotation(rot) == g, 'PAGE reported a genus its embedding does not have'
    return g, dt


def k5_branch_sets(adj, seed=0, tries=200000):
    """K5 branch sets, found by contracting random edges until five vertices
    remain and checking that they are pairwise adjacent."""
    rnd = random.Random(seed)
    for _ in range(tries):
        blocks = {v: {v} for v in adj}
        H = {x: set(s) for x, s in adj.items()}
        while len(H) > 5:
            u = rnd.choice(sorted(H))
            if not H[u]:
                break
            v = rnd.choice(sorted(H[u]))
            for w in H.pop(v):
                H[w].discard(v)
                if w != u:
                    H[w].add(u)
                    H[u].add(w)
            blocks[u] |= blocks.pop(v)
        if len(H) == 5 and all(len(H[x]) == 4 for x in H):
            return [sorted(blocks[x]) for x in sorted(H)]
    return None


def check_k5(adj, sets):
    allv = [v for s in sets for v in s]
    assert len(sets) == 5 and len(allv) == len(set(allv))
    for s in sets:
        assert nx.is_connected(nx.Graph([(a, b) for a in s for b in adj[a] if b in s]) if len(s) > 1
                               else nx.path_graph(1))
    for a, b in combinations(sets, 2):
        assert any(y in adj[x] for x in a for y in b)


def connectivity(adj):
    G = nx.Graph()
    G.add_nodes_from(adj)
    G.add_edges_from(edges(adj))
    return nx.node_connectivity(G)


def labelg():
    for p in (os.environ.get('LABELG'), shutil.which('labelg'), shutil.which('nauty-labelg')):
        if p and os.path.exists(p):
            return p
    raise SystemExit('nauty labelg not found (set LABELG)')


def canonical(g6s):
    r = subprocess.run([labelg(), '-q'], input='\n'.join(g6s) + '\n', capture_output=True, text=True, check=True)
    out = r.stdout.split()
    assert len(out) == len(g6s)
    return out


def invariant(adj):
    deg = {v: len(s) for v, s in adj.items()}
    tri = sum(1 for u in adj for v in adj[u] for w in adj[v] if w in adj[u]) // 6
    return (len(adj), sum(deg.values()) // 2, tuple(sorted(deg.values())),
            tuple(sorted(tuple(sorted(deg[u] for u in adj[v])) for v in adj)), tri)


def mw_graph(line):
    parts = line.split()
    n, bits = int(parts[0]), parts[1]
    adj = {v: set() for v in range(n)}
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits[k] == '1':
                adj[i].add(j)
                adj[j].add(i)
            k += 1
    return adj


def nxg(adj):
    G = nx.Graph()
    G.add_nodes_from(adj)
    G.add_edges_from(edges(adj))
    return G


def main():
    out_root = Path(sys.argv[1])
    # With the option --page, the genus of each graph is also computed by PAGE
    # and recorded. PAGE takes far longer than multi_genus on some of these
    # graphs, so it is not run by default.
    use_page = '--page' in sys.argv[2:]
    o4 = json.load(open(REPO / 'supplemental_data' / 'no_K5_obstructions' / 'O4.json'))
    MW = [mw_graph(l) for l in o4['collection']['lines']]
    mw_canon = set(o4['collection']['canonical graph6 of each entry'])
    by_inv = {}
    for i, adj in enumerate(MW):
        by_inv.setdefault(invariant(adj), []).append(i)
    canon = canonical(GRAPHS)
    assert len(set(canon)) == len(GRAPHS), 'two of the graphs are isomorphic'
    records = []
    for k, (s, cf) in enumerate(zip(GRAPHS, canon), 1):
        t0 = time.time()
        G = from_graph6(s)
        assert nx.is_connected(nxg(G)) and min(len(x) for x in G.values()) >= 3
        es = edges(G)
        minors = [G] + [delete_edge(G, e) for e in es] + [contract_edge(G, e) for e in es]
        emb = multi_genus_embeddings(minors)
        g, rot = emb[0]
        assert g == 2, 'graph %d has genus %d by multi_genus' % (k, g)
        genus = {'multi_genus': g}
        if use_page:
            pg, pt = page_genus(G)
            assert pg == 2, 'graph %d has genus %d by PAGE' % (k, pg)
            genus['PAGE'] = pg
        dele = []
        cont = []
        for i, e in enumerate(es):
            gd, rd = emb[1 + i]
            gc, rc = emb[1 + len(es) + i]
            assert gd <= 1 and gc <= 1, 'graph %d: a one-edge minor at %s is not toroidal' % (k, e)
            dele.append({'edge': list(e), 'genus': gd, 'rotation': rd})
            cont.append({'edge': list(e), 'genus': gc, 'rotation': rc})
        bs = k5_branch_sets(G, seed=k)
        assert bs is not None, 'graph %d: no K5 minor found' % k
        check_k5(G, bs)
        same = by_inv.get(invariant(G), [])
        for i in same:
            assert not nx.is_isomorphic(nxg(G), nxg(MW[i])), 'graph %d is entry %d of the collection' % (k, i + 1)
        assert cf not in mw_canon, 'graph %d: canonical form is in the collection' % k
        records.append({
            'name': 'N1.%d' % k,
            'n': len(G),
            'm': len(es),
            'graph6': s,
            'edges': [list(e) for e in es],
            'connectivity': connectivity(G),
            'genus': genus,
            'genus-two embedding': {'genus': g, 'rotation': rot},
            'G - e': dele,
            'G / e': cont,
            'K5 branch sets': bs,
            'canonical graph6': cf,
            'entries of the collection with the same invariants': [i + 1 for i in same],
        })
        print('N1.%d  n=%d m=%d connectivity %d  genus 2 (%s)  %d one-edge minors toroidal  '
              'K5 minor  not in the collection (%d entries share its invariants)  %.1f s'
              % (k, len(G), len(es), records[-1]['connectivity'], ', '.join(genus), 2 * len(es), len(same),
                 time.time() - t0),
              flush=True)
    item = {
        'item': 'N1',
        'what': ('Thirteen toroidal minor obstructions that are not in the collection of Myrvold and Woodcock '
                 '(item O4), found in searches during this research. Each has a K5 minor, so none is among '
                 'T1,...,T74. Vertices are 0,...,n-1, numbered as in the graph6 string. For each graph: a '
                 'genus-two embedding; a torus embedding of every G-e and of every G/e (for G/e with e=[u,v], v is '
                 'merged into u); K5 branch sets; its connectivity; its canonical graph6; and the entries of the '
                 'collection with the same invariants (numbers of vertices and edges, degrees, degrees of the '
                 'neighbours of each vertex, and number of triangles), none of which is isomorphic to it. An '
                 'embedding is a rotation system (the cyclic order of the neighbours at each vertex) with its '
                 'genus. Since every edge deletion and edge contraction is toroidal, every proper minor is '
                 'toroidal. The genus lower bound (genus two, computed by multi_genus as bundled with graph-genus '
                 '0.1.3) is recorded under "genus"; programs/new_obstructions/make_N1.py recomputes it, and with '
                 'the option --page also computes it with PAGE.'),
        'canonical forms': {'program': 'nauty 2.8.8, labelg (graph6 of the canonically labelled graph)',
                            'all 13 distinct': True,
                            'none among the canonical forms of the collection in O4': True},
        'graphs': records,
    }
    folder = out_root / 'new_obstructions'
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / 'N1.json', 'w') as f:
        f.write(json.dumps(item) + '\n')
    print('wrote', folder / 'N1.json')


if __name__ == '__main__':
    main()
