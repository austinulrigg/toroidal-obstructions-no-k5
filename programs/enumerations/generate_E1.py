"""E1: the bridge H1 of Lemma 4.10 (sole two rails).

The names: Psi_1 has twins p, q on {s1, s2, s3} and r, z on {t1, t2, t3}.
G = H0 + B, where H0 consists of p, q, r, z with their twelve edges and the
edge s3 t3, and the bridge B (H1 in the proof) meets H0 in A = {s1, t1, s2,
t2}. By the proof, B is a part of Construction 3.17 with attachments A and
added edges E = {s1 s2, t1 t2}: A is independent in B, the interior is
nonempty and connected, every attachment has an interior neighbour, and
B + E is nonplanar.

The candidates are the parts of four_root_domain.py with these attachments
and added edges; each records the K3,3 subdivision of B + E it came from.

Classification of each candidate, in order: a vertex of G of degree below
three whose neighbours separate G; an allowed operation of Section 3 after
which B' + E is still nonplanar (a proper nontoroidal minor, by Lemma
4.11); a cut of at most two vertices of G; B without disjoint s1-t1 and
s2-t2 paths through its interior (B must contain P1 and P2); a K5 minor of
B + E or of G (random contractions, small exhaustive search, SAT); and a
model of a graph of the catalogue in G, found after a bounded number of
contractions (minor_search.known_minor), first with a small budget and then
with a larger one; a model found is kept and carried over to the later
candidates whose graph G is isomorphic. The proof finds no obstruction
here.

The classification keeps a file progress.jsonl and continues from it when
it is started again. Run as: python generate_E1.py FOLDER
"""
import itertools as it
import json
import sys
import time
from collections import Counter
from pathlib import Path

import networkx as nx

import common
import four_root_domain as fd
from minor_search import find_k5_minor, known_minor, PositiveModelCache

ROOTS = ("s1", "t1", "s2", "t2")
ADDED_EDGES = (("s1", "s2"), ("t1", "t2"))
edge_list = fd.edge_list


# ---------------------------------------------------------------- the classification

def completion(p):
    q = p.copy()
    q.add_edges_from(ADDED_EDGES)
    return q


def whole_graph(p):
    g = p.copy()
    g.add_edges_from(it.product(('p', 'q'), ('s1', 's2', 's3')))
    g.add_edges_from(it.product(('r', 'z'), ('t1', 't2', 't3')))
    g.add_edge('s3', 't3')
    return g


def valid(p):
    interior = set(p) - set(ROOTS)
    return bool(interior) and nx.is_connected(p.subgraph(interior)) and all(set(p[r]) & interior for r in ROOTS)


def allowed_operations(p):
    # Deletions often leave a visibly retained subdivision, so try them first.
    for edge in edge_list(p):
        q = p.copy()
        q.remove_edge(*edge)
        if valid(q):
            yield {'type': 'delete_edge', 'edge': edge}, q
    for vertex in sorted(set(p) - set(ROOTS)):
        q = p.copy()
        q.remove_node(vertex)
        if valid(q):
            yield {'type': 'delete_vertex', 'vertex': vertex}, q
    for edge in edge_list(p):
        a, b = edge
        if b in ROOTS:
            a, b = b, a
        q = nx.contracted_nodes(p, a, b, self_loops=False)
        q.remove_edges_from(it.combinations(ROOTS, 2))
        if valid(q):
            yield {'type': 'contract_edge', 'edge': edge, 'keep': a}, q


def nonplanarity(g):
    planar, witness = nx.check_planarity(g, counterexample=True)
    assert not planar
    return {'kuratowski_edges': edge_list(witness)}


def cut_outcome(g, cut):
    pieces = [sorted(c) for c in nx.connected_components(g.subgraph(set(g) - set(cut)))]
    if len(pieces) > 1:
        return {'type': 'structural_cut', 'cut': list(cut), 'components': sorted(pieces)}


def has_rails(p):
    """Whether B has disjoint s1-t1 and s2-t2 paths through its interior."""
    interior = sorted(set(p) - set(ROOTS))
    names = sorted(p)
    index = {v: i for i, v in enumerate(names)}
    adjacency = [sum(1 << index[w] for w in p[v]) for v in names]

    def path(a, b, allowed):
        reached = 0
        frontier = 1 << index[a]
        target = 1 << index[b]
        while frontier:
            if frontier & target:
                return True
            reached |= frontier
            neighbors = 0
            while frontier:
                bit = frontier & -frontier
                frontier -= bit
                neighbors |= adjacency[bit.bit_length() - 1]
            frontier = neighbors & allowed & ~reached
        return False
    full = sum(1 << index[v] for v in interior)
    for mask in range(1 << len(interior)):
        first = sum(1 << index[v] for k, v in enumerate(interior) if mask >> k & 1)
        if path('s1', 't1', first | (1 << index['s1']) | (1 << index['t1'])) and \
                path('s2', 't2', (full ^ first) | (1 << index['s2']) | (1 << index['t2'])):
            return True
    return False


def classify(p, catalog, minor_cache):
    g = whole_graph(p)
    for v in sorted(g):
        if g.degree(v) < 3:
            outcome = cut_outcome(g, sorted(g[v]))
            if outcome:
                return outcome
    for operation, q in allowed_operations(p):
        j = completion(q)
        if not nx.check_planarity(j)[0]:
            return {'type': 'proper_allowed_minor', 'operation': operation,
                    'minor_vertices': sorted(q), 'minor_edges': edge_list(q),
                    'nonplanarity': nonplanarity(j)}
    for k in (1, 2):
        for cut in it.combinations(sorted(g), k):
            outcome = cut_outcome(g, cut)
            if outcome:
                return outcome
    if not has_rails(p):
        return {'type': 'no_disjoint_rails', 'interior_partition_count': 1 << (len(p) - 4)}
    for graph_name, h in (('completion', completion(p)), ('whole', g)):
        stats = {}
        model = find_k5_minor(h, trials=64, seed=0, exact_threshold=8,
                              use_sat=True, sat_conflicts=5000, sat_rounds=50, stats=stats)
        if model is not None:
            return {'type': 'k5_minor', 'graph': graph_name, 'branch_sets': model, 'search': stats}
    cached = minor_cache.find(g)
    if cached is not None:
        return {'type': 'proper_obstruction_minor', 'target': cached['metadata']['target'],
                'branch_sets': cached['branch_sets'], 'search': 'checked isomorphism transport'}
    for contractions, states in ((2, 5000), (4, 15000)):
        stats = {}
        found = known_minor(g, catalog, max_contractions=contractions,
                            max_edge_deletions=3, max_states=states, stats=stats)
        if found is not None:
            found['search'] = stats
            target = dict(catalog)[found['target']]
            minor_cache.add(g, target, found['branch_sets'], {'target': found['target']})
            return found
    return {'type': 'unresolved'}


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    records, stats = fd.domain(ROOTS, ADDED_EDGES, 'E1-')
    print('domain:', stats, round(time.time() - t0), 's', flush=True)
    catalog = [(r['name'], nx.from_graph6_bytes(r['graph6'].encode())) for r in common.catalogue_records()]
    progress = out / 'progress.jsonl'
    done = {}
    if progress.exists():
        for line in progress.read_text().splitlines():
            try:
                r = json.loads(line)
            except ValueError:
                continue
            done[r['id']] = r
    minor_cache = PositiveModelCache()
    counts = Counter()
    last = time.time()
    with progress.open('a') as pf:
        for record in records:
            if record['id'] in done:
                record.update(done[record['id']])
                counts[record['outcome']['type']] += 1
                continue
            p = nx.Graph()
            p.add_nodes_from(record['vertices'])
            p.add_edges_from(record['edges'])
            record['nonplanarity'] = {'kuratowski_edges': fd.subdivision_of(record, ROOTS, ADDED_EDGES)}
            record['outcome'] = classify(p, catalog, minor_cache)
            pf.write(json.dumps({'id': record['id'], 'nonplanarity': record['nonplanarity'],
                                 'outcome': record['outcome']}, sort_keys=True, separators=(',', ':')) + '\n')
            pf.flush()
            counts[record['outcome']['type']] += 1
            if time.time() - last > 30:
                print(json.dumps({'classified': sum(counts.values()), 'of': len(records), 'outcomes': dict(counts)}),
                      flush=True)
                last = time.time()
    document = {'item': 'E1', 'candidate_count': len(records), 'outcome_counts': dict(counts),
                'candidates': records}
    (out / 'candidates.json').write_text(json.dumps(document, sort_keys=True, separators=(',', ':')) + '\n')
    summary = dict(stats, item='E1', where='Lemma 4.10', outcome_counts=dict(counts),
                   unresolved=counts.get('unresolved', 0), seconds=round(time.time() - t0, 1))
    (out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
