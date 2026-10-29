"""E3: a bridge with four attachments and G = (Psi_1 + B) + af (Proposition 4.21).

The names: a = s1, b = t1, c = s2, d = t2, e = s3, f = t3, with twins p, q on
{a, c, e} and r, z on {b, d, f}; the paths are the edges ab, cd, ef. As in
the proof, G = (Psi_1 + B) + af, and G is toroidal exactly when
B + ac + cd + bd is planar. So B is a part of Construction 3.17 with
attachments a, b, c, d, added edges E = {ac, cd, bd} and no required edges
(four_root_domain.py); each candidate records the K3,3 subdivision of B + E
it came from.

Classification of each candidate, in order: a vertex of degree three of G
with two adjacent neighbours (recorded with that edge); an interior vertex
of degree two; the interior of B nonplanar (a K5 minor, or two disjoint
K3,3 minors, a model of a disconnected graph of the catalogue); an allowed
operation of Section 3 after which B' + E is still nonplanar; one edge
deletion or contraction after which every choice of the joint test of
Lemma 4.3 fails (joint_test.py); a K5 minor of B + E or of G. The data
folder gives every candidate a proper nontoroidal minor obtained by
deleting or contracting one edge, with the reason that minor is
nontoroidal (package.py). The proof finds no obstruction here.

The classification keeps a file progress.jsonl and continues from it when
it is started again. Run as: python generate_E3.py FOLDER
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
from joint_test import joint_test, CORE_VERTICES, CORE_EDGES
from minor_search import find_k5_minor, first_isomorphism

ROOTS = ("a", "b", "c", "d")
ADDED_EDGES = (("a", "c"), ("b", "d"), ("c", "d"))
edge_list = fd.edge_list


def completion(p):
    q = p.copy()
    q.add_edges_from(ADDED_EDGES)
    return q


def whole_graph(p):
    g = p.copy()
    g.add_edges_from(it.product(('p', 'q'), ('a', 'c', 'e')))
    g.add_edges_from(it.product(('r', 'z'), ('b', 'd', 'f')))
    g.add_edges_from((('a', 'b'), ('c', 'd'), ('e', 'f'), ('a', 'f')))
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


def kuratowski_model(witness):
    """Suppress degree-two paths and keep their vertices in connected bags."""
    branches = sorted(v for v in witness if witness.degree(v) != 2)
    bags = {v: [v] for v in branches}
    kernel = nx.Graph()
    kernel.add_nodes_from(branches)
    seen = set()
    for v in branches:
        for w in sorted(witness[v]):
            if frozenset((v, w)) in seen:
                continue
            previous, current = v, w
            internal = []
            seen.add(frozenset((previous, current)))
            while current not in bags:
                internal.append(current)
                nxt = next(x for x in witness[current] if x != previous)
                previous, current = current, nxt
                seen.add(frozenset((previous, current)))
            kernel.add_edge(v, current)
            bags[v].extend(internal)
    assert len(branches) in (5, 6)
    return kernel, bags


def interior_nonplanarity(p, catalog):
    interior = p.subgraph(set(p) - set(ROOTS))
    if nx.check_planarity(interior)[0]:
        return None
    _, witness = nx.check_planarity(interior, counterexample=True)
    kernel, bags = kuratowski_model(witness)
    if len(kernel) == 5:
        return {'type': 'k5_minor', 'graph': 'whole',
                'branch_sets': [bags[v] for v in sorted(kernel)], 'search': 'interior Kuratowski K5'}
    # The interior K3,3 and the exterior Psi1 are vertex-disjoint: two
    # disjoint K3,3 models, a model of the disconnected graph of the catalogue.
    outside = {'A': ['a'], 'C': ['c'], 'E': ['e'], 'P': ['p'], 'Q': ['q'], 'R': ['r', 'b', 'd', 'f']}
    kernel.add_edges_from(it.product(('A', 'C', 'E'), ('P', 'Q', 'R')))
    bags.update(outside)
    for name, target in catalog:
        if len(target) == 12 and not nx.is_connected(target):
            mapping = first_isomorphism(kernel, target)
            if mapping is not None:
                return {'type': 'proper_obstruction_minor', 'target': name,
                        'branch_sets': {str(mapping[v]): bag for v, bag in bags.items()},
                        'search': 'two vertex-disjoint Kuratowski K3,3 models'}
    raise ValueError('No matching catalogue record for two disjoint K3,3 graphs')


def joint_minor(g):
    possible = [['a', 'f']] + [e for e in edge_list(g) if tuple(e) not in CORE_EDGES and e != ['a', 'f']]
    for kind in ('delete_edge', 'contract_edge'):
        for edge in possible:
            if kind == 'delete_edge':
                h = g.copy()
                h.remove_edge(*edge)
                operation = {'type': kind, 'edge': edge}
            else:
                a, b = edge
                if a in CORE_VERTICES and b in CORE_VERTICES:
                    continue
                if b in CORE_VERTICES:
                    a, b = b, a
                h = nx.contracted_nodes(g, a, b, self_loops=False)
                operation = {'type': kind, 'edge': edge, 'keep': a}
            test = joint_test(h, certify=False, stop_on_pass=True)
            if not test['passes']:
                return {'type': 'proper_joint_minor', 'operation': operation,
                        'minor_vertices': sorted(h), 'minor_edges': edge_list(h),
                        'joint_test': joint_test(h, certify=True)}
    return None


def classify(p, catalog):
    g = whole_graph(p)
    for v in sorted(g):
        if g.degree(v) == 3:
            for a, b in it.combinations(sorted(g[v]), 2):
                if g.has_edge(a, b):
                    h = g.copy()
                    h.remove_edge(a, b)
                    return {'type': 'adjacent_degree_three_neighbors', 'vertex': v,
                            'edge': [a, b], 'minor_vertices': sorted(h), 'minor_edges': edge_list(h)}
    for v in sorted(set(p) - set(ROOTS)):
        if g.degree(v) == 2:
            keep = sorted(g[v])[0]
            h = nx.contracted_nodes(g, keep, v, self_loops=False)
            return {'type': 'degree_two_suppression', 'vertex': v,
                    'operation': {'type': 'contract_edge', 'edge': [keep, v], 'keep': keep},
                    'minor_vertices': sorted(h), 'minor_edges': edge_list(h)}
    witness = interior_nonplanarity(p, catalog)
    if witness is not None:
        return witness
    for operation, q in allowed_operations(p):
        j = completion(q)
        if not nx.check_planarity(j)[0]:
            return {'type': 'proper_allowed_minor', 'operation': operation,
                    'minor_vertices': sorted(q), 'minor_edges': edge_list(q),
                    'nonplanarity': nonplanarity(j)}
    witness = joint_minor(g)
    if witness is not None:
        return witness
    for graph_name, h in (('completion', completion(p)), ('whole', g)):
        stats = {}
        model = find_k5_minor(h, trials=64, seed=0, exact_threshold=8,
                              use_sat=True, sat_conflicts=5000, sat_rounds=50, stats=stats)
        if model is not None:
            return {'type': 'k5_minor', 'graph': graph_name, 'branch_sets': model, 'search': stats}
    return {'type': 'unresolved'}


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    records, stats = fd.domain(ROOTS, ADDED_EDGES, 'E3-')
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
            record['outcome'] = classify(p, catalog)
            pf.write(json.dumps({'id': record['id'], 'nonplanarity': record['nonplanarity'],
                                 'outcome': record['outcome']}, sort_keys=True, separators=(',', ':')) + '\n')
            pf.flush()
            counts[record['outcome']['type']] += 1
            if time.time() - last > 30:
                print(json.dumps({'classified': sum(counts.values()), 'of': len(records), 'outcomes': dict(counts)}),
                      flush=True)
                last = time.time()
    document = {'item': 'E3', 'candidate_count': len(records), 'outcome_counts': dict(counts),
                'candidates': records}
    (out / 'candidates.json').write_text(json.dumps(document, sort_keys=True, separators=(',', ':')) + '\n')
    summary = dict(stats, item='E3', where='Proposition 4.21', outcome_counts=dict(counts),
                   unresolved=counts.get('unresolved', 0), seconds=round(time.time() - t0, 1))
    (out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
