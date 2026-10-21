"""E2: two four-attachment bridges, or one with another bridge on three attachments (Proposition 4.22).

The names: a = s1, b = t1, c = s2, d = t2, e = s3, f = t3, with twins p, q on
{a, c, e} and r, z on {b, d, f}; the paths P1, P2, P3 are the edges ab, cd,
ef. As in the proof, the bridge B with four attachments a, b, c, d has
interior vertices u, p', q' (the program writes C_u, C_d1, C_d2), with u
adjacent to a, c and p', q' adjacent to b, d, u, or the same with A1 and A2
interchanged. The rest of G is one of:

  another component with three attachments, shortened to T(g') with g' an
  end of the edge ab, cd or ef in its attachment triple (8 triples, 2
  choices of g'; interior D_u, D_d1, D_d2);
  a component with four attachments {a, b, e, f} or {c, d, e, f}, of the
  same form as B (2 forms);
  one of the edges af, eb together with one of cf, ed.

With the two forms of B this gives 48 graphs. For each, the program runs
the joint test of Lemma 4.3 (joint_test.py) on the graph itself, and then
on the graphs obtained by deleting or contracting one edge outside Psi_1,
in order, until one fails all twelve choices: a proper nontoroidal minor.
The proof finds such a minor for every graph.

Run as: python generate_E2.py FOLDER
"""
import collections
import itertools as it
import json
import sys
import time
from pathlib import Path

import networkx as nx

from joint_test import (CORE_EDGES, CORE_VERTICES, core_graph, edges,
                                   graph_record, joint_test, preserved_core)

S = set('ace')
T = set('bdf')
TRIPLES = (('abe', 'ab'), ('abf', 'ab'), ('cde', 'cd'), ('cdf', 'cd'),
           ('aef', 'ef'), ('bef', 'ef'), ('cef', 'ef'), ('def', 'ef'))
SYMMETRY_GENERATORS = {
    'swap_S_T': dict(zip('abcdefpqrz', 'badcferzpq')),
    'swap_indices_1_2': dict(zip('abcdefpqrz', 'cdabefpqrz')),
}


def add_four_component(graph, support, orientation, prefix):
    u, d1, d2 = (f'{prefix}_{suffix}' for suffix in ('u', 'd1', 'd2'))
    centers = S if orientation == 'S_center' else T
    assert orientation in ('S_center', 'T_center')
    center_roots = sorted(set(support) & centers)
    opposite_roots = sorted(set(support) - centers)
    assert len(center_roots) == len(opposite_roots) == 2
    graph.add_edges_from((u, v) for v in center_roots + [d1, d2])
    graph.add_edges_from((d, v) for d in (d1, d2) for v in opposite_roots)


def add_triple_component(graph, support, distinguished_root, prefix):
    u, d1, d2 = (f'{prefix}_{suffix}' for suffix in ('u', 'd1', 'd2'))
    others = sorted(set(support) - {distinguished_root})
    graph.add_edge(u, distinguished_root)
    graph.add_edges_from((d, v) for d in (d1, d2) for v in others + [u])


def generate_candidates():
    records = []
    for orientation in ('S_center', 'T_center'):
        base = core_graph()
        add_four_component(base, 'abcd', orientation, 'C')
        alternatives = []
        for support, distinguished_roots in TRIPLES:
            for distinguished_root in distinguished_roots:
                graph = base.copy()
                add_triple_component(graph, support, distinguished_root, 'D')
                alternatives.append((graph, {
                    'alternative_kind': 'triple_component', 'support': list(support),
                    'distinguished_root': distinguished_root}))
        for support in ('abef', 'cdef'):
            for second_orientation in ('S_center', 'T_center'):
                graph = base.copy()
                add_four_component(graph, support, second_orientation, 'D')
                alternatives.append((graph, {
                    'alternative_kind': 'four_attachment_component', 'support': list(support),
                    'D_orientation': second_orientation}))
        for first, second in it.product(('af', 'eb'), ('cf', 'ed')):
            graph = base.copy()
            graph.add_edges_from((tuple(first), tuple(second)))
            alternatives.append((graph, {
                'alternative_kind': 'cross_edges',
                'cross_edges': sorted(tuple(sorted(e)) for e in (first, second))}))
        for graph, recipe in alternatives:
            recipe['C_orientation'] = orientation
            records.append({'id': f'E2-{len(records)+1:03d}',
                            'recipe': recipe, 'graph': graph_record(graph)})
    return records


def graph_key(record):
    return tuple(record['vertices']), tuple(tuple(e) for e in record['edges'])


def source_symmetry_group():
    identity = {v: v for v in CORE_VERTICES}
    st = SYMMETRY_GENERATORS['swap_S_T']
    ij = SYMMETRY_GENERATORS['swap_indices_1_2']
    return [identity, st, ij, {v: st[ij[v]] for v in CORE_VERTICES}]


def symmetry_orbits(records):
    by_key = {graph_key(r['graph']): r for r in records}
    assert len(by_key) == len(records)
    group = source_symmetry_group()
    seen = set()
    orbits = []
    for representative in records:
        if representative['id'] in seen:
            continue
        members = {}
        graph = nx.Graph()
        graph.add_nodes_from(representative['graph']['vertices'])
        graph.add_edges_from(representative['graph']['edges'])
        for core_mapping in group:
            vertex_map = {v: core_mapping.get(v, v) for v in graph}
            transformed = nx.relabel_nodes(graph, vertex_map)
            member = by_key[graph_key(graph_record(transformed))]
            # The four group elements are involutions.  Record member -> rep.
            inverse = {v: u for u, v in vertex_map.items()}
            members.setdefault(member['id'], {'candidate_id': member['id'],
                                              'vertex_map': inverse})
        orbit_id = f'E2-orbit-{len(orbits)+1:02d}'
        for member in members:
            seen.add(member)
        orbits.append({'id': orbit_id, 'representative_id': representative['id'],
                       'vertex_map_direction': 'member_to_representative',
                       'members': [members[k] for k in sorted(members)]})
    for orbit in orbits:
        for member in orbit['members']:
            next(r for r in records if r['id'] == member['candidate_id'])['orbit_id'] = orbit['id']
    return orbits


def proper_operations(graph):
    outside = [e for e in edges(graph) if e not in CORE_EDGES]
    for edge in outside:
        minor = graph.copy()
        minor.remove_edge(*edge)
        assert preserved_core(minor)
        yield {'type': 'delete_edge', 'edge': edge}, minor
    for edge in outside:
        u, v = edge
        if u in CORE_VERTICES and v in CORE_VERTICES:
            continue
        if v in CORE_VERTICES:
            u, v = v, u
        # Explicit vertex map avoids preserving NetworkX contraction metadata.
        mapping = {w: u if w == v else w for w in graph}
        minor = nx.Graph()
        minor.add_nodes_from(mapping.values())
        minor.add_edges_from((mapping[a], mapping[b]) for a, b in graph.edges()
                             if mapping[a] != mapping[b])
        assert preserved_core(minor)
        yield {'type': 'contract_edge', 'edge': list(edge), 'keep': u, 'remove': v,
               'vertex_map': mapping}, minor


def main():
    out = Path(sys.argv[1])
    start = time.monotonic()
    records = generate_candidates()
    orbits = symmetry_orbits(records)
    assert len(records) == 48 and len(orbits) == 13
    unresolved = []
    for record in records:
        graph = nx.Graph()
        graph.add_nodes_from(record['graph']['vertices'])
        graph.add_edges_from(record['graph']['edges'])
        original = joint_test(graph, certify=True)
        record['original_joint_test'] = original
        if original['passes']:
            record['status'] = 'unresolved_original_has_passing_choice'
            unresolved.append(record['id'])
            continue
        attempts = []
        for operation, minor in proper_operations(graph):
            trial = joint_test(minor, certify=False, stop_on_pass=True)
            attempts.append({'operation': operation, 'has_passing_choice': trial['passes']})
            if not trial['passes']:
                certificate = joint_test(minor, certify=True)
                assert certificate['complete'] and not certificate['passes']
                record['proper_minor'] = {
                    'operation': operation, 'graph': graph_record(minor),
                    'core_preserved': preserved_core(minor),
                    'bridges': certificate['bridges'], 'joint_test': certificate}
                record['status'] = 'proper_minor_all_twelve_fail'
                break
        record['minor_search_attempts'] = attempts
        if 'proper_minor' not in record:
            record['status'] = 'unresolved_no_core_preserving_one_edge_witness'
            unresolved.append(record['id'])
        print(f"{record['id']}: {record['status']}, tried {len(attempts)}", flush=True)
    counts = {
        'labeled graphs': len(records), 'symmetry orbits': len(orbits),
        'recipe_families': dict(collections.Counter(r['recipe']['alternative_kind'] for r in records)),
        'orbit_sizes': dict(collections.Counter(len(o['members']) for o in orbits)),
        'originals_all_twelve_fail': sum(not r['original_joint_test']['passes'] for r in records),
        'proper_minor_witnesses': sum('proper_minor' in r for r in records),
        'unresolved': len(unresolved),
    }
    result = {
        'item': 'E2', 'core': graph_record(core_graph()),
        'symmetry_generators': SYMMETRY_GENERATORS,
        'counts': counts, 'orbits': orbits, 'candidates': records,
        'unresolved_candidate_ids': unresolved,
        'minor_search_scope': 'One edge deletion or contraction outside the fixed Psi1; no two core vertices identified',
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / 'candidates.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    summary = {'item': 'E2', 'where': 'Proposition 4.22', 'counts': counts,
               'unresolved_candidate_ids': unresolved, 'seconds': round(time.monotonic() - start, 1)}
    (out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
