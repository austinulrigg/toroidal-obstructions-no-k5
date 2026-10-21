"""E12: a cut {a, u, v} with one component containing b, c, d (Proposition 5.13).

The names: Psi_2 has twins p1, p2 on abc and q1, q2 on abd (p, q, r, z in
the paper), and the bridge B has attachments a, b, c, d. It contains u, v
and the twins p, q with neighbourhood X = {a, u, v} (p', q' in the paper),
and one component D of B - X, containing b, c, d, with N_B(D) = {u, v}. As
in the proof, the components of D - {b, c, d} are contracted to vertices,
and D is reduced to a bipartite tree with them, which is a star whose
center z1 is joined to b, c, d, or a path g1 z1 g2 z2 g3 with
{g1, g2, g3} = {b, c, d} (the paper writes m, or m1, m2). For the star, z1
meets both u, v, or one of them and a vertex of {b, c, d} meets the other.
For the path, z1, z2 meet different vertices of {u, v}, or both meet one
and a vertex of {b, c, d} meets the other. Every such choice is made, and
the bridges are compared up to permutations of the vertices of the tree
outside b, c, d (all eight named vertices fixed): 31 forms.

For each form the program records a plane embedding of B and Kuratowski
subgraphs of the four tests of Lemma 5.3 with (n_c, n_d) = (2,2), and
classifies G = Psi_2 + B: a cut of at most two vertices; an isomorphism
onto some T_i; a K5 minor (found by a random search); or the deletion of
one edge of B after which all four tests are still nonplanar. The proof
finds T40, T53, T56, T73.

Run as: python generate_E12.py FOLDER
"""
import collections
import itertools as it
import json
import sys
from pathlib import Path

import networkx as nx

import common
import graphs as gr

FIXED = ('a', 'b', 'c', 'd', 'u', 'v', 'p', 'q')
TESTS = [('ac', 'bd'), ('ad', 'bc'), ('ab', 'ac', 'ad'), ('ab', 'bc', 'bd')]


def edges(graph):
    return sorted(sorted(e) for e in graph.edges())


def canonical_part(P):
    """Fix all eight named vertices and permute only tree-interior vertices."""
    zs = sorted(set(P) - set(FIXED))
    choices = []
    for perm in it.permutations(zs):
        renaming = dict(zip(zs, perm))
        Q = nx.relabel_nodes(P, renaming)
        choices.append(tuple(tuple(e) for e in edges(Q)))
    return min(choices)


def base_part():
    P = nx.Graph()
    P.add_nodes_from(FIXED)
    P.add_edges_from((p, x) for p in ('p', 'q') for x in ('a', 'u', 'v'))
    return P


def recipes():
    # One z joined to both u,v, or to one and an attachment to the other.
    tree = [('z1', r) for r in ('b', 'c', 'd')]
    yield {'family': 'star', 'joins': [['z1', 'u'], ['z1', 'v']]}, tree
    for x, y in [('u', 'v'), ('v', 'u')]:
        for r in ('b', 'c', 'd'):
            yield {'family': 'star', 'joins': [['z1', x], [r, y]]}, tree
    # Every ordered path; recipes giving the same form are all recorded.
    for r1, r2, r3 in it.permutations(('b', 'c', 'd')):
        tree = [(r1, 'z1'), ('z1', r2), (r2, 'z2'), ('z2', r3)]
        for x, y in [('u', 'v'), ('v', 'u')]:
            yield {'family': 'path', 'order': [r1, r2, r3],
                   'joins': [['z1', x], ['z2', y]]}, tree
            for r in ('b', 'c', 'd'):
                yield {'family': 'path', 'order': [r1, r2, r3],
                       'joins': [['z1', x], ['z2', x], [r, y]]}, tree


def generate():
    table = {}
    for recipe, tree in recipes():
        P = base_part()
        P.add_edges_from(tree)
        P.add_edges_from(recipe['joins'])
        key = canonical_part(P)
        if key not in table:
            Q = nx.Graph()
            Q.add_nodes_from(P)
            Q.add_edges_from(key)
            table[key] = (Q, [])
        table[key][1].append(recipe)
    return [table[k] for k in sorted(table)]


def whole_graph(P):
    G = P.copy()
    G.add_edges_from((p, r) for p in ('p1', 'p2') for r in ('a', 'b', 'c'))
    G.add_edges_from((q, r) for q in ('q1', 'q2') for r in ('a', 'b', 'd'))
    return G


def planar_certificate(G):
    planar, witness = nx.check_planarity(G, counterexample=True)
    if planar:
        return {'planar': True, 'rotation': {v: list(witness.neighbors_cw_order(v)) for v in G}}
    return {'planar': False, 'kuratowski_edges': edges(witness)}


def connectivity_witness(G):
    for size in range(3):
        for cut in it.combinations(sorted(G), size):
            H = G.copy()
            H.remove_nodes_from(cut)
            comps = sorted([sorted(c) for c in nx.connected_components(H)])
            if len(comps) > 1:
                return {'cut': list(cut), 'components': comps}
    return None


def proper_four_test_subgraph(P):
    """An edge of B whose deletion keeps the hypotheses of Lemma 5.3 and all
    four tests nonplanar: a proper nontoroidal minor of G."""
    roots = set('abcd')
    interior = set(P) - roots
    for edge in edges(P):
        Q = P.copy()
        Q.remove_edge(*edge)
        if not nx.is_connected(Q.subgraph(interior)):
            continue
        if any(not (set(Q[r]) & interior) for r in roots):
            continue
        tests = []
        for pairs in TESTS:
            H = Q.copy()
            H.add_edges_from(tuple(e) for e in pairs)
            if nx.check_planarity(H)[0]:
                break
            tests.append({'added_edges': pairs, **planar_certificate(H)})
        if len(tests) == 4:
            return {'type': 'proper_four_test_subgraph', 'deleted_edges': [edge],
                    'part_edges': edges(Q), 'completions': tests}
    return None


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    entries = common.catalogue_records()
    records = []
    for index, (P, recipe_list) in enumerate(generate(), 1):
        G = whole_graph(P)
        completions = []
        for pairs in TESTS:
            H = P.copy()
            H.add_edges_from(tuple(pair) for pair in pairs)
            completions.append({'added_edges': pairs, **planar_certificate(H)})
        matches = []
        for entry in entries:
            H = nx.from_graph6_bytes(entry['graph6'].encode())
            if len(G) != len(H) or G.number_of_edges() != H.number_of_edges():
                continue
            matcher = nx.isomorphism.GraphMatcher(G, H)
            if matcher.is_isomorphic():
                matches.append({'target': entry['name'], 'mapping': matcher.mapping})
        cut = connectivity_witness(G)
        if cut:
            outcome = {'type': 'small_cut', **cut}
        elif matches:
            outcome = {'type': 'survivor', **matches[0]}
        else:
            sets = gr.find_k5_minor(G, trials=1000, seed=0)
            if sets:
                outcome = {'type': 'k5_minor', 'branch_sets': sets}
            else:
                outcome = proper_four_test_subgraph(P) or {'type': 'unresolved'}
        print('E12-%03d: %s %s' % (index, outcome['type'], outcome.get('target', '')), flush=True)
        records.append({'id': 'E12-%03d' % index, 'vertices': sorted(P), 'part_edges': edges(P),
                        'recipes': recipe_list, 'graph_vertices': sorted(G), 'graph_edges': edges(G),
                        'planar_part': planar_certificate(P), 'completions': completions,
                        'small_cut': cut, 'catalog_matches': matches, 'outcome': outcome})
    counts = collections.Counter(r['outcome']['type'] for r in records)
    result = {'item': 'E12', 'raw_recipes': sum(len(r['recipes']) for r in records),
              'forms': len(records), 'candidates': records}
    (out / 'candidates.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    summary = {'item': 'E12', 'where': 'Proposition 5.13', 'recipes': result['raw_recipes'], 'forms': len(records),
               'outcomes': dict(counts),
               'catalogue graphs found': sorted({r['outcome']['target'] for r in records
                                                 if r['outcome']['type'] == 'survivor'}, key=lambda x: int(x[1:])),
               'complete': counts.get('unresolved', 0) == 0}
    (out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
