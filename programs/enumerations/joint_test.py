"""The joint test of Lemma 4.3 for a graph containing Psi_1 whose three paths
are the edges ab, cd, ef (a = s1, b = t1, c = s2, d = t2, e = s3, f = t3),
with twins p, q on {a, c, e} and r, z on {b, d, f}.

For each of the twelve choices (k, x, y) of Lemma 4.2 the four faces are:
J, the pair of simple faces (attachments V(P_i) and V(P_j), test graph with
the auxiliary edges s_i s_j, t_i t_j and the two paths); W, the repeated
face (attachments V(P_k) and s_x, t_y, auxiliary edges s_k s_x, t_k t_y);
and small_S, small_T, the two faces of Lemma 4.2(iii). A choice passes when
some assignment of the Psi_1-bridges to faces containing their attachments
makes every test graph planar. Every assignment of a failing choice is
recorded with the face whose test graph is nonplanar and a Kuratowski
subgraph of it.
"""
import itertools as it

import networkx as nx


S = {1: 'a', 2: 'c', 3: 'e'}
T = {1: 'b', 2: 'd', 3: 'f'}
CORE_VERTICES = tuple('abcdefpqrz')
CORE_EDGES = tuple(sorted(tuple(sorted(e)) for e in (
    [(v, s) for v in ('p', 'q') for s in S.values()]
    + [(v, t) for v in ('r', 'z') for t in T.values()]
    + [(S[i], T[i]) for i in S])))
REGION_NAMES = ('J', 'W', 'small_S', 'small_T')


def edges(graph):
    return sorted(tuple(sorted(e)) for e in graph.edges())


def graph_record(graph):
    return {'vertices': sorted(graph), 'edges': edges(graph)}


def core_graph():
    graph = nx.Graph()
    graph.add_nodes_from(CORE_VERTICES)
    graph.add_edges_from(CORE_EDGES)
    return graph


def preserved_core(graph):
    return (set(CORE_VERTICES) <= set(graph)
            and all(graph.has_edge(*e) for e in CORE_EDGES)
            and all(graph.degree(v) == 3 for v in ('p', 'q', 'r', 'z')))


def bridges(graph):
    """Recompute every component bridge and every non-core edge bridge."""
    if not preserved_core(graph):
        raise ValueError('The fixed Psi1 core and twin degrees must be preserved')
    core = set(CORE_VERTICES)
    result = []
    components = sorted(sorted(c) for c in nx.connected_components(
        graph.subgraph(set(graph) - core)))
    for comp in components:
        inside = set(comp)
        attachments = sorted(set().union(*(set(graph[v]) & core for v in comp)))
        bridge_edges = [e for e in edges(graph) if inside.intersection(e)]
        result.append({'id': f'B{len(result):02d}', 'kind': 'component',
                       'interior': comp, 'attachments': attachments,
                       'vertices': sorted(comp + attachments), 'edges': bridge_edges})
    for edge in edges(graph.subgraph(core)):
        if edge not in CORE_EDGES:
            result.append({'id': f'B{len(result):02d}', 'kind': 'edge',
                           'interior': [], 'attachments': list(edge),
                           'vertices': list(edge), 'edges': [edge]})
    assert not any(set(b['attachments']) & set('pqrz') for b in result)
    return result


def choices():
    for k in (1, 2, 3):
        other = [i for i in (1, 2, 3) if i != k]
        for x, y in it.product(other, repeat=2):
            yield k, x, y


def regions(k, x, y):
    i, j = [i for i in (1, 2, 3) if i != k]
    if x not in (i, j) or y not in (i, j):
        raise ValueError('Invalid joint-test choice')
    xp = next(h for h in (i, j) if h != x)
    yp = next(h for h in (i, j) if h != y)
    raw = [
        ('J', [S[i], T[i], S[j], T[j]],
         [(S[i], T[i]), (S[j], T[j]), (S[i], S[j]), (T[i], T[j])]),
        ('W', [S[k], T[k], S[x], T[y]],
         [(S[k], T[k]), (S[k], S[x]), (T[k], T[y])]),
        ('small_S', [S[k], S[xp]], [(S[k], S[xp])]),
        ('small_T', [T[k], T[yp]], [(T[k], T[yp])]),
    ]
    return [{'name': name, 'support': sorted(support),
             'base_edges': sorted(tuple(sorted(e)) for e in base)}
            for name, support, base in raw]


def region_graph(region, assigned_bridges):
    graph = nx.Graph()
    graph.add_nodes_from(region['support'])
    graph.add_edges_from(region['base_edges'])
    for bridge in assigned_bridges:
        graph.add_nodes_from(bridge['vertices'])
        graph.add_edges_from(bridge['edges'])
    return graph


def joint_test(graph, *, certify=True, stop_on_pass=False):
    """Return every assignment and a failing region, or a planar assignment.

    For an all-failing result every compatible bridge assignment is explicitly
    present unless its choice has an empty eligibility set.  A certificate
    includes one Kuratowski subgraph per failed assignment.  No genus solver
    is used.  Empty-attachment components are eligible for all four regions.
    """
    bridge_list = bridges(graph)
    bridge_ids = [b['id'] for b in bridge_list]
    tests = []
    planar_cache = {}
    rotation_cache = {}
    witness_cache = {}
    any_pass = False
    for k, x, y in choices():
        region_list = regions(k, x, y)
        eligibility = {
            b['id']: [r['name'] for r in region_list
                      if set(b['attachments']) <= set(r['support'])]
            for b in bridge_list}
        blocked = [b for b in bridge_ids if not eligibility[b]]
        choice = {'k': k, 'x': x, 'y': y, 'regions': region_list,
                  'eligible_regions': eligibility, 'blocked_bridges': blocked,
                  'assignments': [], 'passes': False}
        for assignment in it.product(*(eligibility[b] for b in bridge_ids)):
            mapping = dict(zip(bridge_ids, assignment))
            failure = None
            planar_regions = []
            planar_rotations = {}
            for region in region_list:
                part = region_graph(region, [b for b in bridge_list
                    if mapping[b['id']] == region['name']])
                key = (tuple(sorted(part)), tuple(edges(part)))
                if key not in planar_cache:
                    planar, embedding = nx.check_planarity(part)
                    planar_cache[key] = planar
                    if planar and certify:
                        rotation_cache[key] = {
                            v: list(embedding.neighbors_cw_order(v)) for v in sorted(part)}
                if not planar_cache[key]:
                    failure = {'region': region['name'], 'graph': graph_record(part)}
                    if certify:
                        if key not in witness_cache:
                            planar, witness = nx.check_planarity(part, counterexample=True)
                            assert not planar
                            witness_cache[key] = edges(witness)
                        failure['kuratowski_edges'] = witness_cache[key]
                    break
                if certify:
                    planar_regions.append(region['name'])
                    planar_rotations[region['name']] = rotation_cache[key]
            if failure is None:
                choice['passes'] = True
                any_pass = True
                record = {'mapping': mapping, 'passes': True}
                if certify:
                    record['planar_regions'] = planar_regions
                    record['planar_rotations'] = planar_rotations
                choice['assignments'].append(record)
                if stop_on_pass:
                    tests.append(choice)
                    return {'bridges': bridge_list, 'choices': tests,
                            'passes': True, 'complete': False}
            else:
                choice['assignments'].append({'mapping': mapping, 'failure': failure,
                                              'passes': False})
        tests.append(choice)
    return {'bridges': bridge_list, 'choices': tests,
            'passes': any_pass, 'complete': True}
