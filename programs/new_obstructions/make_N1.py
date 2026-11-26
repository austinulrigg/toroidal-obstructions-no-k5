"""Generate N1 from the 120 published website records.

    python programs/new_obstructions/make_N1.py OUTPUT_FOLDER [--recompute | --page]

Retain the original N1.1-N1.13 records, compute canonical forms with nauty
and compare graph invariants with the Myrvold-Woodcock topological collection.
With --recompute, multi_genus computes the genus of each graph, each edge
deletion, every edge contraction of a minor obstruction, and one nontoroidal
contraction of each remaining graph. With --page, also run PAGE on the 120
graphs. All three runs write the published witnesses."""
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import networkx as nx

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'genus'))
from common import decode_rotations, genus_of_rotation, multicode, tool, zero_indexed  # noqa: E402



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
    import gzip
    import hashlib

    out_root = Path(sys.argv[1])
    source_root = REPO / 'supplemental_data' / 'new_obstructions'
    with gzip.open(source_root / 'N1_sources.json.gz', 'rt') as f:
        sources = json.load(f)
    raw = sources['source evidence']
    website = {key: json.loads(value) for key, value in raw.items()}
    legacy = {record['graph6']: record for record in sources['legacy records']}
    first_ids = [next(key for key, value in website.items() if value['graph6'] == record['graph6'])
                 for record in sources['legacy records']]
    order = first_ids + [key for key in sorted(website) if key not in first_ids]
    graph_strings = [website[key]['graph6'] for key in order]
    forms = canonical(graph_strings)
    wanted = {invariant(from_graph6(s)) for s in graph_strings}
    matches = {key: [] for key in wanted}
    with gzip.open(source_root / 'N1_myrvold_topological.txt.gz', 'rt') as f:
        for index, line in enumerate(f, 1):
            G = mw_graph(line)
            key = invariant(G)
            if key in matches:
                matches[key].append(index)
    assert index == 250815

    def embedding(operation, adj):
        component = operation['components'][0]
        assert len(operation['components']) == 1
        nodes = sorted(adj)
        rotation = {str(nodes[i]): [nodes[j] for j in row]
                    for i, row in enumerate(component['rotation'])}
        return {'genus': operation['genus'], 'rotation': rotation}

    records = []
    for index, (key, form) in enumerate(zip(order, forms), 1):
        report = website[key]
        G = from_graph6(report['graph6'])
        is_minor = report['minorObstruction']
        if report['graph6'] in legacy:
            record = dict(legacy[report['graph6']])
        else:
            record = {'name': 'N1.%d' % index, 'n': len(G), 'm': len(edges(G)),
                      'graph6': report['graph6'], 'edges': [list(e) for e in edges(G)],
                      'connectivity': connectivity(G),
                      'genus': {report['checks']['graph']['components'][0]['method'].split()[0]: 2},
                      'genus-two embedding': embedding(report['checks']['graph'], G),
                      'G - e': [], 'G / e': [], 'canonical graph6': form}
            for operation in report['checks']['deletions']:
                record['G - e'].append({'edge': operation['edge'],
                    **embedding(operation, delete_edge(G, operation['edge']))})
            if is_minor:
                for operation in report['checks']['contractions']:
                    record['G / e'].append({'edge': operation['edge'],
                        **embedding(operation, contract_edge(G, operation['edge']))})
            if key in sources['K5 models']:
                record['K5 branch sets'] = sources['K5 models'][key]
            else:
                assert key == 'TO-250816'
                record['isomorphic to'] = 'T74'
                record['isomorphism to T74'] = {str(v): v for v in G}
        record.update({'website id': key, 'minor obstruction': is_minor,
                       'source evidence SHA256': hashlib.sha256(raw[key].encode()).hexdigest(),
                       'source evidence URL': 'https://austinulrigg.com/data/toroidal/2026-11-11/evidence/' + key + '.json',
                       'topological collection entries with the same invariants': matches[invariant(G)]})
        if not is_minor:
            record['nontoroidal contraction'] = next(operation for operation in report['checks']['contractions']
                                                     if operation['genus'] == 2)
        if '--recompute' in sys.argv[2:] or '--page' in sys.argv[2:]:
            operations = [G] + [delete_edge(G, e) for e in edges(G)]
            expected = [2] + [1] * len(edges(G))
            if is_minor:
                operations.extend(contract_edge(G, e) for e in edges(G))
                expected.extend([1] * len(edges(G)))
            else:
                operations.append(contract_edge(G, record['nontoroidal contraction']['edge']))
                expected.append(2)
            results = multi_genus_embeddings(operations)
            assert all(actual == want for (actual, _), want in zip(results, expected)), key
            if '--page' in sys.argv[2:]:
                assert page_genus(G)[0] == 2, key
        records.append(record)
        print('%s %s %s' % (record['name'], key, 'minor' if is_minor else 'topological'), flush=True)
    item = {
        'item': 'N1', 'schema version': 2,
        'what': 'The 120 topological obstructions in the website release 2026-11-11 absent from the '
                '250,815-graph Myrvold collection. Exactly 22 are minor obstructions. N1.1-N1.13 retain '
                'their original records. N1.14 is T74. Every other graph has recorded K5 branch sets. '
                'The original website evidence and completed exact genus search records are bundled in '
                'N1_sources.json.gz. Their input and output checks do not rerun the exhaustive searches.',
        'counts': {'topological obstructions': 120, 'minor obstructions': 22,
                   'topological but not minor obstructions': 98, 'graphs with a K5 minor': 119},
        'website release': sources['website release'], 'website commit': sources['website commit'],
        'comparison collection': {'file': 'N1_myrvold_topological.txt.gz', 'graphs': 250815,
            'sha256': hashlib.sha256((source_root / 'N1_myrvold_topological.txt.gz').read_bytes()).hexdigest()},
        'canonical forms': {'program': 'nauty labelg', 'all distinct': True},
        'graphs': records,
    }
    folder = out_root / 'new_obstructions'
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'N1.json').write_text(json.dumps(item) + '\n')
    print('wrote', folder / 'N1.json')


if __name__ == '__main__':
    main()
