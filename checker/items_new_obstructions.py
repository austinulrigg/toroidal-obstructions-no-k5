"""Check N1's 120 topological obstructions and their recorded computations.

Rotation systems are checked independently. The completed PAGE and multi_genus
records are checked for input identity, search mode and consistent output.
Their exhaustive searches can be rerun with programs/new_obstructions/make_N1.py.
"""
import gzip
import hashlib
import json
import os
import re

import networkx as nx

from witnesses import (check_embedding, check_iso, check_k5, check_no_k5,
                       contract_edge, delete_edge, edges, from_graph6, load,
                       n_edges, need, same_graph)


def read_sources(root):
    with gzip.open(os.path.join(root, 'new_obstructions', 'N1_sources.json.gz'), 'rt') as f:
        return json.load(f)


def int_rotation(rotation):
    return {int(v): nb for v, nb in rotation.items()}


def check_search(record):
    """Check a completed search's recorded input and output, without rerunning it."""
    G = from_graph6(record['graph6'])
    rotation = record['rotation']
    check_embedding(G, {'type': 'embedding', 'genus': record['genus'],
                        'rotation': dict(enumerate(rotation))}, genus=record['genus'])
    need(re.fullmatch('[0-9a-f]{64}', record['binarySha256']) is not None, 'missing binary checksum')
    if record['method'] == 'PAGE completed exact search':
        degree = max(map(len, G.values()))
        lines = ['%d %d' % (len(G), n_edges(G))]
        for v in G:
            lines.append(' '.join(map(str, sorted(G[v]) + [65535] * (degree - len(G[v])))))
        data = ('\n'.join(lines) + '\n').encode()
        need(record['settings'] == {'STDOUT': 1, 'PBN': 1, 'S': 0, 'DEG': degree}, 'PAGE search settings')
        match = re.search(r'Solution with (\d+) cycles \(genus (\d+)\) found', record['stdout'])
        need(match is not None and int(match[2]) == record['genus'], 'PAGE output genus')
        need(int(match[1]) == 2 - 2 * record['genus'] - len(G) + n_edges(G), 'PAGE output face count')
        need('The genus is %d.' % record['genus'] in record['stderr'], 'PAGE did not report a completed exact search')
    else:
        need(record['method'] == 'multi_genus completed exact search', 'unrecognized genus computation')
        need(record['arguments'] == ['w'], 'multi_genus search mode')
        values = [len(G)]
        for v in range(len(G) - 1):
            values.extend(u + 1 for u in sorted(G[v]) if u > v)
            values.append(0)
        data = bytes(values)
        output = bytes.fromhex(record['stdoutHex'])
        need(output[:1] == bytes([len(G)]) and output.endswith(b'\0'), 'multi_genus embedding format')
        rows = [[u - 1 for u in row] for row in output[1:-1].split(b'\0')]
        need(rows == rotation, 'multi_genus output rotation')
        counts = re.findall(r'^graphs with genus (\d+): (\d+)$', record['stderr'], re.MULTILINE)
        need(counts == [(str(record['genus']), '1')] and 'Wrote 1 embeddings.' in record['stderr'],
             'multi_genus did not report one completed exact search')
    need(hashlib.sha256(data).hexdigest() == record['inputSha256'], 'genus computation input checksum')


def check_evidence(G, record):
    """Website records use consecutive labels after contracting an edge."""
    nodes = sorted(G)
    normalized = {i: {nodes.index(u) for u in G[v]} for i, v in enumerate(nodes)}
    need(len(record['components']) == 1, 'N1 expects connected operation graphs')
    component = record['components'][0]
    need(component['labels'] == list(range(len(G))), 'component labels')
    need(same_graph(normalized, from_graph6(component['graph6'])), 'recorded computation is for another graph')
    need(record['genus'] == component['genus'], 'component genus')
    check_search(component)


def check(root):
    from items_obstructions import canonical, connectivity, distinct, invariant, labelg, mw_graph, nxg, to_graph6
    folder = os.path.join(root, 'new_obstructions')
    data = load(os.path.join(folder, 'N1.json'))
    sources = read_sources(root)
    need(data['website release'] == sources['website release'] == sources['source metadata']['release'],
         'website release provenance')
    need(data['website commit'] == sources['website commit'] == sources['evidence provenance']['remote_commit'],
         'website commit provenance')
    baseline = next(item for item in sources['source metadata']['sourceFiles'] if item['file'] == 'all_obs.txt.gz')
    baseline_path = os.path.join(folder, 'N1_myrvold_topological.txt.gz')
    with open(baseline_path, 'rb') as f:
        baseline_bytes = f.read()
    baseline_hash = hashlib.sha256(baseline_bytes).hexdigest()
    need(baseline_hash == baseline['sha256'] and len(baseline_bytes) == baseline['bytes'],
         'Myrvold source checksum')
    need(data['comparison collection'] == {'file': 'N1_myrvold_topological.txt.gz',
                                           'graphs': baseline['records'], 'sha256': baseline_hash},
         'Myrvold comparison provenance')
    need(baseline['records'] == 250815, 'Myrvold source record count')
    legacy = sources['legacy records']
    need([g['name'] for g in legacy] == ['N1.%d' % i for i in range(1, 14)], 'original N1 record names')
    provenance = {item['file']: item for item in sources['evidence provenance']['files']}
    records = data['graphs']
    need([g['name'] for g in records] == ['N1.%d' % i for i in range(1, 121)], 'N1 graph names')
    need(data['counts'] == {'topological obstructions': 120, 'minor obstructions': 22,
                            'topological but not minor obstructions': 98, 'graphs with a K5 minor': 119}, 'N1 counts')
    expected_ids = {'TO-%d' % i for i in range(250816, 250936)}
    need({g['website id'] for g in records} == expected_ids, 'N1 website identifiers')
    need(set(sources['source evidence']) == expected_ids, 'N1 source evidence identifiers')
    for old, new in zip(legacy, records[:13]):
        need(all(new.get(k) == v for k, v in old.items()), 'N1 changed an original record', old['name'])
    graphs = []
    n_embeddings = n_contractions = n_minor = n_k5 = 0
    o1 = load(os.path.join(root, 'no_K5_obstructions', 'O1.json'))
    o2 = load(os.path.join(root, 'no_K5_obstructions', 'O2.json'))
    for record in records:
        name = record['name']
        G = from_graph6(record['graph6'])
        need(len(G) == record['n'] and n_edges(G) == record['m'], name, 'size')
        need({frozenset(e) for e in record['edges']} == {frozenset(e) for e in edges(G)} and len(record['edges']) == n_edges(G), name, 'edge list')
        need(min(map(len, G.values())) >= 3, name, 'minimum degree')
        need(record['connectivity'] == connectivity(G), name, 'connectivity')
        check_embedding(G, {'type': 'embedding', 'genus': 2,
                            'rotation': int_rotation(record['genus-two embedding']['rotation'])}, genus=2)
        for key, operation in [('G - e', delete_edge), ('G / e', contract_edge)]:
            seen = []
            for result in record[key]:
                edge = tuple(result['edge'])
                seen.append(frozenset(edge))
                check_embedding(operation(G, edge), {'type': 'embedding', 'genus': result['genus'],
                                'rotation': int_rotation(result['rotation'])}, at_most=1)
                n_embeddings += 1
            need(len(set(seen)) == len(seen), name, 'repeated operation')
            if key == 'G - e' or record['minor obstruction']:
                need(set(seen) == {frozenset(e) for e in edges(G)}, name, 'missing edge operations', key)
        if record['minor obstruction']:
            n_minor += 1
        else:
            failed = record['nontoroidal contraction']
            check_evidence(contract_edge(G, tuple(failed['edge'])), failed)
            need(failed['genus'] == 2, name, 'nonminor witness genus')
            n_contractions += 1
        if 'K5 branch sets' in record:
            check_k5(G, record['K5 branch sets'])
            n_k5 += 1
        else:
            need(record['isomorphic to'] == 'T74', name, 'unexpected no-K5 graph')
            target = next(g for g in o1['graphs'] if g['name'] == 'T74')
            target_graph = from_graph6(target['graph6'])
            check_iso(G, target_graph, {int(v): u for v, u in record['isomorphism to T74'].items()})
            certificate = next(g for g in o2['graphs'] if g['name'] == 'T74')
            check_no_k5(target_graph, certificate['no K5 certificate'])
        raw = sources['source evidence'][record['website id']]
        source_bytes = raw.encode()
        source_hash = hashlib.sha256(source_bytes).hexdigest()
        filename = record['website id'] + '.json'
        manifest = sources['source metadata']['files']['evidence/' + filename]
        need(source_hash == record['source evidence SHA256'] == manifest['sha256'] == provenance[filename]['sha256'],
             name, 'source checksum')
        need(len(source_bytes) == manifest['bytes'], name, 'source byte count')
        git_blob = hashlib.sha1(b'blob ' + str(len(source_bytes)).encode() + b'\0' + source_bytes).hexdigest()
        need(git_blob == provenance[filename]['git_blob'], name, 'source Git blob')
        website = json.loads(raw)
        need(website['graph6'] == record['graph6'] and website['status'] == 'verified', name, 'website source graph')
        need(website['topologicalObstruction'] is True and website['minorObstruction'] == record['minor obstruction'],
             name, 'website source classification')
        check_evidence(G, website['checks']['graph'])
        need(website['checks']['graph']['genus'] == 2, name, 'recorded graph genus')
        graphs.append((name, G))
    need((n_minor, n_contractions, n_k5) == (22, 98, 119), 'N1 class counts')
    distinct(graphs, 'N1:')
    wanted = {invariant(G) for _, G in graphs}
    by_inv = {key: [] for key in wanted}
    baseline_g6 = []
    with gzip.open(os.path.join(folder, 'N1_myrvold_topological.txt.gz'), 'rt') as source:
        for index, line in enumerate(source, 1):
            G = mw_graph(line)
            key = invariant(G)
            if key in by_inv:
                by_inv[key].append((index, G))
            if labelg():
                baseline_g6.append(to_graph6(G))
    need(index == 250815, 'Myrvold topological collection size')
    for record, (name, G) in zip(records, graphs):
        matches = by_inv[invariant(G)]
        need(record['topological collection entries with the same invariants'] == [i for i, _ in matches],
             name, 'topological collection invariant matches')
        for index, H in matches:
            need(not nx.is_isomorphic(nxg(G), nxg(H)), name, 'isomorphic to Myrvold entry', index)
    if labelg():
        forms = canonical([record['graph6'] for record in records])
        need(forms == [record['canonical graph6'] for record in records], 'N1 canonical forms')
        need(not set(forms).intersection(canonical(baseline_g6)), 'N1 canonical form is in Myrvold collection')
    return ('120 topological obstructions, including 22 minor obstructions; %d torus embeddings; '
            '98 recorded genus-two contractions; 119 K5 models and T74; original N1.1-N1.13 preserved; '
            'pairwise nonisomorphic and absent from all 250,815 Myrvold topological obstructions; '
            '120 recorded exact graph searches checked for input/output consistency' % n_embeddings)
