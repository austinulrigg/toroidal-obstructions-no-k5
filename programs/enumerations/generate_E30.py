"""Compute the 24 attachment-edge completions in Lemma 3.10.

Usage: python generate_E30.py OUTPUT_FOLDER

The rotation systems and Kuratowski subgraphs use the vertex names of the
figure. The checker can verify the certificates without running PAGE.
"""
import argparse
import json
from importlib.metadata import version
from pathlib import Path

import networkx as nx
from graph_genus import embed


# In T and S the common neighbors of p and q are x, y, u.
# In Q the common neighbors are u, v, x.
BASE_EDGES = {
    'T(z)': 'px py pu qx qy qu uz'.split(),
    'S(x;y)': 'px py pu qx qy qu xf fu fz'.split(),
    'Q(x;y,z)': 'pu pv px qu qv qx yu vz'.split(),
}
ATTACHMENT_EDGES = ('xy', 'xz', 'yz')


def sorted_edges(graph):
    return sorted(sorted(edge) for edge in graph.edges())


def page_embedding(graph):
    """Prune pendant edges for PAGE, then put them back in its rotation."""
    core = graph.copy()
    removed = []
    while True:
        leaf = next((v for v in core if core.degree(v) == 1), None)
        if leaf is None:
            break
        parent = next(iter(core[leaf]))
        removed.append((leaf, parent))
        core.remove_node(leaf)
    vertices = sorted(core)
    index = {v: i for i, v in enumerate(vertices)}
    adjacency = [sorted(index[u] for u in core[v]) for v in vertices]
    genus, cycles = embed(adjacency, algorithm='page')
    rotations = {v: [vertices[j] for j in row] for v, row in zip(vertices, cycles, strict=True)}
    for leaf, parent in reversed(removed):
        rotations[parent].append(leaf)
        rotations[leaf] = [parent]
    return {'type': 'embedding', 'genus': genus, 'rotation': rotations}, removed


def generate(output_folder):
    candidates = []
    verdicts = []
    counts = {}
    for name, base in BASE_EDGES.items():
        counts[name] = {'planar': 0, 'nonplanar': 0}
        for mask in range(8):
            added = [e for i, e in enumerate(ATTACHMENT_EDGES) if mask & (1 << i)]
            graph = nx.Graph()
            vertices = sorted(set(''.join(base)))
            graph.add_nodes_from(vertices)
            graph.add_edges_from(tuple(e) for e in base + added)
            record_id = name[0] + str(mask)
            candidates.append({
                'id': record_id, 'base': name, 'added_edges': added,
                'vertices': vertices, 'edges': sorted_edges(graph),
                'graph6': nx.to_graph6_bytes(graph, header=False).decode('ascii').strip(),
            })
            embedding, removed = page_embedding(graph)
            genus = embedding['genus']
            planar, certificate = nx.check_planarity(graph, counterexample=True)
            if planar != (genus == 0):
                raise ValueError(f'{record_id}: PAGE genus {genus} disagrees with planarity test')
            row = {'id': record_id, 'genus': genus, 'planar': planar,
                   'pendant_edges_removed_for_PAGE': removed, 'embedding': embedding}
            if not planar:
                row['kuratowski'] = {'type': 'Kuratowski subgraph', 'edges': sorted_edges(certificate)}
            verdicts.append(row)
            counts[name]['planar' if planar else 'nonplanar'] += 1

    inputs = {
        'item': 'E30', 'where': 'Lemma 3.10 (core:three-graph-completions)',
        'attachment_edges': list(ATTACHMENT_EDGES), 'base_edges': BASE_EDGES,
        'id_rule': 'The letter identifies the base graph. The digit is the sum of 1 for xy, 2 for xz and 4 for yz.',
        'graph6_order': 'The vertices array gives the graph6 vertex order.',
        'calculation': f'PAGE through graph-genus {version("graph-genus")}, checked with NetworkX {nx.__version__}.',
        'pendant_edges': 'A pendant edge can be drawn in a small disk at its neighbor in any embedding. We remove pendant edges before PAGE and restore them in reverse order in its rotation system. Each stored rotation is for the full graph.',
        'certificates': 'A genus-zero rotation certifies planarity. A Kuratowski subdivision certifies nonplanarity. For genus-one cases the stored rotation also gives the upper bound.',
    }
    summary = {'total': len(verdicts), 'planar': sum(row['planar'] for row in verdicts),
               'nonplanar': sum(not row['planar'] for row in verdicts), 'by_base': counts}
    target = output_folder / 'enumerations' / 'E30'
    target.mkdir(parents=True, exist_ok=True)
    for filename, data in [('input', inputs), ('candidates', candidates),
                           ('verdicts', verdicts), ('summary', summary)]:
        (target / (filename + '.json')).write_text(json.dumps(data, indent=1) + '\n', encoding='utf-8')
    print(f'E30: {summary["total"]} completions, {summary["planar"]} planar, {summary["nonplanar"]} nonplanar')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_folder', type=Path)
    generate(parser.parse_args().output_folder)
