"""Check all attachment-edge completions of T(z), S(x;y), and Q(x;y,z)."""

import os

from items_enumerations_a import check_graph6
from witnesses import check_embedding, check_kuratowski, edges, graph, load, need


# Read directly from the labeled figure. Each string names one undirected edge.
FIGURE = {
    'T(z)': 'px py pu qx qy qu uz',
    'S(x;y)': 'px py pu qx qy qu xf uf fz',
    'Q(x;y,z)': 'up pv vq qu px xq yu vz',
}
ATTACHMENT_EDGES = ('xy', 'xz', 'yz')


def edge_pairs(names):
    return {tuple(sorted(edge)) for edge in names}


def lemma_predicts_planar(base, added):
    if base == 'T(z)':
        return not {'xz', 'yz'} <= added
    if base == 'S(x;y)':
        return 'yz' not in added
    return added <= {'xy', 'xz'} or added == {'yz'}


def check_pendant_edges(adj, removals, case):
    core = {u: neighbors.copy() for u, neighbors in adj.items()}
    for leaf, neighbor in removals:
        need(core.get(leaf) == {neighbor}, case, 'vertex was not pendant at removal', leaf)
        core[neighbor].remove(leaf)
        del core[leaf]
    need(all(len(neighbors) >= 2 for neighbors in core.values()), case, 'incomplete pruning')


def check(root, item):
    directory = os.path.join(root, 'enumerations', item)
    inputs = load(os.path.join(directory, 'input.json'))
    candidates = load(os.path.join(directory, 'candidates.json'))
    verdicts = load(os.path.join(directory, 'verdicts.json'))
    summary = load(os.path.join(directory, 'summary.json'))
    need(inputs['item'] == item == 'E30', 'wrong item')
    need(inputs['attachment_edges'] == list(ATTACHMENT_EDGES), 'wrong attachment edge order')
    need(set(inputs['base_edges']) == set(FIGURE), 'wrong set of base graphs')
    expected_ids = {base[0] + str(mask) for base in FIGURE for mask in range(8)}
    for records, name in ((candidates, 'candidates'), (verdicts, 'verdicts')):
        need(len(records) == 24 and {r['id'] for r in records} == expected_ids,
             name, 'must list each of the 24 cases exactly once')
    candidate_by_id = {r['id']: r for r in candidates}
    verdict_by_id = {r['id']: r for r in verdicts}

    by_base = {}
    for base, figure in FIGURE.items():
        base_edges = edge_pairs(figure.split())
        recorded_edges = inputs['base_edges'][base]
        need(len(recorded_edges) == len(base_edges) and edge_pairs(recorded_edges) == base_edges,
             base, 'input edges differ from the figure')
        labels = sorted({v for edge in base_edges for v in edge})
        counts = {'planar': 0, 'nonplanar': 0}
        for mask in range(8):
            case = base[0] + str(mask)
            candidate = candidate_by_id[case]
            verdict = verdict_by_id[case]
            added = {e for i, e in enumerate(ATTACHMENT_EDGES) if mask & (1 << i)}
            recorded = candidate['added_edges']
            need(candidate['base'] == base and len(recorded) == len(added)
                 and set(recorded) == added, case, 'wrong base or attachment edges')
            need(candidate['vertices'] == labels, case, 'incorrect vertex labels or order')
            adj = graph(candidate)
            need(set(edges(adj)) == base_edges | edge_pairs(added), case, 'wrong completed graph')
            check_graph6(adj, labels, candidate['graph6'])
            check_pendant_edges(adj, verdict['pendant_edges_removed_for_PAGE'], case)

            planar = lemma_predicts_planar(base, added)
            need(verdict['planar'] is planar, case, 'planarity disagrees with the lemma')
            genus = 0 if planar else 1
            need(type(verdict['genus']) is int and verdict['genus'] == genus, case, 'wrong genus')
            check_embedding(adj, verdict['embedding'], genus=genus)
            if planar:
                need('kuratowski' not in verdict, case, 'unexpected nonplanarity witness')
                counts['planar'] += 1
            else:
                witness = verdict['kuratowski']
                need(witness['type'] == 'Kuratowski subgraph', case, 'wrong lower-bound witness')
                check_kuratowski(adj, witness['edges'])
                counts['nonplanar'] += 1
        by_base[base] = counts

    planar_total = sum(counts['planar'] for counts in by_base.values())
    need(summary['total'] == 24 and summary['planar'] == planar_total
         and summary['nonplanar'] == 24 - planar_total and summary['by_base'] == by_base,
         'summary differs from independently checked cases')
    return f'24 completions: {planar_total} plane embeddings, {24 - planar_total} genus-one embeddings with Kuratowski witnesses'
