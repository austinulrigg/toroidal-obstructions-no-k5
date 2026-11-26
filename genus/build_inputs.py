"""Rebuild the genus inputs from O1 and E1–E29.

    python build_inputs.py ../supplemental_data

The enumeration graphs are indexed by their canonical graph6 strings from
nauty's labelg. Each index records every candidate that uses the graph.
"""
import gzip
import json
import os
import shutil
import subprocess
import sys
from collections import Counter

from common import HERE, components, from_graph6, n_edges, save_json, to_graph6


def load(p):
    if not os.path.exists(p):
        p += '.gz'
    with (gzip.open(p, 'rt') if p.endswith('.gz') else open(p)) as f:
        return json.load(f)


def adj_from(vertices, edge_list):
    adj = {v: set() for v in vertices}
    for u, v in edge_list:
        assert u != v and u in adj and v in adj
        adj[u].add(v)
        adj[v].add(u)
    return adj


def apply_op(adj, op, e):
    u, v = e
    assert v in adj[u], 'the edge of the witness is not in the graph'
    H = {x: set(s) for x, s in adj.items()}
    if op == 'delete':
        H[u].discard(v)
        H[v].discard(u)
    else:
        for w in H.pop(v):
            H[w].discard(v)
            if w != u:
                H[w].add(u)
                H[u].add(w)
    return H


def covering_graph(base, atoms, counts):
    """The base graph with counts[A] copies of each atom A; the interior vertex
    '@x' of copy k of atom A is named A#k@x (as in the records of E4 and E12)."""
    adj = {}

    def add(x, y):
        adj.setdefault(x, set()).add(y)
        adj.setdefault(y, set()).add(x)
    for x, y in base:
        add(x, y)
    for aid, c in counts.items():
        for k in range(c):
            f = lambda v: '%s#%d@%s' % (aid, k, v[1:]) if v.startswith('@') else v
            for x, y in atoms[aid]['edges']:
                add(f(x), f(y))
    return adj


def bare(w):
    """A proper-minor witness that names only the edge, or whose recorded
    reason is the genus computation (three in E2)."""
    if not (isinstance(w, dict) and w.get('type') == 'proper minor'):
        return False
    why = w.get('why the minor is nontoroidal')
    if why is None:
        return True
    return isinstance(why, dict) and 'genus computations' in why.get('reason', '')


def canonical(g6s, labelg, chunk=200000):
    out = []
    for i in range(0, len(g6s), chunk):
        part = g6s[i:i + chunk]
        r = subprocess.run([labelg, '-q'], input='\n'.join(part) + '\n', capture_output=True, text=True, check=True)
        res = r.stdout.split()
        assert len(res) == len(part)
        out += res
    return out


def main(data):
    labelg = os.environ.get('LABELG') or shutil.which('labelg')
    if not labelg:
        raise SystemExit('labelg (from nauty) is needed; give its path in the environment variable LABELG')
    E = os.path.join(data, 'enumerations')
    graphs, uses = [], []
    count = Counter()

    def add(adj, use, kind):
        assert len(components(adj)) == 1, use
        graphs.append(to_graph6(adj, sorted(adj)))
        uses.append(use)
        count[(use['item'], kind)] += 1

    def minor_use(item, cid, w, extra=None):
        u = {'item': item, 'candidate': cid}
        if extra:
            u.update(extra)
        u['graph'] = '%s of the edge %s%s' % ('deletion' if w['operation'] == 'delete' else 'contraction',
                                              *w['edge'])
        return u

    for number in range(1, 30):
        item = f'E{number}'
        d = os.path.join(E, item)
        cands = {c['id']: c for c in load(os.path.join(d, 'candidates.json'))['candidates']}
        verds = load(os.path.join(d, 'verdicts.json'))['verdicts']
        inp = load(os.path.join(d, 'input.json'))
        rest = None
        if item == 'E2':
            r = inp['the rest of G']
            rest = adj_from(r['vertices'], r['edges'])
        for v in verds:
            if v['verdict'] == 'proper nontoroidal minor' and bare(v['witness']):
                c = cands[v['id']]
                if item == 'E2':
                    b = c['bridge']
                    G = {x: set(s) for x, s in rest.items()}
                    for x, s in adj_from(b['vertices'], b['edges']).items():
                        G.setdefault(x, set()).update(s)
                else:
                    g = c['graph']
                    G = adj_from(g['vertices'], g['edges'])
                w = v['witness']
                add(apply_op(G, w['operation'], w['edge']), minor_use(item, v['id'], w), 'one-edge minor')
            for k, wg in enumerate(v.get('whole graphs', [])):
                if wg['verdict'] == 'proper nontoroidal minor' and bare(wg['witness']):
                    g = wg['graph']
                    w = wg['witness']
                    add(apply_op(adj_from(g['vertices'], g['edges']), w['operation'], w['edge']),
                        minor_use(item, v['id'], w, {'whole graph': k}), 'one-edge minor')
        if item in ('E4', 'E12'):
            atoms = {a['id']: a for a in inp['atoms']}
            base = inp['the base graph']['edges']
            for c in cands.values():
                if c.get('kind') == 'nontoroidal, minimal':
                    G = covering_graph(base, atoms, c['atoms'])
                    if 'graph' in c:
                        assert G == adj_from(c['graph']['vertices'], c['graph']['edges'])
                    add(G, {'item': item, 'candidate': c['id'], 'graph': 'recorded nontoroidal graph of the covering'},
                        'recorded nontoroidal graph')
        if item == 'E19':
            for v in verds:
                if v.get('lattice status') == 'nontoroidal' and v['verdict'] != 'proper nontoroidal minor' \
                        and not v['verdict'].startswith('isomorphic to'):
                    g = cands[v['id']]['graph']
                    add(adj_from(g['vertices'], g['edges']),
                        {'item': item, 'candidate': v['id'], 'graph': 'combination with nontoroidal status'},
                        'nontoroidal combination')
            for a in inp['lattice']['auxiliary']:
                if a['lattice status'] == 'nontoroidal':
                    g = a['graph']
                    add(adj_from(g['vertices'], g['edges']),
                        {'item': item, 'candidate': a['id'], 'graph': 'auxiliary combination with nontoroidal status'},
                        'nontoroidal combination')
    print('graphs with repetitions:', len(graphs), flush=True)
    classes, order = {}, []
    for c, u in zip(canonical(graphs, labelg), uses):
        if c not in classes:
            classes[c] = []
            order.append(c)
        classes[c].append(u)
    order.sort(key=lambda c: (len(from_graph6(c)), n_edges(from_graph6(c)), c))
    with open(HERE / 'enumeration_graphs.txt', 'w') as f:
        f.write('# index graph6 vertices edges uses\n')
        for i, c in enumerate(order):
            g = from_graph6(c)
            f.write('%d %s %d %d %d\n' % (i + 1, c, len(g), n_edges(g), len(classes[c])))
    save_json(HERE / 'enumeration_graph_uses.json.gz', {'graphs': [{'index': i + 1, 'graph6': c, 'uses': classes[c]}
                                                                   for i, c in enumerate(order)]})
    o1 = load(os.path.join(data, 'no_K5_obstructions', 'O1.json'))
    with open(HERE / 'obstructions.txt', 'w') as f:
        f.write('# name graph6 vertices edges\n')
        for r in o1['graphs']:
            g = from_graph6(r['graph6'])
            assert len(g) == r['n'] and n_edges(g) == r['m']
            f.write('%s %s %d %d\n' % (r['name'], r['graph6'], r['n'], r['m']))
    summary = {'graphs with repetitions': len(graphs), 'graphs up to isomorphism': len(order),
               'by item': {'%s, %s' % k: n for k, n in sorted(count.items(), key=lambda x: (int(x[0][0][1:]), x[0][1]))},
               'up to isomorphism, by number of vertices': dict(sorted(Counter(len(from_graph6(c)) for c in order).items()))}
    save_json(HERE / 'inputs_summary.json', summary, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
