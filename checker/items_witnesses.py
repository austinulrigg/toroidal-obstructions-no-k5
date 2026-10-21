"""The witnesses W1-W21: single certificates cited in proofs.

Each item gives one graph (or a list of graphs) with the vertex names of the
text, and a witness of one of the seven kinds. The checker rebuilds every
graph and checks its witness against it: an embedding by tracing its faces,
K5 branch sets, or an isomorphism onto the stated graph of O1. For W1 it also
checks that the graph is the Wagner graph V8 (the claim of Section 2.2) and
that the embedding has genus exactly one.
"""
import os

from witnesses import catalogue, check_iso, check_witness, graph, load, need, v8


def check(root, item):
    d = load(os.path.join(root, 'witnesses', item + '.json'))
    need(d.get('item') == item, item, ': wrong item label inside the file')
    T = catalogue(root)
    entries = d['graphs'] if 'graphs' in d else [dict(d['graph'], witness=d['witness'])]
    kinds = []
    for e in entries:
        G = graph({'vertices': e['vertices'], 'edges': e['edges']})
        w = e['witness']
        check_witness(G, w, T)
        kinds.append(w['type'] + (' onto ' + w['onto'] if w['type'] == 'isomorphism' else ''))
        if item == 'W1':
            check_iso(G, v8(), {v: int(v) for v in G})
            need(w['type'] == 'embedding' and w['genus'] == 1, 'W1 must give an embedding of genus one')
    return '%d graph%s: %s' % (len(entries), '' if len(entries) == 1 else 's', '; '.join(kinds))
