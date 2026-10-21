"""Compare the output of topological_search.py with items TO1 and TO2 of the
data folder:

    python compare.py ../../supplemental_data OUTPUT_FOLDER

The graphs are compared up to isomorphism: the same 381 topological
obstructions (with the same graphs named T1, ..., T74), and the same 7,006
graphs of the search, each with the same verdict. The witnesses may differ
(the search may find other branch sets or embeddings); the checker checks
those of the data folder, and the program's own can be checked by placing its
TO1.json and TO2.json in a copy of the data folder.
"""
import json
import os
import sys

from topological_search import IsoIndex, from_graph6


def load(folder, name):
    for p in (os.path.join(folder, name), os.path.join(folder, 'no_K5_topological_obstructions', name)):
        if os.path.exists(p):
            return json.load(open(p))
    raise SystemExit('cannot find %s in %s' % (name, folder))


def main():
    data, out = sys.argv[1], sys.argv[2]
    ok = True
    a1, b1 = load(data, 'TO1.json'), load(out, 'TO1.json')
    idx = IsoIndex()
    for k, g in enumerate(a1['graphs']):
        idx.add(from_graph6(g['graph6']), k)
    hit = {}
    for g in b1['graphs']:
        k = idx.find(from_graph6(g['graph6']))
        if k is None or k in hit:
            ok = False
            print('TO1: a graph of the output is not in the data folder, or twice:', g['graph6'])
            continue
        hit[k] = g
        if a1['graphs'][k].get('is') != g.get('is'):
            ok = False
            print('TO1: names differ for', g['graph6'])
    if len(hit) != len(a1['graphs']):
        ok = False
        print('TO1: %d graphs of the data folder are not in the output' % (len(a1['graphs']) - len(hit)))
    print('TO1: %d graphs in the data folder, %d in the output, %d matched up to isomorphism' % (
        len(a1['graphs']), len(b1['graphs']), len(hit)))
    a2, b2 = load(data, 'TO2.json'), load(out, 'TO2.json')
    idx = IsoIndex()
    for k, r in enumerate(a2['graphs']):
        idx.add(from_graph6(r['graph6']), k)
    hit, differ = set(), 0
    for r in b2['graphs']:
        k = idx.find(from_graph6(r['graph6']))
        if k is None or k in hit:
            ok = False
            print('TO2: a graph of the output is not in the data folder, or twice:', r['graph6'])
            continue
        hit.add(k)
        if a2['graphs'][k]['verdict'] != r['verdict']:
            differ += 1
            ok = False
            print('TO2: verdicts differ for', r['graph6'], ':', a2['graphs'][k]['verdict'], 'and', r['verdict'])
    if len(hit) != len(a2['graphs']):
        ok = False
        print('TO2: %d graphs of the data folder are not in the output' % (len(a2['graphs']) - len(hit)))
    print('TO2: %d graphs in the data folder, %d in the output, %d matched up to isomorphism, %d with other '
          'verdicts' % (len(a2['graphs']), len(b2['graphs']), len(hit), differ))
    print('counts in the data folder:', a2['counts'])
    print('counts in the output:     ', b2['counts'])
    if a2['counts'] != b2['counts']:
        ok = False
    print('the output agrees with the data folder' if ok else 'THE OUTPUT DIFFERS FROM THE DATA FOLDER')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
