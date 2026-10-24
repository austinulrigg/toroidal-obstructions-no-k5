"""Independently retrace R12 and check every recorded claim.

Usage: python items_rotation_R12.py SUPPLEMENTAL_DATA_FOLDER
Uses only the standard library and does not import the generator.
The check(root, item) entry point can be registered in check_supplement.py.
"""
import itertools
import json
import math
import sys
from collections import Counter
from pathlib import Path


def least_shift(sequence):
    if not sequence:
        return ()
    n = len(sequence)
    return min(tuple(sequence[(j+k) % n] for k in range(n)) for j in range(n))


def text(sequence):
    return ' '.join(sequence) if sequence else '-'


def check(root, item='R12'):
    if item != 'R12':
        raise AssertionError('This module checks only R12')
    table = json.loads((Path(root)/'rotations'/'R12.json').read_text())
    assert table['item'] == 'R12'
    (graph,) = table['graphs']
    vertices = graph['vertices']
    roots = {'s1', 's2', 't1', 't2'}
    assert set(vertices) == {'p','q','r','z','s1','s2','s3','t1','t2','t3'}
    assert len(vertices) == len(set(vertices)) == 10
    expected = {frozenset((twin, 's'+str(i))) for twin in ('p','q') for i in (1,2,3)}
    expected |= {frozenset((twin, 't'+str(i))) for twin in ('r','z') for i in (1,2,3)}
    expected.add(frozenset(('s3','t3')))
    assert {frozenset(e) for e in graph['edges']} == expected
    assert len(graph['edges']) == 13
    assert set(graph['attachments']) == roots and len(graph['attachments']) == 4
    neighbours = {v: sorted(set().union(*(e-{v} for e in expected if v in e))) for v in vertices}
    orders = graph['cyclic orders']
    for v in vertices:
        nb = neighbours[v]
        want = [(nb[0],)+p for p in itertools.permutations(nb[1:])]
        assert [tuple(c) for c in orders[v]] == want
    n = math.prod(math.factorial(len(neighbours[v])-1) for v in vertices)
    assert graph['number of rotation systems'] == n == len(graph['rows']) == 64
    darts = [(v,w) for v in vertices for w in neighbours[v]]
    index = {d:i for i,d in enumerate(darts)}
    genera, counts = Counter(), {}
    full = Counter()
    pair_counts, pair_types, unique = Counter(), Counter(), 0
    signatures = Counter()
    for number in range(n):
        quotient = number
        rotation = {}
        for v in reversed(vertices):
            quotient, digit = divmod(quotient, len(orders[v]))
            rotation[v] = orders[v][digit]
        assert quotient == 0
        permutation = []
        for u,v in darts:
            cyclic = rotation[v]
            following = cyclic[(cyclic.index(u)+1) % len(cyclic)]
            permutation.append(index[(v,following)])
        assert sorted(permutation) == list(range(26))
        visited = [False]*26
        walks = []
        for start in range(26):
            if visited[start]:
                continue
            current = start
            cycle = []
            while not visited[current]:
                visited[current] = True
                cycle.append(darts[current][0])
                current = permutation[current]
            assert current == start
            walks.append(cycle)
        assert all(visited) and sum(map(len,walks)) == 26
        numerator = 5-len(walks)
        assert numerator >= 0 and numerator % 2 == 0
        genus = numerator//2
        sequences = [[x for x in w if x in roots] for w in walks]
        row = [genus, sorted(text(least_shift(s)) for s in sequences)]
        assert row == graph['distinct rows'][graph['rows'][number]], ('row',number)
        genera[genus] += 1
        counts.setdefault(genus,Counter())[len(walks)] += 1
        full[genus] += sum(roots <= set(s) for s in sequences)
        if genus == 0:
            covering = [(i,j) for i in range(len(sequences)) for j in range(i+1,len(sequences))
                        if roots <= set(sequences[i]) | set(sequences[j])]
            pair_counts[len(covering)] += 1
            for i,j in covering:
                s,t = sequences[i],sequences[j]
                pair_types[' | '.join(sorted([text(least_shift(s)),text(least_shift(t))]))] += 1
                if all(s.count(v)+t.count(v)==1 for v in roots):
                    unique += 1
            signatures[' / '.join(row[1])] += 1
    by_genus = {str(g):genera[g] for g in sorted(genera)}
    assert graph['rotation systems by genus'] == by_genus == {'0':16,'1':32,'2':16}
    assert full[0] == full[1] == 0
    assert pair_counts == {1:16} and pair_types == {'s1 s2 | t1 t2':16} and unique == 16
    expected_signatures = set()
    for i,j in itertools.product((1,2),repeat=2):
        sequences = [('s'+str(i),'t'+str(j)),('s1','s2'),('s'+str(3-i),),('t1','t2'),('t'+str(3-j),)]
        expected_signatures.add(' / '.join(sorted(text(least_shift(s)) for s in sequences)))
    assert set(signatures) == expected_signatures and all(v==4 for v in signatures.values())
    computed = [
        {'vertices':10,'edges':13,'rotation systems':64,'by genus':by_genus,
         'facial walk counts by genus':{str(g):{str(k):v for k,v in sorted(c.items())} for g,c in sorted(counts.items())}},
        {'genus-one rotations':32,'walks containing all four attachments':0},
        {'genus-zero rotations':16,'walks containing all four attachments':0,
         'covering pairs per rotation':{'1':16},'covering pair attachment sequences':dict(pair_types),
         'covering pairs with each attachment appearing once':16},
        {'boundary signatures':dict(sorted(signatures.items()))},
    ]
    assert len(table['claims']) == len(computed)
    for claim,want in zip(table['claims'],computed):
        assert claim['holds'] is True and claim['computed'] == want
    return '64 rotation systems retraced; every stored row and all four claims verified'


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python items_rotation_R12.py SUPPLEMENTAL_DATA_FOLDER')
    print(check(sys.argv[1]))
