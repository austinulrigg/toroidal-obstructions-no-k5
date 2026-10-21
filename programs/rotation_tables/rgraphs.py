"""The graphs of the rotation tables R1-R11, with the vertex names of the text."""
from rotlib import RotationTable


def psi1(names=None):
    """Psi_1 with P1, P2, P3 the edges s1t1, s2t2, s3t3. names renames vertices."""
    n = names or {}
    R = lambda x: n.get(x, x)
    V = ['s1', 't1', 's2', 't2', 's3', 't3', 'p', 'q', 'r', 'z']
    E = [(w, 's%d' % i) for w in 'pq' for i in (1, 2, 3)] + [(w, 't%d' % i) for w in 'rz' for i in (1, 2, 3)] + \
        [('s%d' % i, 't%d' % i) for i in (1, 2, 3)]
    return [R(v) for v in V], [(R(u), R(w)) for u, w in E]


def contract(V, E, keep, gone):
    V2 = [v for v in V if v != gone]
    E2 = set()
    for u, w in E:
        u2 = keep if u == gone else u
        w2 = keep if w == gone else w
        if u2 != w2:
            E2.add(tuple(sorted((u2, w2))))
    return V2, sorted(E2)


def R1_graph():
    V, E = psi1()
    return RotationTable('Psi_1 with P1, P2, P3 the edges s1t1, s2t2, s3t3', V, E,
                         ['s1', 't1', 's2', 't2', 's3', 't3'], sep=' ')


def R2_graph():
    # Lemma 4.9: a = s1, b = t1, c = s2, d = t2; t3 contracted into s3
    V, E = psi1({'s1': 'a', 't1': 'b', 's2': 'c', 't2': 'd'})
    V, E = contract(V, E, 's3', 't3')
    return RotationTable('Psi_1/s3t3 (a = s1, b = t1, c = s2, d = t2; the edge s3t3 contracted into s3)',
                         V, E, list('abcd'))


def R2_psi1_graph():
    V, E = psi1({'s1': 'a', 't1': 'b', 's2': 'c', 't2': 'd'})
    return RotationTable('Psi_1 (a = s1, b = t1, c = s2, d = t2), for comparison', V, E, list('abcd'))


def R3_graph():
    # Proposition 4.26: a = s1, b = t1, c = s2, d = t3, e = s3, f = t2; f contracted into c
    V, E = psi1({'s1': 'a', 't1': 'b', 's2': 'c', 't3': 'd', 's3': 'e', 't2': 'f'})
    V, E = contract(V, E, 'c', 'f')
    return RotationTable('Psi_1/cf (a = s1, b = t1, c = s2, d = t3, e = s3, f = t2; f contracted into c)',
                         V, E, list('abcd'))


def psi2(nc=2, nd=2, minus=()):
    C = ['p', 'q', 'p3'][:nc]
    D = ['r', 'z', 'r3'][:nd]
    V = list('abcd') + C + D
    E = [(t, x) for t in C for x in 'abc'] + [(t, x) for t in D for x in 'abd']
    V = [v for v in V if v not in minus]
    E = [e for e in E if not set(e) & set(minus)]
    return V, E


def R4_graphs():
    out = []
    for nc, nd in ((2, 2), (3, 2), (3, 3)):
        V, E = psi2(nc, nd)
        extra = {(2, 2): 'Psi_2', (3, 2): 'Psi_2 + p3 (N(p3) = {a,b,c})',
                 (3, 3): 'Psi_2 + p3 + r3 (N(p3) = {a,b,c}, N(r3) = {a,b,d})'}[(nc, nd)]
        out.append(((nc, nd), RotationTable('%s, (n_c, n_d) = (%d, %d)' % (extra, nc, nd), V, E, list('abcd'))))
    return out


def R5_graph():
    V, E = psi2(2, 2, minus=('z',))
    return RotationTable('Psi_2 - z', V, E, list('abcd'))


def psi3(extra_edges=(), extra_vertices=(), ne=2):
    D = ['r', 'z', 'r3'][:ne]
    V = list('abcde') + ['p', 'q'] + D
    E = [(t, x) for t in ('p', 'q') for x in 'abc'] + [(t, x) for t in D for x in 'ade']
    for name, nb in extra_vertices:
        V.append(name)
        E += [(name, x) for x in nb]
    E += list(extra_edges)
    return V, E


def R6_graph():
    V, E = psi3()
    return RotationTable('Psi_3', V, E, list('abcde'))


def R7_graph():
    V, E = psi3(extra_vertices=[('w', 'acde')])
    return RotationTable('Psi_3 + w, N(w) = {a,c,d,e}', V, E, list('abcde'))


R8_ROWS = [
    ('bd, or {a,b,d}', [((('b', 'd'),), ()), ((), (('v_abd', 'abd'),))]),
    ('bd,cd, or {a,b,d} and cd', [((('b', 'd'), ('c', 'd')), ()), ((('c', 'd'),), (('v_abd', 'abd'),))]),
    ('{a,c,d} and bd, or {a,b,d},{a,c,d}', [((('b', 'd'),), (('v_acd', 'acd'),)),
                                           ((), (('v_abd', 'abd'), ('v_acd', 'acd')))]),
]


def R8_graphs():
    out = []
    for row, alts in R8_ROWS:
        for ne in (2, 3):
            for k, (es, vs) in enumerate(alts):
                V, E = psi3(extra_edges=es, extra_vertices=vs, ne=ne)
                parts = ['Psi_3'] + (['r3 (N(r3) = {a,d,e})'] if ne == 3 else []) + \
                        [x + y for x, y in es] + ['%s (N(%s) = {%s})' % (n, n, ','.join(nb)) for n, nb in vs]
                name = ' + '.join(parts) + ', n_e = %d' % ne
                out.append({'row': row, 'n_e': ne, 'alternative': k + 1,
                            'table': RotationTable(name, V, E, list('acde'))})
    return out


def R9_graph():
    V, E = psi3(extra_edges=[('c', 'e')])
    return RotationTable('Psi_3 + ce', V, E, list('abcde'))


def R10_graphs():
    out = []
    for nb in ('abc', 'ade', 'acd'):
        V, E = psi3(extra_edges=[('c', 'e')], extra_vertices=[('v', nb)])
        out.append((nb, RotationTable('Psi_3 + ce + v, N(v) = {%s}' % ','.join(nb), V, E, list('abcde'))))
    V, E = psi3(extra_vertices=[('v', 'cde')])
    out.append(('cde', RotationTable('Psi_3 + v, N(v) = {c,d,e}', V, E, list('abcde'))))
    return out


def R11_graphs():
    V, E = psi3(extra_edges=[('c', 'e')], extra_vertices=[('v', 'bcd')])
    a = RotationTable('Psi_3 + ce + v, N(v) = {b,c,d}', V, E, list('abcde'))
    V, E = psi3(extra_edges=[('c', 'e'), ('c', 'd')])
    b = RotationTable('Psi_3 + ce + cd', V, E, list('abcde'))
    return [('v on {b,c,d}', a), ('cd', b)]
