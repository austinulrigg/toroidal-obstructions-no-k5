"""The rotation tables R1-R11.

For every graph of every table: the cyclic orders listed at each vertex are
exactly the cyclic orders of its neighbours (each starting with its first
neighbour in alphabetical order); every rotation system is decoded from its
number in the mixed radix and its facial walks are traced again; its genus,
attachment sequences and, where given, regions along a path must equal the
stored row; the counts by genus must agree; and the data quoted in each claim
of the table is recomputed from the retraced faces and must equal the stored
data. R1 also checks the automorphisms of Psi_1 and the orbits of its torus
rotation systems by brute force.

R8 has 21,220,992 rotation systems and takes the longest (several minutes).
"""
import gzip
import itertools
import json
import os
import sys
from collections import Counter, defaultdict

ROOT = None          # the rotations folder, set by check()


def _read(p):
    if not os.path.exists(p):
        p += '.gz'
    op = gzip.open if p.endswith('.gz') else open
    with op(p, 'rt') as f:
        return json.load(f)


def load(item):
    """A table is R<i>.json(.gz), or a folder R<i> holding R<i>.json(.gz) and the
    rows of its largest graphs in parts."""
    folder = os.path.join(ROOT, item)
    if not os.path.isdir(folder):
        return _read(os.path.join(ROOT, item + '.json'))
    d = _read(os.path.join(folder, item + '.json'))
    for g in d['graphs']:
        if isinstance(g['rows'], dict):
            rows = []
            for k, name in enumerate(g['rows']['stored in']):
                part = _read(os.path.join(folder, name[:-3] if name.endswith('.gz') else name))
                assert part['item'] == item and part['graph'] == g['name'] and part['part'] == k + 1
                assert part['first rotation system'] == len(rows) == k * g['rows']['rows per part']
                rows += part['rows']
            g['rows'] = rows
    return d


def min_rot(seq):
    seq = list(seq)
    best = None
    for i in range(len(seq)):
        c = seq[i:] + seq[:i]
        if best is None or c < best:
            best = c
    return best or []


def canon2(seq):
    a = min_rot(seq)
    b = min_rot(list(reversed(list(seq))))
    return min(a, b)


def word(seq, sep):
    return sep.join(seq) if seq else '-'


class G:
    def __init__(self, rec, sep=''):
        self.rec = rec
        self.sep = sep
        self.V = rec['vertices']
        self.adj = {v: set() for v in self.V}
        for u, w in rec['edges']:
            self.adj[u].add(w)
            self.adj[w].add(u)
        self.att = set(rec['attachments'])
        self.E = sum(len(s) for s in self.adj.values()) // 2
        co = rec['cyclic orders']
        for v in self.V:
            nb = sorted(self.adj[v])
            k = len(nb)
            want = 1
            for x in range(2, k):
                want *= x
            lst = co[v]
            assert len(lst) == max(1, want), ('number of cyclic orders', v)
            seen = set()
            for c in lst:
                assert sorted(c) == nb and c[0] == nb[0], ('cyclic order', v, c)
                seen.add(tuple(c))
            assert len(seen) == len(lst)
        self.co = co
        self.radix = [len(co[v]) for v in self.V]
        n = 1
        for r in self.radix:
            n *= r
        assert n == rec['number of rotation systems'] == len(rec['rows'])
        self.n = n

    def rotation(self, i):
        rot = {}
        for v, r in zip(reversed(self.V), reversed(self.radix)):
            i, d = divmod(i, r)
            rot[v] = self.co[v][d]
        assert i == 0
        return rot

    def faces(self, rot):
        succ = {}
        for v, nb in rot.items():
            for j, u in enumerate(nb):
                succ[(v, u)] = nb[(j + 1) % len(nb)]
        seen = set()
        out = []
        for v in self.V:
            for u in sorted(self.adj[v]):
                d0 = (v, u)
                if d0 in seen:
                    continue
                f = []
                d = d0
                while d not in seen:
                    seen.add(d)
                    f.append(d)
                    d = (d[1], succ[(d[1], d[0])])
                out.append(f)
        return out

    def seq(self, f):
        return [d[0] for d in f if d[0] in self.att]

    def genus(self, faces):
        x = 2 - len(self.V) + self.E - len(faces)
        assert x >= 0 and x % 2 == 0
        return x // 2

    def region(self, faces, s, t):
        pos = {}
        for fi, f in enumerate(faces):
            for j, d in enumerate(f):
                pos[d] = (fi, j)
        (a, i), (b, j) = pos[(s, t)], pos[(t, s)]
        if a != b:
            F1, F2 = faces[a], faces[b]
            w = [F1[(i + 1 + k) % len(F1)] for k in range(len(F1) - 1)] + \
                [F2[(j + 1 + k) % len(F2)] for k in range(len(F2) - 1)]
            return word(canon2(self.seq(w)), self.sep)
        F = faces[a]
        n = len(F)
        x = [F[(i + 1 + k) % n] for k in range((j - i - 1) % n)]
        y = [F[(j + 1 + k) % n] for k in range((i - j - 1) % n)]
        p, q = sorted([word(canon2(self.seq(x)), self.sep), word(canon2(self.seq(y)), self.sep)])
        return p + ' | ' + q + ' (annulus)'


def run_graph(rec, sep='', regions=None, hook=None):
    """Retrace every rotation system and compare with the stored rows."""
    g = G(rec, sep)
    rows, distinct = rec['rows'], rec['distinct rows']
    by_genus = Counter()
    for i in range(g.n):
        rot = g.rotation(i)
        faces = g.faces(rot)
        gen = g.genus(faces)
        seqs = [g.seq(f) for f in faces]
        row = distinct[rows[i]]
        assert row[0] == gen, ('genus', rec['name'], i)
        assert row[1] == sorted(word(min_rot(s), sep) for s in seqs), ('faces', rec['name'], i)
        if regions is not None:
            got = regions(g, faces, gen)
            assert row[2] == got, ('regions', rec['name'], i, row[2], got)
        by_genus[gen] += 1
        if hook is not None:
            hook(g, i, faces, gen, seqs)
    assert {str(k): v for k, v in sorted(by_genus.items())} == rec['rotation systems by genus']
    return g, dict(by_genus)


def cdict(c):
    return {k: c[k] for k in sorted(c)}


def restrict(s, keep):
    return [x for x in s if x in keep]


def H_edges(s):
    out = set()
    for i in range(len(s)):
        u, v = s[i], s[(i + 1) % len(s)]
        if u != v:
            out.add(''.join(sorted((u, v))))
    return sorted(out)


def H_kind(s):
    if any(s[i] == s[(i + 1) % len(s)] for i in range(len(s))):
        return None
    es = H_edges(s)
    vs = sorted(set(s))
    deg = sorted(sum(v in e for e in es) for v in vs)
    n, m = len(vs), len(es)
    if m == n and deg == [2] * n and len(s) == n:
        return 'cycle'
    if (n, m, deg, len(s)) == (4, 4, [1, 2, 2, 3], 5):
        return 'triangle with a pendant edge'
    if (n, m, deg, len(s)) == (4, 3, [1, 1, 1, 3], 6):
        return 'three-edge star'
    if m == n - 1 and max(deg) <= 2 and len(s) == 2 * (n - 1):
        return 'path'
    return None


def tests(seqs):
    return {s: {'H': H_kind(s), 'edges of H': H_edges(s)} for s in seqs}


def claims_of(d):
    return [c['computed'] for c in d['claims']]


# ------------------------------------------------------------------ items

def check_R1():
    d = load('R1')
    (rec,) = d['graphs']
    paths = [('s%d' % m, 't%d' % m) for m in (1, 2, 3)]

    def regions(g, faces, gen):
        if gen != 1:
            return None
        return {'P%d' % m: g.region(faces, *paths[m - 1]) for m in (1, 2, 3)}
    torus = []

    def hook(g, i, faces, gen, seqs):
        if gen == 1:
            torus.append((i, faces, seqs))
    g, by = run_graph(rec, ' ', regions, hook)
    c = claims_of(d)
    assert c[0]['by genus'] == {str(k): v for k, v in sorted(by.items())} and \
        c[0]['faces of the genus-one rotations'] == {'5': len(torus)} and all(len(f) == 5 for _, f, _ in torus)
    # automorphisms by brute force over vertex permutations preserving degrees
    V = g.V
    autos = []
    for perm in itertools.permutations(V):
        m = dict(zip(V, perm))
        if all(m[w] in g.adj[m[v]] for v in V for w in g.adj[v]):
            autos.append(m)
    assert len(autos) == c[1]['automorphisms of Psi_1']
    tor = {i for i, _, _ in torus}
    index = {}
    for i in tor:
        index[json.dumps(norm(g.rotation(i)), sort_keys=True)] = i

    def orbit_count(mirror):
        parent = {i: i for i in tor}

        def f(x):
            while parent[x] != x:
                x = parent[x]
            return x
        for i in tor:
            rot = g.rotation(i)
            for m in autos:
                r2 = {m[v]: [m[u] for u in rot[v]] for v in rot}
                for r3 in ([r2, {v: list(reversed(c_)) for v, c_ in r2.items()}] if mirror else [r2]):
                    j = index[json.dumps(norm(r3), sort_keys=True)]
                    a, b = f(i), f(j)
                    if a != b:
                        parent[a] = b
        return len({f(i) for i in tor})
    assert orbit_count(False) == c[1]['orbits under the automorphisms'] == 2
    assert orbit_count(True) == c[1]['orbits under the automorphisms and reversal of all cyclic orders (mirror images)'] == 2
    # Lemma 4.2 by an independent reading of the five faces
    per = Counter()
    for i, faces, seqs in torus:
        ws = sorted(word(canon2(s), ' ') for s in seqs)
        found = []
        for k in (1, 2, 3):
            ii, jj = sorted({1, 2, 3} - {k})
            for x in (ii, jj):
                for y in (ii, jj):
                    xp, yp = ii + jj - x, ii + jj - y
                    S = lambda m: 's%d' % m
                    T = lambda m: 't%d' % m
                    want = sorted([word(canon2([S(ii), T(ii), T(jj), S(jj)]), ' ')] * 2 +
                                  [word(canon2([S(k), T(k), T(y), T(k), S(k), S(x)]), ' '),
                                   word(canon2([S(k), S(xp)]), ' '), word(canon2([T(k), T(yp)]), ' ')])
                    if ws == want:
                        found.append('%d%d%d' % (k, x, y))
        assert len(found) == 1
        per[found[0]] += 1
        info = c[2]['faces and regions by triple (k,x,y)'][found[0]]
        assert info['faces'] == ws
        assert info['regions'] == rec['distinct rows'][rec['rows'][i]][2]
    assert dict(per) == c[2]['torus rotation systems per triple'] and len(per) == 12
    return 'R1 ok'


def norm(rot):
    out = {}
    for v, c in rot.items():
        i = c.index(min(c))
        out[v] = c[i:] + c[:i]
    return out


def faces_and_regions(item, edge, keep3=True):
    d = load(item)
    (rec,) = d['graphs']
    st = {'faces': Counter(), 'regions': Counter(), 'others': Counter()}

    def regions(g, faces, gen):
        return g.region(faces, *edge) if gen == 1 else None

    def hook(g, i, faces, gen, seqs):
        if gen != 1:
            return
        for s in seqs:
            if len(set(s)) >= 3:
                st['faces'][word(canon2(s), '')] += 1
        r = g.region(faces, *edge)
        st['regions'][r] += 1
        letters = set(r.replace(' | ', '').replace(' (annulus)', ''))
        if {'c', 'd'} <= letters:
            others = sorted(word(canon2(g.seq(f)), '') for f in faces
                            if edge not in f and (edge[1], edge[0]) not in f)
            st['others'][(r, ', '.join(others))] += 1
    g, by = run_graph(rec, '', regions, hook)
    return d, by, st


def check_R2():
    d, by, st = faces_and_regions('R2', ('a', 'b'))
    c = claims_of(d)
    assert c[0]['by genus'] == {str(k): v for k, v in sorted(by.items())}
    assert c[1]['Psi_1/s3t3'] == {'faces': cdict(st['faces']), 'regions along ab': cdict(st['regions'])}
    # Psi_1 from R1, renamed a = s1, b = t1, c = s2, d = t2
    r1 = load('R1')['graphs'][0]
    ren = {'s1': 'a', 't1': 'b', 's2': 'c', 't2': 'd'}
    fo, ro = Counter(), Counter()
    for k in r1['rows']:
        gen, faces, regs = r1['distinct rows'][k]
        if gen != 1:
            continue
        for f in faces:
            s = [ren[x] for x in f.split() if x in ren]
            if len(set(s)) >= 3:
                fo[word(canon2(s), '')] += 1
        r = regs['P1']
        parts = r.replace(' (annulus)', '').split(' | ')
        parts = [word(canon2([ren[x] for x in p.split() if x in ren]), '') for p in parts]
        ro[(' | '.join(sorted(parts)) + ' (annulus)') if len(parts) == 2 else parts[0]] += 1
    assert c[1]['Psi_1'] == {'faces': cdict(fo), 'regions along ab': cdict(ro)}
    new_f = sorted(set(st['faces']) - set(fo))
    new_r = sorted(set(st['regions']) - set(ro))
    assert c[2] == {'faces of Psi_1/s3t3 not occurring for Psi_1': new_f,
                    'regions along ab of Psi_1/s3t3 not occurring for Psi_1': new_r}
    for p in new_f + new_r:
        assert 'annulus' not in p and (p == 'abcd' or (len(p) == 3 and len(set(p)) == 3))
    assert c[3] == {'facial walks of Psi_1 with sequence abdc': fo['abdc']}
    assert fo['abdc'] > 0 and ro['acdbdc'] > 0
    return 'R2 ok'


def check_R3():
    d, by, st = faces_and_regions('R3', ('a', 'b'))
    c = claims_of(d)
    assert c[0]['by genus'] == {str(k): v for k, v in sorted(by.items())}
    got = [{'region along ab': r, 'other faces': o, 'torus rotation systems': n} for (r, o), n in sorted(st['others'].items())]
    assert c[1] == got and {(x['region along ab'], x['other faces']) for x in got} == {
        ('ac | bd (annulus)', 'a, bc, cd, cd'), ('acbdc', 'a, bcd, cd'), ('acbcdc', 'a, bd, cd')}
    return 'R3 ok'


def check_R4():
    d = load('R4')
    c = d['claims']
    out = []
    for rec in d['graphs']:
        st = {'disk': Counter(), 'planar': Counter(), 'single': Counter(), 'pairs': Counter(), 'cd': 0, 'alt3': 0,
              'alt': 0}

        def hook(g, i, faces, gen, seqs):
            if gen == 1:
                for s in seqs:
                    if set('abcd') <= set(s):
                        st['disk'][word(canon2(s), '')] += 1
                if any({'c', 'd'} <= set(s) for s in seqs):
                    st['cd'] += 1
                    for s in seqs:
                        if alternation(s):
                            st['alt'] += 1
                            if len(set(s)) >= 3:
                                st['alt3'] += 1
            if gen == 0:
                ws = [word(canon2(s), '') for s in seqs]
                st['planar'][' '.join(sorted(ws))] += 1
                for w in ws:
                    if set('abcd') <= set(w):
                        st['single'][w] += 1
                ps = set()
                for x, y in itertools.combinations(ws, 2):
                    if set('abcd') <= set(x) | set(y):
                        ps.add(' | '.join(sorted([x, y])))
                for p in ps:
                    st['pairs'][p] += 1
        g, by = run_graph(rec, '', None, hook)
        out.append((rec['name'], by, st))
    # compare with the claims
    by_name = {}
    for cl in c:
        by_name.setdefault(cl['statement'][:24], []).append(cl)
    (n22, b22, s22), (n32, b32, s32), (n33, b33, s33) = out
    seqs = [cl for cl in c if 'genus-one attachment sequences meeting all four' in cl['statement']]
    for cl, (_, by, st) in zip(seqs, out):
        assert cl['computed']['number of facial walks with each sequence over all genus-one rotations'] == cdict(st['disk'])
        assert cl['computed']['H'] == tests(sorted(st['disk']))
        assert all(t['H'] in ('cycle', 'triangle with a pendant edge', 'three-edge star')
                   for t in cl['computed']['H'].values())
    assert sorted(s33['disk']) == ['acbd'] and sorted(s32['disk']) == ['abcad', 'abdbc', 'acbd']
    pl = [cl for cl in c if cl['statement'].startswith('Only (n_c, n_d) = (2, 2)')][0]['computed']
    assert pl == {'planar rotations by their face sequences': cdict(s22['planar']),
                  'planar walks containing a, b, c, d': cdict(s22['single']),
                  'pairs of faces of a planar rotation together containing a, b, c, d': cdict(s22['pairs'])}
    assert set(s22['planar']) == {'ab ab ac ad bc bd'} and not s22['single'] and \
        set(s22['pairs']) == {'ac | bd', 'ad | bc'}
    assert b32.get(0, 0) == 0 and b33.get(0, 0) == 0
    l516 = [cl for cl in c if cl['where'] == 'Lemma 5.29'][0]['computed']
    assert l516 == {'rotations with a walk containing c and d': s22['cd'],
                    'alternation in a face with three labels': s22['alt3'], 'alternation in any face': s22['alt']}
    assert s22['cd'] == 80 and s22['alt3'] == 0
    return 'R4 ok'


def alternation(s):
    for u, v in itertools.combinations(sorted(set(s)), 2):
        r = [x for x in s if x in (u, v)]
        cc = []
        for x in r:
            if not cc or cc[-1] != x:
                cc.append(x)
        if len(cc) > 1 and cc[0] == cc[-1]:
            cc.pop()
        if len(cc) >= 4:
            return True
    return False


def merge_repeats(s):
    cc = []
    for x in s:
        if not cc or cc[-1] != x:
            cc.append(x)
    while len(cc) > 1 and cc[0] == cc[-1]:
        cc.pop()
    return cc


def check_R5():
    d = load('R5')
    (rec,) = d['graphs']
    st = {'planar': Counter(), 'd': Counter(), 'raw': Counter()}

    def hook(g, i, faces, gen, seqs):
        red = [merge_repeats(restrict(s, 'acd')) for s in seqs]
        if gen == 0:
            st['planar'][' '.join(sorted(word(canon2(r), '') for r in red))] += 1
        if gen == 1:
            for s, r in zip(seqs, red):
                if 'd' in r:
                    st['d'][word(canon2(r), '')] += 1
                    st['raw'][word(canon2(s), '')] += 1
    g, by = run_graph(rec, '', None, hook)
    c = claims_of(d)[0]
    assert c['genus zero: the four reduced sequences'] == cdict(st['planar']) == {'a ac ad c': by[0]}
    assert c['genus one: the reduced sequence containing d'] == cdict(st['d'])
    assert set(st['d']) <= {'ad', 'acd', 'acad', 'acacd'}
    assert c['genus one: that walk before the reduction (a, b, c, d retained)'] == cdict(st['raw'])
    return 'R5 ok'


def check_R6():
    d = load('R6')
    (rec,) = d['graphs']
    st = {'big': Counter(), 'oth': Counter(), 'planar': Counter(), 'single': 0, 'pairs': Counter(), 'f': True}

    def hook(g, i, faces, gen, seqs):
        st['f'] &= 2 * gen == 5 - len(faces)
        ws = [word(canon2(s), '') for s in seqs]
        if gen == 1:
            big = [w for w in ws if set('bcde') <= set(w)]
            if big:
                st['big'][big[0]] += 1
                rest = list(ws)
                rest.remove(big[0])
                st['oth'][' '.join(sorted(rest))] += 1
        if gen == 0:
            st['planar'][' '.join(sorted(ws))] += 1
            st['single'] += sum(set('bcde') <= set(w) for w in ws)
            for x, y in itertools.combinations(ws, 2):
                if set('bcde') <= set(x) | set(y):
                    st['pairs'][' | '.join(sorted([x, y]))] += 1
    g, by = run_graph(rec, '', None, hook)
    c = claims_of(d)
    assert st['f'] and c[1] == {'that walk': cdict(st['big']), 'the other two walks': cdict(st['oth'])}
    assert set(st['big']) == {'abadacae'} and set(st['oth']) == {'bc de'}
    assert c[2] == {'planar rotations by their face sequences': cdict(st['planar']), 'walks meeting b, c, d, e':
                    st['single'], 'pairs meeting b, c, d, e': cdict(st['pairs'])}
    for p in st['planar']:
        ok = False
        for sw in ({}, {'b': 'c', 'c': 'b'}, {'d': 'e', 'e': 'd'}, {'b': 'c', 'c': 'b', 'd': 'e', 'e': 'd'}):
            q = ' '.join(sorted(word(canon2([sw.get(x, x) for x in w]), '') for w in p.split()))
            ok |= q == 'abad ac ae bc de'
        assert ok
    assert st['single'] == 0 and set(st['pairs']) == {'bc | de'}
    return 'R6 ok'


def abd_check(item, want):
    d = load(item)
    comp = d['claims'][0]['computed']
    for rec in d['graphs']:
        st = {'abd': Counter(), 'raw': Counter()}

        def hook(g, i, faces, gen, seqs):
            if gen == 1:
                for s in seqs:
                    r = restrict(s, 'abd')
                    if set('abd') <= set(r):
                        st['abd'][word(canon2(r), '')] += 1
                        st['raw'][word(canon2(s), '')] += 1
        g, by = run_graph(rec, '', None, hook)
        assert by.get(0, 0) == 0
        mine = comp if item == 'R7' else comp[rec['name']]
        assert mine['walks containing a, b, d, retaining a, b, d' if item == 'R7' else
                    'genus-one walks containing a, b, d, retaining a, b, d'] == cdict(st['abd'])
        assert set(st['abd']) == want(rec['name']), (rec['name'], st['abd'])
    return item + ' ok'


def check_R8():
    d = load('R8')
    for rec, ex in zip(d['graphs'], d['exteriors']):
        assert rec['name'] == ex['exterior']
        st = {'disk': Counter(), 'single': Counter(), 'pairs': Counter()}

        def hook(g, i, faces, gen, seqs):
            rs = [restrict(s, 'acde') for s in seqs]
            if gen == 1:
                for r in rs:
                    if set('acde') <= set(r):
                        st['disk'][word(canon2(r), '')] += 1
            if gen == 0:
                ws = [word(canon2(r), '') for r in rs]
                for w in ws:
                    if set('acde') <= set(w):
                        st['single'][w] += 1
                ps = set()
                for x, y in itertools.combinations(ws, 2):
                    if set('acde') <= set(x) | set(y):
                        ps.add(' | '.join(sorted([x, y])))
                for p in ps:
                    st['pairs'][p] += 1
        g, by = run_graph(rec, '', None, hook)
        assert ex['disk sequences (torus rotations in which they occur)'] == cdict(st['disk'])
        assert ex['H for each disk sequence'] == tests(sorted(st['disk']))
        assert ex['planar walks containing a, c, d, e'] == cdict(st['single']) == {}
        ann = cdict(st['pairs'])
        assert ex['annular pairs (faces of a planar rotation together containing a, c, d, e)'] == ann
        assert set(ann) <= {'ac | de'}
        ts = sorted({tuple(H_edges(s)) for s in st['disk']} | ({('ac', 'de')} if ann else set()))
        mins = [list(t) for t in ts if not any(set(u) < set(t) for u in ts)]
        assert ex['tests'] == [list(t) for t in ts] and ex['tests after omitting supersets'] == mins
    table = d['claims'][-1]['computed']
    want = {
        'bd, or {a,b,d}; n_e = 2': [['ac', 'ad', 'ae'], ['ac', 'de']],
        'bd, or {a,b,d}; n_e = 3': [['ac', 'ad', 'ae', 'de']],
        'bd,cd, or {a,b,d} and cd; n_e = 2': [['ac', 'ad', 'ae', 'cd'], ['ac', 'cd', 'de']],
        'bd,cd, or {a,b,d} and cd; n_e = 3': [['ac', 'ae', 'cd', 'de']],
        '{a,c,d} and bd, or {a,b,d},{a,c,d}; n_e = 2': [['ac', 'ad', 'ae', 'cd'], ['ac', 'ad', 'cd', 'de'],
                                                        ['ac', 'ae', 'cd', 'de']],
        '{a,c,d} and bd, or {a,b,d},{a,c,d}; n_e = 3': [['ac', 'ae', 'cd', 'de']]}
    for k, alts in table.items():
        assert len(alts) == 2 and all(sorted(a['tests after omitting supersets']) == want[k] for a in alts), k
    for ex in d['exteriors']:
        key = '%s; n_e = %d' % (ex['row of the table'], ex['n_e'])
        assert sorted(sorted(t) for t in ex['tests after omitting supersets']) == want[key]
    return 'R8 ok'


def check_R9():
    d = load('R9')
    (rec,) = d['graphs']
    st = {'g1': Counter(), 'planar': Counter(), 'single': 0, 'pairs': Counter()}

    def hook(g, i, faces, gen, seqs):
        ws = [word(canon2(s), '') for s in seqs]
        if gen == 1:
            for s in seqs:
                r = restrict(s, 'bcd')
                if set('bcd') <= set(r):
                    st['g1'][word(canon2(r), '')] += 1
        if gen == 0:
            st['planar'][' '.join(sorted(ws))] += 1
            st['single'] += sum(set('bcd') <= set(w) for w in ws)
            ps = set()
            for x, y in itertools.combinations(ws, 2):
                if set('bcd') <= set(x) | set(y):
                    ps.add(' | '.join(sorted([x, y])))
            for p in ps:
                st['pairs'][p] += 1
    g, by = run_graph(rec, '', None, hook)
    c = claims_of(d)
    assert c[0]['faces containing b, c, d'] == cdict(st['g1']) and sorted(st['g1']) == ['bcd', 'bcdc']
    assert c[0]['H'] == tests(['bcd', 'bcdc'])
    assert c[1] == {'planar rotations by their face sequences': cdict(st['planar']), 'walks containing b, c, d':
                    st['single'], 'pairs': cdict(st['pairs'])}
    assert set(st['planar']) == {'ab ace ace ad bc de'} and st['single'] == 0 and \
        set(st['pairs']) == {'ad | bc', 'bc | de'}
    return 'R9 ok'


CHECKS = {
    'R1': check_R1, 'R2': check_R2, 'R3': check_R3, 'R4': check_R4, 'R5': check_R5, 'R6': check_R6,
    'R7': lambda: abd_check('R7', lambda n: {'abad'}),
    'R8': check_R8, 'R9': check_R9,
    'R10': lambda: abd_check('R10', lambda n: {'abad', 'abd'} if '{a,c,d}' in n else {'abad'}),
    'R11': lambda: abd_check('R11', lambda n: {'abad', 'abd', 'adbd'}),
}

def check(root, item):
    global ROOT
    ROOT = os.path.join(root, 'rotations')
    return CHECKS[item]().replace(item + ' ok', 'every rotation system retraced and every claim recomputed')
