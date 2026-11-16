"""The rotation tables R1-R12 (folder rotations of the data folder).

For each graph of a table, every rotation system is run through: its facial
walks are traced, its genus computed, and the attachment sequence of each
facial walk recorded (rotlib.py). The rows, and every statement the text makes
about the table, are written in the format of the data folder; a statement
that fails stops the program. The passages of the paper quoted in each file
are read from ../paper_text.json.

Usage:
    python make_rotation_tables.py OUTPUT_FOLDER            (all twelve tables)
    python make_rotation_tables.py OUTPUT_FOLDER R1 R5      (some of them)

The files are written to OUTPUT_FOLDER/rotations. R4 takes a few minutes and
R8 about a quarter of an hour. Then compare them with the data folder:
    python ../compare_files.py ../../supplemental_data OUTPUT_FOLDER rotations
"""
import gzip
import itertools
import json
import os
import sys
from array import array
from collections import Counter, defaultdict

from rotlib import *
from rgraphs import *
from make_R12 import generate as build_R12

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = None                         # set by main()
GZ_LIMIT = 25 * 1024 * 1024        # larger files are written gzip-compressed
with open(os.path.join(HERE, '..', 'paper_text.json'), encoding='utf-8') as _f:
    PAPER_TEXT = json.load(_f)


def paper_passages(label, before=None, after=None):
    """The passages of the paper quoted in the table (from ../paper_text.json)."""
    return PAPER_TEXT[label]


def write(relpath, obj):
    """Write JSON under OUT: indented when short, compact otherwise, and
    gzip-compressed (as .json.gz) above GZ_LIMIT bytes."""
    path = os.path.join(OUT, relpath)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    text = json.dumps(obj, separators=(',', ':'), ensure_ascii=False)
    if len(text) < 2_000_000:
        text = json.dumps(obj, indent=1, ensure_ascii=False)
    data = (text + '\n').encode()
    if len(data) > GZ_LIMIT:
        path += '.gz'
        with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', compresslevel=9,
                                                    fileobj=raw, mtime=0) as f:
            f.write(data)
    else:
        with open(path, 'wb') as f:
            f.write(data)
    return path, len(data), os.path.getsize(path)


HOW_TO_READ = (
    "For each graph: 'cyclic orders' lists, for every vertex, all cyclic orders of its neighbours, each starting "
    "with its alphabetically first neighbour (the other neighbours follow in every possible order, in lexicographic "
    "order). Rotation system number i (counting from 0) takes at each vertex the cyclic order whose position is the "
    "corresponding digit of i, written in the mixed radix given by the numbers of cyclic orders in the order of "
    "'vertices', the last vertex varying fastest (the order of Python's itertools.product). Facial walks are traced "
    "by the rule used throughout the supplement: after the dart u->v comes v->w, where w follows u in the cyclic "
    "order at v. The genus is (2 - V + E - F)/2. The attachment sequence of a facial walk lists the attachments it "
    "passes through, in walk order, one entry per passage; it is written from its lexicographically least cyclic "
    "shift, keeping the direction of the walk, and the faces of a row are sorted. Since many rotation systems give "
    "the same row, the distinct rows are listed once under 'distinct rows' and 'rows' gives, for every rotation "
    "system in order, the position of its row in that list. 'Up to cyclic shift and reversal' in the claims means "
    "the least sequence among all cyclic shifts of the sequence and of its reverse.")


def S(m):
    return 's%d' % m


def Tt(m):
    return 't%d' % m


def restrict(seq, keep):
    return tuple(x for x in seq if x in keep)


def tabulate(T, extra=None, collect=None):
    """Run through every rotation system of T. extra(ds, faces, g, seqs) adds
    hashable information to the row; collect(ds, faces, g, seqs) gathers
    statistics for the claims."""
    index = {}
    distinct = []
    rows = array('I')
    by_genus = Counter()
    for ds in T.all_digits():
        faces = T.faces_of_digits(ds)
        g = T.genus(faces)
        seqs = [T.seq(f) for f in faces]
        by_genus[g] += 1
        r = (g, tuple(sorted(fmt(least_shift(s), T.sep) for s in seqs)))
        if extra is not None:
            r = r + (extra(ds, faces, g, seqs),)
        k = index.get(r)
        if k is None:
            k = index[r] = len(distinct)
            distinct.append(r)
        rows.append(k)
        if collect is not None:
            collect(ds, faces, g, seqs)
    assert len(rows) == T.count
    return distinct, rows, dict(sorted(by_genus.items()))


def graph_record(T, distinct, rows, by_genus, row_format, extra_name=None):
    rec = T.header()
    rec['rotation systems by genus'] = {str(g): n for g, n in by_genus.items()}
    rec['row format'] = row_format
    if extra_name is None:
        rec['distinct rows'] = [[g, list(fs)] for g, fs in distinct]
    else:
        rec['distinct rows'] = [[g, list(fs), ex] for g, fs, ex in distinct]
    rec['rows'] = list(rows)
    return rec


ROW = '[genus, attachment sequences of the facial walks]'


def claim(statement, computed, holds, where=None):
    assert holds is True, ('claim fails', statement, computed)
    c = {'statement': statement, 'computed': computed, 'holds': True}
    if where:
        c['where'] = where
    return c


def counter_dict(c):
    return {k: c[k] for k in sorted(c)}


# ------------------------------------------------------------------ R1


def lemma42_triples(seqs):
    hits = []
    got = sorted(canon(s) for s in seqs)
    for k in (1, 2, 3):
        i, j = sorted({1, 2, 3} - {k})
        for x in (i, j):
            for y in (i, j):
                xp = j if x == i else i
                yp = j if y == i else i
                want = sorted([canon((S(i), Tt(i), Tt(j), S(j)))] * 2 +
                              [canon((S(k), Tt(k), Tt(y), Tt(k), S(k), S(x))),
                               canon((S(k), S(xp))), canon((Tt(k), Tt(yp)))])
                if got == want:
                    hits.append((k, x, y))
    return hits


def automorphisms(T):
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(T.vertices)
    G.add_edges_from(T.edges)
    gm = nx.algorithms.isomorphism.GraphMatcher(G, G)
    return [dict(m) for m in gm.isomorphisms_iter()]


def digits_of_rotation(T, rot):
    out = []
    for v in T.vertices:
        cyc = rot[v]
        i = cyc.index(min(cyc))
        cyc = cyc[i:] + cyc[:i]
        out.append(T.orders[v].index(cyc))
    return tuple(out)


def orbits(T, members, autos, mirror):
    members = set(members)
    parent = {m: m for m in members}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for ds in members:
        rot = T.rotation(T_index(T, ds))
        imgs = []
        for pi in autos:
            r2 = {pi[v]: [pi[u] for u in rot[v]] for v in rot}
            imgs.append(r2)
            if mirror:
                imgs.append({v: list(reversed(c)) for v, c in r2.items()})
        for r2 in imgs:
            d2 = digits_of_rotation(T, r2)
            assert d2 in members, 'an automorphism image is not a torus rotation'
            a, b = find(ds), find(d2)
            if a != b:
                parent[a] = b
    groups = defaultdict(list)
    for m in members:
        groups[find(m)].append(m)
    return sorted(sorted(g) for g in groups.values())


def T_index(T, ds):
    i = 0
    for r, d in zip(T.radix, ds):
        i = i * r + d
    return i


def build_R1():
    T = R1_graph()
    paths = [(S(m), Tt(m)) for m in (1, 2, 3)]
    torus = []
    triple_of = {}
    by_triple = defaultdict(Counter)
    five = Counter()

    def extra(ds, faces, g, seqs):
        if g != 1:
            return None
        return {'P%d' % m: region_key(T.region_along(faces, *paths[m - 1]), ' ') for m in (1, 2, 3)}

    def extra_h(ds, faces, g, seqs):
        e = extra(ds, faces, g, seqs)
        return None if e is None else tuple(sorted(e.items()))

    def collect(ds, faces, g, seqs):
        if g != 1:
            return
        torus.append(ds)
        five[len(faces)] += 1
        hits = lemma42_triples(seqs)
        assert len(hits) == 1, hits
        k, x, y = hits[0]
        triple_of[ds] = hits[0]
        regs = extra(ds, faces, g, seqs)
        pattern = (tuple(sorted(fmt(canon(s), ' ') for s in seqs)), tuple(sorted(regs.items())))
        by_triple[hits[0]][pattern] += 1

    distinct, rows, by_genus = tabulate(T, extra_h, collect)
    rec = graph_record(T, [(g, fs, (dict(ex) if ex else None)) for g, fs, ex in distinct], rows, by_genus,
                       '[genus, attachment sequences of the facial walks, for a torus rotation the region along each '
                       'path P_m (the faces on the two sides of the edge s_m t_m glued along it: a disk, or an annulus '
                       'when the same face lies on both sides), up to cyclic shift and reversal]',
                       extra_name='regions')
    claims = []
    claims.append(claim('Psi_1 is nonplanar and every torus embedding has five faces (Lemma 4.2).',
                        {'rotation systems': T.count, 'by genus': by_genus, 'faces of the genus-one rotations':
                         dict(five)}, by_genus.get(0, 0) == 0 and set(five) == {5}))
    autos = automorphisms(T)
    o_kept = orbits(T, torus, autos, mirror=False)
    o_mirror = orbits(T, torus, autos, mirror=True)
    claims.append(claim('Up to automorphisms of Psi_1, there are exactly two torus embeddings (proof of Lemma 4.2).',
                        {'automorphisms of Psi_1': len(autos), 'torus rotation systems': len(torus),
                         'orbits under the automorphisms': len(o_kept),
                         'orbits under the automorphisms and reversal of all cyclic orders (mirror images)':
                             len(o_mirror),
                         'a representative row of each orbit': [T_index(T, o[0]) for o in o_mirror]},
                        len(o_kept) == 2 and len(o_mirror) == 2))
    occ = Counter(triple_of.values())
    choice_table = {}
    for (k, x, y) in sorted(by_triple):
        pats = by_triple[(k, x, y)]
        assert len(pats) == 1, ('several face patterns for one triple', (k, x, y))
        (fs, regs), n = next(iter(pats.items()))
        choice_table['%d%d%d' % (k, x, y)] = {'torus rotation systems': n, 'faces': list(fs), 'regions': dict(regs)}
    claims.append(claim(
        'Lemma 4.2: every torus rotation has, for exactly one triple (k,x,y), two simple faces with sequence '
        's_i t_i t_j s_j, a repeated face s_k t_k t_y t_k s_k s_x, and faces s_k s_x\' and t_k t_y\' '
        '({i,j} = {1,2,3} - {k}, {x,x\'} = {y,y\'} = {i,j}); all twelve triples occur.',
        {'torus rotation systems per triple': {'%d%d%d' % t: occ[t] for t in sorted(occ)},
         'faces and regions by triple (k,x,y)': choice_table},
        len(occ) == 12 and sum(occ.values()) == len(torus)))
    return {'item': 'R1', 'where': 'Lemma 4.2 (proof)',
            'paper text': paper_passages('R1', before=30, after=4),
            'what': 'Every rotation system of Psi_1 with P_1, P_2, P_3 the edges s1t1, s2t2, s3t3 (subdividing an edge '
                    'does not change the faces). Attachments: the path vertices s1, t1, s2, t2, s3, t3; twins '
                    'omitted. Vertex names in a sequence are separated by spaces.',
            'how to read': HOW_TO_READ, 'graphs': [rec], 'claims': claims}


# ------------------------------------------------------------------ R2


def patterns(T, edge):
    """Face patterns (faces with at least three distinct attachments) and regions
    along edge, over the torus rotations."""
    faces3, regions = Counter(), Counter()
    for ds in T.all_digits():
        f = T.faces_of_digits(ds)
        if T.genus(f) != 1:
            continue
        for x in f:
            s = T.seq(x)
            if len(set(s)) >= 3:
                faces3[fmt(canon(s))] += 1
        regions[region_key(T.region_along(f, *edge))] += 1
    return faces3, regions


def is_three_restriction(p, of='abdc'):
    if len(p) != 3 or len(set(p)) != 3:
        return False
    return any(canon(restrict(of, set(of) - {x})) == canon(p) for x in of)


def build_R2():
    T = R2_graph()
    P = R2_psi1_graph()

    def extra(ds, faces, g, seqs):
        return region_key(T.region_along(faces, 'a', 'b')) if g == 1 else None
    distinct, rows, by_genus = tabulate(T, extra)
    rec = graph_record(T, distinct, rows, by_genus,
                       '[genus, attachment sequences of the facial walks, for a torus rotation the region along '
                       'P_1 = ab (the faces on the two sides of ab glued along it: a disk, or an annulus when the same '
                       'face lies on both sides), up to cyclic shift and reversal]', extra_name='region along ab')
    f_new, r_new = patterns(T, ('a', 'b'))
    f_psi, r_psi = patterns(P, ('a', 'b'))
    new_faces = sorted(set(f_new) - set(f_psi))
    new_regions = sorted(set(r_new) - set(r_psi))
    new_all = sorted(set(new_faces) | {r for r in new_regions})
    claims = []
    claims.append(claim('Psi_1/s3t3 is connected and nonplanar, so its torus embeddings are cellular.',
                        {'rotation systems': T.count, 'by genus': by_genus}, by_genus.get(0, 0) == 0))
    claims.append(claim(
        'Comparison with Psi_1 (the table R1, with a = s1, b = t1, c = s2, d = t2): the faces meeting at least three '
        'of a, b, c, d and the regions glued along P_1 = ab, over all torus rotations, up to cyclic shift and '
        'reversal, with the number of torus rotations in which each occurs.',
        {'Psi_1/s3t3': {'faces': counter_dict(f_new), 'regions along ab': counter_dict(r_new)},
         'Psi_1': {'faces': counter_dict(f_psi), 'regions along ab': counter_dict(r_psi)}}, True))
    ok = all(p == canon('abcd') or is_three_restriction(p) for p in map(tuple, new_all)) and \
        all('annulus' not in r for r in new_regions)
    claims.append(claim(
        'The only new disk attachment sequences are a,b,c,d and three-vertex orders that, up to reversal, are '
        'restrictions of a,b,d,c; no new annulus occurs.',
        {'faces of Psi_1/s3t3 not occurring for Psi_1': new_faces,
         'regions along ab of Psi_1/s3t3 not occurring for Psi_1': new_regions}, ok))
    claims.append(claim('In (a) and (b), a simple face of Psi_1 has order a,b,d,c.',
                        {'faces of Psi_1 with sequence abdc (torus rotations)': f_psi.get(fmt(canon('abdc')), 0)},
                        f_psi.get(fmt(canon('abdc')), 0) > 0))
    claims.append(claim('In (c), gluing the two simple facial disks of Psi_1 along P_1 gives a disk with attachment '
                        'sequence a,c,d,b,d,c.',
                        {'regions of Psi_1 along ab with sequence acdbdc (torus rotations)':
                         r_psi.get(fmt(canon('acdbdc')), 0)}, r_psi.get(fmt(canon('acdbdc')), 0) > 0))
    return {'item': 'R2', 'where': 'Lemma 4.8 (proof)', 'paper text': paper_passages('R2', before=22, after=16),
            'what': 'Every rotation system of Psi_1/s3t3, with P_1, P_2 the edges ab, cd (a = s1, b = t1, c = s2, '
                    'd = t2) and the edge s3t3 contracted into s3. Attachments a, b, c, d; twins and s3 omitted.',
            'how to read': HOW_TO_READ, 'graphs': [rec], 'claims': claims}


# ------------------------------------------------------------------ R3


def build_R3():
    T = R3_graph()
    pat = Counter()

    def extra(ds, faces, g, seqs):
        return region_key(T.region_along(faces, 'a', 'b')) if g == 1 else None

    def collect(ds, faces, g, seqs):
        if g != 1:
            return
        reg = T.region_along(faces, 'a', 'b')
        att = set(reg[1]) | (set(reg[2]) if reg[0] == 'annulus' else set())
        if not {'c', 'd'} <= att:
            return
        dab, dba = T.dix[('a', 'b')], T.dix[('b', 'a')]
        others = sorted(fmt(canon(T.seq(x))) for x in faces if dab not in x and dba not in x)
        pat[(region_key(reg), ', '.join(others))] += 1
    distinct, rows, by_genus = tabulate(T, extra, collect)
    rec = graph_record(T, distinct, rows, by_genus,
                       '[genus, attachment sequences of the facial walks, for a torus rotation the region along '
                       'P_1 = ab, up to cyclic shift and reversal]', extra_name='region along ab')
    want = {('ac | bd (annulus)', 'a, bc, cd, cd'), ('acbdc', 'a, bcd, cd'), ('acbcdc', 'a, bd, cd')}
    claims = [claim('Psi_1/cf is connected and nonplanar, so its torus embeddings are cellular.',
                    {'rotation systems': T.count, 'by genus': by_genus}, by_genus.get(0, 0) == 0),
              claim('The torus embeddings of Psi_1/cf give exactly three regions incident with P_1 whose boundary '
                    'contains both c, d: ac | bd (annulus) with a, bc, cd, cd in the other faces; a,c,b,d,c with '
                    'a, bcd, cd; a,c,b,c,d,c with a, bd, cd (up to cyclic shift and reversal).',
                    [{'region along ab': r, 'other faces': o, 'torus rotation systems': n}
                     for (r, o), n in sorted(pat.items())], set(pat) == want)]
    return {'item': 'R3', 'where': 'Proposition 4.25 (proof)', 'paper text': paper_passages('R3', before=6, after=10),
            'what': 'Every rotation system of Psi_1/cf in the notation of Proposition 4.25 (a = s1, b = t1, c = s2, '
                    'd = t3, e = s3, f = t2; the paths are the edges ab, cf, ed, and f is contracted into c). '
                    'Attachments a, b, c, d; twins and e omitted.',
            'how to read': HOW_TO_READ, 'graphs': [rec], 'claims': claims}


# ------------------------------------------------------------------ Section 5


def disk_stats(T, need, keep=None):
    """Collector of torus faces containing all of need (restricted to keep)."""
    st = {'disk': Counter(), 'planar faces': Counter(), 'planar pairs': Counter(), 'planar single': Counter()}
    keep = set(keep or T.attachments)

    def collect(ds, faces, g, seqs):
        rs = [restrict(s, keep) for s in seqs]
        if g == 1:
            for r in rs:
                if set(need) <= set(r):
                    st['disk'][fmt(canon(r))] += 1
        if g == 0:
            st['planar faces'][' '.join(sorted(fmt(canon(r)) for r in rs))] += 1
            for r in rs:
                if set(need) <= set(r):
                    st['planar single'][fmt(canon(r))] += 1
            pairs = set()
            for i, j in itertools.combinations(range(len(rs)), 2):
                if set(need) <= set(rs[i]) | set(rs[j]):
                    pairs.add(' | '.join(sorted([fmt(canon(rs[i])), fmt(canon(rs[j]))])))
            for p in pairs:
                st['planar pairs'][p] += 1
    return st, collect


def tests_of(seqs):
    out = {}
    for s in seqs:
        out[s] = {'H': shape(tuple(s)), 'edges of H': sorted(''.join(e) for e in consecutive_edges(tuple(s)))}
    return out


def build_R4():
    graphs, claims = [], []
    want_disk = {(3, 3): ['acbd'], (3, 2): ['abcad', 'abdbc', 'acbd']}
    for (nc, nd), T in R4_graphs():
        st, col = disk_stats(T, 'abcd')
        cd = {'rotations with a walk containing c and d': 0, 'alternation in a face with three labels': 0,
              'alternation in any face': 0}

        def collect(ds, faces, g, seqs, col=col):
            col(ds, faces, g, seqs)
            if (nc, nd) == (2, 2) and g == 1 and any({'c', 'd'} <= set(s) for s in seqs):
                cd['rotations with a walk containing c and d'] += 1
                for s in seqs:
                    if alternates(s):
                        cd['alternation in any face'] += 1
                        if len(set(s)) >= 3:
                            cd['alternation in a face with three labels'] += 1
        distinct, rows, by_genus = tabulate(T, None, collect)
        graphs.append(graph_record(T, distinct, rows, by_genus, ROW))
        disk = sorted(st['disk'])
        tests = tests_of(disk)
        key = '(n_c, n_d) = (%d, %d)' % (nc, nd)
        claims.append(claim('%s: rotation systems by genus.' % key, {'rotation systems': T.count, 'by genus': by_genus},
                            ((nc, nd) == (2, 2)) == (by_genus.get(0, 0) > 0), where='Lemma 5.3'))
        claims.append(claim(
            '%s: the genus-one attachment sequences meeting all four attachments, up to cyclic shift and reversal, '
            'with H obtained by joining consecutive attachments; each H is a cycle, a triangle with a pendant edge '
            'or a three-edge star, and the sequence is a facial walk of H.' % key,
            {'sequences (torus rotations in which they occur)': counter_dict(st['disk']), 'H': tests},
            all(t['H'] in ('cycle', 'triangle with a pendant edge', 'three-edge star') for t in tests.values())
            and ((nc, nd) not in want_disk or disk == want_disk[(nc, nd)]), where='Lemma 5.3'))
        if (nc, nd) == (2, 2):
            planar = st['planar faces']
            six = 'ab ab ac ad bc bd'
            claims.append(claim(
                'Only (n_c, n_d) = (2, 2) has planar rotations. Their boundary attachment sequences are ab, ab, ac, '
                'bc, ad, bd. No single walk contains all four attachments, and both annular pairs ac | bd and '
                'ad | bc occur.',
                {'planar rotations by their face sequences': counter_dict(planar),
                 'planar walks containing a, b, c, d': counter_dict(st['planar single']),
                 'pairs of faces of a planar rotation together containing a, b, c, d':
                     counter_dict(st['planar pairs'])},
                set(planar) == {six} and not st['planar single'] and
                set(st['planar pairs']) == {'ac | bd', 'ad | bc'}, where='Lemma 5.3'))
            four_five = [s for s in disk if len(s) in (4, 5)]
            six_l = [s for s in disk if len(s) == 6]
            ok = all({'ac', 'bd'} <= set(tests[s]['edges of H']) or {'ad', 'bc'} <= set(tests[s]['edges of H'])
                     for s in four_five) and sorted(tuple(tests[s]['edges of H']) for s in six_l) == \
                [('ab', 'ac', 'ad'), ('ab', 'bc', 'bd')]
            claims.append(claim(
                'For (2, 2) in genus one, the sequences of length four or five give graphs H containing ac, bd or '
                'ad, bc; the remaining sequences give the tests B+ab+ac+ad and B+ab+bc+bd, and each occurs.',
                {'length four or five': four_five, 'length six': six_l}, ok, where='Lemma 5.3'))
            # Lemmas 5.30 and 5.32: genus-zero boundary sequences and the pairs containing c and d
            pairs_cd = set()
            for p in planar:
                fs = p.split()
                for i, j in itertools.combinations(range(len(fs)), 2):
                    if {'c', 'd'} <= set(fs[i]) | set(fs[j]):
                        pairs_cd.add(' | '.join(sorted([fs[i], fs[j]])))
            single_bad = [f for p in planar for f in p.split()
                          if len(set(f)) >= 3 or {'c', 'd'} <= set(f)]
            claims.append(claim(
                'When a closed regular neighbourhood of Psi_2 has genus zero, its boundary attachment sequences are '
                'ab, ab, ac, bc, ad, bd. No single boundary walk contains three attachments or both c and d, and the '
                'pairs of walks together containing c and d are ac | bd, ad | bc, ac | ad, bc | bd.',
                {'pairs containing c and d': sorted(pairs_cd), 'walks with three attachments or both c, d':
                 single_bad},
                set(planar) == {six} and not single_bad and
                pairs_cd == {'ac | bd', 'ad | bc', 'ac | ad', 'bc | bd'}, where='Lemmas 5.30 and 5.32'))
            claims.append(claim(
                'There are 576 rotations of Psi_2; 80 of genus one have a facial walk containing both c and d, and in '
                'every face of these rotations containing at least three attachment labels, no two labels alternate '
                'as u, v, u, v.', dict(cd),
                T.count == 576 and cd['rotations with a walk containing c and d'] == 80 and
                cd['alternation in a face with three labels'] == 0, where='Lemma 5.29'))
        if (nc, nd) == (3, 2):
            claims.append(claim(
                'For (3, 2), the sequences abcad and abdbc are facial walks of the triangle abc with pendant edge ad '
                'or bd, giving B+ab+ac+bc+ad and B+ab+ac+bc+bd; acbd gives B+ac+bc+bd+ad.',
                {s: tests[s] for s in disk},
                tests['abcad']['edges of H'] == ['ab', 'ac', 'ad', 'bc'] and
                tests['abdbc']['edges of H'] == ['ab', 'ac', 'bc', 'bd'] and
                tests['acbd']['edges of H'] == ['ac', 'ad', 'bc', 'bd'], where='Lemma 5.3'))
    return {'item': 'R4', 'where': 'Lemma 5.3 (proof), Lemma 5.29 (proof), Lemma 5.30 (proof) and Lemma 5.32 (proof)',
            'paper text': paper_passages('R4', before=12, after=10),
            'what': 'Every rotation system of Psi_2 with n_c vertices of neighbourhood {a,b,c} (p, q, and p3 when '
                    'n_c = 3) and n_d of neighbourhood {a,b,d} (r, z, and r3 when n_d = 3), for (n_c, n_d) = (2, 2), '
                    '(3, 2) and (3, 3). Attachments a, b, c, d; twins omitted.',
            'how to read': HOW_TO_READ, 'graphs': graphs, 'claims': claims}


def build_R5():
    T = R5_graph()
    st = {'planar': Counter(), 'with d': Counter(), 'with d raw': Counter()}

    def collect(ds, faces, g, seqs):
        red = [combine_repeats(restrict(s, 'acd')) for s in seqs]
        if g == 0:
            st['planar'][' '.join(sorted(fmt(canon(r)) for r in red))] += 1
        if g == 1:
            for s, r in zip(seqs, red):
                if 'd' in r:
                    st['with d'][fmt(canon(r))] += 1
                    st['with d raw'][fmt(canon(s))] += 1
    distinct, rows, by_genus = tabulate(T, None, collect)
    rec = graph_record(T, distinct, rows, by_genus, ROW)
    claims = [claim(
        'Retaining only a, c, d in each boundary walk and combining consecutive equal appearances: up to cyclic shift '
        'and reversal the boundary attachment sequences are a, ac, ad, c when the neighbourhood has genus zero, and '
        'the sequence containing d is ad, acd, acad or acacd when it has genus one.',
        {'rotation systems': T.count, 'by genus': by_genus,
         'genus zero: the four reduced sequences': counter_dict(st['planar']),
         'genus one: the reduced sequence containing d': counter_dict(st['with d']),
         'genus one: that walk before the reduction (a, b, c, d retained)': counter_dict(st['with d raw'])},
        set(st['planar']) == {'a ac ad c'} and set(st['with d']) <= {'ad', 'acd', 'acad', 'acacd'})]
    return {'item': 'R5', 'where': 'Lemma 5.5 (proof)', 'paper text': paper_passages('R5', before=10, after=4),
            'what': 'Every rotation system of Psi_2 - z. Attachments a, b, c, d; twins omitted. The claim retains '
                    'only a, c, d, as the text does.',
            'how to read': HOW_TO_READ, 'graphs': [rec], 'claims': claims}


# ------------------------------------------------------------------ Section 6


def swaps(seqs_str):
    """The multiset of sequences under interchanging b,c and/or d,e."""
    out = set()
    for sw in ({}, {'b': 'c', 'c': 'b'}, {'d': 'e', 'e': 'd'}, {'b': 'c', 'c': 'b', 'd': 'e', 'e': 'd'}):
        out.add(' '.join(sorted(fmt(canon(tuple(sw.get(x, x) for x in s))) for s in seqs_str.split())))
    return out


def build_R6():
    T = R6_graph()
    st = {'faces': Counter(), 'g1 bcde': Counter(), 'g1 others': Counter(), 'planar': Counter(),
          'planar single': 0, 'planar pairs': Counter(), 'formula': True}

    def collect(ds, faces, g, seqs):
        st['formula'] &= (2 * g == 5 - len(faces))
        if g == 1:
            big = [s for s in seqs if set('bcde') <= set(s)]
            if big:
                assert len(big) == 1
                st['g1 bcde'][fmt(canon(big[0]))] += 1
                st['g1 others'][' '.join(sorted(fmt(canon(s)) for s in seqs if s is not big[0]))] += 1
        if g == 0:
            st['planar'][' '.join(sorted(fmt(canon(s)) for s in seqs))] += 1
            st['planar single'] += sum(set('bcde') <= set(s) for s in seqs)
            for i, j in itertools.combinations(range(len(seqs)), 2):
                if set('bcde') <= set(seqs[i]) | set(seqs[j]):
                    st['planar pairs'][' | '.join(sorted([fmt(canon(seqs[i])), fmt(canon(seqs[j]))]))] += 1
    distinct, rows, by_genus = tabulate(T, None, collect)
    rec = graph_record(T, distinct, rows, by_genus, ROW)
    claims = [
        claim('There are 6 * 2^4 = 96 rotations; a rotation with f facial walks has genus (5 - f)/2.',
              {'rotation systems': T.count, 'by genus': by_genus}, T.count == 96 and st['formula']),
        claim('In genus one, a walk containing b, c, d, e has sequence a,b,a,d,a,c,a,e up to cyclic shift and '
              'reversal, and the other two walks have sequences bc, de.',
              {'that walk': counter_dict(st['g1 bcde']), 'the other two walks': counter_dict(st['g1 others'])},
              set(st['g1 bcde']) == {'abadacae'} and set(st['g1 others']) == {'bc de'}),
        claim('In genus zero, the sequences are abad, ac, ae, bc, de, up to interchanging b, c and interchanging d, e; '
              'no single walk meets b, c, d, e, and the only pair that does is bc, de.',
              {'planar rotations by their face sequences': counter_dict(st['planar']),
               'walks meeting b, c, d, e': st['planar single'],
               'pairs meeting b, c, d, e': counter_dict(st['planar pairs'])},
              all('abad ac ae bc de' in swaps(p) for p in st['planar']) and
              st['planar single'] == 0 and set(st['planar pairs']) == {'bc | de'})]
    return {'item': 'R6', 'where': 'Lemma 6.3 (proof)', 'paper text': paper_passages('R6', before=20, after=12),
            'what': 'Every rotation system of Psi_3. Attachments a, b, c, d, e; twins omitted.',
            'how to read': HOW_TO_READ, 'graphs': [rec], 'claims': claims}


def abd_stats(T):
    st = {'abd': Counter(), 'raw': Counter()}

    def collect(ds, faces, g, seqs):
        if g == 1:
            for s in seqs:
                r = restrict(s, 'abd')
                if set('abd') <= set(r):
                    st['abd'][fmt(canon(r))] += 1
                    st['raw'][fmt(canon(s))] += 1
    return st, collect


def build_R7():
    T = R7_graph()
    st, col = abd_stats(T)
    distinct, rows, by_genus = tabulate(T, None, col)
    rec = graph_record(T, distinct, rows, by_genus, ROW)
    claims = [claim('The exterior contains K_{3,3}, so it has no planar rotation; in every genus-one rotation, each '
                    'facial walk containing a, b, d has attachment sequence abad, up to cyclic shift and reversal '
                    '(retaining a, b, d).',
                    {'rotation systems': T.count, 'by genus': by_genus,
                     'walks containing a, b, d, retaining a, b, d': counter_dict(st['abd']),
                     'the same walks retaining a, b, c, d, e': counter_dict(st['raw'])},
                    by_genus.get(0, 0) == 0 and set(st['abd']) == {'abad'})]
    return {'item': 'R7', 'where': 'Lemma 6.14 (proof)', 'paper text': paper_passages('R7', before=8, after=4),
            'what': 'Every rotation system of Psi_3 + w with N(w) = {a,c,d,e}. Attachments a, b, c, d, e; twins and '
                    'w omitted. The claim retains a, b, d, the attachments of B_1.',
            'how to read': HOW_TO_READ, 'graphs': [rec], 'claims': claims}


R8_TABLE = {
    ('bd, or {a,b,d}', 2): [['ac', 'de'], ['ac', 'ad', 'ae']],
    ('bd, or {a,b,d}', 3): [['ac', 'ad', 'ae', 'de']],
    ('bd,cd, or {a,b,d} and cd', 2): [['ac', 'cd', 'de'], ['ac', 'ad', 'ae', 'cd']],
    ('bd,cd, or {a,b,d} and cd', 3): [['ac', 'ae', 'cd', 'de']],
    ('{a,c,d} and bd, or {a,b,d},{a,c,d}', 2): [['ac', 'ad', 'ae', 'cd'], ['ac', 'ad', 'cd', 'de'],
                                               ['ac', 'cd', 'de', 'ae']],
    ('{a,c,d} and bd, or {a,b,d},{a,c,d}', 3): [['ac', 'ae', 'cd', 'de']],
}


def build_R8():
    graphs, per = [], []
    for x in R8_graphs():
        T = x['table']
        st, col = disk_stats(T, 'acde', keep='acde')
        distinct, rows, by_genus = tabulate(T, None, col)
        graphs.append(graph_record(T, distinct, rows, by_genus, ROW))
        disk = sorted(st['disk'])
        tests = tests_of(disk)
        ann = sorted(st['planar pairs'])
        all_tests = [tuple(tests[s]['edges of H']) for s in disk]
        if by_genus.get(0, 0):
            assert ann == ['ac | de'], ann
            all_tests.append(('ac', 'de'))
        all_tests = sorted(set(all_tests))
        minimal = [t for t in all_tests if not any(set(u) < set(t) for u in all_tests)]
        per.append({'exterior': T.name, 'row of the table': x['row'], 'n_e': x['n_e'],
                    'alternative': x['alternative'], 'rotation systems': T.count, 'by genus': by_genus,
                    'disk sequences (torus rotations in which they occur)': counter_dict(st['disk']),
                    'H for each disk sequence': tests,
                    'planar walks containing a, c, d, e': counter_dict(st['planar single']),
                    'annular pairs (faces of a planar rotation together containing a, c, d, e)':
                        counter_dict(st['planar pairs']),
                    'tests': [list(t) for t in all_tests], 'tests after omitting supersets': [list(t) for t in minimal]})
    claims = []
    ok_shape = all(t['H'] in ('cycle', 'triangle with a pendant edge', 'three-edge star', 'path')
                   for p in per for t in p['H for each disk sequence'].values())
    claims.append(claim(
        'If the regular neighbourhood has genus one, Lemma 2.2 applies to each disk sequence: H is a cycle, a '
        'triangle with a pendant edge, a three-edge star or a path, and the sequence is a facial walk of H.',
        {p['exterior']: p['H for each disk sequence'] for p in per}, ok_shape))
    ok_g0 = all(not p['planar walks containing a, c, d, e'] and
                set(p['annular pairs (faces of a planar rotation together containing a, c, d, e)']) <= {'ac | de'}
                for p in per)
    claims.append(claim(
        'In genus zero, no single walk meets all four attachments and the only possible annulus is ac | de.',
        {p['exterior']: {'planar rotations': p['by genus'].get(0, 0),
                         'annular pairs': p['annular pairs (faces of a planar rotation together containing a, c, d, '
                                            'e)']} for p in per}, ok_g0))
    ex = per[0]
    claims.append(claim(
        'For example, Psi_3 + bd has disk sequences acade, acadae and the annular pair ac | de; their tests add '
        '{ac,ad,ae,de}, {ac,ad,ae} and {ac,de}, and the first contains the others.',
        {'disk sequences': sorted(ex['disk sequences (torus rotations in which they occur)']),
         'annular pairs': sorted(ex['annular pairs (faces of a planar rotation together containing a, c, d, e)']),
         'tests': ex['tests']},
        ex['exterior'] == 'Psi_3 + bd, n_e = 2' and
        sorted(ex['disk sequences (torus rotations in which they occur)']) == ['acadae', 'acade'] and
        sorted(ex['annular pairs (faces of a planar rotation together containing a, c, d, e)']) == ['ac | de']))
    table_ok = True
    table = {}
    for p in per:
        want = sorted(sorted(t) for t in R8_TABLE[(p['row of the table'], p['n_e'])])
        got = sorted(sorted(t) for t in p['tests after omitting supersets'])
        table_ok &= want == got
        table.setdefault('%s; n_e = %d' % (p['row of the table'], p['n_e']), []).append(
            {'alternative': p['alternative'], 'tests after omitting supersets': got})
    claims.append(claim(
        'Omitting supersets gives the table of Lemma 6.15; the alternatives in each row have the same tests, and each '
        'retained test is realized by a rotation.', table, table_ok))
    return {'item': 'R8', 'where': 'Lemma 6.15 (proof)', 'paper text': paper_passages('R8', before=40, after=14),
            'what': 'Every rotation system of each of the twelve exteriors of Lemma 6.15: Psi_3, with a third vertex '
                    'r3 of neighbourhood {a,d,e} when n_e = 3, plus the additional edges or vertices of a row of the '
                    'table (v_abd and v_acd are the additional vertices with neighbourhoods {a,b,d} and {a,c,d}). '
                    'Attachments a, c, d, e: the text retains the appearances of a, c, d, e in the facial walks.',
            'how to read': HOW_TO_READ, 'graphs': graphs, 'exteriors': per, 'claims': claims}


def build_R9():
    T = R9_graph()
    st = {'g1': Counter(), 'planar': Counter(), 'single': 0, 'pairs': Counter()}

    def collect(ds, faces, g, seqs):
        if g == 1:
            for s in seqs:
                r = restrict(s, 'bcd')
                if set('bcd') <= set(r):
                    st['g1'][fmt(canon(r))] += 1
        if g == 0:
            st['planar'][' '.join(sorted(fmt(canon(s)) for s in seqs))] += 1
            st['single'] += sum(set('bcd') <= set(s) for s in seqs)
            ps = set()
            for i, j in itertools.combinations(range(len(seqs)), 2):
                if set('bcd') <= set(seqs[i]) | set(seqs[j]):
                    ps.add(' | '.join(sorted([fmt(canon(seqs[i])), fmt(canon(seqs[j]))])))
            for p in ps:
                st['pairs'][p] += 1
    distinct, rows, by_genus = tabulate(T, None, collect)
    rec = graph_record(T, distinct, rows, by_genus, ROW)
    g1 = sorted(st['g1'])
    claims = [
        claim('There are 384 rotations. The genus-one rotations give only bcd, bcdc for a face containing all three '
              'attachments b, c, d (retaining b, c, d), requiring B+bc+bd+cd or B+bc+cd planar.',
              {'rotation systems': T.count, 'by genus': by_genus, 'faces containing b, c, d': counter_dict(st['g1']),
               'H': tests_of(g1)},
              T.count == 384 and g1 == ['bcd', 'bcdc'] and tests_of(g1)['bcd']['edges of H'] == ['bc', 'bd', 'cd']
              and tests_of(g1)['bcdc']['edges of H'] == ['bc', 'cd']),
        claim('The planar rotations have sequences ab, ace, ace, ad, bc, de. No single boundary contains b, c, d, and '
              'the only pairs containing them are bc | ad and bc | de; a planar rotation has both a bc face and an ad '
              'face.',
              {'planar rotations by their face sequences': counter_dict(st['planar']),
               'walks containing b, c, d': st['single'], 'pairs': counter_dict(st['pairs'])},
              set(st['planar']) == {'ab ace ace ad bc de'} and st['single'] == 0 and
              set(st['pairs']) == {'ad | bc', 'bc | de'})]
    return {'item': 'R9', 'where': 'Lemma 6.18 (proof)', 'paper text': paper_passages('R9', before=6, after=16),
            'what': 'Every rotation system of Psi_3 + ce. Attachments a, b, c, d, e; twins omitted.',
            'how to read': HOW_TO_READ, 'graphs': [rec], 'claims': claims}


def build_R10():
    graphs, comp = [], {}
    ok = True
    for nb, T in R10_graphs():
        st, col = abd_stats(T)
        distinct, rows, by_genus = tabulate(T, None, col)
        graphs.append(graph_record(T, distinct, rows, by_genus, ROW))
        want = {'abad', 'abd'} if nb == 'acd' else {'abad'}
        ok &= by_genus.get(0, 0) == 0 and set(st['abd']) == want
        comp[T.name] = {'rotation systems': T.count, 'by genus': by_genus,
                        'genus-one walks containing a, b, d, retaining a, b, d': counter_dict(st['abd']),
                        'the same walks retaining a, b, c, d, e': counter_dict(st['raw'])}
    claims = [claim('Enumerating all rotations of each exterior gives no planar rotations. The genus-one facial walks '
                    'containing a, b, d have sequence abad, with abd also possible for the {a,c,d} exterior.',
                    comp, ok)]
    return {'item': 'R10', 'where': 'Lemma 6.19 (proof)', 'paper text': paper_passages('R10', before=8, after=6),
            'what': 'Every rotation system of the exteriors Psi_3 + ce + v with N(v) = {a,b,c}, {a,d,e}, {a,c,d}, '
                    'and Psi_3 + v with N(v) = {c,d,e}. Attachments a, b, c, d, e; twins and v omitted.',
            'how to read': HOW_TO_READ, 'graphs': graphs, 'claims': claims}


def build_R11():
    graphs, comp = [], {}
    ok = True
    for lab, T in R11_graphs():
        st, col = abd_stats(T)
        distinct, rows, by_genus = tabulate(T, None, col)
        graphs.append(graph_record(T, distinct, rows, by_genus, ROW))
        ok &= by_genus.get(0, 0) == 0 and set(st['abd']) == {'abad', 'abd', 'adbd'}
        comp[T.name] = {'rotation systems': T.count, 'by genus': by_genus,
                        'genus-one walks containing a, b, d, retaining a, b, d': counter_dict(st['abd']),
                        'the same walks retaining a, b, c, d, e': counter_dict(st['raw'])}
    claims = [claim('Both exteriors have no planar rotation (they contain K_{3,3}). Their genus-one facial walks '
                    'containing a, b, d give exactly abad, abd, adbd, and both four-letter sequences occur (for each '
                    'exterior).', comp, ok)]
    return {'item': 'R11', 'where': 'Lemma 6.19 (proof)', 'paper text': paper_passages('R11', before=10, after=10),
            'what': 'Every rotation system of Psi_3 + ce + v with N(v) = {b,c,d}, and of Psi_3 + ce + cd. '
                    'Attachments a, b, c, d, e; twins and v omitted.',
            'how to read': HOW_TO_READ, 'graphs': graphs, 'claims': claims}


BUILDERS = {'R1': build_R1, 'R2': build_R2, 'R3': build_R3, 'R4': build_R4, 'R5': build_R5, 'R6': build_R6,
            'R7': build_R7, 'R8': build_R8, 'R9': build_R9, 'R10': build_R10, 'R11': build_R11,
            'R12': build_R12}

PART = 9_000_000


def write_table(item, d):
    """Write a table. A table with a graph of more than PART rotation systems
    becomes a folder: the table itself, and the rows of each such graph in parts
    of PART rows, so that no file is larger than about 25 MB."""
    big = [g for g in d['graphs'] if len(g['rows']) > PART]
    if not big:
        return [write('rotations/%s.json' % item, d)]
    out = []
    for g in big:
        rows = g['rows']
        names = []
        for k, start in enumerate(range(0, len(rows), PART)):
            name = 'rows_%d.json' % (k + 1)
            res = write('rotations/%s/%s' % (item, name), {
                'item': item, 'graph': g['name'], 'part': k + 1, 'first rotation system': start,
                'rows': rows[start:start + PART]})
            names.append(os.path.basename(res[0]))
            out.append(res)
        g['rows'] = {'stored in': names, 'rows per part': PART,
                     'meaning': 'The rows of this graph are stored in these files of this folder, in order; part k '
                                'starts with rotation system number (k - 1) * %d.' % PART}
    out.append(write('rotations/%s/%s.json' % (item, item), d))
    return out


def main():
    global OUT
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    OUT = os.path.abspath(sys.argv[1])
    import time
    for item in (sys.argv[2:] or list(BUILDERS)):
        t0 = time.time()
        d = BUILDERS[item]()
        res = write_table(item, d)
        print(item, ', '.join(os.path.relpath(r[0], OUT) for r in res), '%.0f s' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
