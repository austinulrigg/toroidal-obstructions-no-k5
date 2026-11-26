"""Write an enumeration in the layout of the data folder.

A generator writes its records to a working folder, with its own names for
some vertices (the twins of Psi_2 are p0, p1, r0, r1 there, for instance).
The functions here translate the records into the names used in the paper and
in the data folder, attach the place of use and the quoted passage of the
paper, and write the files input.json, candidates.json, verdicts.json and
summary.json of the item.

A file is written with an indentation of one space when its compact form has
fewer than 2,000,000 characters and compactly otherwise; a file of more than
25 MB is written gzip-compressed, with the suffix .json.gz.
"""
import gzip
import json
import os
from collections import Counter

import common

GZ_LIMIT = 25 * 1024 * 1024

# Where the paper uses each enumeration.
WHERE = {
    'E1': 'Lemma 4.10', 'E2': 'Proposition 4.21', 'E3': 'Theorem 4.22',
    'E4': 'Proposition 4.24', 'E5': 'Proposition 4.25', 'E6': 'Proposition 4.28',
    'E7': 'Lemma 4.29 and Proposition 4.30', 'E8': 'Proposition 5.9', 'E9': 'Lemma 5.8',
    'E10': 'Proposition 5.16', 'E11': 'Lemma 5.17', 'E12': 'Proposition 5.24',
    'E13': 'Proposition 5.24, the case of two components', 'E14': 'Proposition 5.27 (earlier proof)', 'E15': 'Proposition 5.27',
    'E16': 'Lemma 5.26', 'E17': 'Lemma 5.28', 'E18': 'Lemma 5.30', 'E19': 'Proposition 5.33',
    'E20': 'Proposition 6.8', 'E21': 'Lemma 6.11', 'E22': 'Proposition 6.13', 'E23': 'Lemma 6.14 (earlier proof)',
    'E24': 'Lemma 6.15', 'E25': 'Proposition 6.16', 'E26': 'Lemma 6.18', 'E27': 'Lemma 6.19',
    'E28': 'Lemma 6.19', 'E29': 'Proposition 6.25',
}


# ---------------------------------------------------------------- files

def write(folder, name, obj):
    """Write obj as folder/name (name ends in .json), following the rule of the
    module docstring. Returns the path written."""
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, name)
    text = json.dumps(obj, separators=(',', ':'), ensure_ascii=False)
    if len(text) < 2_000_000:
        text = json.dumps(obj, indent=1, ensure_ascii=False)
    data = (text + '\n').encode()
    for stale in (path, path + '.gz'):
        if os.path.exists(stale):
            os.remove(stale)
    if len(data) > GZ_LIMIT:
        path += '.gz'
        with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', compresslevel=9,
                                                    fileobj=raw, mtime=0) as f:
            f.write(data)
    else:
        with open(path, 'wb') as f:
            f.write(data)
    return path


def read(path):
    """A JSON file, or the .json.gz file of the same name."""
    if not os.path.exists(path) and os.path.exists(path + '.gz'):
        path += '.gz'
    with (gzip.open(path, 'rt', encoding='utf-8') if path.endswith('.gz') else open(path, encoding='utf-8')) as f:
        return json.load(f)


def write_item(out, label, inp, cands, verdicts, summary):
    folder = os.path.join(out, 'enumerations', label)
    write(folder, 'input.json', inp)
    write(folder, 'candidates.json', {'item': label, 'candidates': cands})
    write(folder, 'verdicts.json', {'item': label, 'verdicts': verdicts})
    write(folder, 'summary.json', summary)


def head(label):
    """The first entries of input.json: the item, where the paper uses it and
    the passages of the paper that cite it."""
    return {'item': label, 'where': WHERE[label], 'paper text': common.paper_text(label)}


# ---------------------------------------------------------------- graphs

def adj_from(vertices, edges):
    adj = {v: set() for v in vertices}
    for u, v in edges:
        assert u != v and u in adj and v in adj, ('bad edge', u, v)
        adj[u].add(v)
        adj[v].add(u)
    return adj


def to_graph6(adj, order):
    """graph6 string of adj with the vertices numbered in the given order."""
    order = list(order)
    n = len(order)
    assert n < 63
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if order[i] in adj[order[j]] else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [n + 63]
    for k in range(0, len(bits), 6):
        x = 0
        for b in bits[k:k + 6]:
            x = (x << 1) | b
        out.append(x + 63)
    return ''.join(chr(c) for c in out)


def components(adj, vertices=None):
    vs = set(adj) if vertices is None else set(vertices)
    seen, comps = set(), []
    for s in vs:
        if s in seen:
            continue
        comp, stack = {s}, [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y in vs and y not in seen:
                    seen.add(y)
                    comp.add(y)
                    stack.append(y)
        comps.append(comp)
    return comps


class Relabel:
    """Renames vertices; names not in the map are kept."""

    def __init__(self, m):
        self.m = dict(m)

    def v(self, x):
        return self.m.get(x, x)

    def vs(self, xs):
        return [self.v(x) for x in xs]

    def es(self, es):
        return [[self.v(u), self.v(w)] for u, w in es]

    def graph(self, g):
        vs = self.vs(g['vertices'])
        assert len(set(vs)) == len(vs), ('two vertices get the same name', self.m, g['vertices'])
        return {'vertices': vs, 'edges': self.es(g['edges'])}


def rel_spec(x, R):
    """Rename the vertices named inside a specification (lists, values of
    dictionaries and single names; keys are kept)."""
    if isinstance(x, list):
        return [rel_spec(y, R) for y in x]
    if isinstance(x, dict):
        return {k: rel_spec(v, R) for k, v in x.items()}
    if isinstance(x, str):
        return R.v(x)
    return x


# ---------------------------------------------------------------- verdicts

def merge_embedding(parts, R):
    """One rotation system from the embeddings of the components."""
    rot, g = {}, 0
    for p in parts:
        for k, nb in p['rotation'].items():
            kk = R.v(k)
            assert kk not in rot
            rot[kk] = R.vs(nb)
        g += p['genus']
    return {'type': 'embedding', 'genus': g, 'rotation': rot}


def convert_outcome(o, R):
    """The verdict and witness of the data folder for an outcome recorded by
    enumeration.classify."""
    k = o['outcome']
    if k == 'K5 minor':
        return 'K5 minor', {'type': 'K5 branch sets', 'branch sets': [R.vs(s) for s in o['branch_sets']]}
    if k == 'cut of order at most two':
        return 'cut of at most two vertices', {'type': 'cut', 'vertices': R.vs(o['cut'])}
    if k == 'catalogue graph':
        return 'isomorphic to ' + o['catalogue'], {'type': 'isomorphism', 'onto': o['catalogue'],
                                                   'map': {R.v(a): b for a, b in o['map'].items()}}
    if k == 'toroidal':
        return 'toroidal', merge_embedding(o['embedding'], R)
    if k == 'proper nontoroidal minor':
        return 'proper nontoroidal minor', {'type': 'proper minor', 'operation': o['operation'],
                                            'edge': R.vs(o['edge'])}
    raise AssertionError('outcome %s: the enumeration is not complete' % k)


_CATALOGUE = {}


def catalogue_graph(name):
    """T_name of item O1, with vertices 0, 1, ... in the order of its graph6."""
    import networkx as nx
    if not _CATALOGUE:
        for r in common.catalogue_records():
            _CATALOGUE[r['name']] = nx.from_graph6_bytes(r['graph6'].encode())
    return _CATALOGUE[name]


def canonical_map(g, onto):
    """The isomorphism of the graph g (vertices and edges) onto T_onto whose
    list of images of the vertices of g, taken in the order of their names,
    comes first; the map lists the vertices in that order. A search may find
    any of the isomorphisms, in an order that can change from run to run;
    this choice does not depend on the search."""
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(g['vertices'])
    G.add_edges_from(tuple(e) for e in g['edges'])
    order = sorted(g['vertices'])
    best = None
    for iso in nx.algorithms.isomorphism.GraphMatcher(G, catalogue_graph(onto)).isomorphisms_iter():
        key = [iso[v] for v in order]
        if best is None or key < best:
            best = key
    assert best is not None, 'not isomorphic to ' + onto
    return dict(zip(order, best))


def verdict_counts(verdicts):
    return dict(sorted(Counter(v['verdict'] for v in verdicts).items()))


# The facts used to exclude a graph by a hypothesis of its case (Sections 5 and 6).
FACT = {
    'disjoint twin triples': 'two twin triples are disjoint',
    'two twin triples share two vertices': 'two twin triples share two vertices',
    'three twin triples share a pair': 'three twin triples share a pair of vertices',
    'component adjacent to all four attachments': 'a choice of twin pairs leaves a component adjacent to all four of their attachments',
    'component adjacent to all five attachments': 'a choice of twin pairs leaves a component adjacent to all five of their attachments',
    'twin triples meet other than once': 'two twin triples meet in zero or two vertices',
    'earlier row': 'a choice of twin pairs leaves a component with at least four neighbours',
    'adjacent neighbors': 'a vertex of degree three has two adjacent neighbours',
}

NONPLANAR = {'P nonplanar': 'B is nonplanar, where B is G minus p, q, r, z and the edges among a, b, c, d (Section 5.4 assumes B planar)'}


def convert_fact(se, R, conditions):
    key = se.get('condition') or se.get('kind')
    out = {}
    if se.get('kind') == 'nonplanar subgraph':
        spec = conditions[se['condition']]
        out['fact'] = NONPLANAR[se['condition']]
        out['subgraph'] = {k: (R.vs(v) if k in ('remove vertices', 'remove edges within', 'keep vertices') else
                               (R.es(v) if k == 'add edges' else v))
                           for k, v in spec.items() if k not in ('kind', 'meaning')}
        out['Kuratowski subgraph'] = R.es(se['kuratowski'])
        return out
    out['fact'] = FACT[key]
    for k in ('triples', 'twins'):
        if k in se:
            out[k] = [R.vs(t) for t in se[k]]
    if 'component' in se:
        out['component'] = R.vs(se['component'])
    if 'neighbors' in se:
        out['neighbours'] = R.vs(se['neighbors'])
    if 'vertex' in se:
        out['vertex'] = R.v(se['vertex'])
    return out


def convert_bridge_tests(tests, R):
    out = []
    for t in tests:
        c = t['certificate']
        if t['answer'] == 'planar':
            cert = {'type': 'embedding', 'genus': 0, 'rotation': {R.v(k): R.vs(v) for k, v in c.items()}}
        else:
            cert = {'type': 'Kuratowski subgraph', 'edges': R.es(c)}
        out.append({'completion': t['name'], 'added edges': R.es(t['added edges']), 'answer': t['answer'],
                    'witness': cert})
    return out


# ---------------------------------------------------------------- enumerations of whole graphs

S6 = {'p1': 'p', 'p2': 'q', 'q1': 'r', 'q2': 'z', 'q3': 'r3', 'q4': 'r4'}
S5a = {'p0': 'p', 'p1': 'q', 'p2': 'p3', 'r0': 'r', 'r1': 'z', 'r2': 'r3'}
S5b = {'p': "p'", 'q': "q'", 't0': 'p', 't1': 'q', 't2': 'r', 't3': 'z'}
S5c = {'p': "p'", 'q': "q'", 'x1': 'p', 'x2': 'q', 'y1': 'r', 'y2': 'z'}
S5d = {'abc0': 'p', 'abc1': 'q', 'abd0': 'r', 'abd1': 'z'}
S5e = {'x1': 'p', 'x2': 'q', 'y1': 'r', 'y2': 'z'}

PSI6 = [('pq', 'abc'), ('rz', 'ade')]
PSI5 = [('pq', 'abc'), ('rz', 'abd')]
PSI5w = [('pq', 'abc'), ('rz', 'abd'), (["p'", "q'"], 'acw')]
PSI4 = [('pq', 'ace'), ('rz', 'bdf')]

NAMES6 = ('Psi_3: twins p, q with neighbourhood {a,b,c} and r, z with neighbourhood {a,d,e}, as in Section 6. '
          'Further vertices with neighbourhood {a,b,c} are p3, p4, and with {a,d,e} r3, r4. '
          'Interior vertices of bridges are u0, u1, ...')
NAMES5 = ('Psi_2: twins p, q with neighbourhood {a,b,c} and r, z with neighbourhood {a,b,d}, as in Section 5. '
          'Further vertices with these neighbourhoods are p3, ... and r3, .... Interior vertices of bridges are u0, u1, ...')
NAMES5w = NAMES5 + " Twins p', q' have neighbourhood {a,c,w}, as in Section 5.4.3."
NAMES4 = ('Section 4.7, type (c): a = s1, b = t1, c = s2, d = t3, e = s3, f = t2; '
          'twins p, q with neighbourhood {a,c,e} and r, z with {b,d,f}.')

# label: (names used by the generator -> names of the paper, the twins, which
# records to keep, the text on vertex names)
WHOLE = {
    'E5': ({}, PSI4, None, NAMES4),
    'E8': (S5a, PSI5, None, NAMES5),
    'E9': (S5a, PSI5, None, NAMES5),
    'E14': (S5b, PSI5w, None, NAMES5w),
    'E15': (S5c, PSI5w, None, NAMES5w),
    'E16': (S5c, PSI5w, None, NAMES5w),
    'E17': (S5d, PSI5, None, NAMES5 + ' Here the other vertices with a given neighbourhood are named by it, '
                                      'as abc2, abd2, acd0, bcd0, ....'),
    'E19': (S5e, PSI5, None, NAMES5 + ' Vertices of the k-th form are f<k>_<i>.'),
    'E20': (S6, PSI6, None, NAMES6),
    'E21': (S6, PSI6, None, NAMES6),
    'E22': (S6, PSI6, None, NAMES6),
    'E23': (S6, PSI6, lambda c: c['case'] == 'C_on_abd', NAMES6),
    'E24': (S6, PSI6, None, NAMES6),
    'E25': (S6, PSI6, None, NAMES6),
    'E26': (S6, PSI6, None, NAMES6),
    'E27': (S6, PSI6, None, NAMES6),
    'E28': (S6, PSI6, None, NAMES6),
    'E29': (S6, PSI6, None, NAMES6 + ' Vertices of the forms are B<k>_<x>; added single vertices are v_<neighbourhood>.'),
}


def raw_folders(work):
    """The folders holding the records of a run: the run folder itself, or,
    for a run written in parts, its subfolders in order of name."""
    if os.path.exists(os.path.join(work, 'candidates.json')):
        return [work]
    subs = sorted(d for d in os.listdir(work) if os.path.isdir(os.path.join(work, d))
                  and (os.path.exists(os.path.join(work, d, 'candidates.json'))
                       or os.path.exists(os.path.join(work, d, 'candidates.json.gz'))))
    assert subs, 'no records in %s' % work
    return [os.path.join(work, d) for d in subs]


def check_psi(adj, psi):
    for tw, nb in psi:
        for t in tw:
            assert t in adj and adj[t] == set(nb), ('a vertex of Psi has the wrong neighbours', t, nb)


def package_whole(label, work, out):
    """Package an enumeration of whole graphs classified with
    enumeration.classify (or the variants of psi2_runs and psi3_runs)."""
    Rmap, psi, keep, names = WHOLE[label]
    R = Relabel(Rmap)
    conditions, lattice, cands, verdicts = {}, None, [], []
    also = Counter()
    for f in raw_folders(work):
        d = read(os.path.join(f, 'candidates.json'))
        conditions.update(d.get('structural conditions', {}))
        if 'lattice' in d:
            lattice = d['lattice']
        for c in d['candidates']:
            if keep and not keep(c):
                continue
            g = R.graph(c['graph'])
            adj = adj_from(g['vertices'], g['edges'])
            check_psi(adj, psi)
            cand = {'id': c['id']}
            for k in ('case', 'variant'):
                if c.get(k) not in (None, ''):
                    cand[k] = c[k]
            cand['graph'] = g
            cand['graph6'] = to_graph6(adj, g['vertices'])
            if c.get('part') and c['part'].get('vertices'):
                cand['bridge'] = R.graph(c['part'])
            for k in ('form', 'bd added', 'config', 'base'):
                if k in c:
                    cand[k] = c[k]
            cands.append(cand)
            verdict, w = convert_outcome(c['outcome'], R)
            if w['type'] == 'isomorphism':
                w['map'] = canonical_map(g, w['onto'])
            vrec = {'id': c['id'], 'verdict': verdict, 'witness': w}
            se = c.get('structural exclusion')
            if se:
                fact = convert_fact(se, R, conditions)
                vrec['also excluded by a hypothesis of the case'] = fact
                also[fact['fact']] += 1
            if c.get('bridge tests'):
                vrec['planarity of the bridge completions'] = convert_bridge_tests(c['bridge tests'], R)
            for k in ('lattice status', 'case flags'):
                if k in c:
                    vrec[k] = c[k]
            verdicts.append(vrec)
    ids = [c['id'] for c in cands]
    assert len(ids) == len(set(ids)), 'two candidates with the same id'
    cases = {k: rel_spec(v, R) for k, v in read(os.path.join(work, 'summary.json')).get('cases', {}).items()}

    def last(v):
        return v['verdict'].split()[-1]
    by_num = lambda s: int(s[1:])
    isos = Counter(last(v) for v in verdicts if v['verdict'].startswith('isomorphic to'))
    remaining = sorted({last(v) for v in verdicts if v['verdict'].startswith('isomorphic to')
                        and 'also excluded by a hypothesis of the case' not in v}, key=by_num)
    excluded = sorted({last(v) for v in verdicts if v['verdict'].startswith('isomorphic to')
                       and 'also excluded by a hypothesis of the case' in v}, key=by_num)
    inp = head(label)
    inp['vertex names'] = names
    inp['cases of the construction'] = cases
    if conditions:
        inp['hypotheses of the case used for exclusions'] = {
            k: rel_spec({kk: vv for kk, vv in v.items() if kk != 'kind'}, R) for k, v in conditions.items()}
    if lattice is not None:
        inp['lattice'] = rel_spec(lattice, R)
        # The auxiliary rotations and isomorphisms also use vertices as dictionary keys.
        for raw, named in zip(lattice['auxiliary'], inp['lattice']['auxiliary']):
            vertex_names = {v: R.v(v) for v in raw['graph']['vertices']}
            named['outcome'] = remap_outcome(raw['outcome'], vertex_names)
    summary = {'item': label, 'where': WHERE[label], 'candidates': len(cands), 'verdicts': verdict_counts(verdicts),
               'isomorphic to T_i (copies)': dict(sorted(isos.items(), key=lambda kv: by_num(kv[0]))),
               'graphs T_i not excluded by a hypothesis of the case': remaining}
    if also:
        summary['also excluded by a hypothesis of the case'] = dict(sorted(also.items()))
    if excluded:
        summary['graphs T_i excluded by a hypothesis of the case'] = excluded
    case_of = {c['id']: c.get('case') for c in cands}
    by_case = {}
    for v in verdicts:
        k = case_of[v['id']]
        if k is None:
            continue
        e = by_case.setdefault(k, {'candidates': 0, 'left': set()})
        e['candidates'] += 1
        if v['verdict'].startswith('isomorphic to') and 'also excluded by a hypothesis of the case' not in v:
            e['left'].add(last(v))
    if len(by_case) > 1:
        summary['by case'] = {k: {'candidates': e['candidates'],
                                  'graphs T_i not excluded by a hypothesis of the case': sorted(e['left'], key=by_num)}
                              for k, e in sorted(by_case.items())}
    write_item(out, label, inp, cands, verdicts, summary)
    return summary


# ---------------------------------------------------------------- E6, E7: bridges with tests on the bridge

TESTS_TEXT = {
    'E6': {'roots': 'the attachments s_x, t_x, s_y, t_y of B, written sx, tx, sy, ty',
           'B + C': 'B + s_x t_x + t_x t_y + t_y s_y + s_y s_x, the cycle completion B ∪ C',
           'B + path': 'B + s_y s_x + s_x t_x + t_x t_y, the path completion B ∪ s_y s_x t_x t_y',
           'B + K4[R]': 'B plus all six edges among s_x, t_x, s_y, t_y',
           'B + Rx + Ry': 'B + s_x t_x + s_y t_y, a subgraph of G once P_x, P_y are edges'},
    'E7': {'roots': 'the attachments a = s_m, b = t_m, c = s_x, d = t_x of B',
           'U + ac + cd + db': 'B + ac + cd + db',
           'U + ac + bd': 'B + ac + bd',
           'U + K4[abcd]': 'B plus all six edges among a, b, c, d',
           'U + cd': 'B + cd, a subgraph of G since P_x = cd',
           'Psi1 + U': 'the graph consisting of Psi_1 and B, with P_x = cd and P_w shortened to ef: '
                       'twins p, q on {a, c, e} and r, z on {b, d, f}, where e = s_w, f = t_w'},
}

# The test graphs under the names of the data folder.
RENAME = {
    'E6': {'B + C': 'B + C', 'B + path': 'B + sy sx + sx tx + tx ty', 'B + K4[R]': 'B + K4[sx,tx,sy,ty]',
           'B + Rx + Ry': 'B + sx tx + sy ty'},
    'E7': {'U + ac + cd + db': 'B + ac + cd + db', 'U + ac + bd': 'B + ac + bd', 'U + K4[abcd]': 'B + K4[a,b,c,d]',
           'U + cd': 'B + cd', 'Psi1 + U': 'Psi1 + B'},
}

# The bridges given in the paper. E6: twins p', q' on {sx, sy, v}, with v
# also adjacent to tx, ty. E7: the bridges B1, ..., B5 of Lemma 4.29.
PAPER_FORMS = {
    'E6': {'the form of Proposition 4.28': [("p'", 'sx'), ("p'", 'sy'), ("p'", 'v'), ("q'", 'sx'), ("q'", 'sy'),
                                           ("q'", 'v'), ('v', 'tx'), ('v', 'ty')]},
    'E7': {'B1': [("p'", 'a'), ("p'", 'c'), ("p'", 'v'), ("q'", 'a'), ("q'", 'c'), ("q'", 'v'),
                  ('a', 'u'), ('b', 'u'), ('d', 'u'), ('u', 'v')],
           'B2': [("p'", 'a'), ("p'", 'c'), ("p'", 'v'), ("q'", 'a'), ("q'", 'c'), ("q'", 'v'),
                  ('a', 'u'), ('b', 'u'), ('d', 'v'), ('u', 'v')],
           'B3': [("p'", 'a'), ("p'", 'u'), ("p'", 'v'), ("q'", 'a'), ("q'", 'u'), ("q'", 'v'),
                  ('c', 'u'), ('b', 'v'), ('d', 'v')],
           'B4': [("p'", 'd'), ("p'", 'u'), ("p'", 'v'), ("q'", 'd'), ("q'", 'u'), ("q'", 'v'),
                  ('c', 'u'), ('a', 'v'), ('b', 'v')],
           'B5': [("p'", 'a'), ("p'", 'b'), ("p'", 'v'), ("q'", 'a'), ("q'", 'b'), ("q'", 'v'),
                  ('c', 'v'), ('d', 'v')]},
}
# Interchanging A1 and A2.
SWAP = {'E6': {'sx': 'tx', 'tx': 'sx', 'sy': 'ty', 'ty': 'sy'},
        'E7': {'a': 'b', 'b': 'a', 'c': 'd', 'd': 'c'}}

FACTS = {
    'path completion nonplanar': (
        'the path completion B + s_y s_x + s_x t_x + t_x t_y is nonplanar',
        'The labels are chosen so that the path completion B ∪ s_y s_x t_x t_y is planar.'),
    'disjoint sx-ty and sy-tx paths': (
        'B has disjoint s_x-t_y and s_y-t_x paths',
        'Lemma 4.7 excludes vertex-disjoint s_x–t_y and s_y–t_x paths in B.'),
    'degree-three vertex with adjacent neighbors (Lemma 3.4(a))': (
        'an interior vertex of degree three has two adjacent neighbours',
        'Lemma 3.4(a); the edges added in the test graph are edges of G.'),
    'B + K4[R] has a cut of order at most two': (
        'B plus all edges among its four attachments has a cut of at most two vertices',
        'Three-connectivity of B ∪ K_4[attachments], used with Construction 3.18.'),
    'U + ac + bd nonplanar': (
        'B + ac + bd is nonplanar',
        'Lemma 4.29 shows that B + ac + bd is planar.'),
    'U + K4[abcd] has a cut of order at most two': (
        'B plus all edges among a, b, c, d has a cut of at most two vertices',
        'Three-connectivity of B ∪ K_4[{a,b,c,d}], used with Construction 3.18.'),
}

WHOLE_GRAPHS_E7 = ('For each bridge that passes, the two graphs of Lemma 4.29(i): G = (Psi_1 ∪ B) + eb and '
                   'G = (Psi_1 ∪ B) + af, with P_x = cd and P_w = ef. Proposition 4.30 cites this part for the '
                   'claim that these give only T66. The exterior formed by a bridge uv on {c,d,e,f} is treated '
                   'in the paper by argument, so it is not included here.')


def allowed_reduction_ok(adj, roots, v, x):
    """Whether deleting the edge vx of the bridge adj is an allowed operation:
    the interior stays nonempty and connected and every attachment keeps an
    interior neighbour."""
    H = {y: set(s) for y, s in adj.items()}
    H[v].discard(x)
    H[x].discard(v)
    inner = [y for y in H if y not in roots]
    return bool(inner) and len(components(H, inner)) == 1 and \
        all(any(u not in roots for u in H[r]) for r in roots)


def paper_form(label, roots, bridge_edges):
    """The bridge given in the paper that a passing bridge is: a bijection of
    the interior vertices, after possibly interchanging A1 and A2, carrying
    its edges onto the edges of a form of PAPER_FORMS."""
    from itertools import permutations
    E = {frozenset(e) for e in bridge_edges}
    inner = sorted({x for e in bridge_edges for x in e} - set(roots))
    for name, fe in PAPER_FORMS[label].items():
        F = {frozenset(e) for e in fe}
        finner = sorted({x for e in fe for x in e} - set(roots))
        if len(finner) != len(inner) or len(F) != len(E):
            continue
        for swap in (False, True):
            sw = SWAP[label] if swap else {}
            for perm in permutations(finner):
                m = dict(zip(inner, perm))
                img = {frozenset(sw.get(m.get(x, x), m.get(x, x)) if x in roots else m[x] for x in e) for e in E}
                if img == F:
                    return {'form': name, 'A1 and A2 interchanged': swap, 'names of the interior vertices': m}
    raise AssertionError('a bridge that passes is not a bridge given in the paper')


def package_bridges(label, work, out):
    """Package E6 or E7 (the bridges B with the tests of the case)."""
    d = read(os.path.join(work, 'candidates.json'))
    roots = d['roots']
    RN = RENAME[label]
    tests = {RN[k]: [list(e) for e in v] for k, v in d['test graphs'].items()}
    meaning = {RN[k]: TESTS_TEXT[label][k] for k in d['test graphs']}
    kinds = {}
    for k, v in d['outcome tests'].items():
        v = list(v)
        for i in range(1, len(v)):
            if isinstance(v[i], str):
                v[i] = RN[v[i]]
        kinds[k] = v
    cands, verdicts = [], []
    whole_counts, found, forms = Counter(), Counter(), Counter()
    for c in d['candidates']:
        B = adj_from(sorted(set(c['part']['vertices']) | set(roots)), c['part']['edges'])
        cands.append({'id': c['id'], 'bridge': {'vertices': c['part']['vertices'], 'edges': c['part']['edges']}})
        o = c['outcome']
        kind = o['outcome']
        how = kinds[kind]
        rec = {'id': c['id']}
        if how[0] == 'nonplanar':
            rec['verdict'] = 'excluded by a hypothesis of the case'
            rec['witness'] = {'fact': FACTS[kind][0], 'test graph': how[1], 'Kuratowski subgraph': o['kuratowski']}
        elif how[0] == 'K5 minor':
            rec['verdict'] = 'K5 minor'
            rec['witness'] = {'type': 'K5 branch sets', 'in': how[1], 'branch sets': o['branch_sets']}
        elif how[0] == 'cut':
            rec['verdict'] = 'excluded by a hypothesis of the case'
            rec['witness'] = {'fact': FACTS[kind][0], 'test graph': how[1], 'cut': o['cut']}
        elif how[0] == 'degree three':
            rec['verdict'] = 'excluded by a hypothesis of the case'
            rec['witness'] = {'fact': FACTS[kind][0], 'in': how[1], 'vertex': o['vertex'],
                              'neighbours': o['adjacent neighbors']}
        elif how[0] == 'disjoint paths':
            rec['verdict'] = 'excluded by a hypothesis of the case'
            rec['witness'] = {'fact': FACTS[kind][0], 'paths': o['paths']}
        elif how[0] == 'allowed minor nonplanar':
            rec['verdict'] = 'proper nontoroidal minor'
            w = {'type': 'allowed reduction', 'operation': o['operation'], 'arguments': o['arguments'],
                 'failing test after the reduction': how[1], 'Kuratowski subgraph': o['kuratowski']}
            if o['operation'] == 'delete vertex':
                # Record the deletion of an edge at the vertex instead: the first
                # edge, attachments first, whose deletion is an allowed operation.
                (v,) = o['arguments']
                assert v not in {x for e in o['kuratowski'] for x in e}
                x = next(x for x in sorted(B[v], key=lambda y: (y not in roots, y))
                         if allowed_reduction_ok(B, set(roots), v, x))
                w['operation'] = 'delete edge'
                w['arguments'] = [v, x]
                w['note'] = ('Deleting the edge %s%s is an allowed operation whose result contains B - %s, and the '
                             'Kuratowski subgraph lies in B - %s.' % (v, x, v, v))
            rec['witness'] = w
        elif how[0] == 'survivor':
            passes = {'%s is nonplanar' % how[1]: {'type': 'Kuratowski subgraph', 'edges': o['B + C kuratowski']},
                      '%s is planar' % how[2]: {'type': 'embedding', 'genus': 0,
                                                'rotation': o['path completion rotation']},
                      'every allowed one-step reduction makes %s planar' % how[1]: [
                          {'operation': m['operation'], 'arguments': m['arguments'],
                           'planar embedding of the test after it': m['rotation']} for m in o['allowed minors']]}
            if c.get('edge deletions keeping B + C nonplanar'):
                passes['edge deletions keeping %s nonplanar' % how[1]] = [
                    {'deleted edge': x['deleted edge'], 'Kuratowski subgraph': x['kuratowski']}
                    for x in c['edge deletions keeping B + C nonplanar']]
            pf = paper_form(label, roots, c['part']['edges'])
            sw = SWAP[label] if pf['A1 and A2 interchanged'] else {}
            m = {x: (sw.get(x, x) if x in roots else pf['names of the interior vertices'][x]) for x in B}
            rec['verdict'] = 'isomorphic to ' + pf['form']
            rec['in the paper'] = pf
            forms[pf['form']] += 1
            wg = []
            for x in c.get('whole graphs', []):
                g = x['graph']
                adj = adj_from(g['vertices'], g['edges'])
                verdict, wit = convert_outcome(x['outcome'], Relabel({}))
                if wit['type'] == 'isomorphism':
                    wit['map'] = canonical_map(g, wit['onto'])
                wg.append({'exterior': x['exterior'], 'graph': g, 'graph6': to_graph6(adj, g['vertices']),
                           'verdict': verdict, 'witness': wit})
                whole_counts[verdict] += 1
                if verdict.startswith('isomorphic'):
                    found[verdict.split()[-1]] += 1
            if wg:
                rec['whole graphs'] = wg
            rec['witness'] = {'type': 'isomorphism', 'onto': pf['form'], 'map': dict(sorted(m.items()))}
            rec['also: B passes every test of the case'] = passes
        else:
            raise AssertionError(how)
        verdicts.append(rec)
    inp = head(label)
    inp.update({
        'attachments': roots,
        'test graphs': {k: {'added edges': v, 'meaning': meaning[k]} for k, v in tests.items()},
        'vertex names': TESTS_TEXT[label]['roots'] + '; interior vertices of B are u0, u1, ...; '
                        'for the bridges that pass, the paper\'s names are given with the verdict',
        'verdicts': {
            'excluded by a hypothesis of the case': [{'fact': v[0], 'why it excludes the bridge': v[1]}
                                                     for k, v in FACTS.items() if k in kinds],
            'K5 minor': 'a K5 minor in the named graph, which is a minor of G',
            'proper nontoroidal minor': 'an allowed reduction of B (Section 3) after which the named test graph '
                                        'is still nonplanar; the proof shows that this gives a proper '
                                        'nontoroidal minor of G',
            'isomorphic to a bridge given in the text': (
                'B passes every test of the case and is isomorphic to a bridge given in the text (listed under '
                '"bridges given in the text"); the record that it passes the tests is kept with the verdict')}})
    if label == 'E7':
        inp['whole graphs'] = WHOLE_GRAPHS_E7
    s = read(os.path.join(work, 'summary.json'))
    for k in ('specification', 'generation'):
        if k in s:
            inp[k] = s[k]
    inp['bridges given in the text'] = {
        form: {'vertices': sorted({x for e in es for x in e}), 'edges': sorted(sorted(e) for e in es)}
        for form, es in PAPER_FORMS[label].items()}
    summary = {'item': label, 'where': WHERE[label], 'bridges': len(cands), 'verdicts': verdict_counts(verdicts),
               'bridges that pass, by form in the paper (labelled, so a form and its interchange count separately)':
                   dict(sorted(forms.items()))}
    if whole_counts:
        summary['whole graphs from the bridges that pass'] = dict(whole_counts)
        summary['isomorphic to T_i (copies)'] = dict(found)
    write_item(out, label, inp, cands, verdicts, summary)
    return summary


# ---------------------------------------------------------------- E18

def planarity_witness(vertices, edges):
    """A plane rotation system or a Kuratowski subgraph of the graph, found
    with networkx (vertices and edges taken in sorted order)."""
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(sorted(vertices))
    G.add_edges_from(sorted(tuple(sorted(e)) for e in edges))
    ok, P = nx.check_planarity(G, counterexample=True)
    if ok:
        return True, {'type': 'embedding', 'genus': 0, 'rotation': {v: list(P.neighbors_cw_order(v)) for v in G}}
    return False, {'type': 'Kuratowski subgraph', 'edges': sorted(sorted(e) for e in P.edges())}


def package_E18(work, out):
    R = Relabel({'x1': 'p', 'x2': 'q', 'y1': 'r', 'y2': 'z'})
    d = read(os.path.join(work, 'candidates.json'))
    s = read(os.path.join(work, 'summary.json'))
    forms = read(os.path.join(work, 'forms.json'))
    cands, verdicts, by_case = [], [], {}
    n_kind = Counter()
    for c in d['candidates']:
        g = R.graph(c['graph'])
        adj = adj_from(g['vertices'], g['edges'])
        _, T, kind, _, _, S = c['case'].split('_')
        check_psi(adj, PSI5)
        assert adj['v'] == set(S)
        inner = [x for x in g['vertices'] if x not in set('abcdpqrzv')]
        bedges = [e for e in g['edges'] if e[0] in inner or e[1] in inner]
        planar, pw = planarity_witness(list(T) + inner, bedges)
        assert planar == (kind == 'planar')
        cand = {'id': c['id'], 'case': c['case'], 'U': 'on {%s}, %s' % (','.join(T), kind), 'N(v)': list(S),
                'graph': g, 'graph6': to_graph6(adj, g['vertices']),
                'U is planar' if planar else 'U is nonplanar': pw}
        if planar:
            x, y, t = T
            comps = []
            for E in (((x, y), (x, t)), ((x, y), (y, t)), ((x, t), (y, t))):
                ok, kw = planarity_witness(list(T) + inner, bedges + [list(e) for e in E])
                assert not ok
                comps.append({'added edges': [list(e) for e in E], 'Kuratowski subgraph': kw['edges']})
            cand['the three two-edge completions of U are nonplanar'] = comps
            cand['form'] = int(c['id'].rsplit('-', 1)[1])
        cands.append(cand)
        n_kind[kind] += 1
        verdict, w = convert_outcome(c['outcome'], R)
        verdicts.append({'id': c['id'], 'verdict': verdict, 'witness': w})
        by_case.setdefault(c['case'], Counter())[verdict] += 1
    counts = verdict_counts(verdicts)
    assert not any(v.startswith('isomorphic') or v == 'toroidal' for v in counts), counts
    gen = s['cases']['U_acd_nonplanar_v_on_abc']
    inp = head('E18')
    inp.update({
        'vertex names': 'Psi_2: twins p, q with neighbourhood {a,b,c} and r, z with {a,b,d}; v is the contracted '
                        'other component. Interior vertices of U are u0, u1, ....',
        'cases': 'U on T = {a,c,d} or {b,c,d}, and N(v) each of the three other triples: six cases.',
        'nonplanar U': {'construction': 'Construction 3.18 with attachments T and no added edges',
                        'bridges per case': gen['parts up to rooted isomorphism'],
                        'generation': {k: v for k, v in gen.items() if k != 'outcomes'}},
        'planar U': {'construction': 'Construction 3.18 with attachments T and one added edge xy, required, for '
                                     'each pair xy of T; kept when planar with all three two-edge completions '
                                     'nonplanar; merged over the three pairs up to isomorphism fixing T',
                     'forms per triple': forms['acd']['count'], 'forms': forms}})
    summary = {'item': 'E18', 'where': WHERE['E18'], 'graphs': len(cands),
               'nonplanar U': n_kind['nonplanar'], 'planar U': n_kind['planar'], 'verdicts': counts,
               'by case': {k: dict(sorted(v.items())) for k, v in sorted(by_case.items())},
               'remaining': 0}
    write_item(out, 'E18', inp, cands, verdicts, summary)
    return summary


# ---------------------------------------------------------------- E29: the bases

def package_E29_bases(work, out):
    """bases.json of E29: every base graph (two vertices on abc and two on
    ade), discarded with the fact used, or kept."""
    R = Relabel(S6)
    b = read(os.path.join(work, 'bases.json'))['candidates']
    used = {c['base'] for c in read(os.path.join(work, 'candidates.json'))['candidates']}
    n_cands = len(read(os.path.join(work, 'candidates.json'))['candidates'])
    bases, cnt = [], Counter()
    for x in b:
        g = R.graph(x['graph'])
        adj = adj_from(g['vertices'], g['edges'])
        check_psi(adj, PSI6)
        rec = {'id': x['id'], 'description': x['description'], 'graph': g, 'graph6': to_graph6(adj, g['vertices'])}
        se = x.get('structural exclusion')
        if se:
            fact = convert_fact(se, R, {})
            rec['verdict'] = 'excluded by a hypothesis of the case'
            rec['witness'] = fact
            assert x['id'] not in used
            cnt[fact['fact']] += 1
        else:
            rec['verdict'] = 'kept'
            assert x['id'] in used
            cnt['kept'] += 1
        bases.append(rec)
    what = ('The {:,} graphs of Construction 6.24 with n_c = n_e = 2, up to isomorphism fixing a,...,e. '
            'Each is discarded, with the fact used, or kept; the {} kept ones give the {} candidates '
            '(field "base" of each candidate) after the five choices of (n_c, n_e).').format(
                len(bases), cnt['kept'], n_cands)
    return write(os.path.join(out, 'enumerations', 'E29'), 'bases.json',
                 {'item': 'E29', 'what': what, 'counts': dict(cnt), 'bases': bases})


# ---------------------------------------------------------------- E11

def package_E11(work, out):
    d = read(os.path.join(work, 'candidates.json'))
    tests = [['ac', 'bd'], ['ad', 'bc'], ['ab', 'ac', 'ad'], ['ab', 'bc', 'bd']]
    cands, verdicts, found = [], [], Counter()
    for c in d['candidates']:
        fam = c['recipes'][0]['family']
        m = {'p': "p'", 'q': "q'", 'p1': 'p', 'p2': 'q', 'q1': 'r', 'q2': 'z'}
        m.update({'z1': 'm'} if fam == 'star' else {'z1': 'm1', 'z2': 'm2'})
        R = Relabel(m)
        g = R.graph({'vertices': c['graph_vertices'], 'edges': c['graph_edges']})
        G = adj_from(g['vertices'], g['edges'])
        Bv = R.vs(c['vertices'])
        Be = R.es(c['part_edges'])
        for w, N in (('p', 'abc'), ('q', 'abc'), ('r', 'abd'), ('z', 'abd')):
            assert G[w] == set(N)
        assert G["p'"] == G["q'"] == {'a', 'u', 'v'}
        comp = []
        for t, cc in zip(tests, c['completions']):
            assert sorted(cc['added_edges']) == sorted(t) and cc['planar'] is False
            comp.append({'test': 'B + ' + ' + '.join(t), 'Kuratowski subgraph': R.es(cc['kuratowski_edges'])})
        rec = c['recipes'][0]
        joins = ', '.join(R.v(x) + R.v(y) for x, y in rec['joins'])
        if fam == 'star':
            form = 'star: m adjacent to b, c, d; joins ' + joins
        else:
            form = 'path %s m1 %s m2 %s; joins %s' % (rec['order'][0], rec['order'][1], rec['order'][2], joins)
        cands.append({'id': c['id'], 'graph': g, 'graph6': to_graph6(G, g['vertices']), 'form of D': form,
                      'bridge B': {'vertices': sorted(Bv), 'edges': sorted(sorted(e) for e in Be)},
                      'B is planar': {'type': 'embedding', 'genus': 0,
                                      'rotation': {R.v(k): R.vs(x) for k, x in c['planar_part']['rotation'].items()}},
                      'the four tests are nonplanar': comp})
        o = c['outcome']
        t = o['type']
        r = {'id': c['id']}
        if t == 'survivor':
            tn = o['target']
            r['verdict'] = 'isomorphic to ' + tn
            r['witness'] = {'type': 'isomorphism', 'onto': tn, 'map': canonical_map(g, tn)}
            found[tn] += 1
        elif t == 'k5_minor':
            r['verdict'] = 'K5 minor'
            r['witness'] = {'type': 'K5 branch sets', 'branch sets': [R.vs(s) for s in o['branch_sets']]}
        elif t == 'small_cut':
            r['verdict'] = 'cut of at most two vertices'
            r['witness'] = {'type': 'cut', 'vertices': R.vs(o['cut'])}
        elif t == 'proper_four_test_subgraph':
            de = R.es(o['deleted_edges'])
            assert len(de) == 1
            fails = []
            for t2, cc in zip(tests, o['completions']):
                assert sorted(cc['added_edges']) == sorted(t2) and cc['planar'] is False
                fails.append({'test': "B' + " + ' + '.join(t2), 'Kuratowski subgraph': R.es(cc['kuratowski_edges'])})
            r['verdict'] = 'proper nontoroidal minor'
            r['witness'] = {'type': 'proper minor', 'operation': 'delete', 'edge': de[0],
                            'why the minor is nontoroidal': {
                                'reason': "B' = B - %s still satisfies the hypotheses of Lemma 5.3 and all four "
                                          "tests for (n_c, n_d) = (2, 2) are nonplanar" % ''.join(de[0]),
                                'tests': fails}}
        else:
            raise AssertionError('outcome %s: the enumeration is not complete' % t)
        verdicts.append(r)
    inp = head('E11')
    inp.update({
        'vertex names': "Psi_2 with twins p, q on {a, b, c} and r, z on {a, b, d}, as in Section 5. The bridge B "
                        "has attachments a, b, c, d and contains u, v and the twins p', q' with neighbourhood "
                        "X = {a, u, v}. The vertices of M are m (star) or m1, m2 (path g1 m1 g2 m2 g3).",
        'the forms': 'One component D of B - X contains b, c, d, with N_B(D) = {u, v}. Its bipartite tree is a '
                     'star whose centre m is joined to b, c, d, or a path g1 m1 g2 m2 g3 with {g1, g2, g3} = '
                     '{b, c, d}. For the star, m meets both u, v, or one of them and a vertex of {b, c, d} meets the '
                     'other. For the path, m1, m2 meet different vertices of {u, v}, or both meet one and a vertex '
                     'of {b, c, d} meets the other. Up to symmetry this gives %d forms.' % len(cands),
        'the four tests': 'Lemma 5.3 with (n_c, n_d) = (2, 2): G is toroidal if and only if one of B + ac + bd, '
                          'B + ad + bc, B + ab + ac + ad, B + ab + bc + bd is planar.'})
    summary = {'item': 'E11', 'where': WHERE['E11'], 'forms': len(cands), 'verdicts': verdict_counts(verdicts),
               'isomorphic to T_i (copies)': dict(sorted(found.items(), key=lambda x: int(x[0][1:])))}
    write_item(out, 'E11', inp, cands, verdicts, summary)
    return summary


# ---------------------------------------------------------------- the joint test of Lemma 4.3 (E2, E3)

def convert_joint_test(jt, R=lambda v: v):
    """The record of joint_test.joint_test for a graph failing every choice,
    in the format of the data folder: the bridges B1, B2, ... and, for each
    of the twelve choices, the bridges with no face or every assignment with
    its nonplanar test graph and a Kuratowski subgraph. R renames vertices."""
    names = {}
    out_b = []
    for b in jt['bridges']:
        names[b['id']] = 'B%d' % (len(names) + 1)
        out_b.append({'name': names[b['id']], 'attachments': sorted(R(v) for v in b['attachments']),
                      'edges': [[R(u), R(v)] for u, v in b['edges']]})
    out_c = []
    for ch in jt['choices']:
        c = {'choice': [ch['k'], ch['x'], ch['y']]}
        assert ch['passes'] is False
        if ch['blocked_bridges']:
            c['bridges with no face'] = [names[n] for n in ch['blocked_bridges']]
        else:
            c['assignments'] = []
            for a in ch['assignments']:
                assert a['passes'] is False
                fl = a['failure']
                c['assignments'].append({'faces': {names[n]: f for n, f in a['mapping'].items()},
                                         'nonplanar test': fl['region'],
                                         'Kuratowski subgraph': [[R(u), R(v)] for u, v in fl['kuratowski_edges']]})
        out_c.append(c)
    return {'bridges': out_b, 'choices': out_c}


# ---------------------------------------------------------------- E1

def is_minor_model(adj, target, blocks):
    """Whether blocks (vertex of target -> list of vertices of adj) is a minor
    model of the graph target (adjacency sets) in adj."""
    allv = [v for b in blocks.values() for v in b]
    if set(blocks) != set(target) or len(allv) != len(set(allv)) or not set(allv) <= set(adj):
        return False
    if not all(b and len(components(adj, b)) == 1 for b in blocks.values()):
        return False
    return all(any(y in adj[x] for x in blocks[u] for y in blocks[v]) for u in target for v in target[u])


def catalogue_adj(name):
    T = catalogue_graph(name)
    return {v: set(T[v]) for v in T}


def one_edge_minor_with_model(G, tn, blocks):
    """An edge of G whose deletion keeps the minor model blocks of T_tn, or,
    failing that, an edge inside a branch set whose contraction keeps it.
    Returns (operation, edge, model, where)."""
    target = catalogue_adj(tn)
    for e in sorted(sorted(e) for e in [(u, v) for u in G for v in G[u] if u < v]):
        H = {x: set(s) for x, s in G.items()}
        H[e[0]].discard(e[1])
        H[e[1]].discard(e[0])
        if is_minor_model(H, target, blocks):
            return 'delete', e, blocks, 'G - %s%s' % tuple(e)
    for k, blk in sorted(blocks.items()):
        for u in blk:
            for v in sorted(G[u]):
                if v in blk and u < v:
                    H = {x: set(s) for x, s in G.items()}
                    for w in H.pop(v):
                        H[w].discard(v)
                        if w != u:
                            H[w].add(u)
                            H[u].add(w)
                    nb = {kk: [x for x in vv if x != v] for kk, vv in blocks.items()}
                    assert is_minor_model(H, target, nb)
                    return 'contract', [u, v], nb, 'G/%s%s (%s merged into %s)' % (u, v, v, u)
    raise AssertionError('no one-edge minor keeps the model')


def package_E1(work, out):
    d = read(os.path.join(work, 'candidates.json'))
    s = read(os.path.join(work, 'summary.json'))
    roots = ['s1', 't1', 's2', 't2']
    H2 = adj_from(['p', 'q', 'r', 'z', 's1', 's2', 's3', 't1', 't2', 't3'],
                  [[w, x] for w in 'pq' for x in ('s1', 's2', 's3')] +
                  [[w, t] for w in 'rz' for t in ('t1', 't2', 't3')] + [['s3', 't3']])
    cands, verdicts = [], []
    found, ops = Counter(), Counter()
    for c in d['candidates']:
        assert sorted(c['roots']) == sorted(roots)
        G = {x: set(nb) for x, nb in H2.items()}
        for x, nb in adj_from(c['vertices'], c['edges']).items():
            G.setdefault(x, set()).update(nb)
        kur = c['nonplanarity']['kuratowski_edges']
        cands.append({'id': c['id'], 'bridge': {'vertices': sorted(c['vertices']),
                                                'edges': sorted(sorted(e) for e in c['edges'])},
                      'B + s1s2 + t1t2 is nonplanar': {'type': 'Kuratowski subgraph', 'edges': kur}})
        o = c['outcome']
        t = o['type']
        rec = {'id': c['id']}
        if t == 'structural_cut':
            rec['verdict'] = 'cut of at most two vertices'
            rec['witness'] = {'type': 'cut', 'vertices': o['cut'], 'in': 'G = H2 ∪ B'}
        elif t == 'k5_minor':
            if o.get('graph') == 'completion':
                rec['witness'] = {'type': 'K5 branch sets', 'in': 'B + s1s2 + t1t2, a minor of G',
                                  'branch sets': o['branch_sets']}
            else:
                assert o.get('graph') == 'whole'
                rec['witness'] = {'type': 'K5 branch sets', 'in': 'G = H2 ∪ B', 'branch sets': o['branch_sets']}
            rec['verdict'] = 'K5 minor'
        elif t == 'proper_allowed_minor':
            op = o['operation']
            if op['type'] == 'delete_edge':
                operation, args = 'delete edge', list(op['edge'])
            else:
                assert op['type'] == 'contract_edge'
                u, v = op['edge']
                operation, args = 'contract edge', [op['keep'], v if op['keep'] == u else u]
            rec['verdict'] = 'proper nontoroidal minor'
            rec['witness'] = {'type': 'allowed reduction', 'operation': operation, 'arguments': args,
                              "B' + s1s2 + t1t2 is nonplanar": {'type': 'Kuratowski subgraph',
                                                                'edges': o['nonplanarity']['kuratowski_edges']}}
        elif t == 'proper_obstruction_minor':
            tn = o['target']
            blocks = {int(k): v for k, v in o['branch_sets'].items()}
            operation, e, nb, where = one_edge_minor_with_model(G, tn, blocks)
            rec['verdict'] = 'proper nontoroidal minor'
            rec['witness'] = {'type': 'proper minor', 'operation': operation, 'edge': e,
                              'why the minor is nontoroidal': {'reason': 'it has a minor model of %s' % tn,
                                                               'type': 'minor model', 'of': tn, 'in': where,
                                                               'branch sets': {str(k): v for k, v in sorted(nb.items())}}}
            found[tn] += 1
            ops[operation] += 1
        elif t == 'no_disjoint_rails':
            rec['verdict'] = 'excluded by a hypothesis of the case'
            rec['witness'] = {'fact': 'B has no disjoint s1-t1 and s2-t2 paths, but B = H1 contains P1 and P2',
                              'checked by': 'trying every s1-t1 path of B'}
        else:
            raise AssertionError('outcome %s: the enumeration is not complete' % t)
        verdicts.append(rec)
    inp = head('E1')
    inp.update({
        'vertex names': 'Psi_1 as in Section 4: twins p, q on {s1, s2, s3} and r, z on {t1, t2, t3}. '
                        'Interior vertices of B are i0, i1, ...',
        'the case': 'A = {s1, t1, s2, t2} and E = {s1s2, t1t2}. G = H2 ∪ B, where H2 consists of p, q, r, z, '
                    'their twelve edges and s3t3, and B (called H1 in the proof) meets H2 in A. The candidates are '
                    'the graphs B of Construction 3.18 with attachments A and added edges E: A independent in B, '
                    'interior nonempty and connected, every attachment with an interior neighbour, and '
                    'B + E nonplanar (each candidate gives a Kuratowski subgraph).',
        'H2': {'vertices': sorted(H2), 'edges': sorted(sorted((u, v)) for u in H2 for v in H2[u] if u < v)},
        'verdicts': {
            'K5 minor': 'branch sets in G, or in B + E, which is a minor of G (paths s1 p s2 and t1 r t2)',
            'cut of at most two vertices': 'a cut of G',
            'proper nontoroidal minor': "either an allowed reduction B' of B with B' + E nonplanar, which makes "
                                        "H2 ∪ B' nontoroidal by the implication proved in Lemma 4.10, or an edge "
                                        "of G = H2 ∪ B whose deletion or contraction leaves a minor model of some T_i",
            'excluded by a hypothesis of the case': 'B has no disjoint s1-t1 and s2-t2 paths; by minimality B is '
                                                    'H1, which contains P1 and P2'},
        'generation': {k: s[k] for k in ('raw_placements', 'placement_symmetry_classes', 'expanded_recipes',
                                         'rooted_candidates')}})
    summary = {'item': 'E1', 'where': WHERE['E1'], 'candidates': len(cands), 'verdicts': verdict_counts(verdicts),
               'proper minors that are T_i (copies)': dict(sorted(found.items(), key=lambda x: int(x[0][1:]))),
               'remaining': 0, 'one-edge minors of G with a minor model of some T_i': dict(sorted(ops.items()))}
    write_item(out, 'E1', inp, cands, verdicts, summary)
    return summary


# ---------------------------------------------------------------- E2

PSI1_S = {1: 'a', 2: 'c', 3: 'e'}
PSI1_T = {1: 'b', 2: 'd', 3: 'f'}
PSI1_V = set('abcdefpqrz')
PSI1_E = {frozenset(e) for e in [(w, PSI1_S[i]) for w in 'pq' for i in PSI1_S] +
          [(w, PSI1_T[i]) for w in 'rz' for i in PSI1_T] + [(PSI1_S[i], PSI1_T[i]) for i in PSI1_S]}
CHOICES = [(k, x, y) for k in (1, 2, 3) for x in sorted({1, 2, 3} - {k}) for y in sorted({1, 2, 3} - {k})]


def sorted_edges(adj):
    return sorted(sorted((u, v)) for u in adj for v in adj[u] if u < v)


def planar_or_kuratowski(adj):
    """(True, rotation system of genus 0) or (False, edges of a Kuratowski
    subgraph) for the graph adj, found with networkx with the vertices and
    edges in sorted order."""
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(sorted(adj))
    G.add_edges_from(tuple(e) for e in sorted_edges(adj))
    ok, P = nx.check_planarity(G, counterexample=True)
    if ok:
        return True, {v: list(P.neighbors_cw_order(v)) for v in G}
    return False, sorted(sorted(e) for e in P.edges())


def faces_of_choice(k, x, y):
    """Lemma 4.2: the attachments and the auxiliary edges of the four test
    graphs of the choice (k, x, y)."""
    S, T = PSI1_S, PSI1_T
    i, j = sorted({1, 2, 3} - {k})
    xp = j if x == i else i
    yp = j if y == i else i
    return {
        'J': ({S[i], T[i], S[j], T[j]}, [(S[i], T[i]), (S[j], T[j]), (S[i], S[j]), (T[i], T[j])]),
        'W': ({S[k], T[k], S[x], T[y]}, [(S[k], T[k]), (S[k], S[x]), (T[k], T[y])]),
        'small_S': ({S[k], S[xp]}, [(S[k], S[xp])]),
        'small_T': ({T[k], T[yp]}, [(T[k], T[yp])]),
    }


def psi1_bridges(H):
    """The Psi_1-bridges of H: the edges outside Psi_1 with both ends on it,
    then the components of H - V(Psi_1) (in the order of their least vertex)
    with their attaching edges."""
    out = []
    for u in sorted(H):
        for v in sorted(H[u]):
            if u < v and u in PSI1_V and v in PSI1_V and frozenset((u, v)) not in PSI1_E:
                out.append({'attachments': sorted({u, v}), 'edges': [[u, v]]})
    rest = [v for v in H if v not in PSI1_V]
    for C in sorted(components(H, rest), key=min):
        es = sorted({tuple(sorted((u, v))) for u in C for v in H[u]})
        att = sorted({v for u in C for v in H[u] if v not in C})
        out.append({'attachments': att, 'edges': [list(e) for e in es]})
    return out


def has_psi1(H):
    return all(next(iter(e - {u}), None) in H.get(u, ()) for e in PSI1_E for u in e) and \
        all(len(H.get(w, ())) == 3 for w in 'pqrz')


def evaluate_joint_test(H):
    """The joint test of Lemma 4.3 on H: (True, a passing assignment) if
    some choice passes, else (False, a record of the failure of every choice)."""
    import itertools
    import networkx as nx
    bl = psi1_bridges(H)
    named = [{'name': 'B%d' % (i + 1), 'attachments': b['attachments'], 'edges': b['edges']} for i, b in enumerate(bl)]
    names = [b['name'] for b in named]
    byname = {b['name']: b for b in named}
    choices = []
    for (k, x, y) in CHOICES:
        F = faces_of_choice(k, x, y)
        elig = {b['name']: [f for f, (sup, _) in F.items() if set(b['attachments']) <= sup] for b in named}
        blocked = [n for n, e in elig.items() if not e]
        if blocked:
            choices.append({'choice': [k, x, y], 'bridges with no face': blocked})
            continue
        asgs = []
        for combo in itertools.product(*[elig[n] for n in names]):
            asg = dict(zip(names, combo))
            bad = None
            for f, (sup, aux) in F.items():
                G = nx.Graph()
                G.add_edges_from(aux)
                for n in names:
                    if asg[n] == f:
                        G.add_edges_from(tuple(e) for e in byname[n]['edges'])
                ok, P = nx.check_planarity(G, counterexample=True)
                if not ok:
                    bad = (f, sorted(sorted(e) for e in P.edges()))
                    break
            if bad is None:
                return True, {'choice': [k, x, y], 'faces': asg}
            asgs.append({'faces': asg, 'nonplanar test': bad[0], 'Kuratowski subgraph': bad[1]})
        choices.append({'choice': [k, x, y], 'assignments': asgs})
    return False, {'bridges': named, 'choices': choices}


def bridge_ok(B, roots):
    """The hypotheses of Section 3 on a bridge: independent attachments, a
    nonempty connected interior, an interior neighbour at every attachment."""
    roots = set(roots)
    inner = [v for v in B if v not in roots]
    return (all(not (B[r] & roots) and any(u not in roots for u in B[r]) for r in roots)
            and bool(inner) and len(components(B, inner)) == 1)


def allowed_reduction(B, roots, op, args):
    """The bridge after an allowed operation (delete an edge, or contract an
    edge into its first end, which must not be an interior vertex merged into
    an attachment the wrong way), or None if the result breaks the
    hypotheses of Section 3."""
    roots = set(roots)
    H = {v: set(x) for v, x in B.items()}
    if op == 'delete edge':
        u, v = args
        H[u].discard(v)
        H[v].discard(u)
    else:
        u, v = args          # v (interior) is merged into u
        if v in roots:
            return None
        for w in H.pop(v):
            H[w].discard(v)
            if w != u:
                H[w].add(u)
                H[u].add(w)
        if u in roots:
            for w in list(H[u]):
                if w in roots:
                    H[u].discard(w)
                    H[w].discard(u)
    return H if bridge_ok(H, roots) else None


def delete_edge_adj(G, e):
    H = {x: set(s) for x, s in G.items()}
    H[e[0]].discard(e[1])
    H[e[1]].discard(e[0])
    return H


def contract_edge_adj(G, e):
    """Contract uv: v merges into u; loops deleted, parallel edges merged."""
    u, v = e
    H = {x: set(s) for x, s in G.items()}
    for w in H.pop(v):
        H[w].discard(v)
        if w != u:
            H[w].add(u)
            H[u].add(w)
    return H


def e2_one_edge_minor(G, B, roots, E):
    """For a candidate of E2 whose G has a vertex of degree three with two
    adjacent neighbours: a one-edge proper minor of G with the reason it is
    nontoroidal, an allowed reduction B' of B with B' + E nonplanar, or a
    one-edge minor keeping Psi_1 that fails every choice of Lemma 4.3.
    Returns (witness, reason) or (None, None)."""
    for (x, y) in sorted_edges(B):
        for op, args in (('delete edge', [x, y]), ('contract edge', [x, y]), ('contract edge', [y, x])):
            Bp = allowed_reduction(B, roots, op, args)
            if Bp is None:
                continue
            H = {v: set(s) for v, s in Bp.items()}
            for u, w in E:
                H.setdefault(u, set()).add(w)
                H.setdefault(w, set()).add(u)
            ok, kw = planar_or_kuratowski(H)
            if ok:
                continue
            oe = {'operation': 'delete' if op == 'delete edge' else 'contract', 'edge': args}
            return ({'type': 'proper minor', **oe,
                     'why the minor is nontoroidal': {
                         'reason': "it contains (Psi_1 ∪ B') + af for the allowed reduction B' of B given by the "
                                   "same operation, and B' + ac + cd + bd is nonplanar",
                         'Kuratowski subgraph': kw}},
                    'an allowed reduction keeping B + ac + cd + bd nonplanar')
    for (x, y) in sorted_edges(G):
        if frozenset((x, y)) in PSI1_E:
            continue
        for kind, (u, w) in (('delete', (x, y)), ('contract', (x, y)), ('contract', (y, x))):
            if kind == 'contract' and w in PSI1_V:
                continue
            M = delete_edge_adj(G, (u, w)) if kind == 'delete' else contract_edge_adj(G, (u, w))
            M = {a: nb for a, nb in M.items() if nb}
            if not has_psi1(M):
                continue
            passes, rec = evaluate_joint_test(M)
            if passes:
                continue
            return ({'type': 'proper minor', 'operation': kind, 'edge': [u, w],
                     'why the minor is nontoroidal': {'reason': 'it fails all twelve choices of Lemma 4.3',
                                                      'the twelve choices': rec}},
                    'the minor fails all twelve choices of Lemma 4.3')
    return None, None


def first_edge_in_bag(G, bag):
    return min(tuple(sorted((u, v))) for u in bag for v in G[u] if v in bag)


def package_E2(work, out):
    d = read(os.path.join(work, 'candidates.json'))
    roots = ['a', 'b', 'c', 'd']
    E = [['a', 'c'], ['c', 'd'], ['b', 'd']]
    X = adj_from(list('abcdefpqrz'), [[w, x] for w in 'pq' for x in 'ace'] + [[w, t] for w in 'rz' for t in 'bdf'] +
                 [['a', 'b'], ['c', 'd'], ['e', 'f'], ['a', 'f']])
    cands, verdicts = [], []
    reasons, ops = Counter(), Counter()
    for c in d['candidates']:
        assert sorted(c['roots']) == roots
        B = adj_from(c['vertices'], c['edges'])
        G = {x: set(nb) for x, nb in X.items()}
        for x, nb in B.items():
            G.setdefault(x, set()).update(nb)
        cands.append({'id': c['id'], 'bridge': {'vertices': sorted(c['vertices']),
                                                'edges': sorted(sorted(e) for e in c['edges'])},
                      'B + ac + cd + bd is nonplanar': {'type': 'Kuratowski subgraph',
                                                        'edges': c['nonplanarity']['kuratowski_edges']}})
        o = c['outcome']
        t = o['type']
        if t == 'adjacent_degree_three_neighbors':
            x, y = o['edge']
            v = o['vertex']
            w, how = e2_one_edge_minor(G, B, roots, E)
            if w is None:
                w = {'type': 'proper minor', 'operation': 'delete', 'edge': [x, y],
                     'why the minor is nontoroidal': {
                         'reason': 'genus at least two, by the genus computations of Appendix A',
                         'also': 'the vertex %s has degree three and its neighbours %s, %s are adjacent, so a torus '
                                 'embedding of G - %s%s would extend to G (proof of Lemma 3.4(a)); G is nontoroidal '
                                 'because B + ac + cd + bd is nonplanar' % (v, x, y, x, y)}}
                how = 'a genus computation'
            reasons[how] += 1
            ops[w['operation']] += 1
        elif t == 'proper_allowed_minor':
            op = o['operation']
            if op['type'] == 'delete_edge':
                oe = {'operation': 'delete', 'edge': list(op['edge'])}
            else:
                assert op['type'] == 'contract_edge'
                u, v = op['edge']
                oe = {'operation': 'contract', 'edge': [op['keep'], v if op['keep'] == u else u]}
            w = {'type': 'proper minor', **oe,
                 'why the minor is nontoroidal': {
                     'reason': "it contains (Psi_1 ∪ B') + af for the allowed reduction B' of B given by the same "
                               "operation, and B' + ac + cd + bd is nonplanar",
                     'Kuratowski subgraph': o['nonplanarity']['kuratowski_edges']}}
            reasons['an allowed reduction keeping B + ac + cd + bd nonplanar'] += 1
            ops[oe['operation']] += 1
        elif t == 'proper_obstruction_minor':
            tn = o['target']
            bags = {str(k): sorted(map(str, v)) for k, v in o['branch_sets'].items()}
            big = [k for k in sorted(bags) if len(bags[k]) > 1]
            assert big
            u, v = first_edge_in_bag(G, set(bags[big[0]]))
            moved = {k: sorted({u if x == v else x for x in g}) for k, g in bags.items()}
            w = {'type': 'proper minor', 'operation': 'contract', 'edge': [u, v],
                 'why the minor is nontoroidal': {'reason': 'it has a minor model of %s' % tn, 'type': 'minor model',
                                                  'of': tn,
                                                  'branch sets': {str(k): g for k, g in
                                                                  sorted(((int(k), g) for k, g in moved.items()))}}}
            reasons['a minor model of %s in the minor' % tn] += 1
            ops['contract'] += 1
        elif t == 'proper_joint_minor':
            assert o['operation']['type'] == 'delete_edge'
            x, y = o['operation']['edge']
            w = {'type': 'proper minor', 'operation': 'delete', 'edge': [x, y],
                 'why the minor is nontoroidal': {'reason': 'it fails all twelve choices of Lemma 4.3',
                                                  'the twelve choices': convert_joint_test(o['joint_test'])}}
            reasons['the minor fails all twelve choices of Lemma 4.3'] += 1
            ops['delete'] += 1
        else:
            raise AssertionError('outcome %s for %s' % (t, c['id']))
        verdicts.append({'id': c['id'], 'verdict': 'proper nontoroidal minor', 'witness': w})
    inp = head('E2')
    inp.update({
        'vertex names': 'a = s1, b = t1, c = s2, d = t2, e = s3, f = t3; twins p, q on {a, c, e} and r, z on '
                        '{b, d, f}; the paths are the edges ab, cd, ef. Interior vertices of B are i0, i1, ...',
        'the case': 'G = (Psi_1 ∪ B) + af, where B has attachments A = {a, b, c, d}. The candidates are the '
                    'graphs B of Construction 3.18 with attachments A, added edges E = {ac, cd, bd} and no '
                    'required edges: A independent in B, interior nonempty and connected, every attachment with '
                    'an interior neighbour, and B + E nonplanar (each candidate gives a Kuratowski subgraph).',
        'the rest of G': {'vertices': sorted(X), 'edges': sorted_edges(X)},
        'verdicts': 'Every candidate gets a proper nontoroidal minor obtained by deleting or contracting one '
                    'edge of G (a contraction merges the second vertex into the first). The witness says why '
                    'that minor is nontoroidal.'})
    summary = {'item': 'E2', 'where': WHERE['E2'], 'candidates': len(cands), 'verdicts': verdict_counts(verdicts),
               'one-edge minors': dict(ops), 'why nontoroidal': dict(reasons),
               'remaining': 0}
    write_item(out, 'E2', inp, cands, verdicts, summary)
    return summary


# ---------------------------------------------------------------- E10

def remap_outcome(o, phi):
    """The same outcome for an isomorphic copy; phi maps the vertices of the
    graph the outcome is about onto those of the copy."""
    o = dict(o)
    k = o['outcome']
    f = phi.__getitem__
    if k == 'K5 minor':
        o['branch_sets'] = [[f(v) for v in s] for s in o['branch_sets']]
    elif k == 'cut of order at most two':
        o['cut'] = [f(v) for v in o['cut']]
    elif k == 'catalogue graph':
        o['map'] = {f(v): t for v, t in o['map'].items()}
    elif k == 'toroidal':
        o['embedding'] = [{'vertices': [f(v) for v in p['vertices']], 'genus': p['genus'],
                           'rotation': {f(v): [f(x) for x in nb] for v, nb in p['rotation'].items()}}
                          for p in o['embedding']]
    elif k == 'proper nontoroidal minor':
        o['edge'] = [f(v) for v in o['edge']]
    else:
        raise AssertionError(k)
    return o


def package_E10(work, out):
    import pickle
    import networkx as nx
    import generate_E10 as gen
    R = Relabel({'t0': 'p', 't1': 'q', 't2': 'r', 't3': 'z'})
    d = read(os.path.join(work, 'candidates.json'))
    forms = read(os.path.join(work, 'forms.json'))
    with open(os.path.join(work, 'forms.pkl'), 'rb') as fh:
        cforms, _, rforms, _ = pickle.load(fh)
    key = lambda info: json.dumps(info, sort_keys=True)
    built = {key(i): G for G, i in gen.build(cforms, rforms)}
    # Each candidate stands for a class of isomorphic assembled graphs. When the
    # graph recorded for a class has B nonplanar but another member has B
    # planar, record that member, with the outcome carried over by an
    # isomorphism (Section 5.4 retains the graphs with B planar). Of the
    # isomorphisms, the one whose list of images of the vertices, taken in the
    # order of their names, comes last is used, so that the choice does not
    # depend on the order of a search.
    for c in d['candidates']:
        infos = [c['recipe']['first']] + c['recipe']['also']
        rep = built[key(infos[0])]
        assert sorted(sorted(map(str, e)) for e in rep.edges()) == c['graph']['edges']
        if gen.b_planar(rep):
            continue
        planar = [i for i in infos if gen.b_planar(built[key(i)])]
        if not planar:
            continue
        new = built[key(planar[0])]
        rs, ns = nx.relabel_nodes(rep, str), nx.relabel_nodes(new, str)
        order = sorted(rs)
        phi = max(nx.algorithms.isomorphism.GraphMatcher(rs, ns).isomorphisms_iter(),
                  key=lambda m: [m[v] for v in order])
        c['outcome'] = remap_outcome(c['outcome'], phi)
        c['graph'] = {'vertices': sorted(map(str, new)), 'edges': sorted(sorted(map(str, e)) for e in new.edges())}
        c['recipe'] = {'first': planar[0], 'also': [i for i in infos if i is not planar[0]]}
    SUBCASE = {'a': 'one component D_0 with bridge B_0, and the components D_c, D_d',
               'b': 'two singleton components with neighbourhood {a,b,w}, and D_c, D_d',
               'c': 'two singleton components with neighbourhood {a,b,w}, and the common component D_1'}
    cands, verdicts = [], []
    planar_found, nonplanar_found, n_recipes = Counter(), Counter(), Counter()
    n_planar = 0
    for c in d['candidates']:
        recs = [c['recipe']['first']] + c['recipe'].get('also', [])
        for r in recs:
            n_recipes[r['case']] += 1
        g = R.graph(c['graph'])
        adj = adj_from(g['vertices'], g['edges'])
        for t, nb in (('p', 'abc'), ('q', 'abc'), ('r', 'abd'), ('z', 'abd')):
            assert adj[t] == set(nb)
        X = set('abcd')
        B = {v: {u for u in nb if u not in 'pqrz' and not (u in X and v in X)}
             for v, nb in adj.items() if v not in ('p', 'q', 'r', 'z')}
        planar, w0 = planar_or_kuratowski(B)
        pw = {'type': 'embedding', 'genus': 0, 'rotation': w0} if planar else {'type': 'Kuratowski subgraph',
                                                                               'edges': w0}
        n_planar += planar
        cases = sorted({r['case'] for r in recs})
        cands.append({'id': c['id'], 'subcases': ['(%s)' % x for x in cases], 'recipes': recs, 'graph': g,
                      'graph6': to_graph6(adj, g['vertices']),
                      'B planar' if planar else 'B nonplanar': pw})
        verdict, w = convert_outcome(c['outcome'], R)
        if w['type'] == 'isomorphism':
            w['map'] = canonical_map(g, w['onto'])
        verdicts.append({'id': c['id'], 'verdict': verdict, 'witness': w})
        if verdict.startswith('isomorphic'):
            for x in cases:
                (planar_found if planar else nonplanar_found)[(verdict.split()[-1], x)] += 1
    by_num = lambda x: int(x[1:])
    inp = head('E10')
    inp.update({
        'vertex names': 'Psi_2: twins p, q with neighbourhood {a,b,c} and r, z with {a,b,d}. X = {a,b,w}. '
                        'Interior vertices: U* in B_0, C* in D_c, D* in D_d, R* in D_1, k1, k2 the two '
                        'singleton components.',
        'subcases': SUBCASE,
        'construction of D_c': 'with all four graphs H_c + E planar, at most two single vertices on c, w and '
                               'some of a, b, with or without cw; otherwise Construction 3.18 with attachments '
                               'a, b, c, w, added edges ab, aw, bw and E (E = {ac}, {bc} or {ac, bc}, required), '
                               'original pair cw as an edge or a once-subdivided segment, or on the quotient '
                               '(H_c + ac + bc)/ab with every lift of the edges at the identified vertex; every '
                               'interior component given its neighbours at c and w, and at most two further '
                               'single vertices for missing neighbours at a, b',
        'construction of D_1': 'Construction 3.18 with attachments a, c, d, w, added edges ac, ad, aw, the '
                               'subdivision using ac or ad, original pairs cw, dw as an edge or a once-subdivided '
                               'segment; the interior pieces divided into components in every way and joined by '
                               'every forest; then at most two further single vertices, each adjacent to w and to '
                               'c or d (and possibly a), and every inclusion-minimal set of attaching edges from '
                               'the retained components that makes every interior component meet w and c or d, '
                               'joins c and d through the interior, and, for N(D_1) = {a,w}, gives a an interior '
                               'neighbour; the whole configuration minimal. Both N(D_1) = {w} and {a,w}.',
        'reduction': 'forms with a one-step proper minor of the same signature (neighbourhood and planar '
                     'completions) that meets the same conditions are removed before combination',
        'B': 'B is G minus p, q, r, z and the edges among a, b, c, d. Each candidate is a class of isomorphic '
             'assembled graphs, recorded by a member whose B is planar whenever the class has one, and gives '
             'a plane rotation system of B or a Kuratowski subgraph of it. The candidates with B planar are '
             'those retained in Section 5.4.',
        'forms': {'forms of B_0 (the eleven of the proof)': forms['U forms'],
                  'forms of D_c (D_d is the same with d for c)': forms['C forms'],
                  'forms of D_1': forms['R forms'],
                  'how the forms of D_c were found': forms['C form statistics'],
                  'how the forms of D_1 were found': forms['R form statistics']},
        'assembled graphs (with repetitions)': sum(n_recipes.values())})
    summary = {'item': 'E10', 'where': WHERE['E10'],
               'forms': {'B_0': len(forms['U forms']), 'D_c': len(forms['C forms']), 'D_1': len(forms['R forms'])},
               'assembled graphs (with repetitions)': sum(n_recipes.values()),
               'candidates (up to isomorphism)': len(cands), 'with B planar': n_planar,
               'verdicts': verdict_counts(verdicts),
               'T_i with B planar': sorted({k[0] for k in planar_found}, key=by_num),
               'subcase of those': 'all (a), one component D_0' if all(k[1] == 'a' for k in planar_found)
               else 'subcases ' + ', '.join(sorted({k[1] for k in planar_found})),
               'T_i with B nonplanar (handled by Proposition 5.9, E8)':
                   sorted({k[0] for k in nonplanar_found}, key=by_num)}
    write_item(out, 'E10', inp, cands, verdicts, summary)
    return summary


# ---------------------------------------------------------------- E13

def package_E13(work, out):
    R = Relabel({'t0': 'p', 't1': 'q', 't2': 'r', 't3': 'z', 'p': "p'", 'q': "q'"})
    ROWS = {
        '1': "D_1 contains b, c; N(D_1) in {u}, {a,u}; N(D_2) = {u,v}",
        '2': "D_1 contains b, c; N(D_1) in {u}, {a,u}; N(D_2) in {v}, {a,v}",
        '3': "D_1 contains b, c; N(D_1) = N(D_2) = {u,v}",
        '4': "D_1 contains b, c; N(D_1) = {u,v}; N(D_2) in {v}, {a,v}",
        '5': "D_1 contains c, d; one neighbourhood is {u,v}",
        '6': "D_1 contains c, d; N(D_1) in {u}, {a,u}; N(D_2) in {v}, {a,v}",
    }
    d = read(os.path.join(work, 'candidates.json'))
    forms = read(os.path.join(work, 'forms.json'))
    s = read(os.path.join(work, 'summary.json'))
    assert s['complete'] and not s['catalogue graphs found']
    cands, verdicts = [], []
    by_row = Counter()
    for c in d['candidates']:
        recs = [c['recipe']['first']] + c['recipe'].get('also', [])
        g = R.graph(c['graph'])
        adj = adj_from(g['vertices'], g['edges'])
        for t, nb in (('p', 'abc'), ('q', 'abc'), ('r', 'abd'), ('z', 'abd'), ("p'", 'auv'), ("q'", 'auv')):
            assert adj[t] == set(nb), t
        rows = sorted({r['row'] for r in recs})
        for x in rows:
            by_row[x] += 1
        cands.append({'id': c['id'], 'graph': g, 'graph6': to_graph6(adj, g['vertices']),
                      'rows of the table': rows,
                      'recipes': [{'partition': r['partition'], 'D_1 form': r['R form'], 'D_2 form': r['D form'],
                                   'N(D_1)': r['N(R)'], 'N(D_2)': r['N(D)'], 'row': r['row']} for r in recs]})
        verdict, w = convert_outcome(c['outcome'], R)
        verdicts.append({'id': c['id'], 'verdict': verdict, 'witness': w})
    counts = verdict_counts(verdicts)
    assert not any(v.startswith('isomorphic') for v in counts)

    def form_list(lst):
        return [{'index': x['index'], 'signature': x['signature'], 'construction': x['origin'],
                 'mirror image under u <-> v': x.get('mirrored', False), 'graph': x['graph']} for x in lst]
    inp = head('E13')
    inp.update({
        'vertex names': "Psi_2: twins p, q with neighbourhood {a,b,c} and r, z with {a,b,d}; p', q' have "
                        "neighbourhood X = {a,u,v}. Interior vertices of D_1 start with R and those of D_2 with D.",
        'partitions': {'P1': 'D_1 contains b, c and D_2 contains d',
                       'P2': 'D_1 contains c, d and D_2 contains b',
                       'note': 'up to interchanging c, d; each neighbourhood lies among {u}, {v}, {a,u}, {a,v}, '
                               '{u,v}, and together they contain u and v'},
        'rows of the table': ROWS,
        'construction of D_1': {
            'no subdivision': 'a path x t y with t in D_1 and some of tu, tv, ta, original edges xu, xv, yu, yv, '
                              'and possibly a second single vertex (a superset of the seven forms of the text)',
            'B_1 + uv + xy nonplanar (P1, N(D_1) = {u,v})':
                'Construction 3.18 with attachments b, c, u, v, added edges bc, uv, required edge bc, original '
                'pairs bu, bv, cu, cv as an edge or a once-subdivided segment',
            'the hard tests of row two (P1)':
                'Construction 3.18 with attachments a, b, c, u, added edges bc, au, required edge bc, original '
                'pairs bu, cu; the interior pieces and u joined into one connected graph by every forest',
            'B_1 + au + ac + ad nonplanar (P2, N(D_1) in {u}, {a,u})':
                'Construction 3.18 with attachments a, c, d, u, added edges au, ac, ad, the subdivision using '
                'ac or ad, original pairs cu, du as an edge or a once-subdivided segment',
            'completion (Construction 3.18, Lemma 5.11)':
                'the interior pieces are divided into components in every way and joined by every forest; then '
                'at most two further single vertices are added, with every admissible neighbourhood, and every '
                'inclusion-minimal set of attaching edges from the retained components (and, for B_1 + uv + xy, '
                'of original edges xu, xv, yu, yv) that makes every interior component meet u or v (u, for the '
                'other two constructions) and x or y, makes D_1 connected, gives u and v a neighbour in D_1 where '
                'N(D_1) = {u,v}, and gives a a neighbour in a retained component where N(D_1) = {a,u}; the whole '
                'configuration is minimal. Mirror images under u <-> v are included.',
            'reduction': 'forms with a one-step proper minor of the same signature (neighbourhood and planar '
                         'completions) that meets the same conditions are removed'},
        'forms of D_1': {p: form_list(forms[p]['R forms']) for p in ('P1', 'P2')},
        'forms of D_2': {p: form_list(forms[p]['D forms']) for p in ('P1', 'P2')},
        'how the forms were found': {p: {'D_1': forms[p]['R statistics'], 'D_2': forms[p]['D statistics']}
                                     for p in ('P1', 'P2')},
        'assembly': 'every pair (D_1, D_2) with N(D_1) and N(D_2) together containing u and v, for both partitions',
        'assembled graphs (with repetitions)': s['cases']['assemblies']['assembled'],
        'verdicts': {'toroidal': 'an embedding of genus at most one',
                     'K5 minor': 'branch sets', 'cut of at most two vertices': 'the cut',
                     'proper nontoroidal minor': 'the edge deleted or contracted; the minor has genus at least two '
                                                 '(the genus computations of Appendix A)'}})
    summary = {'item': 'E13', 'where': WHERE['E13'],
               'forms of D_1': {p: len(forms[p]['R forms']) for p in ('P1', 'P2')},
               'forms of D_2': {p: len(forms[p]['D forms']) for p in ('P1', 'P2')},
               'assembled graphs (with repetitions)': s['cases']['assemblies']['assembled'],
               'candidates (up to isomorphism)': len(cands), 'verdicts': counts,
               'candidates by row of the table': dict(sorted(by_row.items())),
               'obstructions': 'none'}
    write_item(out, 'E13', inp, cands, verdicts, summary)
    return summary


# ---------------------------------------------------------------- E3, E4, E12: coverings

FACES_TEXT = ('Faces of a choice (k,x,y) of Lemma 4.2: J is the pair of simple faces (attachments V(P_i) ∪ V(P_j); '
              'test graph with auxiliary edges s_i s_j, t_i t_j), W the repeated face (attachments V(P_k) ∪ {s_x, '
              't_y}; auxiliary edges s_k s_x, t_k t_y), small_S and small_T the two faces of (iii) (attachments '
              '{s_k, s_x\'} and {t_k, t_y\'}; a bridge set fits such a disk exactly when it plus the edge between '
              'the two attachments is planar).')
DRAT_TEXT = ('drat-trim verifies that cover.cnf is unsatisfiable, with the proof cover.drat; the checker of this '
             'repository runs it')
ONE_CLAUSE_TEXT = ('For a toroidal graph with n_i copies of A<i>: the clause that some atom has more copies than in '
                   'it. For a nontoroidal graph: the clause that some atom in it has fewer copies.')


def read_cnf(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        lines = f.read().split('\n')
    head = lines[0].split()
    assert head[:2] == ['p', 'cnf']
    cl, cur = [], []
    for line in lines[1:]:
        if not line.strip() or line.startswith('c'):
            continue
        for x in map(int, line.split()):
            if x == 0:
                cl.append(tuple(sorted(cur)))
                cur = []
            else:
                cur.append(x)
    assert not cur and len(cl) == int(head[3])
    return int(head[2]), cl


def copy_formula(label, out, src_cnf, src_drat):
    """cover.cnf and cover.drat of the item (gzip-compressed above 25 MB)."""
    import shutil
    folder = os.path.join(out, 'enumerations', label)
    os.makedirs(folder, exist_ok=True)
    files = {}
    for src, name in ((src_cnf, 'cover.cnf'), (src_drat, 'cover.drat')):
        for stale in (os.path.join(folder, name), os.path.join(folder, name + '.gz')):
            if os.path.exists(stale):
                os.remove(stale)
        if os.path.getsize(src) > GZ_LIMIT:
            with open(src, 'rb') as f, open(os.path.join(folder, name + '.gz'), 'wb') as raw, \
                    gzip.GzipFile(filename='', mode='wb', compresslevel=9, fileobj=raw, mtime=0) as g:
                shutil.copyfileobj(f, g)
            files[name] = name + '.gz'
        else:
            shutil.copyfile(src, os.path.join(folder, name))
            files[name] = name
    return files


def genus_of_rotation(adj, rot):
    """The genus of the embedding given by the rotation system, added over components."""
    succ = {}
    for v, nb in rot.items():
        for i, u in enumerate(nb):
            succ[(v, u)] = nb[(i + 1) % len(nb)]
    total = 0
    for comp in components(adj):
        darts = {(u, v) for u in comp for v in adj[u]}
        m = len(darts) // 2
        if not m:
            continue
        faces = 0
        while darts:
            start = d = next(iter(darts))
            faces += 1
            while True:
                darts.discard(d)
                d = (d[1], succ[(d[1], d[0])])
                if d == start:
                    break
        total += (2 - len(comp) + m - faces) // 2
    return total


def torus_embedding(adj):
    """An embedding of genus at most one, found with PAGE (package graph-genus)."""
    import graph_genus
    order = sorted(adj)
    idx = {x: i for i, x in enumerate(order)}
    genus, rot = graph_genus.embed([[idx[u] for u in sorted(adj[x])] for x in order])
    assert genus <= 1
    rotation = {order[i]: [order[j] for j in nb] for i, nb in enumerate(rot)}
    return {'type': 'embedding', 'genus': genus_of_rotation(adj, rotation), 'rotation': rotation}


class AtomNames:
    """Vertices of Psi_1 by a fixed map; interior vertices renamed A<i>#<copy>@<role>."""

    def __init__(self, fixed, pattern, idx):
        import re
        self.fixed, self.pattern, self.idx, self.cache = fixed, re.compile(pattern), idx, {}

    def v(self, x):
        if x not in self.cache:
            m = self.pattern.fullmatch(x)
            if m:
                self.cache[x] = 'A%d#%s@%s' % (self.idx[m.group('atom')] + 1, m.group('copy'), m.group('role'))
            else:
                self.cache[x] = self.fixed.get(x, x)
        return self.cache[x]

    def vs(self, xs):
        return [self.v(x) for x in xs]

    def es(self, es):
        return [[self.v(a), self.v(b)] for a, b in es]

    def graph(self, g):
        vs = self.vs(g['vertices'])
        assert len(set(vs)) == len(vs)
        return {'vertices': vs, 'edges': self.es(g['edges'])}

    def rot(self, r):
        return {self.v(k): self.vs(x) for k, x in r.items()}


PSI1_NAMES = {'a': 's1', 'b': 't1', 'c': 's2', 'd': 't2', 'e': 's3', 'f': 't3'}


def package_E3(work, out):
    import re
    import endpoint_domain as dom
    A = list(dom.ATOMS)
    mx = [a['max_count'] for a in A]
    thr = [a['threshold_variables'] for a in A]
    N = AtomNames(PSI1_NAMES, r'v(?P<atom>\d+)_(?P<copy>\d+)_(?P<role>\w+)', {'%02d' % i: i for i in range(len(A))})
    core = [tuple(e) for e in dom.CORE_EDGES]
    with open(os.path.join(work, 'coverage.jsonl')) as f:
        recs = [json.loads(line) for line in f if line.strip()]
    with open(os.path.join(work, 'candidate_outcomes.jsonl')) as f:
        outs = {o['id']: o['outcome'] for o in (json.loads(line) for line in f if line.strip())}
    nv, cnf = read_cnf(os.path.join(work, 'coverage.cnf'))
    domain_cnf = dom.domain_cnf()
    paper = {}
    for a in A:
        if a['kind'] == 'triple':
            f = re.sub(r'[a-f]', lambda m: PSI1_NAMES[m.group(0)], a['form'])
            f = 'one vertex' if a['form'] in ('single', 'singleton') else f
            paper[a['id']] = '%s on {%s}' % (f, ', '.join(PSI1_NAMES[x] for x in a['support']))
        elif a['kind'] == 'quad':
            paper[a['id']] = 'one vertex on {%s}' % ', '.join(PSI1_NAMES[x] for x in a['support'])
        else:
            paper[a['id']] = 'the edge ' + ''.join(PSI1_NAMES[x] for x in a['support'] if x in PSI1_NAMES)
    quads = [i for i, a in enumerate(A) if a['kind'] == 'quad']
    cands, verdicts = [], []
    counts, found, ops, model_ops = Counter(), Counter(), Counter(), Counter()
    quad_copies = 0
    for r in recs:
        n = r['counts']
        jt = r['joint_test']
        names = {b['id']: 'B%d' % (j + 1) for j, b in enumerate(jt['bridges'])}
        bridges = [{'name': names[b['id']], 'attachments': sorted(b['attachments']), 'edges': b['edges']}
                   for b in jt['bridges']]
        g = N.graph(r['graph'])
        adj = adj_from(g['vertices'], g['edges'])
        cands.append({'id': r['id'], 'kind': 'toroidal, maximal' if r['direction'] == 'down' else 'nontoroidal, minimal',
                      'atoms': {'A%d' % (i + 1): c for i, c in enumerate(n) if c}, 'graph': g,
                      'graph6': to_graph6(adj, g['vertices'])})
        if r['direction'] == 'down':
            (ch,) = [c for c in jt['choices'] if c['passes']][:1]
            (a,) = [a for a in ch['assignments'] if a['passes']][:1]
            w = {'type': 'passing choice of Lemma 4.3',
                 'bridges': [{'name': b['name'], 'attachments': N.vs(b['attachments']), 'edges': N.es(b['edges'])}
                             for b in bridges],
                 'choice': [ch['k'], ch['x'], ch['y']], 'faces': {names[b]: f for b, f in a['mapping'].items()},
                 'plane embeddings of the test graphs': {f: N.rot(x) for f, x in a['planar_rotations'].items()}}
            verdicts.append({'id': r['id'], 'verdict': 'toroidal', 'witness': torus_embedding(adj),
                             'also: a passing choice of Lemma 4.3': w})
            counts['toroidal'] += 1
            continue
        rec = {'id': r['id'], 'nontoroidal': {'reason': 'it fails all twelve choices of Lemma 4.3',
                                               'the twelve choices': convert_joint_test(jt, N.v)}}
        o = outs[r['id']]
        t = o['type']
        if t == 'exact_survivor':
            tn = o['target']
            rec['verdict'] = 'isomorphic to ' + tn
            rec['witness'] = {'type': 'isomorphism', 'onto': tn, 'map': canonical_map(g, tn)}
            found[tn] += 1
            quad_copies += sum(1 for i in quads if n[i])
        elif t == 'proper_obstruction_minor':
            tn = o['target']
            blocks = {int(k): N.vs(v) for k, v in o['branch_sets'].items()}
            op, e, nb, where = one_edge_minor_with_model(adj, tn, blocks)
            rec['verdict'] = 'proper nontoroidal minor'
            rec['witness'] = {'type': 'proper minor', 'operation': op, 'edge': e,
                              'why the minor is nontoroidal': {'reason': 'it has a minor model of %s' % tn,
                                                               'type': 'minor model', 'of': tn, 'in': where,
                                                               'branch sets': {str(k): v for k, v in sorted(nb.items())}}}
            model_ops[op] += 1
        elif t == 'proper_joint_minor':
            op = o['operation']
            if op['type'] == 'delete_edge':
                u, v = op['edge']
                opname = 'delete'
            else:
                assert op['type'] == 'contract_edge'
                u = op['keep']
                (v,) = [x for x in op['edge'] if x != u]
                opname = 'contract'
            ops[opname] += 1
            rec['verdict'] = 'proper nontoroidal minor'
            rec['witness'] = {'type': 'proper minor', 'operation': opname, 'edge': N.vs([u, v]),
                              'why the minor is nontoroidal': {
                                  'reason': 'it fails all twelve choices of Lemma 4.3',
                                  'the twelve choices': convert_joint_test(o['joint_test'], N.v)}}
        else:
            raise AssertionError('outcome %s for %s' % (t, r['id']))
        counts[rec['verdict']] += 1
        verdicts.append(rec)
    files = copy_formula('E3', out, os.path.join(work, 'coverage.cnf'), os.path.join(work, 'coverage.drat'))
    inp = head('E3')
    inp.update({
        'vertex names': 'Psi_1 with twins p, q on {s1, s2, s3} and r, z on {t1, t2, t3}; the paths are the edges '
                        's1t1, s2t2, s3t3. Interior vertices of copy k of atom A<i> are named A<i>#k@<role>.',
        'the base graph': {'edges': sorted(sorted(N.vs(e)) for e in core)},
        'atoms': [{'id': 'A%d' % (i + 1), 'is': paper[a['id']],
                   'attachments': sorted(N.v(x) for x in a['support']),
                   'edges': [[x if x.startswith('@') else N.v(x), y if y.startswith('@') else N.v(y)]
                             for x, y in a['edges']],
                   'most copies': mx[i], 'variables': thr[i]} for i, a in enumerate(A)],
        'the formula': {
            'variables': 'Variable v of atom A<i> with list [v1, ..., vm] means "at least j copies of A<i>" for '
                         'v = vj. Variables after %d are auxiliary (a counter for the bound of five '
                         'three-attachment bridges).' % sum(mx),
            'base clauses': {g['name']: g['clause_count'] for g in domain_cnf.groups},
            'base clause count': len(domain_cnf.clauses),
            'one clause per recorded graph': ONE_CLAUSE_TEXT,
            'clause count': len(cnf), 'files': files, 'DRAT check': DRAT_TEXT},
        'what this shows': 'Every assembly allowed by the base clauses either lies inside a recorded toroidal '
                           'graph, so it is toroidal, or contains a recorded nontoroidal graph. An obstruction of '
                           'this case therefore equals a recorded nontoroidal graph, and the verdicts of those '
                           'leave only T48 and T58.',
        'how a choice fails or passes': FACES_TEXT})
    summary = {'item': 'E3', 'where': WHERE['E3'], 'recorded graphs': len(cands),
               'toroidal, maximal': sum(1 for c in cands if c['kind'] == 'toroidal, maximal'),
               'nontoroidal, minimal': sum(1 for c in cands if c['kind'] != 'toroidal, maximal'),
               'verdicts': dict(sorted(counts.items())),
               'isomorphic to T_i (copies)': dict(found),
               'four-attachment bridges in the copies of T48 and T58': quad_copies,
               'one-edge minors failing all twelve choices': dict(ops),
               'formula': {'variables': nv, 'clauses': len(cnf), 'DRAT': 'verified by drat-trim'},
               'one-edge minors with a minor model of some T_i': dict(model_ops)}
    write_item(out, 'E3', inp, cands, verdicts, summary)
    return summary


COVER_INPUT = {
    'what this shows': 'Each allowed assembly is a subgraph of a recorded toroidal graph or contains a recorded '
                       'nontoroidal graph. By minor minimality, each obstruction in the case is one of the '
                       'recorded nontoroidal graphs.',
    'not in this item': 'The genus computations in genus/ certify the recorded nontoroidal graphs (Appendix A).'}


def form_name(form, roles, R):
    x, y, z = [R(v) for v in roles]
    return {'single': 'one vertex', 'Tx': 'T(%s)' % x, 'Ty': 'T(%s)' % y, 'Tz': 'T(%s)' % z,
            'Sxy': 'S(%s;%s)' % (x, y), 'Syx': 'S(%s;%s)' % (y, x), 'Q': 'Q(%s;%s,%s)' % (z, x, y),
            'FC1': 'the first four-cycle form of Lemma 3.14 with (x,y,z) = (%s,%s,%s)' % (x, y, z),
            'FC2': 'the second four-cycle form of Lemma 3.14 with (x,y,z) = (%s,%s,%s)' % (x, y, z)}[form]


def e4_names(atoms, N):
    out = []
    for a in atoms:
        k = a['kind']
        if k in ('matching', 'cross'):
            out.append('the edge ' + ''.join(N.v(x) for x in a['attachments']))
        elif k == 'triple':
            out.append('%s on {%s}' % (form_name(a['form'], a['roles'], N.v), ', '.join(N.vs(a['support']))))
        elif k == 'quad':
            out.append('one vertex on {%s}' % ', '.join(N.vs(a['support'])))
        else:
            out.append('one vertex on {%s}, type (c), only as the bridge containing P_%d' %
                       (', '.join(N.vs(a['support'])), a['path m']))
    return out


def e4_part(shape):
    if shape['all paths edges']:
        return 'all three edges s_m t_m present'
    if shape['singletons only']:
        return '(a)'
    tsq, fc = shape['uses T, S or Q'], shape['uses a four-cycle form']
    return '(b) and (c)' if tsq and fc else ('(b)' if tsq else '(c)')


def e12_names(atoms, N):
    out = []
    for a in atoms:
        if a['kind'] == 'edge':
            out.append('the edge %s%s' % (a['r'], a['end']))
        else:
            s, t = a['support']
            if a['form'] == 'single':
                out.append('one vertex on {%s, %s, %s}' % (a['r'], s, t))
            else:
                out.append('%s on {%s, %s, %s}' % (form_name(a['form'], (s, t, a['r']), N.v), a['r'], s, t))
    return out


E4_INPUT = {
    'vertex names': 'Psi_1 with twins p, q on {s1, s2, s3} and r, z on {t1, t2, t3}. Interior vertices of copy k of '
                    'atom A<i> are named A<i>#k@<role>. In the forms on {c1, c2, h}, the roles x, y, z are c1, c2, h; '
                    'for a support containing s_m, t_m they are g, e, h as in Lemma 4.15.',
    'constraints of the case (the base clauses)': [
        'Lemma 3.4(a): an atom whose degree-three interior vertex sees both ends of an edge s_i t_j excludes that '
        'edge.',
        'Lemma 4.7: not both s_i t_j and s_j t_i.',
        'At most three bridges on one attachment triple, or four when it contains s_m, t_m and s_m t_m is absent '
        '(Lemma 4.16); at most three single vertices on one set {s_i, t_i, s_j, t_j} (most copies 3).',
        'The Q, four-cycle and type (c) forms occur only as the bridge containing P_m: at most one of them for each '
        'm, and only when s_m t_m is absent.',
        'For each m, s_m t_m is present or some bridge attaches to both s_m and t_m.'],
    'parts': {'(a)': 'every bridge has at most one interior vertex', '(b)': 'some bridge has a T, S or Q form',
              '(c)': 'some bridge has a four-cycle form',
              'all three edges s_m t_m present': 'outside the case of the proposition (some s_m t_m is absent); the '
                                                 'text notes that these recover T48 and T58'}}

E12_MAP = {'t0': 'p', 't1': 'q', 't2': 'r', 't3': 'z', 'p': "p'", 'q': "q'"}
E12_INPUT = {
    'vertex names': "Psi_2 with twins p, q on {a, b, c} and r, z on {a, b, d}; the twins p', q' of B have "
                    "neighbourhood X = {a, u, v}. Interior vertices of copy k of atom A<i> are named A<i>#k@<role>.",
    'constraints of the case (the base clauses)': [
        'For each g in {b, c, d}: at most one form, on one pair st in {au, av, uv} (two singleton interiors are two '
        'copies of the one-vertex atom).',
        'A form on {g, s, t} allows only the original edges gs, gt; the edge ag is never an atom.',
        'Some atom is present for each g: a form, or for an empty interior an edge gu or gv.']}


def package_cover(label, work, out, fixed, names_fn, extra_input, parts_fn=None):
    """E4 and E12: the records of a covering (generate_E4.py, generate_E12.py)."""
    from pysat.solvers import Solver
    dom = read(os.path.join(work, 'domain.json'))
    atoms, mx = dom['atoms'], dom['maxima']
    thr, v = [], 0
    for m in mx:
        thr.append(list(range(v + 1, v + m + 1)))
        v += m
    ids = {a['id']: i for i, a in enumerate(atoms)}
    nv, cnf = read_cnf(os.path.join(work, 'cover.cnf'))
    passing = read(os.path.join(work, 'passing.json'))['passing']
    failing = read(os.path.join(work, 'failing.json'))['failing']
    base = cnf[:len(cnf) - len(passing) - len(failing)]

    def vec(d):
        n = [0] * len(atoms)
        for k, c in d.items():
            n[ids[k]] = c
        return n
    N = AtomNames(fixed, r'(?P<atom>[^#@]+)#(?P<copy>\d+)@(?P<role>\w+)', ids)
    solver = Solver(name='cadical153', bootstrap_with=[list(c) for c in base])
    cands, verdicts = [], []
    counts, found, two_single = Counter(), Counter(), Counter()
    by_part = {}
    outside = 0
    for p in passing:
        n = vec(p['vector'])
        c = {'id': p['id'], 'kind': 'toroidal, maximal', 'atoms': {'A%d' % (i + 1): x for i, x in enumerate(n) if x}}
        if p.get('image of'):
            c['image of'] = p['image of']
        cands.append(c)
        (e,) = p['embedding']
        verdicts.append({'id': p['id'], 'verdict': 'toroidal',
                         'witness': {'type': 'embedding', 'genus': e['genus'], 'rotation': N.rot(e['rotation'])}})
        counts['toroidal'] += 1
    for f in failing:
        n = vec(f['vector'])
        assum = [(thr[i][k - 1] if k <= n[i] else -thr[i][k - 1]) for i in range(len(n)) for k in range(1, mx[i] + 1)]
        inside = solver.solve(assumptions=assum)
        assert inside == f['in domain'], f['id']
        c = {'id': f['id'], 'kind': 'nontoroidal, minimal', 'atoms': {'A%d' % (i + 1): x for i, x in enumerate(n) if x},
             'in the family': inside}
        if f.get('image of'):
            c['image of'] = f['image of']
        cands.append(c)
        if not inside:
            outside += 1
            continue
        g = N.graph(f['graph'])
        c['graph'] = g
        c['graph6'] = to_graph6(adj_from(g['vertices'], g['edges']), g['vertices'])
        if parts_fn:
            c['part'] = parts_fn(f['shape'])
        verdict, w = convert_outcome(f['outcome'], N)
        if w['type'] == 'isomorphism':
            w['map'] = canonical_map(g, w['onto'])
        verdicts.append({'id': f['id'], 'verdict': verdict, 'witness': w})
        counts[verdict] += 1
        if verdict.startswith('isomorphic'):
            tn = verdict.split()[-1]
            found[tn] += 1
            if parts_fn:
                by_part.setdefault(c['part'], Counter())[tn] += 1
            if any(x == 2 for i, x in enumerate(n) if atoms[i].get('form') == 'single'):
                two_single[tn] += 1
    files = copy_formula(label, out, os.path.join(work, 'cover.cnf'), os.path.join(work, 'cover.drat'))
    inp = head(label)
    inp.update({
        'the base graph': {'edges': sorted(sorted(N.vs(e)) for e in dom['base edges'])},
        'atoms': [{'id': 'A%d' % (i + 1), 'is': name, 'attachments': sorted(N.v(x) for x in a['attachments']),
                   'edges': [[x if x.startswith('@') else N.v(x), y if y.startswith('@') else N.v(y)]
                             for x, y in a['edges']],
                   'most copies': mx[i], 'variables': thr[i]}
                  for i, (a, name) in enumerate(zip(atoms, names_fn(atoms, N)))],
        'the formula': {
            'variables': 'Variable vj of atom A<i> (its list [v1, ..., vm]) means "at least j copies of A<i>". '
                         'Variables after %d are auxiliary (cardinality counters).' % v,
            'base clauses': '%d clauses: v(j+1) implies vj for each atom, then the constraints of the case listed '
                            'below.' % len(base),
            'one clause per recorded graph': ONE_CLAUSE_TEXT, 'clause count': len(cnf), 'files': files,
            'DRAT check': DRAT_TEXT}})
    inp.update(COVER_INPUT)
    inp.update(extra_input)
    nt = sum(1 for c in cands if c['kind'] != 'toroidal, maximal')
    by_num = lambda kv: int(kv[0][1:])
    summary = {'item': label, 'where': WHERE[label], 'recorded graphs': len(cands),
               'toroidal, maximal': len(passing), 'nontoroidal, minimal': nt,
               'nontoroidal, minimal, in the family': sum(1 for c in cands if c.get('in the family')),
               'nontoroidal, minimal, outside the family (used only for the covering)': outside,
               'verdicts': dict(sorted(counts.items())),
               'isomorphic to T_i (copies)': dict(sorted(found.items(), key=by_num)),
               'formula': {'variables': nv, 'clauses': len(cnf), 'base clauses': len(base),
                           'DRAT': 'verified by drat-trim'}}
    if label == 'E12':
        summary['copies of T_i with a two-singleton interior'] = sum(two_single.values())
    if by_part:
        summary['T_i by part'] = {k: dict(sorted(c.items(), key=by_num)) for k, c in sorted(by_part.items())}
    write_item(out, label, inp, cands, verdicts, summary)
    return summary


def package_E4(work, out):
    return package_cover('E4', work, out, PSI1_NAMES, e4_names, E4_INPUT, parts_fn=e4_part)


def package_E12(work, out):
    return package_cover('E12', work, out, E12_MAP, e12_names, E12_INPUT)
