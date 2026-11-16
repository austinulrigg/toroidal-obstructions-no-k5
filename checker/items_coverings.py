"""The covering enumerations E4 (Proposition 4.24) and E12 (Proposition 5.24).

A recorded graph of these items is the base graph of the case plus a stated
number of copies of each atom (an edge, or a bridge given by its edges); the
interior vertex '@x' of copy k of atom A is named A#k@x. The candidates are
the recorded graphs: "toroidal, maximal" ones, whose verdict is a torus
embedding, and "nontoroidal, minimal" ones. Those in the family of the case
have a verdict and a witness; those outside it are used only for the covering
and are recorded by their atom counts. The formula cover.cnf says that some
atom-count vector allowed by the constraints of the case lies below no
toroidal record and above no nontoroidal record, and cover.drat refutes it.

What is checked, using only the files of the item and the graphs T1-T74 of O1:

The case
  * the base graph, the vertex names, the constraint list, the description of
    the record clauses and (E4) the parts are the texts whose reading is
    implemented here;
  * each atom's description ('is') is parsed, and its edges are exactly the
    form it names, up to renaming its interior vertices, by the definitions
    of the paper: one vertex joined to the attachments; T(w), two vertices
    joined to u and to the other two attachments, and the edge uw; S(x;y),
    two vertices joined to x, y, u and a vertex f joined to x, z, u;
    Q(x;y,z), the four-cycle upvq with the edges xp, xq, uy, vz; the two
    four-cycle forms of Lemma 3.14, the cycle a1 b1 a2 b2 with z a1, z a2 and
    either y b1, y b2, x a2 or x b1, x b2, y a2. Its attachments are the ends
    of its edges outside the interior;
  * the atoms are exactly the edge and bridge types that the paper lists for
    the case, each once. E4: the nine edges s_i t_j; on each of the eighteen
    triples {c1, c2, h} (c1, c2 on one side) one vertex, T(c1), T(c2), T(h),
    S(c1;c2), S(c2;c1), and on the twelve triples containing some s_m, t_m
    also Q(h;g,e) and the two four-cycle forms with (x,y,z) = (g,e,h), where
    e is the third attachment and g the end of P_m on its side; one vertex on
    each set {s_i,t_i,s_j,t_j}, and one vertex of type (c) for P_m on each set
    {s_m,t_m,s_x,t_y}. E12: for each g in {b,c,d} and st in {au,av,uv}, one
    vertex, T(g), T(s), T(t), S(s;t), S(t;s), Q(g;s,t) and the two four-cycle
    forms with (x,y,z) = (s,t,g), all on {g,s,t}, and the edges gu, gv;
  * the most copies of each atom is the bound that the constraints give
    (E4: four on a triple containing some s_m, t_m, three on the other triples
    and on the sets {s_i,t_i,s_j,t_j}, one for the edges and for the Q,
    four-cycle and type (c) forms; E12: two for one vertex, one otherwise),
    and the variables of the atoms are 1, ..., N, one for each copy.

The records
  * every record is rebuilt from its atom counts (each between 1 and the most
    copies); where a graph is recorded it equals the rebuilt graph, and its
    graph6 string encodes it with the vertices in the order listed;
  * each toroidal record has an embedding of the rebuilt graph of genus at
    most one, checked by tracing its faces;
  * each nontoroidal record in the family has exactly one verdict, checked by
    its kind: K5 branch sets, a cut of at most two vertices, an isomorphism
    onto the T_i it names, or a proper minor whose edge is in the graph (why
    the minor is nontoroidal is left to the genus computations in the folder
    genus unless a minor model or isomorphism of some T_i is recorded for it);
    toroidal records have the verdict "toroidal" and records outside the
    family have no verdict;
  * 'image of' (E4): the record is the image of the named record of the same
    kind under a relabelling of the attachments that maps the base graph and
    the set of atoms onto themselves;
  * 'part' (E4): recomputed from the atoms present.

The formula
  * cover.cnf is, as a multiset of clauses with sorted literals, exactly
      (a) the chain clauses v(j+1) -> vj of each atom,
      (b) the clauses of the constraints of the case, rebuilt here from the
          constraint list and the atoms. In E4 the bound on the number of
          bridges on one attachment triple is encoded, for each triple in the
          order of its first atom, by a sequential counter for at most four,
          then one for at most three (every clause of this second counter
          also containing "not s_m t_m" when the triple contains s_m, t_m);
          the inputs are the variables of the atoms on the triple, in order,
          and the auxiliary variables are numbered from N + 1 in order of
          first use,
      (c) one clause per record: for a toroidal record with n_i copies of A_i
          the clause "some atom has more than n_i copies", for a nontoroidal
          record the clause "some atom in it has fewer copies";
    and its header gives the number of variables used and of clauses;
  * each sequential counter is shown, by dynamic programming over all
    assignments of its inputs, to have an assignment of its auxiliary
    variables exactly when at most its bound of inputs are true; since the
    auxiliary variables of different counters are distinct and occur nowhere
    else, (b) is equivalent to the constraint list as read here;
  * the 'in the family' flag of every nontoroidal record agrees with two
    evaluations of its count vector: the constraint list as read here, and
    the clauses (b) (the counters by the bound they were shown to express);
  * in E12, which has no auxiliary variables and whose constraints separate
    by g, the models of (b) are enumerated for each g and are exactly the
    count vectors that the constraint list allows;
  * cover.drat is checked with drat-trim, given by the environment variable
    DRAT_TRIM, found on the PATH, or compiled from tools/drat-trim.c of the
    repository; without it the proof is reported as not checked.

The summary
  * every count of summary.json is recomputed from the records and the
    formula.

Not checked here: that the nontoroidal records, and the proper minors named
in the verdicts, are nontoroidal (their genus computations are in the tables
of Appendix A). The minimality and maximality of the records are not checked
either; the covering argument does not use them.
"""
import gzip
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from itertools import combinations, permutations

from witnesses import (CheckFailed, need, load, graph, same_graph, from_graph6, check_witness,
                       check_embedding, catalogue, delete_edge, contract_edge)

TOR, NONTOR = 'toroidal, maximal', 'nontoroidal, minimal'
IN_FAMILY = 'nontoroidal, minimal, in the family'
OUTSIDE = 'nontoroidal, minimal, outside the family (used only for the covering)'

# The statements of input.json whose reading is implemented below. If one of
# them changes, this module must be read again against the new text.
RECORD_CLAUSES = ('For a toroidal graph with n_i copies of A<i>: the clause that some atom has more copies than in '
                  'it. For a nontoroidal graph: the clause that some atom in it has fewer copies.')
VARIABLES_TEXT = re.compile(r'Variable vj of atom A<i> \(its list \[v1, \.\.\., vm\]\) means "at least j copies of '
                            r'A<i>"\. Variables after (\d+) are auxiliary \(cardinality counters\)\.')
BASE_TEXT = re.compile(r'(\d+) clauses: v\(j\+1\) implies vj for each atom, then the constraints of the case listed '
                       r'below\.')
STATEMENTS = {
    'E4': {
        'where': 'Proposition 4.24',
        'vertex names': 'Psi_1 with twins p, q on {s1, s2, s3} and r, z on {t1, t2, t3}. Interior vertices of copy k '
                        'of atom A<i> are named A<i>#k@<role>. In the forms on {c1, c2, h}, the roles x, y, z are c1, '
                        'c2, h; for a support containing s_m, t_m they are g, e, h as in Lemma 4.15.',
        'constraints': [
            'Lemma 3.4(a): an atom whose degree-three interior vertex sees both ends of an edge s_i t_j excludes '
            'that edge.',
            'Lemma 4.7: not both s_i t_j and s_j t_i.',
            'At most three bridges on one attachment triple, or four when it contains s_m, t_m and s_m t_m is absent '
            '(Lemma 4.16); at most three single vertices on one set {s_i, t_i, s_j, t_j} (most copies 3).',
            'The Q, four-cycle and type (c) forms occur only as the bridge containing P_m: at most one of them for '
            'each m, and only when s_m t_m is absent.',
            'For each m, s_m t_m is present or some bridge attaches to both s_m and t_m.'],
        'parts': {'(a)': 'every bridge has at most one interior vertex',
                  '(b)': 'some bridge has a T, S or Q form',
                  '(c)': 'some bridge has a four-cycle form',
                  'all three edges s_m t_m present': 'outside the case of the proposition (some s_m t_m is absent); '
                                                     'the text notes that these recover T48 and T58'},
    },
    'E12': {
        'where': 'Proposition 5.24',
        'vertex names': "Psi_2 with twins p, q on {a, b, c} and r, z on {a, b, d}; the twins p', q' of B have "
                        "neighbourhood X = {a, u, v}. Interior vertices of copy k of atom A<i> are named "
                        "A<i>#k@<role>.",
        'constraints': [
            'For each g in {b, c, d}: at most one form, on one pair st in {au, av, uv} (two singleton interiors are '
            'two copies of the one-vertex atom).',
            'A form on {g, s, t} allows only the original edges gs, gt; the edge ag is never an atom.',
            'Some atom is present for each g: a form, or for an empty interior an edge gu or gv.'],
    },
}


# ---------------------------------------------------------------- the forms of the paper

def _support(text):
    return frozenset(x.strip() for x in text.split(','))


_DESCRIPTIONS = [
    ('edge', re.compile(r'the edge (\S+)')),
    ('one', re.compile(r'one vertex on \{([^{}]*)\}')),
    ('type c', re.compile(r'one vertex on \{([^{}]*)\}, type \(c\), only as the bridge containing P_(\d+)')),
    ('T', re.compile(r'T\(([^(),;]+)\) on \{([^{}]*)\}')),
    ('S', re.compile(r'S\(([^(),;]+);([^(),;]+)\) on \{([^{}]*)\}')),
    ('Q', re.compile(r'Q\(([^(),;]+);([^(),;]+),([^(),;]+)\) on \{([^{}]*)\}')),
    ('FC', re.compile(r'the (first|second) four-cycle form of Lemma 3\.14 with \(x,y,z\) = '
                      r'\(([^(),;]+),([^(),;]+),([^(),;]+)\) on \{([^{}]*)\}')),
]


def parse_description(text, edges):
    """The form that an atom's description names, as a key:
    ('edge', {x, y}), ('one', support, None), ('one', support, {s_m, t_m}) for
    type (c) with the path P_m, ('T', w, support), ('S', x, y, support) for
    S(x;y), ('Q', x, y, z) for Q(x;y,z), ('FC1' or 'FC2', x, y, z)."""
    for name, rx in _DESCRIPTIONS:
        m = rx.fullmatch(text)
        if m:
            break
    else:
        raise CheckFailed('unknown atom description: %s' % text)
    g = m.groups()
    if name == 'edge':
        need(len(edges) == 1 and not any(v.startswith('@') for v in edges[0]), 'the edge atom', text, 'is not one edge')
        x, y = edges[0]
        need(g[0] in (x + y, y + x), 'the description', text, 'does not name the edge', x, y)
        return ('edge', frozenset((x, y)))
    if name == 'one':
        return ('one', _support(g[0]), None)
    if name == 'type c':
        k = int(g[1])
        return ('one', _support(g[0]), frozenset(('s%d' % k, 't%d' % k)))
    if name == 'T':
        sup = _support(g[1])
        need(len(sup) == 3 and g[0] in sup, 'bad description', text)
        return ('T', g[0], sup)
    if name == 'S':
        sup = _support(g[2])
        need(len(sup) == 3 and g[0] != g[1] and {g[0], g[1]} <= sup, 'bad description', text)
        return ('S', g[0], g[1], sup)
    if name == 'Q':
        need(_support(g[3]) == {g[0], g[1], g[2]} and len({g[0], g[1], g[2]}) == 3, 'bad description', text)
        return ('Q', g[0], g[1], g[2])
    need(_support(g[4]) == {g[1], g[2], g[3]} and len({g[1], g[2], g[3]}) == 3, 'bad description', text)
    return ('FC1' if g[0] == 'first' else 'FC2', g[1], g[2], g[3])


def key_support(key):
    if key[0] in ('edge', 'one', 'T', 'S'):
        return key[-1] if key[0] in ('T', 'S') else key[1]
    return frozenset(key[1:])


def form_edges(key):
    """The edges of the form named by key, by the definitions of the paper
    (interior vertices start with '@')."""
    kind = key[0]
    if kind == 'edge':
        return [tuple(sorted(key[1]))]
    if kind == 'one':
        return [('@h', v) for v in sorted(key[1])]
    if kind == 'T':                     # two vertices joined to u and the other two attachments, and uw
        w = key[1]
        o1, o2 = sorted(key[2] - {w})
        return [(d, v) for d in ('@d1', '@d2') for v in (o1, o2, '@u')] + [('@u', w)]
    if kind == 'S':                     # S(x;y): two vertices joined to x, y, u; f joined to x, z, u
        x, y = key[1], key[2]
        (z,) = key[3] - {x, y}
        return [(d, v) for d in ('@d1', '@d2') for v in (x, y, '@u')] + [('@f', v) for v in (x, z, '@u')]
    if kind == 'Q':                     # Q(x;y,z): the four-cycle upvq and xp, xq, uy, vz
        x, y, z = key[1:]
        return [('@u', '@p'), ('@p', '@v'), ('@v', '@q'), ('@q', '@u'), (x, '@p'), (x, '@q'), ('@u', y), ('@v', z)]
    x, y, z = key[1:]                   # Lemma 3.14: the cycle a1 b1 a2 b2, z a1, z a2, then ...
    cyc = [('@a1', '@b1'), ('@b1', '@a2'), ('@a2', '@b2'), ('@b2', '@a1'), (z, '@a1'), (z, '@a2')]
    if kind == 'FC1':
        return cyc + [(y, '@b1'), (y, '@b2'), (x, '@a2')]
    return cyc + [(x, '@b1'), (x, '@b2'), (y, '@a2')]


def same_up_to_interior(es, fs):
    """Some renaming of the interior vertices ('@...') of es gives exactly fs."""
    ia = sorted({v for e in es for v in e if v.startswith('@')})
    ib = sorted({v for e in fs for v in e if v.startswith('@')})
    target = {frozenset(e) for e in fs}
    if len(ia) != len(ib) or len(es) != len(target):
        return False
    for perm in permutations(ib):
        m = dict(zip(ia, perm))
        if {frozenset((m.get(x, x), m.get(y, y))) for x, y in es} == target:
            return True
    return False


def map_key(key, sg):
    """The image of an atom key under a relabelling sg of the attachments."""
    f = lambda v: sg.get(v, v)
    fs = lambda s: frozenset(f(v) for v in s)
    kind = key[0]
    if kind == 'edge':
        return ('edge', fs(key[1]))
    if kind == 'one':
        return ('one', fs(key[1]), None if key[2] is None else fs(key[2]))
    if kind == 'T':
        return ('T', f(key[1]), fs(key[2]))
    if kind == 'S':
        return ('S', f(key[1]), f(key[2]), fs(key[3]))
    return (kind,) + tuple(f(v) for v in key[1:])


class Atom:
    __slots__ = ('idx', 'id', 'key', 'form', 'support', 'edges', 'mx', 'vars', 'deg3', 'path', 'g', 'end')

    def __init__(self, idx, a):
        self.idx = idx
        self.id = a.get('id')
        need(self.id == 'A%d' % (idx + 1), 'atom', idx + 1, 'has id', self.id)
        need(set(a) == {'id', 'is', 'attachments', 'edges', 'most copies', 'variables'},
             'atom', self.id, 'has fields', sorted(a))
        es = [tuple(e) for e in a['edges']]
        need(es and all(len(e) == 2 and e[0] != e[1] and all(isinstance(v, str) for v in e) for e in es),
             'atom', self.id, 'has a bad edge')
        need(len({frozenset(e) for e in es}) == len(es), 'atom', self.id, 'repeats an edge')
        self.edges = es
        att = {v for e in es for v in e if not v.startswith('@')}
        need(len(a['attachments']) == len(set(a['attachments'])) and set(a['attachments']) == att,
             'atom', self.id, ': the attachments are not the ends of its edges outside the interior')
        self.key = parse_description(a['is'], es)
        self.form = self.key[0]
        self.support = frozenset(att)
        need(key_support(self.key) == self.support, 'atom', self.id, ': the description names other attachments')
        need(same_up_to_interior(es, form_edges(self.key)), 'atom', self.id, ': its edges are not the form', a['is'])
        self.mx = a['most copies']
        self.vars = list(a['variables'])
        need(isinstance(self.mx, int) and self.mx >= 1 and len(self.vars) == self.mx
             and all(isinstance(v, int) for v in self.vars), 'atom', self.id, ': bad most copies or variables')
        # pairs of attachments that are both neighbours of an interior vertex of degree three
        nb = defaultdict(set)
        for x, y in es:
            nb[x].add(y)
            nb[y].add(x)
        self.deg3 = {frozenset(p) for w in nb if w.startswith('@') and len(nb[w]) == 3
                     for p in combinations(sorted(v for v in nb[w] if not v.startswith('@')), 2)}
        self.path = self.g = self.end = None


# ---------------------------------------------------------------- sequential counters

def counter_symbolic(n, k):
    """The clauses of the sequential counter for "at most k of x1, ..., xn".
    A literal is (sign, symbol), the symbols being ('x', i) and ('s', i, j),
    where s(i, j) stands for "at least j of x1, ..., xi"; it exists for
    1 <= j <= k and j <= i <= j + n - k - 1. For t = 1, ..., n - k, in turn:
        x_t -> s(t,1), and for j = 1, ..., k, with i = t + j - 1,
        s(i,j) -> s(i+1,j)                  (if t < n - k)
        x_{i+1} and s(i,j) -> s(i+1,j+1)    (if j < k)
        not both x_{i+1} and s(i,k)         (if j = k).
    For k >= n there is no clause, for k = n - 1 the one clause "not all of
    x1, ..., xn", for k = 0 the clauses "not x_i". (This is the layout of the
    sequential counter of PySAT, so the auxiliary variables numbered in order
    of first use are those of the formula.)"""
    need(k >= 0, 'negative bound')
    X = lambda i: ('x', i)
    if k >= n:
        return []
    if k == n - 1:
        return [[(-1, X(i)) for i in range(1, n + 1)]]
    if k == 0:
        return [[(-1, X(i))] for i in range(1, n + 1)]
    out = []
    for t in range(1, n - k + 1):
        out.append([(-1, X(t)), (1, ('s', t, 1))])
        for j in range(1, k + 1):
            i = t + j - 1
            if t < n - k:
                out.append([(-1, ('s', i, j)), (1, ('s', i + 1, j))])
            if j < k:
                out.append([(-1, X(i + 1)), (-1, ('s', i, j)), (1, ('s', i + 1, j + 1))])
            else:
                out.append([(-1, X(i + 1)), (-1, ('s', i, j))])
    return out


def counter_clauses(lits, k, guard, top):
    """The counter for "at most k of lits" with its auxiliary variables numbered
    top + 1, ... in order of first use; with a guard g every clause also gets
    the literal -g. Returns the clauses (sorted tuples) and the new top."""
    num = {}
    out = []
    for cl in counter_symbolic(len(lits), k):
        c = []
        for sign, s in cl:
            if s[0] == 'x':
                v = lits[s[1] - 1]
            else:
                if s not in num:
                    top += 1
                    num[s] = top
                v = num[s]
            c.append(sign * v)
        if guard is not None:
            c.append(-guard)
        out.append(tuple(sorted(c)))
    return out, top


_COUNTER_CHECKED = set()


def check_counter(n, k):
    """By dynamic programming over x1, ..., xn: the clauses of
    counter_symbolic(n, k) have an assignment of their s-symbols exactly when
    at most k of x1, ..., xn are true, for each of the 2^n assignments of the
    x. Every clause involves only x_i and the s-symbols of layers i - 1 and i
    for one i, so the assignments of layer i reachable from a satisfying
    assignment of the earlier layers can be carried along, together with the
    number of true x (capped at k + 1)."""
    if (n, k) in _COUNTER_CHECKED:
        return
    sym = counter_symbolic(n, k)
    if k >= n:
        need(not sym, 'counter')
    elif k == n - 1:
        need(sym == [[(-1, ('x', i)) for i in range(1, n + 1)]], 'counter')
    else:
        layers, regs = defaultdict(list), defaultdict(set)
        for cl in sym:
            idx = [s[1] for _, s in cl]
            top = max(idx)
            xs = [s for _, s in cl if s[0] == 'x']
            need(min(idx) >= top - 1 and len(xs) <= 1 and all(s[1] == top for s in xs), 'counter clause not layered')
            layers[top].append(cl)
            for _, s in cl:
                if s[0] == 's':
                    regs[s[1]].add(s)
        need(not regs[0] and not regs[n], 'counter registers out of range')
        pos = {i: {s: b for b, s in enumerate(sorted(regs[i]))} for i in range(n + 1)}
        states = {(0, frozenset([0]))}
        for i in range(1, n + 1):
            conds = []                      # each clause as the condition that falsifies it
            for cl in layers[i]:
                xv, pm, pb, cm, cb = None, 0, 0, 0, 0
                for sign, s in cl:
                    bad = 0 if sign > 0 else 1          # the value making this literal false
                    if s[0] == 'x':
                        xv = bad
                    elif s[1] == i - 1:
                        b = 1 << pos[i - 1][s]
                        need(not pm & b, 'counter clause repeats a variable')
                        pm, pb = pm | b, pb | (b if bad else 0)
                    else:
                        b = 1 << pos[i][s]
                        need(not cm & b, 'counter clause repeats a variable')
                        cm, cb = cm | b, cb | (b if bad else 0)
                conds.append((xv, pm, pb, cm, cb))
            width = len(regs[i])
            trans, new = {}, set()
            for c, R in states:
                for x in (0, 1):
                    if (R, x) not in trans:
                        trans[(R, x)] = frozenset(
                            rc for rp in R for rc in range(1 << width)
                            if not any((xv is None or x == xv) and rp & pm == pb and rc & cm == cb
                                       for xv, pm, pb, cm, cb in conds))
                    new.add((min(c + x, k + 1), trans[(R, x)]))
            states = new
        for c, R in states:
            need(bool(R) == (c <= k), 'the counter clauses for at most', k, 'of', n, 'inputs do not express that bound')
    _COUNTER_CHECKED.add((n, k))


# ---------------------------------------------------------------- formulas

def read_cnf(path):
    """The header (number of variables) and the clauses (sorted tuples) of a DIMACS file."""
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        text = f.read()
    nv = nc = None
    clauses, cur = [], []
    for line in text.splitlines():
        s = line.strip()
        if not s or s[0] == 'c':
            continue
        if s[0] == 'p':
            parts = s.split()
            need(nv is None and len(parts) == 4 and parts[1] == 'cnf', 'cover.cnf: bad header', s)
            nv, nc = int(parts[2]), int(parts[3])
            continue
        need(nv is not None, 'cover.cnf: a clause before the header')
        for tok in s.split():
            x = int(tok)
            if x == 0:
                clauses.append(tuple(sorted(cur)))
                cur = []
            else:
                need(abs(x) <= nv, 'cover.cnf: literal', x, 'beyond the', nv, 'variables of the header')
                cur.append(x)
    need(nv is not None and not cur, 'cover.cnf: no header, or the last clause is not ended by 0')
    need(nc == len(clauses), 'cover.cnf: the header gives', nc, 'clauses, the file has', len(clauses))
    return nv, clauses


def models(variables, clauses):
    """All assignments of the variables satisfying the clauses (whose variables
    are among them), as frozensets of the true variables, by backtracking."""
    order = sorted(variables)
    at = {v: p for p, v in enumerate(order)}
    by_last = defaultdict(list)
    for cl in clauses:
        by_last[max(at[abs(l)] for l in cl)].append(cl)
    val, out = {}, []

    def rec(p):
        if p == len(order):
            out.append(frozenset(v for v in order if val[v]))
            return
        for b in (False, True):
            val[order[p]] = b
            if all(any(val[abs(l)] == (l > 0) for l in cl) for cl in by_last[p]):
                rec(p + 1)
        del val[order[p]]
    rec(0)
    return out


class Clauses:
    """Evaluate a set of clauses without auxiliary variables, together with
    counters (lits, bound, guard) taken as the bounds they express, on the
    assignment whose true variables are a given set."""

    def __init__(self, clauses, counters):
        self.by_neg = defaultdict(list)
        self.positive = []
        for cl in clauses:
            negs = [-l for l in cl if l < 0]
            (self.by_neg[min(negs)] if negs else self.positive).append(cl)
        self.counters = [(frozenset(lits), k, guard) for lits, k, guard in counters]
        self.counters_of = defaultdict(list)
        for b, (lits, k, guard) in enumerate(self.counters):
            for v in lits:
                self.counters_of[v].append(b)

    def satisfied(self, true):
        for cl in self.positive:
            if not any(l in true for l in cl):
                return False
        for v in true:
            for cl in self.by_neg.get(v, ()):
                if all(((-l) in true) if l < 0 else (l not in true) for l in cl):
                    return False
        for b in {b for v in true for b in self.counters_of.get(v, ())}:
            lits, k, guard = self.counters[b]
            if (guard is None or guard in true) and len(lits & true) > k:
                return False
        return True


TOOLS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tools')


def find_drat_trim():
    """drat-trim given by the environment variable DRAT_TRIM, or on the PATH, or
    tools/drat-trim of the repository, which is compiled from tools/drat-trim.c
    the first time it is needed if a C compiler (cc) is available."""
    exe = os.environ.get('DRAT_TRIM')
    if exe:
        need(os.path.isfile(exe) and os.access(exe, os.X_OK), 'DRAT_TRIM =', exe, 'is not an executable file')
        return exe
    exe = shutil.which('drat-trim')
    if exe:
        return exe
    local = os.path.join(TOOLS, 'drat-trim')
    if os.path.isfile(local) and os.access(local, os.X_OK):
        return local
    src = os.path.join(TOOLS, 'drat-trim.c')
    cc = shutil.which('cc') or shutil.which('clang') or shutil.which('gcc')
    if os.path.isfile(src) and cc:
        r = subprocess.run([cc, '-O2', '-o', local, src], capture_output=True, text=True)
        if r.returncode == 0 and os.path.isfile(local):
            print('       (compiled drat-trim from tools/drat-trim.c)', flush=True)
            return local
    return None


def check_drat(cnf, drat):
    """Run drat-trim on the formula and its proof; False if drat-trim is missing."""
    exe = find_drat_trim()
    if not exe:
        return False
    with tempfile.TemporaryDirectory() as tmp:
        files = []
        for path, name in ((cnf, 'cover.cnf'), (drat, 'cover.drat')):
            if path.endswith('.gz'):
                plain = os.path.join(tmp, name)
                with gzip.open(path, 'rb') as f, open(plain, 'wb') as g:
                    shutil.copyfileobj(f, g)
                path = plain
            files.append(path)
        r = subprocess.run([exe] + files, capture_output=True, text=True)
    lines = [l.strip() for l in r.stdout.splitlines()]
    need('s VERIFIED' in lines, 'drat-trim does not verify the proof:', ' | '.join(lines[-5:]))
    return True


# ---------------------------------------------------------------- the two cases

class Case:
    """The atoms, the base graph and the constraints of one item."""

    def __init__(self, item, inp):
        self.item = item
        st = STATEMENTS[item]
        need(inp.get('item') == item and inp.get('where') == st['where'], 'input.json: item or where')
        need(inp.get('vertex names') == st['vertex names'],
             'input.json: the vertex names differ from the text this check reads')
        need(inp.get('constraints of the case (the base clauses)') == st['constraints'],
             'input.json: the constraints of the case differ from the statements this check reads')
        if 'parts' in st:
            need(inp.get('parts') == st['parts'], 'input.json: the parts differ from the text this check reads')
        F = inp['the formula']
        need(F.get('one clause per recorded graph') == RECORD_CLAUSES,
             'input.json: the record clauses are described differently from the text this check reads')
        m = VARIABLES_TEXT.fullmatch(F.get('variables', ''))
        need(m, 'input.json: the description of the variables differs from the text this check reads')
        self.n_thr_stated = int(m.group(1))
        m = BASE_TEXT.fullmatch(F.get('base clauses', ''))
        need(m, 'input.json: the description of the base clauses differs from the text this check reads')
        self.n_base_stated = int(m.group(1))
        self.n_clauses_stated = F.get('clause count')
        self.files = F.get('files', {})
        need(set(self.files) == {'cover.cnf', 'cover.drat'}, 'input.json: the files of the formula')
        # the base graph
        es = [tuple(e) for e in inp['the base graph']['edges']]
        need(set(inp['the base graph']) == {'edges'}, 'input.json: the base graph')
        self.base = graph({'vertices': sorted({v for e in es for v in e}), 'edges': es})
        need({frozenset(e) for e in es} == {frozenset(e) for e in self.expected_base()},
             'the base graph is not the one described under vertex names')
        need(not any(v.startswith('@') or '#' in v for v in self.base), 'a base vertex name looks like an interior name')
        # the atoms
        self.atoms = [Atom(i, a) for i, a in enumerate(inp['atoms'])]
        self.by_id = {a.id: a for a in self.atoms}
        allv = [v for a in self.atoms for v in a.vars]
        self.n_thr = len(allv)
        need(sorted(allv) == list(range(1, self.n_thr + 1)), 'the variables of the atoms are not 1, ...,', self.n_thr)
        need(self.n_thr == self.n_thr_stated, 'input.json says auxiliary variables start after', self.n_thr_stated,
             'but the atoms have', self.n_thr, 'variables')
        keys = Counter(a.key for a in self.atoms)
        need(all(c == 1 for c in keys.values()), 'two atoms have the same form on the same attachments')
        expected = self.expected_atoms()
        need(set(keys) == set(expected), 'the atoms are not the edge and bridge types of the case:',
             len(set(keys) - set(expected)), 'unexpected,', len(set(expected) - set(keys)), 'missing')
        for a in self.atoms:
            need(a.mx == expected[a.key], 'atom', a.id, 'has most copies', a.mx, 'where the constraints give',
                 expected[a.key])
        self.setup()

    # rebuilding a recorded graph
    def vector(self, rid, counts):
        need(isinstance(counts, dict) and counts, rid, ': no atom counts')
        vec = {}
        for aid, c in counts.items():
            a = self.by_id.get(aid)
            need(a is not None, rid, ': unknown atom', aid)
            need(isinstance(c, int) and 1 <= c <= a.mx, rid, ':', c, 'copies of', aid, 'where at most', a.mx)
            vec[a.idx] = c
        return vec

    def rebuild(self, vec):
        adj = {v: set(s) for v, s in self.base.items()}
        for i, c in vec.items():
            a = self.atoms[i]
            for k in range(c):
                pre = '%s#%d' % (a.id, k)
                for x, y in a.edges:
                    X = pre + x if x.startswith('@') else x
                    Y = pre + y if y.startswith('@') else y
                    sx = adj.setdefault(X, set())
                    need(Y not in sx, 'the atoms give the edge', X, Y, 'twice')
                    sx.add(Y)
                    adj.setdefault(Y, set()).add(X)
        return adj

    def true_vars(self, vec):
        return {self.atoms[i].vars[j] for i, c in vec.items() for j in range(c)}

    def record_clause(self, kind, vec):
        """Toroidal: some atom has more copies than in the record. Nontoroidal:
        some atom of the record has fewer copies."""
        if kind == TOR:
            return tuple(sorted(a.vars[vec.get(a.idx, 0)] for a in self.atoms if vec.get(a.idx, 0) < a.mx))
        return tuple(sorted(-self.atoms[i].vars[c - 1] for i, c in vec.items()))

    def chain(self):
        return [(-a.vars[j], a.vars[j - 1]) for a in self.atoms for j in range(1, a.mx)]


class E4(Case):
    """Proposition 4.24: Psi_1 without its paths, and the 162 edge and bridge types."""
    S, T = ('s1', 's2', 's3'), ('t1', 't2', 't3')

    @staticmethod
    def partner(v):
        return ('t' if v[0] == 's' else 's') + v[1:]

    def expected_base(self):
        return [(w, s) for w in 'pq' for s in self.S] + [(w, t) for w in 'rz' for t in self.T]

    def expected_atoms(self):
        out = {('edge', frozenset((s, t))): 1 for s in self.S for t in self.T}
        for side, other in ((self.S, self.T), (self.T, self.S)):
            for c1, c2 in combinations(side, 2):
                for h in other:
                    sup = frozenset((c1, c2, h))
                    path = self.partner(h) in (c1, c2)
                    mx = 4 if path else 3
                    for key in (('one', sup, None), ('T', c1, sup), ('T', c2, sup), ('T', h, sup),
                                ('S', c1, c2, sup), ('S', c2, c1, sup)):
                        out[key] = mx
                    if path:
                        g = self.partner(h)
                        e = c2 if g == c1 else c1
                        for key in (('Q', h, g, e), ('FC1', g, e, h), ('FC2', g, e, h)):
                            out[key] = 1
        for i, j in combinations((1, 2, 3), 2):
            out[('one', frozenset(('s%d' % i, 't%d' % i, 's%d' % j, 't%d' % j)), None)] = 3
        for m in (1, 2, 3):
            x, y = [k for k in (1, 2, 3) if k != m]
            for a, b in ((x, y), (y, x)):
                out[('one', frozenset(('s%d' % m, 't%d' % m, 's%d' % a, 't%d' % b)),
                     frozenset(('s%d' % m, 't%d' % m)))] = 1
        return out

    def setup(self):
        A = self.atoms
        self.edge_atom = {a.support: a.idx for a in A if a.form == 'edge'}
        self.M = {m: self.edge_atom[frozenset(('s%d' % m, 't%d' % m))] for m in (1, 2, 3)}
        ends = {m: frozenset(('s%d' % m, 't%d' % m)) for m in (1, 2, 3)}
        for a in A:
            if a.form != 'edge':
                ms = [m for m in (1, 2, 3) if ends[m] <= a.support]
                if len(a.support) == 3:
                    need(len(ms) <= 1, 'atom', a.id)
                    a.path = ms[0] if ms else None
        # the Q, four-cycle and type (c) forms, by the m of their path
        self.path_only = {m: [] for m in (1, 2, 3)}
        for a in A:
            if a.form in ('Q', 'FC1', 'FC2'):
                need(a.path is not None, 'atom', a.id, 'is a Q or four-cycle form on a triple without s_m, t_m')
                self.path_only[a.path].append(a.idx)
            elif a.form == 'one' and a.key[2] is not None:
                (m,) = [m for m in (1, 2, 3) if ends[m] == a.key[2]]
                need(ends[m] <= a.support, 'atom', a.id, ': type (c) for a path whose ends it does not attach to')
                self.path_only[m].append(a.idx)
        self.triples = defaultdict(list)                       # attachment triple -> its atoms, in order
        for a in A:
            if a.form != 'edge' and len(a.support) == 3:
                self.triples[a.support].append(a.idx)
        self.excl3 = [(a.idx, self.edge_atom[p]) for a in A if a.form != 'edge'
                      for p in sorted(a.deg3, key=sorted) if p in self.edge_atom]
        self.opposite = [(self.edge_atom[frozenset(('s%d' % i, 't%d' % j))],
                          self.edge_atom[frozenset(('s%d' % j, 't%d' % i))]) for i, j in combinations((1, 2, 3), 2)]
        self.covers = {m: [a.idx for a in A if a.form != 'edge' and ends[m] <= a.support] for m in (1, 2, 3)}

    def base_formula(self):
        """The clauses of the constraints, as listed in input.json, and the counters."""
        A = self.atoms
        v1 = lambda i: A[i].vars[0]
        cons = []
        # Lemma 3.4(a): an atom whose degree-three interior vertex sees both ends of an edge s_i t_j excludes it
        for i, e in self.excl3:
            cons.append(tuple(sorted((-v1(i), -v1(e)))))
        # Lemma 4.7: not both s_i t_j and s_j t_i
        for e, f in self.opposite:
            cons.append(tuple(sorted((-v1(e), -v1(f)))))
        # at most three bridges on one attachment triple, or four when it contains s_m, t_m and s_m t_m is absent
        counters, blocks, top = [], [], self.n_thr
        for sup, ids in self.triples.items():
            lits = [v for i in ids for v in A[i].vars]
            m = A[ids[0]].path
            for k, guard in ((4, None), (3, v1(self.M[m]) if m is not None else None)):
                cl, top = counter_clauses(lits, k, guard, top)
                counters += cl
                blocks.append((lits, k, guard))
        # the Q, four-cycle and type (c) forms: at most one for each m, and only without s_m t_m
        for m in (1, 2, 3):
            P = self.path_only[m]
            for i in P:
                cons.append(tuple(sorted((-v1(i), -v1(self.M[m])))))
            for i, j in combinations(P, 2):
                cons.append(tuple(sorted((-v1(i), -v1(j)))))
        # for each m, s_m t_m is present or some bridge attaches to both s_m and t_m
        for m in (1, 2, 3):
            cons.append(tuple(sorted([v1(self.M[m])] + [v1(i) for i in self.covers[m]])))
        return cons, counters, blocks, top

    def allowed(self, vec):
        """The constraint list, read directly on a count vector."""
        A = self.atoms
        for i, e in self.excl3:
            if i in vec and e in vec:
                return False
        for e, f in self.opposite:
            if e in vec and f in vec:
                return False
        on_triple, path_only = Counter(), Counter()
        for i, c in vec.items():
            a = A[i]
            if a.form != 'edge' and len(a.support) == 3:
                on_triple[a.support] += c
        for m in (1, 2, 3):
            path_only[m] = sum(vec.get(i, 0) for i in self.path_only[m])
        for sup, n in on_triple.items():
            m = A[self.triples[sup][0]].path
            if n > (4 if m is not None and self.M[m] not in vec else 3):
                return False
        for m in (1, 2, 3):
            if path_only[m] > 1 or (path_only[m] and self.M[m] in vec):
                return False
            if self.M[m] not in vec and not any(i in vec for i in self.covers[m]):
                return False
        return True

    def part(self, vec):
        if all(self.M[m] in vec for m in (1, 2, 3)):
            return 'all three edges s_m t_m present'
        forms = {self.atoms[i].form for i in vec} - {'edge'}
        if forms <= {'one'}:
            return '(a)'
        tsq, fc = bool(forms & {'T', 'S', 'Q'}), bool(forms & {'FC1', 'FC2'})
        need(tsq or fc, 'no part for a vector')
        return '(b) and (c)' if tsq and fc else ('(b)' if tsq else '(c)')


class E12(Case):
    """Proposition 5.24: Psi_2 with the twins p', q' on {a, u, v}, and the bridges
    at b, c, d."""
    G, PAIRS = ('b', 'c', 'd'), (('a', 'u'), ('a', 'v'), ('u', 'v'))

    def expected_base(self):
        return ([(w, x) for w in 'pq' for x in 'abc'] + [(w, x) for w in 'rz' for x in 'abd']
                + [(w, x) for w in ("p'", "q'") for x in 'auv'])

    def expected_atoms(self):
        out = {}
        for g in self.G:
            for s, t in self.PAIRS:
                sup = frozenset((g, s, t))
                out[('one', sup, None)] = 2
                for key in (('T', g, sup), ('T', s, sup), ('T', t, sup), ('S', s, t, sup), ('S', t, s, sup),
                            ('Q', g, s, t), ('FC1', s, t, g), ('FC2', s, t, g)):
                    out[key] = 1
            for e in 'uv':
                out[('edge', frozenset((g, e)))] = 1
        return out

    def setup(self):
        for a in self.atoms:
            (a.g,) = a.support & set(self.G)
            if a.form == 'edge':
                (a.end,) = a.support - {a.g}
        self.forms = {g: [a.idx for a in self.atoms if a.g == g and a.form != 'edge'] for g in self.G}
        self.edges = {g: [a.idx for a in self.atoms if a.g == g and a.form == 'edge'] for g in self.G}

    def base_formula(self):
        A = self.atoms
        v1 = lambda i: A[i].vars[0]
        cons = []
        for g in self.G:
            # at most one form (one atom; two singleton interiors are two copies of the one-vertex atom)
            for i, j in combinations(self.forms[g], 2):
                cons.append(tuple(sorted((-v1(i), -v1(j)))))
            # a form on {g, s, t} allows only the original edges gs, gt
            for i in self.forms[g]:
                for e in self.edges[g]:
                    if A[e].end not in A[i].support:
                        cons.append(tuple(sorted((-v1(i), -v1(e)))))
            # some atom is present for g
            cons.append(tuple(sorted(v1(i) for i in self.forms[g] + self.edges[g])))
        return cons, [], [], self.n_thr

    def allowed_at(self, g, vec):
        forms = [i for i in self.forms[g] if i in vec]
        edges = [i for i in self.edges[g] if i in vec]
        if len(forms) > 1:
            return False
        if forms and any(self.atoms[e].end not in self.atoms[forms[0]].support for e in edges):
            return False
        return bool(forms or edges)

    def allowed(self, vec):
        """The constraint list, read directly on a count vector."""
        return all(self.allowed_at(g, vec) for g in self.G)

    def check_models(self, base):
        """The constraints separate by g: for each g, the models of the clauses on
        the variables of its atoms are exactly the count vectors that the
        constraint list allows (those have at most one form atom, by the first
        constraint, so it suffices to test the vectors with at most one)."""
        seen = 0
        for g in self.G:
            ids = self.forms[g] + self.edges[g]
            mine = {v for i in ids for v in self.atoms[i].vars}
            cls = [cl for cl in base if abs(cl[0]) in mine]
            need(all({abs(l) for l in cl} <= mine for cl in cls), 'a base clause joins the atoms of two vertices g')
            seen += len(cls)
            found = set()
            for true in models(mine, cls):
                vec = {}
                for i in ids:
                    c = sum(1 for v in self.atoms[i].vars if v in true)
                    need(set(self.atoms[i].vars[:c]) <= true, 'a model of the chain clauses that is not a count')
                    if c:
                        vec[i] = c
                found.add(frozenset(vec.items()))
            allowed = set()
            for f in [None] + self.forms[g]:
                for c in (range(1, self.atoms[f].mx + 1) if f is not None else [0]):
                    for r in range(len(self.edges[g]) + 1):
                        for es in combinations(self.edges[g], r):
                            vec = {e: 1 for e in es}
                            if f is not None:
                                vec[f] = c
                            if self.allowed_at(g, vec):
                                allowed.add(frozenset(vec.items()))
            need(found == allowed, 'the base clauses at', g, 'have', len(found), 'models, the constraints allow',
                 len(allowed), 'vectors, and the two sets differ')
        need(seen == len(base), 'some base clause involves no atom of b, c, d')

    def two_singletons(self, vec):
        return any(c == 2 for i, c in vec.items() if self.atoms[i].form == 'one')


def symmetries(case):
    """The permutations of the atoms induced by the relabellings of the
    attachments that map the base graph and the set of atoms onto themselves
    (keeping the most copies). Every base vertex outside the attachments has
    only attachments as neighbours, so such a relabelling extends to the base
    graph exactly when it maps the neighbourhoods of these vertices, with
    their multiplicities, onto themselves."""
    att = sorted(set().union(*(a.support for a in case.atoms)))
    others = [v for v in case.base if v not in att]
    need(all(case.base[v] <= set(att) for v in others) and all(not (case.base.get(v, set()) & set(att)) for v in att),
         'the base graph is not bipartite between attachments and other vertices')
    nbhd = Counter(frozenset(case.base[v]) for v in others)
    by_key = {a.key: a for a in case.atoms}
    out = []
    for perm in permutations(att):
        sg = dict(zip(att, perm))
        if Counter(frozenset(sg[x] for x in N) for N in nbhd.elements()) != nbhd:
            continue
        img = []
        for a in case.atoms:
            b = by_key.get(map_key(a.key, sg))
            if b is None or b.mx != a.mx:
                break
            img.append(b.idx)
        else:
            out.append(img)
    return out


# ---------------------------------------------------------------- the verdicts

def check_verdict(adj, v, cat, tally):
    """The verdict and witness of a nontoroidal record in the family."""
    need(set(v) == {'id', 'verdict', 'witness'}, v.get('id'), ': verdict fields', sorted(v))
    verdict, w = v['verdict'], v['witness']
    need(isinstance(w, dict), v['id'], ': no witness')
    t = w.get('type')
    if verdict == 'K5 minor':
        need(t == 'K5 branch sets' and set(w) == {'type', 'branch sets'}, v['id'], ': a K5 minor needs branch sets')
        check_witness(adj, w, cat)
    elif verdict == 'cut of at most two vertices':
        need(t == 'cut' and set(w) == {'type', 'vertices'}, v['id'], ': a cut is needed')
        check_witness(adj, w, cat)
    elif verdict.startswith('isomorphic to '):
        name = verdict[len('isomorphic to '):]
        need(t == 'isomorphism' and w.get('onto') == name and name in cat and set(w) == {'type', 'onto', 'map'},
             v['id'], ': an isomorphism onto', name, 'is needed')
        check_witness(adj, w, cat)
    elif verdict == 'proper nontoroidal minor':
        need(t == 'proper minor' and w.get('operation') in ('delete', 'contract')
             and set(w) <= {'type', 'operation', 'edge', 'why the minor is nontoroidal'}, v['id'], ': bad proper minor')
        e = w['edge']
        need(isinstance(e, list) and len(e) == 2 and e[0] in adj and e[1] in adj[e[0]],
             v['id'], ': the edge of the proper minor is not in the graph')
        why = w.get('why the minor is nontoroidal')
        if why is None:
            tally['left to the genus computations'] += 1
        else:
            minor = (delete_edge if w['operation'] == 'delete' else contract_edge)(adj, e)
            need(isinstance(why, dict) and why.get('type') in ('minor model', 'isomorphism'), v['id'],
                 ': only a minor model or an isomorphism of some T_i is accepted for why the minor is nontoroidal')
            check_witness(minor, why, cat)
            tally['minor with a model of some T_i'] += 1
    else:
        raise CheckFailed('%s: unexpected verdict %r' % (v['id'], verdict))


# ---------------------------------------------------------------- the item

_CATALOGUE = {}


def check(root, item):
    need(item in ('E4', 'E12'), 'not a covering item:', item)
    folder = os.path.join(root, 'enumerations', item)
    inp = load(os.path.join(folder, 'input.json'))
    case = (E4 if item == 'E4' else E12)(item, inp)
    if root not in _CATALOGUE:
        _CATALOGUE[root] = catalogue(root)
    cat = _CATALOGUE[root]

    # ---- the formula: chain, constraints and counters
    chain = case.chain()
    cons, counters, blocks, top = case.base_formula()
    for lits, k, guard in blocks:
        check_counter(len(lits), k)
    base = chain + cons + counters
    need(len(base) == case.n_base_stated, 'input.json gives', case.n_base_stated, 'base clauses, the constraints give',
         len(base))
    family = Clauses(chain + cons, blocks)
    if item == 'E12':
        need(not counters, 'E12 has no counters')
        case.check_models(chain + cons)

    # ---- the records
    verdicts = load(os.path.join(folder, 'verdicts.json'))
    need(set(verdicts) == {'item', 'verdicts'} and verdicts['item'] == item, 'verdicts.json: fields or item')
    vby = {}
    for v in verdicts['verdicts']:
        need(isinstance(v, dict) and 'id' in v and v['id'] not in vby, 'verdicts.json: a verdict without id, or twice:',
             v.get('id') if isinstance(v, dict) else v)
        vby[v['id']] = v
    cands = load(os.path.join(folder, 'candidates.json'))
    need(set(cands) == {'item', 'candidates'} and cands['item'] == item, 'candidates.json: fields or item')
    cands = cands['candidates']

    nv, cnf = read_cnf(os.path.join(folder, case.files['cover.cnf']))
    need(case.n_clauses_stated == len(cnf), 'input.json gives', case.n_clauses_stated, 'clauses, cover.cnf has', len(cnf))

    expected = Counter(base)
    counts = Counter()
    verdict_counts, found, by_part, tally = Counter(), Counter(), defaultdict(Counter), Counter()
    two_single = 0
    vectors, kinds, used = {}, {}, set()
    for rec in cands:
        rid = rec.get('id')
        need(isinstance(rid, str) and rid not in vectors, 'candidates.json: a record without id, or twice:', rid)
        kind = rec.get('kind')
        need(kind in (TOR, NONTOR), rid, ': unknown kind', kind)
        allowed_keys = {'id', 'kind', 'atoms', 'image of'}
        vec = case.vector(rid, rec.get('atoms'))
        vectors[rid], kinds[rid] = vec, kind
        expected[case.record_clause(kind, vec)] += 1
        adj = case.rebuild(vec)
        if kind == TOR:
            counts[TOR] += 1
            v = vby.get(rid)
            need(v is not None and set(v) == {'id', 'verdict', 'witness'} and v['verdict'] == 'toroidal'
                 and isinstance(v['witness'], dict) and set(v['witness']) == {'type', 'genus', 'rotation'},
                 rid, ': a toroidal record needs the verdict toroidal with an embedding')
            check_embedding(adj, v['witness'], at_most=1)
            used.add(rid)
            verdict_counts['toroidal'] += 1
        else:
            counts[NONTOR] += 1
            flag = rec.get('in the family')
            need(isinstance(flag, bool), rid, ': no in-the-family flag')
            read, clauses_say = case.allowed(vec), family.satisfied(case.true_vars(vec))
            need(read == clauses_say, rid, ': the constraints as read and the base clauses disagree on this vector')
            need(flag == read, rid, ': flagged', 'in' if flag else 'outside', 'the family, but the constraints',
                 'allow' if read else 'exclude', 'its counts')
            allowed_keys |= {'in the family'}
            if not flag:
                counts[OUTSIDE] += 1
                need(rid not in vby, rid, ': a record outside the family has a verdict')
            else:
                counts[IN_FAMILY] += 1
                allowed_keys |= {'graph', 'graph6'} | ({'part'} if item == 'E4' else set())
                g = rec.get('graph')
                need(isinstance(g, dict) and set(g) == {'vertices', 'edges'}, rid, ': no recorded graph')
                rg = graph(g)
                need(same_graph(rg, adj), rid, ': the recorded graph is not the base graph plus its atoms')
                V = g['vertices']
                pos = {x: k for k, x in enumerate(V)}
                h = from_graph6(rec.get('graph6', ''))
                need(len(h) == len(V) and all(h[k] == {pos[y] for y in rg[x]} for k, x in enumerate(V)),
                     rid, ': graph6 does not encode the recorded graph in the order of its vertices')
                v = vby.get(rid)
                need(v is not None, rid, ': a record in the family has no verdict')
                check_verdict(rg, v, cat, tally)
                used.add(rid)
                verdict_counts[v['verdict']] += 1
                if item == 'E4':
                    part = case.part(vec)
                    need(rec.get('part') == part, rid, ': part', rec.get('part'), 'where the atoms give', part)
                if v['verdict'].startswith('isomorphic to '):
                    tn = v['verdict'].split()[-1]
                    found[tn] += 1
                    if item == 'E4':
                        by_part[part][tn] += 1
                    else:
                        two_single += case.two_singletons(vec)
        need(set(rec) <= allowed_keys and {'id', 'kind', 'atoms'} <= set(rec), rid, ': unexpected fields',
             sorted(set(rec) - allowed_keys))
    need(used == set(vby), 'verdicts.json has verdicts for ids that are not toroidal records or records in the family:',
         sorted(set(vby) - used)[:5])

    # ---- 'image of': the image of another record under a symmetry of the case
    images = [rec for rec in cands if 'image of' in rec]
    if images:
        syms = symmetries(case)
        for rec in images:
            src = rec['image of']
            need(src in vectors and src != rec['id'] and kinds[src] == kinds[rec['id']],
                 rec['id'], ': image of', src, 'which is not another record of the same kind')
            sv, tv = vectors[src], vectors[rec['id']]
            need(any({p[i]: c for i, c in sv.items()} == tv for p in syms),
                 rec['id'], ': not the image of', src, 'under a symmetry of the case')

    # ---- the formula, compared clause by clause
    have = Counter(cnf)
    if have != expected:
        extra, missing = have - expected, expected - have
        raise CheckFailed('cover.cnf is not the chain, constraint and record clauses: %d clauses not expected '
                          '(e.g. %s), %d expected clauses missing (e.g. %s)'
                          % (sum(extra.values()), list(extra)[:2], sum(missing.values()), list(missing)[:2]))
    need(nv == top, 'cover.cnf declares', nv, 'variables, the clauses use', top)

    # ---- the summary
    summary = load(os.path.join(folder, 'summary.json'))
    want = {'item': item, 'where': inp['where'], 'recorded graphs': len(cands),
            TOR: counts[TOR], NONTOR: counts[NONTOR], IN_FAMILY: counts[IN_FAMILY], OUTSIDE: counts[OUTSIDE],
            'verdicts': dict(verdict_counts), 'isomorphic to T_i (copies)': dict(found),
            'formula': {'variables': nv, 'clauses': len(cnf), 'base clauses': len(base), 'DRAT': 'verified by drat-trim'}}
    if item == 'E4':
        want['T_i by part'] = {p: dict(c) for p, c in by_part.items()}
    else:
        want['copies of T_i with a two-singleton interior'] = two_single
    drat_ok = check_drat(os.path.join(folder, case.files['cover.cnf']), os.path.join(folder, case.files['cover.drat']))
    if not drat_ok:
        summary = dict(summary, formula={k: x for k, x in summary.get('formula', {}).items() if k != 'DRAT'})
        del want['formula']['DRAT']
    for k in sorted(set(want) | set(summary)):
        need(summary.get(k) == want.get(k), 'summary.json:', k, 'is', summary.get(k), 'but the records give', want.get(k))

    left = tally['left to the genus computations']
    iso = ', '.join('%s x%d' % (t, n) for t, n in sorted(found.items(), key=lambda x: int(x[0][1:])))
    return ('%d records rebuilt; %d torus embeddings; %d in the family: %d K5 minor, %d cut, %d isomorphic (%s), '
            '%d proper minors (%d left to the genus computations); cover.cnf = %d chain + %d constraint + %d counter '
            '+ %d record clauses exactly, %d variables; %d flags agree with the constraints%s; summary agrees; %s'
            % (len(cands), counts[TOR], counts[IN_FAMILY], verdict_counts['K5 minor'],
               verdict_counts['cut of at most two vertices'], sum(found.values()), iso or 'none',
               verdict_counts['proper nontoroidal minor'], left, len(chain), len(cons), len(counters), len(cands), nv,
               counts[NONTOR], '; the base clauses have exactly the models the constraints allow'
               if item == 'E12' else '', 'DRAT proof verified by drat-trim' if drat_ok else 'DRAT proof not checked'))
