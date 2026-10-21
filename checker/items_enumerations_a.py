"""The enumerations E6, E9, E10, E11, E14, E15, E16, E17, E18 and E19.

Each check reads only the item's folder (input.json, candidates.json,
verdicts.json, summary.json), the graphs T1, ..., T74 of O1, and, for one
cross-reference of E11, the summary of E9. For every item it verifies:

  * every candidate has its own id, and verdicts.json gives exactly one
    verdict for each candidate and none for anything else; records carry no
    field that this program does not check;
  * each graph is rebuilt from its lists of vertices and edges (no loops,
    repeated vertices or repeated edges), and its graph6 string encodes the
    same graph with the vertices numbered in the order listed;
  * the twins of Psi have exactly the neighbourhoods given in the item's
    vertex names;
  * the graph is the one the case describes (see each item below);
  * the witness of the verdict, by its kind:
      toroidal                     an embedding; its faces are traced and its
                                   genus is at most one
      K5 minor                     five disjoint connected branch sets with an
                                   edge between every two
      cut of at most two vertices  deleting them disconnects the graph
      isomorphic to Ti             a bijection onto Ti of O1 carrying edges
                                   exactly onto edges
      proper nontoroidal minor     the operation deletes or contracts an edge
                                   of the graph; where a minor of Ti is
                                   recorded, its steps are applied and the
                                   result is checked to be isomorphic to Ti;
                                   otherwise the minor is left to the genus
                                   computations of Appendix A, and the report
                                   counts these
      excluded by a hypothesis     the fact is checked again (in these items:
                                   B is nonplanar, with B formed as stated and
                                   a Kuratowski subgraph in it)
  * summary.json is recomputed from the records and must equal the stored
    one, and the counts stated in input.json agree with the records.

Two general tools are defined here: an isomorphism search (colour refinement
of both graphs together, then individualisation with backtracking), and a
planarity test (the algorithm of Demoucron, Malgrange and Pertuiset). Every
answer of the planarity test is turned into a certificate, a plane rotation
system or a Kuratowski subgraph (found by deleting edges while the graph stays
nonplanar), and the certificate is checked by witnesses.py before it is used.
"""
import itertools
import os
import re
from collections import Counter, deque

from witnesses import (CheckFailed, need, load, graph, add_edges, edges, n_edges, components,
                       delete_edge, contract_edge, from_graph6, genus_of_rotation, check_embedding,
                       check_k5, check_cut, check_iso, check_kuratowski, as_vertex_of, catalogue)


# ================================================================ general tools

def keys(rec, allowed, required=None, what='record'):
    """The record has only fields this program checks, and the required ones."""
    extra = set(rec) - set(allowed)
    need(not extra, what, 'has a field not checked by this program:', sorted(extra))
    missing = set(allowed if required is None else required) - set(rec)
    need(not missing, what, 'lacks the field(s)', sorted(missing))


def adjacency(vertices, edge_list):
    adj = {v: set() for v in vertices}
    for u, v in edge_list:
        need(u != v and u in adj and v in adj, 'bad edge', u, v)
        adj[u].add(v)
        adj[v].add(u)
    return adj


def edge_set(adj):
    return frozenset(frozenset((u, v)) for u in adj for v in adj[u])


def graph_key(adj):
    """A string determined by the labelled graph: its vertices and edges."""
    return '|'.join(sorted(map(str, adj))) + '#' + ';'.join(
        sorted('-'.join(sorted(map(str, e))) for e in edge_set(adj)))


def check_graph6(adj, order, s):
    """The graph6 string s encodes adj with vertex i = order[i]."""
    g = from_graph6(s)
    need(len(g) == len(order) == len(adj), 'graph6 string with', len(g), 'vertices')
    for i, x in enumerate(order):
        need({order[j] for j in g[i]} == adj[x], 'the graph6 string differs from the graph at', x)


def tsort(names):
    """T_i names sorted by their index."""
    return sorted(names, key=lambda s: int(s[1:]))


# ---------------------------------------------------------------- isomorphism search

def _refine(col, nb):
    k = len(set(col))
    while True:
        sig = [(col[v], tuple(sorted(col[u] for u in nb[v]))) for v in range(len(col))]
        ids = {s: i for i, s in enumerate(sorted(set(sig)))}
        new = [ids[s] for s in sig]
        if len(ids) == k:
            return new
        col, k = new, len(ids)


def _iso_search(col, nb, nbs, n):
    col = _refine(col, nb)
    ca = Counter(col[:n])
    if ca != Counter(col[n:]):
        return None
    big = [c for c, k in ca.items() if k > 1]
    if not big:
        where = {col[j]: j for j in range(n, 2 * n)}
        m = [where[col[i]] for i in range(n)]
        for i in range(n):
            for j in nb[i]:
                if m[j] not in nbs[m[i]]:
                    return None
        return m
    c = min(big, key=lambda c: (ca[c], c))
    x = col.index(c)
    ys = [j for j in range(n, 2 * n) if col[j] == c]
    y0 = ys[0]
    if all(nbs[y] - {y0} == nbs[y0] - {y} for y in ys):
        # every candidate image is a twin of y0: exchanging the two is an
        # automorphism of the second graph keeping the colours, so y0 suffices
        ys = [y0]
    fresh = max(col) + 1
    for y in ys:
        c2 = list(col)
        c2[x] = c2[y] = fresh
        m = _iso_search(c2, nb, nbs, n)
        if m is not None:
            return m
    return None


def find_isomorphism(A, B, colour=None):
    """An isomorphism from the graph A onto the graph B (a dictionary), or None
    if there is none. colour(v), if given, is a colour of each vertex that the
    isomorphism must keep. The search is complete; a map it returns is checked
    once more by check_iso."""
    if len(A) != len(B) or n_edges(A) != n_edges(B):
        return None
    va, vb = list(A), list(B)
    n = len(va)
    ia = {v: i for i, v in enumerate(va)}
    ib = {v: n + i for i, v in enumerate(vb)}
    nb = [[ia[u] for u in A[v]] for v in va] + [[ib[u] for u in B[v]] for v in vb]
    nbs = [frozenset(x) for x in nb]
    if colour is None:
        col = [0] * (2 * n)
    else:
        init = [colour(v) for v in va] + [colour(v) for v in vb]
        names = {c: i for i, c in enumerate(sorted(set(init), key=repr))}
        col = [names[c] for c in init]
    m = _iso_search(col, nb, nbs, n)
    if m is None:
        return None
    phi = {va[i]: vb[m[i] - n] for i in range(n)}
    check_iso(A, B, phi)
    return phi


# ---------------------------------------------------------------- planarity with certificates

def _blocks(adj):
    """The blocks of adj (maximal two-connected subgraphs, and bridges), as
    lists of edges."""
    index, low, stack, out = {}, {}, [], []
    counter = 0
    for root in sorted(adj, key=str):
        if root in index:
            continue
        index[root] = low[root] = counter
        counter += 1
        todo = [(root, None, iter(sorted(adj[root], key=str)))]
        while todo:
            v, parent, it = todo[-1]
            for w in it:
                if w == parent:
                    continue
                if w not in index:
                    stack.append((v, w))
                    index[w] = low[w] = counter
                    counter += 1
                    todo.append((w, v, iter(sorted(adj[w], key=str))))
                    break
                if index[w] < index[v]:
                    stack.append((v, w))
                    low[v] = min(low[v], index[w])
            else:
                todo.pop()
                if parent is not None:
                    low[parent] = min(low[parent], low[v])
                    if low[v] >= index[parent]:
                        blk = []
                        while True:
                            e = stack.pop()
                            blk.append(e)
                            if e == (parent, v):
                                break
                        out.append(blk)
    return out


def _a_cycle(adj):
    """A cycle of a graph of minimum degree at least two."""
    start = min(adj, key=str)
    path, pos = [start], {start: 0}
    while True:
        v = path[-1]
        prev = path[-2] if len(path) > 1 else None
        w = min((u for u in adj[v] if u != prev), key=str)
        if w in pos:
            return path[pos[w]:]
        pos[w] = len(path)
        path.append(w)


def _fragments(adj, Hv, He):
    out = []
    for u in Hv:
        for v in adj[u]:
            e = frozenset((u, v))
            if v in Hv and e not in He and str(u) < str(v):
                out.append(({u, v}, set()))
    seen = set()
    for s in adj:
        if s in Hv or s in seen:
            continue
        comp, contacts, todo = {s}, set(), [s]
        seen.add(s)
        while todo:
            x = todo.pop()
            for y in adj[x]:
                if y in Hv:
                    contacts.add(y)
                elif y not in seen:
                    seen.add(y)
                    comp.add(y)
                    todo.append(y)
        out.append((contacts, comp))
    return out


def _fragment_path(adj, contacts, inner, Hv):
    if not inner:
        return sorted(contacts, key=str)
    x = min(contacts, key=str)
    prev, queue = {}, deque()
    for y in sorted(adj[x], key=str):
        if y in inner and y not in prev:
            prev[y] = x
            queue.append(y)
    while queue:
        y = queue.popleft()
        for z in sorted(adj[y], key=str):
            if z in Hv and z != x:
                path = [z, y]
                while path[-1] != x:
                    path.append(prev[path[-1]])
                return path[::-1]
            if z in inner and z not in prev:
                prev[z] = y
                queue.append(z)
    raise CheckFailed('planarity test: a fragment with a single contact (graph not two-connected)')


def _plane_faces(adj):
    """For a two-connected graph: the faces of a plane embedding, each a cyclic
    list of vertices, oriented so that every edge is used once in each
    direction; None if the graph is nonplanar."""
    m = n_edges(adj)
    cyc = _a_cycle(adj)
    Hv = set(cyc)
    He = {frozenset((cyc[i - 1], cyc[i])) for i in range(len(cyc))}
    faces = [cyc, cyc[::-1]]
    while len(He) < m:
        fsets = [set(f) for f in faces]
        best = None
        for contacts, inner in _fragments(adj, Hv, He):
            adm = [i for i, s in enumerate(fsets) if contacts <= s]
            if not adm:
                return None
            if best is None or len(adm) < len(best[2]):
                best = (contacts, inner, adm)
                if len(adm) == 1:
                    break
        contacts, inner, adm = best
        path = _fragment_path(adj, contacts, inner, Hv)
        f = faces[adm[0]]
        i, j = f.index(path[0]), f.index(path[-1])
        if i < j:
            forward, back = f[i:j + 1], f[j:] + f[:i + 1]
        else:
            forward, back = f[i:] + f[:j + 1], f[j:i + 1]
        mid = path[1:-1]
        faces[adm[0]] = forward + mid[::-1]
        faces.append(back + mid)
        Hv.update(path)
        He.update(frozenset(path[k:k + 2]) for k in range(len(path) - 1))
    return faces


def _nonplanar_block(adj):
    """The edges of a nonplanar block of adj, or None if every block is planar."""
    for blk in _blocks(adj):
        if len(blk) >= 9 and _plane_faces(adjacency({x for e in blk for x in e}, blk)) is None:
            return blk
    return None


def planarity_certificate(adj):
    """(True, a rotation system of genus 0) or (False, the edges of a
    Kuratowski subgraph); the certificate is checked before it is returned."""
    blk = _nonplanar_block(adj)
    if blk is not None:
        # Delete edges while the graph stays nonplanar. The edges kept are each
        # needed for nonplanarity, so the result is a minimal nonplanar graph:
        # a subdivision of K5 or K3,3 by Kuratowski's theorem.
        kept, todo = [], [tuple(e) for e in blk]
        while todo:
            e = todo.pop()
            trial = kept + todo
            smaller = _nonplanar_block(adjacency({x for f in trial for x in f}, trial))
            if smaller is None:
                kept.append(e)
            else:
                s = {frozenset(f) for f in smaller}
                need(all(frozenset(f) in s for f in kept), 'planarity test: inconsistent answers')
                todo = [f for f in todo if frozenset(f) in s]
        ks = [list(e) for e in kept]
        check_kuratowski(adj, ks)
        return False, ks
    rot = {v: [] for v in adj}
    for blk in _blocks(adj):
        if len(blk) == 1:
            (u, v), = blk
            rot[u].append(v)
            rot[v].append(u)
            continue
        B = adjacency({x for e in blk for x in e}, blk)
        faces = _plane_faces(B)
        need(faces is not None, 'planarity test: inconsistent answers')
        succ = {}
        for f in faces:
            for t in range(len(f)):
                succ[(f[t], f[t - 1])] = f[(t + 1) % len(f)]
        for v in B:
            start = min(B[v], key=str)
            order, x = [start], succ[(v, start)]
            while x != start:
                order.append(x)
                x = succ[(v, x)]
            rot[v].extend(order)
    need(genus_of_rotation(adj, rot) == 0, 'planarity test: the rotation system found is not plane')
    return True, rot


# ================================================================ verdicts and facts

TORUS, K5, CUT, MINOR = 'toroidal', 'K5 minor', 'cut of at most two vertices', 'proper nontoroidal minor'
EXCLUDED = 'excluded by a hypothesis of the case'
ISO = 'isomorphic to '
B_FACT = ('B is nonplanar, where B is G minus p, q, r, z and the edges among a, b, c, d '
          '(Section 5.4 assumes B planar)')
B_SPEC = {'remove vertices': ['p', 'q', 'r', 'z'], 'remove edges within': ['a', 'b', 'c', 'd']}


def subgraph_by_spec(adj, spec):
    """G minus the vertices 'remove vertices' and the edges with both ends in
    'remove edges within'."""
    drop, within = set(spec['remove vertices']), set(spec['remove edges within'])
    return {v: {u for u in nb if u not in drop and not (u in within and v in within)}
            for v, nb in adj.items() if v not in drop}


def b_graph(adj):
    """B of Section 5.4: G minus p, q, r, z and the edges among a, b, c, d."""
    return subgraph_by_spec(adj, B_SPEC)


def same_spec(spec):
    return set(spec) == set(B_SPEC) and all(sorted(spec[k]) == B_SPEC[k] for k in B_SPEC)


def check_fact(adj, f):
    """A hypothesis of the case that the graph violates, with its certificate.
    The only fact used in these items is that B is nonplanar."""
    need(f.get('fact') == B_FACT, 'a fact not checked by this program:', f.get('fact'))
    keys(f, {'fact', 'subgraph', 'Kuratowski subgraph'}, what='the fact')
    need(same_spec(f['subgraph']), 'B is not formed as stated:', f['subgraph'])
    check_kuratowski(b_graph(adj), f['Kuratowski subgraph'])
    return f['fact']


def check_steps_to_T(adj, op, e, m, cat):
    """A minor of the proper minor that is isomorphic to some Ti: the steps
    begin with the operation of the verdict, every step deletes or contracts an
    edge, isolated vertices are then deleted, and the map is an isomorphism."""
    keys(m, {'steps', 'onto', 'map'}, what='the minor of T_i')
    steps = m['steps']
    need(steps and steps[0][0] == op and list(steps[0][1]) == list(e),
         'the steps do not begin with the operation of the verdict')
    H = adj
    for o, f in steps:
        need(o in ('delete', 'contract'), 'unknown step', o)
        H = delete_edge(H, f) if o == 'delete' else contract_edge(H, f)
    H = {x: s for x, s in H.items() if s}
    need(m['onto'] in cat, 'unknown graph', m['onto'])
    T = cat[m['onto']]
    check_iso(H, T, {x: as_vertex_of(T, y) for x, y in m['map'].items()})


def check_verdict(adj, rec, cat):
    """Check the witness of one verdict on the graph adj. Returns 'genus
    tables' for a proper minor whose nontoroidality is left to Appendix A,
    'minor of Ti' when a minor of Ti is recorded, and None otherwise."""
    v, w = rec['verdict'], rec['witness']
    t = w.get('type')
    if v == TORUS:
        keys(w, {'type', 'genus', 'rotation'}, what='the embedding')
        need(t == 'embedding', 'a toroidal verdict needs an embedding')
        check_embedding(adj, w, at_most=1)
    elif v == K5:
        keys(w, {'type', 'branch sets'}, what='the K5 witness')
        need(t == 'K5 branch sets', 'a K5 verdict needs branch sets')
        check_k5(adj, w['branch sets'])
    elif v == CUT:
        keys(w, {'type', 'vertices'}, what='the cut')
        need(t == 'cut', 'a cut verdict needs a cut')
        check_cut(adj, w['vertices'])
    elif v.startswith(ISO):
        name = v[len(ISO):]
        keys(w, {'type', 'onto', 'map'}, what='the isomorphism')
        need(t == 'isomorphism' and w['onto'] == name and name in cat, 'the isomorphism is not onto', name)
        T = cat[name]
        check_iso(adj, T, {x: as_vertex_of(T, y) for x, y in w['map'].items()})
    elif v == MINOR:
        keys(w, {'type', 'operation', 'edge', 'minor of T_i'}, {'type', 'operation', 'edge'}, 'the proper minor')
        need(t == 'proper minor' and w['operation'] in ('delete', 'contract'), 'not an edge deletion or contraction')
        e = w['edge']
        need(len(e) == 2 and e[0] in adj and e[1] in adj[e[0]], 'the edge', e, 'is not in the graph')
        if 'minor of T_i' in w:
            check_steps_to_T(adj, w['operation'], e, w['minor of T_i'], cat)
            return 'minor of Ti'
        return 'genus tables'
    elif v == EXCLUDED:
        check_fact(adj, w)
    else:
        raise CheckFailed('unknown verdict: %s' % v)
    return None


class Tally:
    """Counts of verdicts and of the minors left to the genus tables."""

    def __init__(self):
        self.verdicts = Counter()
        self.left = 0
        self.to_T = 0

    def add(self, rec, kind):
        self.verdicts[rec['verdict']] += 1
        self.left += kind == 'genus tables'
        self.to_T += kind == 'minor of Ti'

    def text(self):
        v = self.verdicts
        iso = sum(k for x, k in v.items() if x.startswith(ISO))
        parts = []
        for name, short in ((K5, 'K5 minor'), (CUT, 'cut'), (TORUS, 'toroidal'), (EXCLUDED, 'excluded')):
            if v[name]:
                parts.append('%d %s' % (v[name], short))
        if iso:
            parts.append('%d isomorphic to T_i (%s)' % (iso, ', '.join(tsort({x[len(ISO):] for x in v
                                                                           if x.startswith(ISO)}))))
        if v[MINOR]:
            s = '%d proper minor' % v[MINOR]
            if self.to_T:
                s += ' (%d with a minor of T_i checked)' % self.to_T
            s += ', %d left to the genus tables' % self.left
            parts.append(s)
        return ', '.join(parts)


# ================================================================ twins of Psi

_RULE_CACHE = {}


def check_twins(adj, rules, required):
    """rules: (regular expression for vertex names, neighbourhood). Every
    vertex whose name matches has exactly that neighbourhood, and the vertices
    in required are present."""
    for x in required:
        need(x in adj, 'the vertex', x, 'of Psi is missing')
    rules = tuple(rules)
    for v in adj:
        k = (rules, v)
        if k not in _RULE_CACHE:
            _RULE_CACHE[k] = [set(nb) for pat, nb in rules if re.fullmatch(pat, str(v))]
        for nb in _RULE_CACHE[k]:
            need(adj[v] == nb, 'the twin', v, 'has neighbourhood', sorted(adj[v]), 'not', sorted(nb))


PSI2 = [(r'p|q|p\d+', 'abc'), (r'r|z|r\d+', 'abd')]
PRIMED = [(r"p'|q'", 'acw')]


def psi_edges(n_c=2, n_d=2, tc='abc', td='abd'):
    """The twins of Psi: p, q, p3, ... on tc and r, z, r3, ... on td."""
    out = []
    for t in ['p', 'q'] + ['p%d' % k for k in range(3, n_c + 1)]:
        out += [(t, x) for x in tc]
    for t in ['r', 'z'] + ['r%d' % k for k in range(3, n_d + 1)]:
        out += [(t, x) for x in td]
    return out


# ================================================================ items given by a bridge and an exterior

def exterior(edge_list):
    vs = {x for e in edge_list for x in e}
    return adjacency(vs, edge_list)


def construction_counts(case_info, by_case, n_bridges, label):
    """The counts of 'cases of the construction' in input.json agree with the
    records: the outcomes of each case are the verdicts of its candidates (a
    candidate that also violates a hypothesis being counted as excluded), and
    the parts up to rooted isomorphism are its distinct bridges."""
    name = {'cut of order at most two': CUT, 'catalogue graph': 'isomorphic', 'K5 minor': K5,
            'toroidal': TORUS, 'proper nontoroidal minor': MINOR, 'excluded: P nonplanar': 'excluded'}
    need(set(case_info) == set(by_case), label, 'cases of the construction differ from the cases of the candidates')
    for case, info in case_info.items():
        got = Counter(by_case[case])
        want = Counter()
        for k, n in info['outcomes'].items():
            need(k in name, label, 'unknown outcome', k)
            want[name[k]] += n
        need(got == want, label, 'case', case, ': the outcomes', dict(want), 'differ from the records', dict(got))
        if 'parts up to rooted isomorphism' in info:
            need(info['parts up to rooted isomorphism'] == n_bridges[case], label, 'case', case,
                 ': parts up to rooted isomorphism', info['parts up to rooted isomorphism'], 'but',
                 n_bridges[case], 'distinct bridges')


def bridge_item(root, label, exteriors, twin_rules, twin_required, flags=None):
    """The items whose candidates are a bridge of the case together with the
    exterior of the case (E6, E9, E10, E15, E16, E17). For each candidate:
    the bridge is a subgraph of the graph; its vertices are the attachments of
    the case (the roots, or lifted roots, of its specification in input.json),
    each with a neighbour inside, and interior vertices; every edge at an
    interior vertex is in the bridge; an edge of the bridge joining two
    attachments is a permitted original edge of the case; and the graph minus
    the bridge is one of the exteriors of the case. Within a case every
    bridge occurs exactly once with each exterior."""
    item = Item(root, label)
    cat = catalogue(root)
    info = item.input['cases of the construction']
    hyp = item.input.get('hypotheses of the case used for exclusions')
    if hyp is not None:
        need(len(hyp) == 1 and all(same_spec({k: v for k, v in s.items() if k != 'meaning'}) for s in hyp.values()),
             label, 'the hypothesis used for exclusions is not B as stated in the fact')
    tally = Tally()
    also = Counter()
    pairs = {}
    by_case = {}
    t_found = {}
    cache = {}

    def ext_of(case):
        if case not in cache:
            cache[case] = [(set(X), edge_set(X)) for X in exteriors(case)]
        return cache[case]
    for c in item.cands:
        keys(c, {'id', 'case', 'graph', 'graph6', 'bridge'}, what=c.get('id'))
        case = c['case']
        need(case in info, label, c['id'], 'has an unknown case', case)
        adj = item.rebuild(c)
        check_twins(adj, twin_rules, twin_required)
        spec = info[case]['specification']
        att = set(spec['roots'] if 'roots' in spec else spec['lifted roots'])
        permitted = {frozenset(e) for e in spec.get('permitted original edges', [])}
        br = graph(c['bridge'])
        inner = set(br) - att
        need(att <= set(br), c['id'], 'the bridge lacks an attachment')
        rest = {v: set(s) for v, s in adj.items()}
        for u in br:
            need(u in adj, c['id'], 'bridge vertex', u, 'not in the graph')
            for v in br[u]:
                need(v in adj[u], c['id'], 'bridge edge', u, v, 'not in the graph')
                rest[u].discard(v)
                if u in att and v in att:
                    need(frozenset((u, v)) in permitted, c['id'], 'the bridge edge', u, v, 'joins two attachments')
        for a in att:
            need(br[a] & inner, c['id'], 'the attachment', a, 'has no neighbour inside the bridge')
        for v in inner:
            need(not rest[v], c['id'], 'an edge at', v, 'is not in the bridge')
            del rest[v]
        ext = edge_set(rest)
        found = [i for i, (V, E) in enumerate(ext_of(case)) if set(rest) == V and ext == E]
        need(len(found) == 1, label, c['id'], ': the graph minus the bridge is not the exterior of case', case)
        k = (case, found[0], graph_key(br))
        need(k not in pairs, label, c['id'], 'repeats the bridge and exterior of', pairs.get(k))
        pairs[k] = c['id']
        rec = item.verdict[c['id']]
        keys(rec, {'id', 'verdict', 'witness', 'also excluded by a hypothesis of the case'} |
             ({'case flags'} if flags else set()), {'id', 'verdict', 'witness'} | ({'case flags'} if flags else set()),
             what='the verdict of ' + c['id'])
        tally.add(rec, check_verdict(adj, rec, cat))
        outcome = rec['verdict']
        if 'also excluded by a hypothesis of the case' in rec:
            also[check_fact(adj, rec['also excluded by a hypothesis of the case'])] += 1
            outcome = 'excluded'
        elif outcome.startswith(ISO):
            outcome = 'isomorphic'
        by_case.setdefault(case, []).append(outcome)
        if rec['verdict'].startswith(ISO):
            t_found[c['id']] = (case, rec['verdict'][len(ISO):], 'also excluded by a hypothesis of the case' in rec)
        if flags:
            flags(c, rec['case flags'], br)
    # every bridge with every exterior of its case
    n_bridges = {}
    for case in by_case:
        sets = [{k[2] for k in pairs if k[0] == case and k[1] == i} for i in range(len(ext_of(case)))]
        need(all(s == sets[0] for s in sets), label, 'case', case, ': the bridges differ between its exteriors')
        n_bridges[case] = len(sets[0])
    construction_counts(info, by_case, n_bridges, label)
    # the summary
    isos = Counter(t for _, t, _ in t_found.values())
    want = {'item': label, 'where': item.input['where'], 'candidates': len(item.cands),
            'verdicts': dict(tally.verdicts),
            'isomorphic to T_i (copies)': isos,
            'graphs T_i not excluded by a hypothesis of the case': tsort({t for _, t, x in t_found.values() if not x})}
    if also:
        want['also excluded by a hypothesis of the case'] = dict(also)
    excl = tsort({t for _, t, x in t_found.values() if x})
    if excl:
        want['graphs T_i excluded by a hypothesis of the case'] = excl
    if len(by_case) > 1:
        want['by case'] = {case: {'candidates': len(by_case[case]),
                                  'graphs T_i not excluded by a hypothesis of the case':
                                      tsort({t for cs, t, x in t_found.values() if cs == case and not x})}
                           for case in by_case}
    item.compare_summary(want)
    extra = '%d bridges' % sum(n_bridges.values())
    if also:
        extra += ', %d also with B nonplanar' % sum(also.values())
    return item, tally, extra


class Item:
    """The four files of an enumeration, with the ids checked."""

    def __init__(self, root, label):
        self.label = label
        folder = os.path.join(root, 'enumerations', label)
        self.input = load(os.path.join(folder, 'input.json'))
        c = load(os.path.join(folder, 'candidates.json'))
        v = load(os.path.join(folder, 'verdicts.json'))
        self.summary = load(os.path.join(folder, 'summary.json'))
        need(set(c) == {'item', 'candidates'} and c['item'] == label, label, 'candidates.json: unexpected layout')
        need(set(v) == {'item', 'verdicts'} and v['item'] == label, label, 'verdicts.json: unexpected layout')
        need(self.input.get('item') == label and self.summary.get('item') == label, label, 'wrong item name')
        self.cands = c['candidates']
        ids = [x.get('id') for x in self.cands]
        need(all(isinstance(x, str) and x for x in ids), label, 'a candidate without an id')
        need(len(ids) == len(set(ids)), label, 'a candidate id is repeated')
        self.verdict = {}
        for r in v['verdicts']:
            need(r.get('id') not in self.verdict, label, 'two verdicts for', r.get('id'))
            self.verdict[r.get('id')] = r
        need(set(self.verdict) == set(ids), label, 'the verdicts do not correspond one to one to the candidates')

    def rebuild(self, c):
        adj = graph(c['graph'])
        check_graph6(adj, c['graph']['vertices'], c['graph6'])
        return adj

    def compare_summary(self, want):
        got = self.summary
        for k in set(got) | set(want):
            need(k in got, self.label, 'summary.json lacks', k)
            need(k in want, self.label, 'summary.json has a field not recomputed here:', k)
            need(got[k] == want[k], self.label, 'summary.json,', k, ':', got[k], 'but the records give', want[k])


# ---------------------------------------------------------------- E6, E9, E10, E15, E16, E17

def check_E6(root):
    """E6 (Proposition 4.26, Section 4.7 type (c), a = s1, b = t1, c = s2,
    d = t3, e = s3, f = t2). The exterior is Psi (p, q on a, c, e and r, z on
    b, d, f) with the edges cf and de; the case names add the edges ad and bc
    (cases with _ad, _bc) and the vertices hA, hB with the neighbourhoods given
    under 'outside' in input.json (cases with _hA, _hB). The bridge has the
    attachments a, b, c, d."""
    info = load(os.path.join(root, 'enumerations', 'E6', 'input.json'))['cases of the construction']

    def exteriors(case):
        es = psi_edges(tc='ace', td='bdf') + [('c', 'f'), ('d', 'e')]
        parts = case.split('_')[1:]
        need(case.startswith('Z') and set(parts) <= {'ad', 'bc', 'hA', 'hB'}, 'E6: unknown case', case)
        es += [tuple(x) for x in parts if x in ('ad', 'bc')]
        outside = info[case]['specification']['outside']
        need(set(outside) == {x for x in parts if x.startswith('h')}, 'E6: case', case, 'and its outside vertices differ')
        for h, nb in outside.items():
            es += [(h, x) for x in nb]
        return [exterior(es)]
    item, tally, extra = bridge_item(root, 'E6', exteriors, [(r'p|q', 'ace'), (r'r|z', 'bdf')],
                                     ['p', 'q', 'r', 'z', 'e', 'f'])
    return '%d candidates in %d cases, %s: %s' % (len(item.cands), len(info), extra, tally.text())


def check_E9(root):
    """E9 (Proposition 5.9). The bridge has the attachments a, b, c, d. Case A
    (B nonplanar): every bridge with Psi_2 having two twins on each of
    {a,b,c}, {a,b,d}, with three and two, and with three and three; case B
    (B + ab + ac + bc nonplanar): with three and two. The T_i found in case A
    all have two and two twins, the setting of the proposition (G = B u Psi_2)."""
    ext = {'A': [exterior(psi_edges(2, 2)), exterior(psi_edges(3, 2)), exterior(psi_edges(3, 3))],
           'B': [exterior(psi_edges(3, 2))]}
    item, tally, extra = bridge_item(root, 'E9', lambda case: ext[case], PSI2, ['p', 'q', 'r', 'z'])
    for c in item.cands:
        if c['case'] == 'A' and item.verdict[c['id']]['verdict'].startswith(ISO):
            need(len(c['graph']['vertices']) - len(c['bridge']['vertices']) == 4, 'E9', c['id'],
                 ': a T_i of case A with more than two twins on a triple')
    return '%d candidates, %s: %s' % (len(item.cands), extra, tally.text())


def check_E10(root):
    """E10 (Proposition 5.10): every bridge (attachments a, b, c, d) with Psi_2
    having three twins on {a,b,c} and two on {a,b,d}."""
    item, tally, extra = bridge_item(root, 'E10', lambda case: [exterior(psi_edges(3, 2))], PSI2,
                                     ['p', 'q', 'p3', 'r', 'z'])
    return '%d candidates, %s: %s' % (len(item.cands), extra, tally.text())


def primed_exterior(*more):
    return exterior(psi_edges() + [(t, x) for t in ("p'", "q'") for x in 'acw'] + list(more))


def check_E15(root):
    """E15 (Proposition 5.14, D_d without a neighbour at a): the bridge H_b on
    b, c, w; the exterior is Psi_2, the twins p', q' on {a,c,w} and the edge dw
    (the bridge of D_d on d, w). The case flags of each verdict are computed
    again with certificates: U is the bridge, and the flags say whether U + cw,
    U + cw + bw and U + cw + bc are planar; the candidate is in the item's case
    exactly when the first is planar and the other two are not."""
    planar_answers = Counter()

    def flags(c, f, U):
        keys(f, {'U+cw planar', 'U+cw+bw planar', 'U+cw+bc planar', "in the item's case"}, what='the case flags')
        got = {}
        for k, extra in (('U+cw planar', [('c', 'w')]), ('U+cw+bw planar', [('c', 'w'), ('b', 'w')]),
                         ('U+cw+bc planar', [('c', 'w'), ('b', 'c')])):
            got[k] = planarity_certificate(add_edges(U, extra))[0]
            planar_answers[got[k]] += 1
        got["in the item's case"] = got['U+cw planar'] and not got['U+cw+bw planar'] and not got['U+cw+bc planar']
        need(got == f, 'E15', c['id'], ': the case flags', f, 'but the certificates give', got)
    item, tally, extra = bridge_item(root, 'E15', lambda case: [primed_exterior(('d', 'w'))], PSI2 + PRIMED,
                                     ['p', 'q', 'r', 'z', "p'", "q'"], flags=flags)
    inside = sum(item.verdict[c['id']]['case flags']["in the item's case"] for c in item.cands)
    return '%d candidates (%d in the case by their flags, all flags certified: %d planar, %d nonplanar), %s: %s' % (
        len(item.cands), inside, planar_answers[True], planar_answers[False], extra, tally.text())


def check_E16(root):
    """E16 (Proposition 5.14, D_d with a neighbour at a): case altB, the bridge
    on b, c, w with the other side contracted; case altD, the bridge on a, d, w.
    The exterior is Psi_2, p', q' on {a,c,w}, and the contracted vertex with
    the neighbourhood its 'contracted side' gives in input.json. The records
    that also have B nonplanar carry a Kuratowski subgraph of B."""
    info = load(os.path.join(root, 'enumerations', 'E16', 'input.json'))['cases of the construction']

    def exteriors(case):
        m = re.fullmatch(r'(\S+) adjacent to (.+)', info[case]['specification']['contracted side'])
        need(m, 'E16: unreadable contracted side')
        return [primed_exterior(*[(m.group(1), x) for x in m.group(2).split(', ')])]
    item, tally, extra = bridge_item(root, 'E16', exteriors, PSI2 + PRIMED, ['p', 'q', 'r', 'z', "p'", "q'"])
    return '%d candidates, %s: %s' % (len(item.cands), extra, tally.text())


def check_E17(root):
    """E17 (Proposition 5.14, b, d in one component): the bridge H_1 on
    b, c, d, w (permitted original edges bw, dw); the exterior is Psi_2 and
    p', q' on {a,c,w}. The records that also have B nonplanar carry a
    Kuratowski subgraph of B."""
    item, tally, extra = bridge_item(root, 'E17', lambda case: [primed_exterior()], PSI2 + PRIMED,
                                     ['p', 'q', 'r', 'z', "p'", "q'"])
    return '%d candidates, %s: %s' % (len(item.cands), extra, tally.text())


# ---------------------------------------------------------------- E18

def check_E18(root):
    """E18 (Lemma 5.15): the candidates are exactly the 252 graphs of the case,
    each once: a, b, c, d, the given numbers of vertices on each of {a,b,c}
    (p, q, abc2, ...), {a,b,d} (r, z, abd2, ...), {a,c,d} (acd0, ...),
    {b,c,d} (bcd0, ...), with n_abc, n_abd at least two, at most nine in all,
    and cd absent or present. The bridge recorded is the whole graph. The
    graphs isomorphic to some Ti have no edge cd, as the text says."""
    item = Item(root, 'E18')
    cat = catalogue(root)
    (case, info), = item.input['cases of the construction'].items()
    want = {}
    for n1, n2, n3, n4 in itertools.product(range(2, 10), range(2, 10), range(10), range(10)):
        if n1 + n2 + n3 + n4 > 9:
            continue
        names = {'abc': ['p', 'q'] + ['abc%d' % i for i in range(2, n1)],
                 'abd': ['r', 'z'] + ['abd%d' % i for i in range(2, n2)],
                 'acd': ['acd%d' % i for i in range(n3)], 'bcd': ['bcd%d' % i for i in range(n4)]}
        es = [(t, x) for T, ts in names.items() for t in ts for x in T]
        for cd in (False, True):
            adj = adjacency(set('abcd') | {t for ts in names.values() for t in ts}, es + ([('c', 'd')] if cd else []))
            want[graph_key(adj)] = (n1, n2, n3, n4, cd)
    need(len(want) == info['count vectors with cd absent or present'], 'E18: the case has', len(want), 'graphs')
    tally = Tally()
    seen = {}
    outcomes = []
    t_found = Counter()
    rules = [(r'p|q', 'abc'), (r'r|z', 'abd')] + [(T + r'\d+', T) for T in ('abc', 'abd', 'acd', 'bcd')]
    for c in item.cands:
        keys(c, {'id', 'case', 'graph', 'graph6', 'bridge'}, what=c.get('id'))
        need(c['case'] == case, 'E18', c['id'], 'unknown case')
        adj = item.rebuild(c)
        need(graph(c['bridge']) == adj, 'E18', c['id'], 'the bridge is not the whole graph')
        check_twins(adj, rules, ['p', 'q', 'r', 'z'])
        k = graph_key(adj)
        need(k in want, 'E18', c['id'], 'is not a graph of the case')
        need(k not in seen, 'E18', c['id'], 'repeats', seen.get(k))
        seen[k] = c['id']
        rec = item.verdict[c['id']]
        keys(rec, {'id', 'verdict', 'witness'}, what='the verdict of ' + c['id'])
        tally.add(rec, check_verdict(adj, rec, cat))
        v = rec['verdict']
        outcomes.append('isomorphic' if v.startswith(ISO) else v)
        if v.startswith(ISO):
            t_found[v[len(ISO):]] += 1
            need('d' not in adj['c'], 'E18', c['id'], 'is isomorphic to a T_i and has the edge cd')
    need(len(seen) == len(want), 'E18: graphs of the case without a candidate')
    construction_counts({case: info}, {case: outcomes}, {case: len(seen)}, 'E18')
    item.compare_summary({'item': 'E18', 'where': item.input['where'], 'candidates': len(item.cands),
                          'verdicts': dict(tally.verdicts), 'isomorphic to T_i (copies)': dict(t_found),
                          'graphs T_i not excluded by a hypothesis of the case': tsort(t_found)})
    return 'all %d graphs of the case, each once: %s; those isomorphic to T_i have no cd' % (
        len(item.cands), tally.text())


# ---------------------------------------------------------------- E11 and E14: graphs assembled from forms

def renamed_edges(g, keep, prefix, swap=None):
    """The edges of a form, its vertices in keep named as they are (or by
    swap), the others prefixed."""
    def name(x):
        return (swap or {}).get(x, x) if x in keep else prefix + x
    return [(name(u), name(v)) for u, v in g['edges']]


def side_of(g, X):
    """A form as a graph, its vertices outside X, and their neighbours in X
    (as a word, letters in order)."""
    adj = graph(g)
    D = set(adj) - set(X)
    N = {u for v in D for u in adj[v]} & set(X)
    return adj, D, ''.join(sorted(N))


def check_E11(root):
    """E11 (Proposition 5.12, the cut X = {a, b, w}).

    Forms: the eleven forms of B_0, each with a connected interior; the forms
    of D_c, each with {c} and its interior connected, and the forms of D_1,
    with {c, d} and its interior connected, numbered in order, each having
    the neighbourhood in X its signature gives.

    Recipes: the recipes of all candidates are exactly every combination of
    the paper, each once: (a) a form of B_0, a form of D_c, a form of D_d (a
    form of D_c with d for c) and a subset of {aw, bw}; (b) the two singletons
    k1, k2 on {a,b,w} with forms of D_c and D_d; (c) k1, k2 with a form of D_1
    (its recipe giving the neighbourhood of the form). Their number is the
    number of assembled graphs. Each recipe is assembled from Psi_2 and its
    forms (interior vertices prefixed U, C, D, R); the first gives exactly the
    recorded graph, and every other is isomorphic to it. The subcases of a
    candidate are those of its recipes.

    B: each candidate carries a plane rotation system of B, or a Kuratowski
    subgraph of B, for its recorded graph (an assembled graph, by the above).
    For a candidate with B nonplanar, B is nonplanar for every graph of its
    recipes: the Kuratowski subgraph is carried over by an isomorphism (one
    keeping {a,b,c,d} and {p,q,r,z} is sought if the first does not serve), or
    else a certificate is computed. So the candidates with B planar are
    exactly the classes with a member whose B is planar. The T_i among them
    have all recipes in subcase (a), and the T_i with B nonplanar are among the
    graphs Proposition 5.9 lists for nonplanar B (case A of the summary of E9).
    """
    item = Item(root, 'E11')
    cat = catalogue(root)
    F = item.input['forms']
    keys(F, {'forms of B_0 (the eleven of the proof)', 'forms of D_c (D_d is the same with d for c)',
             'forms of D_1', 'how the forms of D_c were found', 'how the forms of D_1 were found'}, what='E11 forms')
    X = 'abw'
    B0 = {}
    for f in F['forms of B_0 (the eleven of the proof)']:
        keys(f, {'name', 'graph'}, what='a form of B_0')
        adj, D, N = side_of(f['graph'], X)
        need(N == 'abw' and D and not D & set('cd') and len(components(adj, D)) == 1,
             'E11: the form', f['name'], 'of B_0 is not a bridge on a, b, w with connected interior')
        B0[f['name']] = f['graph']
    need(len(B0) == 11, 'E11: there are not eleven forms of B_0')
    DC, D1 = F['forms of D_c (D_d is the same with d for c)'], F['forms of D_1']
    for lst, attach, core, other in ((DC, 'abcw', {'c'}, 'd'), (D1, 'acdw', {'c', 'd'}, 'b')):
        for i, f in enumerate(lst):
            keys(f, {'index', 'signature', 'origin', 'graph'}, what='a form')
            adj, D, N = side_of(f['graph'], X)
            need(f['index'] == i and core <= D and set(attach) <= set(adj) and other not in adj
                 and len(components(adj, D)) == 1,
                 'E11: the form', i, 'is not a connected side with its attachments')
            need(N == f['signature'][0], 'E11: the form', i, 'has neighbourhood', N, 'not', f['signature'][0])
    ATT = set('abcdw')

    def assemble(r):
        es = psi_edges()
        if r['case'] in ('b', 'c'):
            es += [(k, x) for k in ('k1', 'k2') for x in X]
        if r['case'] in ('a', 'b'):
            es += renamed_edges(DC[r['C form']]['graph'], ATT, 'C')
            es += renamed_edges(DC[r['D form']]['graph'], ATT, 'D', {'c': 'd'})
        if r['case'] == 'a':
            es += renamed_edges(B0[r['U']], ATT, 'U') + [tuple(e) for e in r['edges in X']]
        if r['case'] == 'c':
            es += renamed_edges(D1[r['R form']]['graph'], ATT, 'R')
        return exterior(es)
    combos = set()
    for U in B0:
        for i, j in itertools.product(range(len(DC)), repeat=2):
            for E in ((), ('aw',), ('bw',), ('aw', 'bw')):
                combos.add(('a', U, i, j, E))
    for i, j in itertools.product(range(len(DC)), repeat=2):
        combos.add(('b', i, j))
    for i in range(len(D1)):
        combos.add(('c', i))
    seen = {}
    tally = Tally()
    n_planar = 0
    t_planar, t_nonplanar = {}, {}
    carried = Counter()
    for c in item.cands:
        keys(c, {'id', 'subcases', 'recipes', 'graph', 'graph6', 'B planar', 'B nonplanar'},
             {'id', 'subcases', 'recipes', 'graph', 'graph6'}, what=c.get('id'))
        need(('B planar' in c) != ('B nonplanar' in c), 'E11', c['id'], 'needs exactly one of B planar, B nonplanar')
        adj = item.rebuild(c)
        check_twins(adj, PSI2 + [(r'k1|k2', 'abw')], ['p', 'q', 'r', 'z'])
        B = b_graph(adj)
        if 'B planar' in c:
            w = c['B planar']
            keys(w, {'type', 'genus', 'rotation'}, what='B planar')
            need(w['type'] == 'embedding', 'E11', c['id'], 'B planar needs an embedding')
            check_embedding(B, w, genus=0)
            n_planar += 1
        else:
            w = c['B nonplanar']
            keys(w, {'type', 'edges'}, what='B nonplanar')
            need(w['type'] == 'Kuratowski subgraph', 'E11', c['id'], 'B nonplanar needs a Kuratowski subgraph')
            check_kuratowski(B, w['edges'])
        need(c['recipes'], 'E11', c['id'], 'has no recipe')
        for n, r in enumerate(c['recipes']):
            case = r.get('case')
            if case == 'a':
                keys(r, {'case', 'U', 'C form', 'D form', 'edges in X'}, what='a recipe')
                E = tuple(r['edges in X'])
                k = ('a', r['U'], r['C form'], r['D form'], E)
            elif case == 'b':
                keys(r, {'case', 'C form', 'D form'}, what='a recipe')
                k = ('b', r['C form'], r['D form'])
            else:
                keys(r, {'case', 'R form', 'N(R)'}, what='a recipe')
                need(case == 'c' and 0 <= r['R form'] < len(D1), 'E11', c['id'], 'unknown recipe', r)
                need(r['N(R)'] == D1[r['R form']]['signature'][0], 'E11', c['id'], 'N(R) of', r)
                k = ('c', r['R form'])
            need(k in combos, 'E11', c['id'], 'a recipe outside the combinations of the paper:', r)
            need(k not in seen, 'E11', c['id'], 'repeats the recipe of', seen.get(k))
            seen[k] = c['id']
            G = assemble(r)
            if n == 0:
                need(G == adj, 'E11', c['id'], 'the first recipe does not give the recorded graph')
                continue
            phi = find_isomorphism(G, adj)
            need(phi is not None, 'E11', c['id'], 'the recipe', r, 'gives a graph not isomorphic to the recorded one')
            if 'B nonplanar' in c:
                inv = {y: x for x, y in phi.items()}
                BG = b_graph(G)
                ks = [[inv[x], inv[y]] for x, y in c['B nonplanar']['edges']]
                if all(y in BG.get(x, ()) for x, y in ks):
                    check_kuratowski(BG, ks)
                    carried['by the isomorphism'] += 1
                    continue
                cls = lambda v: 1 if v in 'abcd' else (2 if v in ('p', 'q', 'r', 'z') else 0)
                psi = find_isomorphism(G, adj, cls)
                if psi is not None:
                    inv = {y: x for x, y in psi.items()}
                    check_kuratowski(BG, [[inv[x], inv[y]] for x, y in c['B nonplanar']['edges']])
                    carried['by an isomorphism keeping the classes'] += 1
                    continue
                planar, cert = planarity_certificate(BG)
                need(not planar, 'E11', c['id'], ': the recipe', r, 'gives a graph whose B is planar, '
                     'but the candidate is recorded with B nonplanar')
                carried['by a computed certificate'] += 1
        cases = sorted({r['case'] for r in c['recipes']})
        need(c['subcases'] == ['(%s)' % x for x in cases], 'E11', c['id'], 'subcases', c['subcases'], 'but recipes', cases)
        rec = item.verdict[c['id']]
        keys(rec, {'id', 'verdict', 'witness'}, what='the verdict of ' + c['id'])
        tally.add(rec, check_verdict(adj, rec, cat))
        if rec['verdict'].startswith(ISO):
            (t_planar if 'B planar' in c else t_nonplanar).setdefault(rec['verdict'][len(ISO):], []).append(cases)
    need(set(seen) == combos, 'E11:', len(combos - set(seen)), 'combinations of the paper have no recipe')
    need(item.input['assembled graphs (with repetitions)'] == len(combos), 'E11: input.json, assembled graphs')
    for t, lst in t_planar.items():
        need(all(cases == ['a'] for cases in lst), 'E11:', t, 'with B planar has a recipe outside subcase (a)')
    e9 = load(os.path.join(root, 'enumerations', 'E9', 'summary.json'))
    listed = e9['by case']['A']['graphs T_i not excluded by a hypothesis of the case']
    need(set(t_nonplanar) <= set(listed), 'E11: T_i with B nonplanar not listed by Proposition 5.9:',
         sorted(set(t_nonplanar) - set(listed)))
    item.compare_summary({
        'item': 'E11', 'where': item.input['where'],
        'forms': {'B_0': len(B0), 'D_c': len(DC), 'D_1': len(D1)},
        'assembled graphs (with repetitions)': len(combos),
        'candidates (up to isomorphism)': len(item.cands), 'with B planar': n_planar,
        'verdicts': dict(tally.verdicts),
        'T_i with B planar': tsort(t_planar),
        'subcase of those': 'all (a), one component D_0',
        'T_i with B nonplanar (handled by Proposition 5.9, E9)': tsort(t_nonplanar)})
    return ('%d candidates from all %d combinations of the forms (each assembled and matched); %d with B planar '
            '(T_i: %s, all recipes in (a)); B nonplanar carried to the other recipes %s: %s' % (
                len(item.cands), len(combos), n_planar, ', '.join(tsort(t_planar)),
                ', '.join('%d %s' % (n, k) for k, n in sorted(carried.items())), tally.text()))


def e14_row(P, n1, n2):
    """The row of the table of Proposition 5.13 (up to exchanging u, v)."""
    su, sv = ('u', 'au'), ('v', 'av')
    single = lambda x, y: (x in su and y in sv) or (x in sv and y in su)
    if P == 'P1':
        if n1 in su + sv and n2 == 'uv':
            return '1'
        if single(n1, n2):
            return '2'
        if n1 == n2 == 'uv':
            return '3'
        if n1 == 'uv' and n2 in su + sv:
            return '4'
    elif P == 'P2':
        if 'uv' in (n1, n2) and {n1, n2} <= set(su + sv + ('uv',)):
            return '5'
        if single(n1, n2):
            return '6'
    return None


def check_E14(root):
    """E14 (Proposition 5.13, X = {a, u, v}, two components containing b, c, d).

    Forms: for each partition (P1: D_1 contains b, c and D_2 contains d; P2:
    D_1 contains c, d and D_2 contains b) the forms are numbered in order, each
    is its component (connected) with attachments among a, u, v, and its
    neighbourhood in X is the first entry of its signature.

    Recipes: the recipes of all candidates are exactly the pairs (D_1, D_2) of
    a partition whose neighbourhoods together contain u and v, each once, and
    their number is that of the assembled graphs. Each recipe states the
    neighbourhoods of its forms, and its row of the table is recomputed from
    the partition and the two neighbourhoods; the rows of a candidate are
    those of its recipes, and the count of candidates by row is recomputed.
    Each recipe is assembled from Psi_2, p', q' on {a,u,v} and its two forms
    (interior vertices prefixed R and D); the first gives exactly the recorded
    graph and every other is isomorphic to it. No candidate is a T_i."""
    item = Item(root, 'E14')
    cat = catalogue(root)
    inp = item.input
    X = 'auv'
    sides = {('P1', 1): 'bc', ('P1', 2): 'd', ('P2', 1): 'cd', ('P2', 2): 'b'}
    forms = {}
    for (P, k), core in sides.items():
        lst = inp['forms of D_1' if k == 1 else 'forms of D_2'][P]
        for i, f in enumerate(lst):
            keys(f, {'index', 'signature', 'construction', 'mirror image under u <-> v', 'graph'}, what='a form')
            adj, D, N = side_of(f['graph'], X)
            need(f['index'] == i and set(core) <= D and not (D & (set('abcd') - set(core)))
                 and set(X) <= set(adj) and len(components(adj, D)) == 1,
                 'E14: the form', P, k, i, 'is not a connected side')
            need(N == f['signature'][0], 'E14: the form', P, k, i, 'has neighbourhood', N, 'not', f['signature'][0])
            forms[(P, k, i)] = (f['graph'], N)
    combos = {(P, i, j) for P in ('P1', 'P2') for i in range(len(inp['forms of D_1'][P]))
              for j in range(len(inp['forms of D_2'][P]))
              if set('uv') <= set(forms[(P, 1, i)][1]) | set(forms[(P, 2, j)][1])}
    keep = set('abcduv')
    seen = {}
    by_row = Counter()
    tally = Tally()
    for c in item.cands:
        keys(c, {'id', 'graph', 'graph6', 'rows of the table', 'recipes'}, what=c.get('id'))
        adj = item.rebuild(c)
        check_twins(adj, PSI2 + [(r"p'|q'", 'auv')], ['p', 'q', 'r', 'z', "p'", "q'"])
        rows = set()
        for n, r in enumerate(c['recipes']):
            keys(r, {'partition', 'D_1 form', 'D_2 form', 'N(D_1)', 'N(D_2)', 'row'}, what='a recipe')
            P, i, j = r['partition'], r['D_1 form'], r['D_2 form']
            k = (P, i, j)
            need(k in combos, 'E14', c['id'], 'a recipe outside the assembly:', r)
            need(k not in seen, 'E14', c['id'], 'repeats the recipe of', seen.get(k))
            seen[k] = c['id']
            (g1, n1), (g2, n2) = forms[(P, 1, i)], forms[(P, 2, j)]
            need((r['N(D_1)'], r['N(D_2)']) == (n1, n2), 'E14', c['id'], 'the neighbourhoods of', r)
            row = e14_row(P, n1, n2)
            need(row is not None and row == r['row'], 'E14', c['id'], 'the recipe', r, 'is in row', row)
            rows.add(row)
            G = exterior(psi_edges() + [(t, x) for t in ("p'", "q'") for x in X] +
                         renamed_edges(g1, keep, 'R') + renamed_edges(g2, keep, 'D'))
            if n == 0:
                need(G == adj, 'E14', c['id'], 'the first recipe does not give the recorded graph')
            else:
                need(find_isomorphism(G, adj) is not None, 'E14', c['id'], 'the recipe', r,
                     'gives a graph not isomorphic to the recorded one')
        need(c['rows of the table'] == sorted(rows), 'E14', c['id'], 'rows of the table')
        for x in rows:
            by_row[x] += 1
        rec = item.verdict[c['id']]
        keys(rec, {'id', 'verdict', 'witness'}, what='the verdict of ' + c['id'])
        tally.add(rec, check_verdict(adj, rec, cat))
    need(set(seen) == combos, 'E14:', len(combos - set(seen)), 'pairs of the assembly have no recipe')
    need(inp['assembled graphs (with repetitions)'] == len(combos), 'E14: input.json, assembled graphs')
    need(set(inp['rows of the table']) == set(by_row) == set('123456'), 'E14: the rows of the table')
    obstructions = [v for v in tally.verdicts if v.startswith(ISO)]
    item.compare_summary({
        'item': 'E14', 'where': inp['where'],
        'forms of D_1': {P: len(inp['forms of D_1'][P]) for P in ('P1', 'P2')},
        'forms of D_2': {P: len(inp['forms of D_2'][P]) for P in ('P1', 'P2')},
        'assembled graphs (with repetitions)': len(combos), 'candidates (up to isomorphism)': len(item.cands),
        'verdicts': dict(tally.verdicts), 'candidates by row of the table': dict(by_row),
        'obstructions': 'none' if not obstructions else obstructions})
    return ('%d candidates from all %d pairs of forms (each assembled and matched, rows recomputed: %s): %s' % (
        len(item.cands), len(combos), ', '.join('%s: %d' % kv for kv in sorted(by_row.items())), tally.text()))


# ---------------------------------------------------------------- E19

def check_E19(root):
    """E19 (Lemma 5.17). Each case U_T_kind_v_on_S has the bridge U on the
    triple T = {a,c,d} or {b,c,d} and the vertex v with neighbourhood S, one of
    the three other triples; the record's U and N(v) agree with the case. U
    is the graph of the vertices other than a, b, c, d, p, q, r, z, v with
    their edges: its interior is connected, every vertex of T has a neighbour
    in it, and the rest of the graph is exactly Psi_2 with v.

    Nonplanar U: a Kuratowski subgraph of U. Planar U: a plane rotation system
    of U, and a Kuratowski subgraph of each of the three two-edge completions
    U + xy + xt, U + xy + yt, U + xt + yt; U is exactly the planar form of
    input.json named by the record. For each T and kind, the same bridges U
    occur in the three cases, each once: 869 for nonplanar U (bridges per
    case) and the 120 planar forms of T for planar U. Each planar form of
    input.json is on T, and the subgraph K' of its recipe (its vertices i0,
    i1, ... being the form's u0, u1, ...) is a Kuratowski subgraph of the form
    plus its added edge xy, which it uses."""
    item = Item(root, 'E19')
    cat = catalogue(root)
    inp = item.input
    triples = ('abc', 'abd', 'acd', 'bcd')
    planar_forms = inp['planar U']['forms']
    need(set(planar_forms) == {'acd', 'bcd'}, 'E19: planar forms of other triples')
    for T, d in planar_forms.items():
        keys(d, {'count', 'by added edge', 'forms'}, what='the planar forms of ' + T)
        need(d['count'] == len(d['forms']) == inp['planar U']['forms per triple'], 'E19: number of forms of', T)
        need(set(d['by added edge']) == {x + y for x, y in itertools.combinations(T, 2)}, 'E19: added edges of', T)
        for i, f in enumerate(d['forms']):
            keys(f, {'index', 'graph', 'recipe'}, what='a planar form')
            need(f['index'] == i, 'E19: form', T, i, 'out of order')
            g = graph(f['graph'])
            need(set(T) <= set(g) and not set(g) & (set('abcd') - set(T)), 'E19: form', T, i, 'is not on', T)
            r = f['recipe']
            xy = r['added edge']
            need(len(xy) == 2 and set(xy) <= set(T), 'E19: form', T, i, 'added edge', xy)
            ks = [[('u' + x[1:]) if x.startswith('i') else x for x in e] for e in r["K'"]]
            need(any(set(e) == set(xy) for e in ks), 'E19: K\' of form', T, i, 'does not use its added edge')
            check_kuratowski(add_edges(g, [xy]), ks)
    tally = Tally()
    by_case = {}
    bridges = {}
    for c in item.cands:
        planar = 'U is planar' in c
        keys(c, {'id', 'case', 'U', 'N(v)', 'graph', 'graph6'} |
             ({'U is planar', 'the three two-edge completions of U are nonplanar', 'form'} if planar
              else {'U is nonplanar'}), what=c.get('id'))
        m = re.fullmatch(r'U_(acd|bcd)_(nonplanar|planar)_v_on_(abc|abd|acd|bcd)', c['case'])
        need(m and m.group(1) != m.group(3), 'E19', c['id'], 'unknown case', c['case'])
        T, kind, S = m.groups()
        need((kind == 'planar') == planar and c['U'] == 'on {%s}, %s' % (','.join(T), kind) and c['N(v)'] == list(S),
             'E19', c['id'], 'U and N(v) do not agree with the case')
        adj = item.rebuild(c)
        check_twins(adj, PSI2[:1] + [(r'r|z', 'abd'), (r'v', S)], ['p', 'q', 'r', 'z', 'v'])
        inner = set(adj) - set('abcdpqrzv')
        U = {x: set() for x in set(T) | inner}
        rest = {x: set(nb) for x, nb in adj.items() if x not in inner}
        for x in inner:
            for y in adj[x]:
                need(y in U, 'E19', c['id'], 'an interior vertex of U is adjacent to', y)
                U[x].add(y)
                U[y].add(x)
            for y in set(T) & adj[x]:
                rest[y].discard(x)
        need(edge_set(rest) == edge_set(exterior(psi_edges() + [('v', x) for x in S])) and
             set(rest) == set('abcdpqrzv'), 'E19', c['id'], 'the graph minus U is not Psi_2 with v')
        need(inner and len(components(U, inner)) == 1 and all(U[x] & inner for x in T),
             'E19', c['id'], 'U is not a bridge on', T)
        if planar:
            w = c['U is planar']
            keys(w, {'type', 'genus', 'rotation'}, what='U is planar')
            need(w['type'] == 'embedding', 'E19', c['id'], 'U is planar needs an embedding')
            check_embedding(U, w, genus=0)
            comp = c['the three two-edge completions of U are nonplanar']
            x, y, t = T
            want = [{(x, y), (x, t)}, {(x, y), (y, t)}, {(x, t), (y, t)}]
            need(len(comp) == 3 and [{tuple(e) for e in k['added edges']} for k in comp] == want,
                 'E19', c['id'], 'not the three two-edge completions')
            for k in comp:
                keys(k, {'added edges', 'Kuratowski subgraph'}, what='a completion')
                check_kuratowski(add_edges(U, k['added edges']), k['Kuratowski subgraph'])
            f = planar_forms[T]['forms'][c['form']]
            need(f['index'] == c['form'] and graph(f['graph']) == U, 'E19', c['id'], 'U is not its form', c['form'])
        else:
            w = c['U is nonplanar']
            keys(w, {'type', 'edges'}, what='U is nonplanar')
            need(w['type'] == 'Kuratowski subgraph', 'E19', c['id'], 'U is nonplanar needs a Kuratowski subgraph')
            check_kuratowski(U, w['edges'])
        key = graph_key(U)
        b = bridges.setdefault((T, kind), {})
        need(key not in b.get(S, set()), 'E19', c['id'], 'repeats its bridge in case', c['case'])
        b.setdefault(S, set()).add(key)
        rec = item.verdict[c['id']]
        keys(rec, {'id', 'verdict', 'witness'}, what='the verdict of ' + c['id'])
        tally.add(rec, check_verdict(adj, rec, cat))
        by_case.setdefault(c['case'], Counter())[rec['verdict']] += 1
    per = {'nonplanar': inp['nonplanar U']['bridges per case'], 'planar': inp['planar U']['forms per triple']}
    need(inp['nonplanar U']['generation']['parts up to rooted isomorphism'] == per['nonplanar'], 'E19: bridges per case')
    for T in ('acd', 'bcd'):
        for kind in ('nonplanar', 'planar'):
            b = bridges.get((T, kind), {})
            need(set(b) == set(triples) - {T}, 'E19: the cases of U on', T, kind)
            need(all(s == b[S] and len(s) == per[kind] for S in b for s in b.values()), 'E19: U on', T, kind,
                 'does not have the same', per[kind], 'bridges in its three cases')
    remaining = sum(n for v, n in tally.verdicts.items() if v == TORUS or v.startswith(ISO))
    n_planar = sum(1 for c in item.cands if 'U is planar' in c)
    item.compare_summary({'item': 'E19', 'where': inp['where'], 'graphs': len(item.cands),
                          'nonplanar U': len(item.cands) - n_planar, 'planar U': n_planar,
                          'verdicts': dict(tally.verdicts),
                          'by case': {k: dict(v) for k, v in by_case.items()}, 'remaining': remaining})
    return ('%d graphs in 12 cases (%d with nonplanar U, %d with planar U and its three nonplanar two-edge '
            'completions; 240 planar forms with their K\'): %s' % (len(item.cands), len(item.cands) - n_planar,
                                                                    n_planar, tally.text()))


CHECKS = {'E6': check_E6, 'E9': check_E9, 'E10': check_E10, 'E11': check_E11, 'E14': check_E14,
          'E15': check_E15, 'E16': check_E16, 'E17': check_E17, 'E18': check_E18, 'E19': check_E19}


def check(root, item):
    return CHECKS[item](root)
