"""The enumerations E1, E2, E3, E4, E7, E8 and E12 (Sections 4 and 5).

Each folder holds input.json (the case: attachments, added edges, test graphs,
the parts of G outside the enumerated graph), candidates.json, verdicts.json
and summary.json. For every item it is checked that every candidate has a
distinct id and exactly one verdict, that every graph is rebuilt from its
record (and its graph6 string, where one is given, encodes it with the
vertices in the order listed), that every witness is right for its kind, and
that every count of summary.json agrees with the records.

Words used below. A bridge B with attachments A satisfies the hypotheses of
its case when A is independent in B, the interior B - A is nonempty and
connected, and every attachment has an interior neighbour. An allowed
operation (Section 3) deletes an interior vertex or an edge, or contracts an
edge with an interior end v into its other end u ('contract edge' [u, v]),
dropping the edges this creates between attachments; the result must again
satisfy the hypotheses. Lemma 4.3 is checked as in the docstring of the class
Psi1 below. A contraction 'contract' [u, v] of a whole graph merges v into u.

E1 (Lemma 4.11). The rest of G, H0, is the graph stated (twins p, q on
  s1, s2, s3 and r, z on t1, t2, t3, and s3t3); each bridge B meets H0 in
  A = {s1, t1, s2, t2}, satisfies the hypotheses, and B + s1s2 + t1t2 is
  nonplanar (a Kuratowski subgraph). The bridges are pairwise distinct up to
  relabelling the interior. Verdicts: a cut of G; K5 branch sets in G or in
  B + s1s2 + t1t2; an allowed operation B' of B with B' + s1s2 + t1t2
  nonplanar; an edge of G whose deletion or contraction leaves a minor model
  of T_i; or 'B has no disjoint s1-t1 and s2-t2 paths', checked by trying
  every s1-t1 path of B.
E2 (Proposition 4.22). The 48 graphs are exactly those described in
  input.json: Psi_1 (paths ab, cd, ef), the bridge B in either form, and one
  of the 16 shortened three-attachment components T(g') on the triples
  {s_i, t_i, v} meeting {e, f} (g' an end of the path s_i t_i), the 4
  four-attachment components on {a,b,e,f} or {c,d,e,f} of the form of B
  (centre on either side), or the 4 pairs of an edge af or eb and an edge cf
  or ed. Each verdict deletes an edge and the minor fails all twelve choices
  of Lemma 4.3.
E3 (Proposition 4.22). G = (Psi_1 + af) together with B on A = {a, b, c, d}
  satisfying the hypotheses, B + ac + cd + bd nonplanar, the bridges
  pairwise distinct up to relabelling the interior. Every verdict is an edge
  of G deleted or contracted, and the minor is nontoroidal because: it
  contains (Psi_1 + af) together with the allowed operation B' of B given by
  the same edge, and B' + ac + cd + bd is nonplanar; or it fails all twelve
  choices of Lemma 4.3; or it has a minor model of T25; or (three graphs) by
  the genus computation of Appendix A, where the stated degree-three vertex
  and its adjacent neighbours are checked.
E4 (Theorem 4.23). The 81 atoms are exactly the bridge and edge types the
  case allows, rebuilt here: on each of the 18 attachment triples with
  vertices on both sides, one vertex, T(r) and S(x;y), S(y;x) (x, y on one
  side), keeping those with no degree-three interior vertex whose neighbours
  are adjacent in Psi_1 together with the atom; one vertex on each set
  {s_i, t_i, s_j, t_j}; the six edges s_i t_j. Each recorded graph is the
  base graph plus the stated copies of the atoms. A toroidal record has a
  torus embedding (and the recorded passing choice of Lemma 4.3 is checked).
  A nontoroidal record fails all twelve choices, and its verdict is an
  isomorphism onto T48 or T58 (with no four-attachment atom) or an edge whose
  deletion or contraction fails all twelve choices or has a minor model of
  T_i. The base clauses of cover.cnf are rebuilt from the constraints listed
  in input.json; the sequential counter is shown, by dynamic programming, to
  allow exactly the count vectors with at most five triple atoms, and five
  only with no other atom; cover.cnf is exactly the base clauses and one
  clause per record; and cover.drat is checked with drat-trim.
E7 (Proposition 4.29), E8 (Lemma 4.30, Proposition 4.31). Each bridge B has
  the four attachments of the case and satisfies the hypotheses; the bridges
  are pairwise distinct up to relabelling the interior; the test graphs are
  those of the paper. Verdicts: a fact excluding the bridge (a Kuratowski
  subgraph of the named test, disjoint paths, a degree-three interior vertex
  with two adjacent neighbours in a test graph whose added edges are edges of
  G, or a cut of B plus all edges among the attachments); K5 branch sets in a
  test graph that is a minor of G; an allowed operation after which the main
  test is still nonplanar; or an isomorphism onto a bridge given in the text
  (attachments fixed, or the two sides interchanged), together with the
  record that B passes every test: the main test is nonplanar, the other is
  planar (an embedding), and every allowed one-step operation makes the main
  test planar (the list is exactly the allowed operations, each with a plane
  embedding). E8 also gives, for each bridge that passes, the graphs
  (Psi_1 + B) + eb and (Psi_1 + B) + af with their verdicts.
E12 (Proposition 5.13). The 31 bridges are exactly the forms described in
  input.json, up to relabelling the twins p', q' and the vertices of M;
  G = Psi_2 with (n_c, n_d) = (2, 2) together with B, where N(p') = N(q') =
  {a, u, v}; B satisfies the hypotheses of Lemma 5.3 and is planar (an
  embedding), and the four tests are nonplanar. Verdicts: an isomorphism onto
  T_i, K5 branch sets, a cut, or an edge of B whose deletion keeps the
  hypotheses of Lemma 5.3 with all four tests nonplanar.

Where a proper minor is only named by its edge, its nontoroidality is left to
the genus computation of Appendix A and is counted as such. That the
candidates are all the graphs of the construction is the business of the
programs that generate them; the counts of the generation in input.json are
not checked here.
"""
import os
import re
from collections import Counter
from itertools import combinations, product

import networkx as nx

from witnesses import (catalogue, check_cut, check_embedding, check_iso, check_k5, check_kuratowski,
                       check_minor_model, components, contract_edge, copy, delete_edge, edges, from_graph6,
                       genus_of_rotation, graph, load, n_edges, need, same_graph, vertex_of)
from items_coverings import check_drat, read_cnf

ISO = re.compile(r'isomorphic to (T\d+)')


# ---------------------------------------------------------------- common pieces

def item_files(root, item):
    d = os.path.join(root, 'enumerations', item)
    inp = load(os.path.join(d, 'input.json'))
    C = load(os.path.join(d, 'candidates.json'))
    V = load(os.path.join(d, 'verdicts.json'))
    S = load(os.path.join(d, 'summary.json'))
    for x, name in ((inp, 'input'), (C, 'candidates'), (V, 'verdicts'), (S, 'summary')):
        need(x.get('item') == item, item, ':', name + '.json is labelled', x.get('item'))
    cands = C['candidates']
    ids = [c['id'] for c in cands]
    need(len(set(ids)) == len(ids), item, ': a candidate id is repeated')
    by = {}
    for v in V['verdicts']:
        need(v['id'] not in by, item, ': two verdicts for', v['id'])
        by[v['id']] = v
    need(set(by) == set(ids), item, ': the verdicts and the candidates have different ids')
    return d, inp, cands, by, S


def plus(adj, es):
    H = copy(adj)
    for u, v in es:
        need(u != v, 'loop at', u)
        H.setdefault(u, set()).add(v)
        H.setdefault(v, set()).add(u)
    return H


def union(*gs):
    H = {}
    for g in gs:
        for v, s in g.items():
            H.setdefault(v, set()).update(s)
    return H


def from_pairs(pairs):
    H = {}
    for u, v in pairs:
        H.setdefault(u, set()).add(v)
        H.setdefault(v, set()).add(u)
    return H


def edge_set(adj):
    return {frozenset(e) for e in edges(adj)}


def subgraph_of(A, B):
    return set(A) <= set(B) and all(A[v] <= B[v] for v in A)


def check_graph6(adj, order, s, what):
    """s encodes adj with its vertices in the order listed."""
    H = from_graph6(s)
    need(len(H) == len(order) == len(adj) and set(order) == set(adj), what, ': graph6 has the wrong vertices')
    need(n_edges(H) == n_edges(adj), what, ': graph6 has the wrong number of edges')
    for i in H:
        for j in H[i]:
            need(order[j] in adj[order[i]], what, ': graph6 does not encode the graph')


def bridge_hypotheses(B, roots):
    R = set(roots)
    need(R <= set(B), 'an attachment is missing from the bridge')
    for r in R:
        need(not (B[r] & R), 'the attachments', r, 'and', sorted(B[r] & R), 'are adjacent')
        need(B[r] - R, 'the attachment', r, 'has no interior neighbour')
    inner = set(B) - R
    need(inner and len(components(B, inner)) == 1, 'the interior is empty or not connected')


def allowed_op(B, roots, op, args):
    """The result of an allowed operation of Section 3, or CheckFailed."""
    R = set(roots)
    H = copy(B)
    if op == 'delete vertex':
        (v,) = args
        need(v in H and v not in R, 'deleting', v, 'is not the deletion of an interior vertex')
        for w in H.pop(v):
            H[w].discard(v)
    elif op == 'delete edge':
        u, v = args
        need(u in H and v in H[u], 'the edge', u, v, 'is not in the bridge')
        H[u].discard(v)
        H[v].discard(u)
    elif op == 'contract edge':
        u, v = args
        need(u in H and v in H[u] and v not in R, 'contracting', u, v, 'does not merge an interior vertex into a neighbour')
        for w in H.pop(v):
            H[w].discard(v)
            if w != u:
                H[w].add(u)
                H[u].add(w)
        if u in R:
            for w in list(H[u] & R):
                H[u].discard(w)
                H[w].discard(u)
    else:
        need(False, 'unknown operation', op)
    bridge_hypotheses(H, roots)
    return H


def all_allowed_ops(B, roots):
    """Every allowed one-step operation, a contraction of two interior
    vertices taken once (either direction gives the same graph)."""
    R = set(roots)
    out = []
    cands = [('delete vertex', (v,)) for v in sorted(set(B) - R)]
    for u, v in edges(B):
        cands.append(('delete edge', (u, v)))
        if v not in R and u in R:
            cands.append(('contract edge', (u, v)))
        elif u not in R and v in R:
            cands.append(('contract edge', (v, u)))
        elif u not in R and v not in R:
            cands.append(('contract edge', (u, v)))
    for op, args in cands:
        try:
            allowed_op(B, roots, op, args)
        except AssertionError:
            continue
        out.append((op, args))
    return out


def op_key(op, args, roots):
    if op == 'contract edge' and args[0] not in roots and args[1] not in roots:
        return (op, frozenset(args))
    if op == 'delete edge':
        return (op, frozenset(args))
    return (op, tuple(args))


def has_linkage(B, st1, st2):
    """Are there vertex-disjoint s1-t1 and s2-t2 paths in B? Every s1-t1 path
    avoiding s2, t2 is tried."""
    G = nx.Graph()
    G.add_nodes_from(B)
    G.add_edges_from(edges(B))
    (s1, t1), (s2, t2) = st1, st2
    G1 = G.copy()
    G1.remove_nodes_from([s2, t2])
    for path in nx.all_simple_paths(G1, s1, t1):
        G2 = G.copy()
        G2.remove_nodes_from(path)
        if nx.has_path(G2, s2, t2):
            return True
    return False


def check_path(B, path, s, t, roots):
    need(path[0] == s and path[-1] == t, 'the path', path, 'does not join', s, t)
    need(len(set(path)) == len(path), 'the path', path, 'repeats a vertex')
    need(all(y in B[x] for x, y in zip(path, path[1:])), 'the path', path, 'is not a path of B')
    need(not set(path[1:-1]) & set(roots), 'the path', path, 'passes through an attachment')


def distinct_rooted(items, roots, what):
    """The graphs are pairwise distinct up to relabelling the vertices other
    than the roots (Weisfeiler-Lehman hashes, then isomorphism tests)."""
    R = set(roots)
    buckets = {}
    for key, B in items:
        G = nx.Graph()
        G.add_nodes_from((v, {'r': v if v in R else '.'}) for v in B)
        G.add_edges_from(edges(B))
        h = nx.weisfeiler_lehman_graph_hash(G, node_attr='r', iterations=3)
        for k2, G2 in buckets.get(h, []):
            need(not nx.is_isomorphic(G, G2, node_match=lambda p, q: p['r'] == q['r']), what, key, 'and', k2,
                 'are the same up to relabelling the interior')
        buckets.setdefault(h, []).append((key, G))


def apply_op(G, op, e):
    need(op in ('delete', 'contract'), 'unknown operation', op)
    M = delete_edge(G, tuple(e)) if op == 'delete' else contract_edge(G, tuple(e))
    return {x: s for x, s in M.items() if s}


def model_of(M, T, branch_sets):
    check_minor_model(M, T, {vertex_of(T, k): v for k, v in branch_sets.items()})


def compare(summary, mine, item):
    for k, v in mine.items():
        need(summary.get(k) == v, item, ': summary.json gives', k, '=', summary.get(k), 'but the records give', v)


def nums(c):
    return dict(sorted(c.items(), key=lambda kv: (len(kv[0]), kv[0])))


# ---------------------------------------------------------------- Lemma 4.3

class Psi1:
    """Psi_1 with paths the edges s_i t_i and twins of degree three on
    {s1, s2, s3} and on {t1, t2, t3}. A choice (k, x, y) of Lemma 4.2 has the
    faces J (attachments V(P_i) and V(P_j), auxiliary edges s_i t_i, s_j t_j,
    s_i s_j, t_i t_j), W (attachments V(P_k), s_x, t_y; auxiliary edges
    s_k t_k, s_k s_x, t_k t_y), small_S ({s_k, s_x'}; the edge s_k s_x') and
    small_T ({t_k, t_y'}; the edge t_k t_y'), where {i, j} = {1, 2, 3} - {k}
    and x', y' are the other elements of {i, j}. The Psi_1-bridges are the
    chords and the components of G - V(Psi_1) with their edges; a bridge fits
    a face whose attachments contain its own. A choice fails when some bridge
    fits no face, or when every assignment of the bridges to faces they fit
    makes the test graph of some face (its auxiliary edges and the edges of
    the bridges assigned to it) nonplanar."""
    CHOICES = [(k, x, y) for k in (1, 2, 3) for x in sorted({1, 2, 3} - {k}) for y in sorted({1, 2, 3} - {k})]

    def __init__(self, S, T, twins_s, twins_t):
        self.S, self.T = S, T
        self.twins = tuple(twins_s) + tuple(twins_t)
        self.V = set(S.values()) | set(T.values()) | set(self.twins)
        self.E = {frozenset(e) for e in [(w, S[i]) for w in twins_s for i in S] +
                  [(w, T[i]) for w in twins_t for i in T] + [(S[i], T[i]) for i in S]}

    def faces(self, k, x, y):
        S, T = self.S, self.T
        i, j = sorted({1, 2, 3} - {k})
        xp = j if x == i else i
        yp = j if y == i else i
        return {'J': ({S[i], T[i], S[j], T[j]}, [(S[i], T[i]), (S[j], T[j]), (S[i], S[j]), (T[i], T[j])]),
                'W': ({S[k], T[k], S[x], T[y]}, [(S[k], T[k]), (S[k], S[x]), (T[k], T[y])]),
                'small_S': ({S[k], S[xp]}, [(S[k], S[xp])]),
                'small_T': ({T[k], T[yp]}, [(T[k], T[yp])])}

    def check_contains(self, H):
        for e in self.E:
            u, v = tuple(e)
            need(u in H and v in H[u], 'the graph does not contain Psi_1: the edge', u, v, 'is missing')
        for w in self.twins:
            need(len(H[w]) == 3, 'the twin', w, 'does not have degree three')

    def bridges(self, H):
        out = []
        for u, v in edges(H):
            if u in self.V and v in self.V and frozenset((u, v)) not in self.E:
                out.append((frozenset((u, v)), frozenset([frozenset((u, v))])))
        rest = [v for v in H if v not in self.V]
        for C in components(H, rest):
            es = frozenset(frozenset((u, v)) for u in C for v in H[u])
            att = frozenset(v for u in C for v in H[u] if v not in C)
            out.append((att, es))
        return sorted(out, key=lambda b: (sorted(b[0]), sorted(sorted(e) for e in b[1])))

    def recorded_bridges(self, H, rec):
        self.check_contains(H)
        mine = self.bridges(H)
        theirs = sorted(((frozenset(b['attachments']), frozenset(frozenset(e) for e in b['edges'])) for b in rec),
                        key=lambda b: (sorted(b[0]), sorted(sorted(e) for e in b[1])))
        need(mine == theirs, 'the recorded Psi_1-bridges are not those of the graph')
        names = [b['name'] for b in rec]
        need(len(set(names)) == len(names), 'two bridges have the same name')
        return {b['name']: b for b in rec}

    def test_graph(self, aux, byname, asg, f):
        es = list(aux) + [tuple(e) for n in byname if asg[n] == f for e in byname[n]['edges']]
        return from_pairs(es)

    def all_choices_fail(self, H, rec):
        byname = self.recorded_bridges(H, rec['bridges'])
        names = sorted(byname)
        seen = set()
        for ch in rec['choices']:
            k, x, y = ch['choice']
            need((k, x, y) in self.CHOICES and (k, x, y) not in seen, 'choice', ch['choice'], 'unknown or repeated')
            seen.add((k, x, y))
            F = self.faces(k, x, y)
            fits = {n: [f for f, (sup, _) in F.items() if set(byname[n]['attachments']) <= sup] for n in names}
            if 'bridges with no face' in ch:
                need(ch['bridges with no face'] and all(n in fits and not fits[n] for n in ch['bridges with no face']),
                     'choice', ch['choice'], ': a bridge said to fit no face fits one')
                continue
            want = {tuple(zip(names, combo)) for combo in product(*[fits[n] for n in names])}
            got = set()
            for a in ch['assignments']:
                asg = a['faces']
                need(set(asg) == set(names) and all(asg[n] in fits[n] for n in names),
                     'choice', ch['choice'], ': an assignment puts a bridge in a face it does not fit')
                got.add(tuple((n, asg[n]) for n in names))
                f = a['nonplanar test']
                check_kuratowski(self.test_graph(F[f][1], byname, asg, f), a['Kuratowski subgraph'])
            need(got == want, 'choice', ch['choice'], ': not every assignment of the bridges is shown to fail')
        need(seen == set(self.CHOICES), 'not all twelve choices are given')

    def choice_passes(self, H, w):
        byname = self.recorded_bridges(H, w['bridges'])
        k, x, y = w['choice']
        need((k, x, y) in self.CHOICES, 'unknown choice', w['choice'])
        F = self.faces(k, x, y)
        asg = w['faces']
        need(set(asg) == set(byname), 'the passing choice does not place every bridge')
        for n, f in asg.items():
            need(f in F and set(byname[n]['attachments']) <= F[f][0], 'the bridge', n, 'does not fit the face', f)
        for f, (sup, aux) in F.items():
            G = self.test_graph(aux, byname, asg, f)
            need(genus_of_rotation(G, w['plane embeddings of the test graphs'][f]) == 0,
                 'the test graph of the face', f, 'is not embedded in the plane')


PSI1_AF = Psi1({1: 'a', 2: 'c', 3: 'e'}, {1: 'b', 2: 'd', 3: 'f'}, 'pq', 'rz')
PSI1_ST = Psi1({1: 's1', 2: 's2', 3: 's3'}, {1: 't1', 2: 't2', 3: 't3'}, 'pq', 'rz')


def psi1_edges(names):
    S, T = names.S, names.T
    return sorted(tuple(sorted(e)) for e in names.E)


# ---------------------------------------------------------------- E1

def check_E1(root):
    d, inp, cands, V, summary = item_files(root, 'E1')
    T = catalogue(root)
    A = ['s1', 't1', 's2', 't2']
    E = [('s1', 's2'), ('t1', 't2')]
    H0 = graph(inp['H0'])
    want = from_pairs([(w, s) for w in 'pq' for s in ('s1', 's2', 's3')] +
                      [(w, t) for w in 'rz' for t in ('t1', 't2', 't3')] + [('s3', 't3')])
    need(same_graph(H0, want), 'E1: H0 is not the graph of the case')
    counts, targets, ops = Counter(), Counter(), Counter()
    bridges = []
    for c in cands:
        cid = c['id']
        b = c['bridge']
        B = graph(b)
        need(set(B) & set(H0) == set(A), 'E1', cid, ': the bridge does not meet H0 in s1, t1, s2, t2')
        bridge_hypotheses(B, A)
        check_kuratowski(plus(B, E), c['B + s1s2 + t1t2 is nonplanar']['edges'])
        bridges.append((cid, B))
        G = union(H0, B)
        r = V[cid]
        v, w = r['verdict'], r['witness']
        if v == 'cut of at most two vertices':
            need(w['type'] == 'cut' and w['in'] == 'G = H0 ∪ B', 'E1', cid, ': bad cut')
            check_cut(G, w['vertices'])
        elif v == 'K5 minor':
            need(w['type'] == 'K5 branch sets', 'E1', cid, ': bad K5 witness')
            if w['in'] == 'G = H0 ∪ B':
                check_k5(G, w['branch sets'])
            else:
                need(w['in'] == 'B + s1s2 + t1t2, a minor of G', 'E1', cid, ': K5 in', w['in'])
                check_k5(plus(B, E), w['branch sets'])
        elif v == 'proper nontoroidal minor':
            if w['type'] == 'allowed reduction':
                Bp = allowed_op(B, A, w['operation'], w['arguments'])
                k = w["B' + s1s2 + t1t2 is nonplanar"]
                need(k['type'] == 'Kuratowski subgraph', 'E1', cid, ': bad certificate')
                check_kuratowski(plus(Bp, E), k['edges'])
            else:
                need(w['type'] == 'proper minor', 'E1', cid, ': unknown witness', w['type'])
                M = apply_op(G, w['operation'], w['edge'])
                why = w['why the minor is nontoroidal']
                need(why['type'] == 'minor model' and why['reason'] == 'it has a minor model of %s' % why['of'],
                     'E1', cid, ': bad minor model')
                model_of(M, T[why['of']], why['branch sets'])
                targets[why['of']] += 1
                ops[w['operation']] += 1
        elif v == 'excluded by a hypothesis of the case':
            need(w['fact'] == 'B has no disjoint s1-t1 and s2-t2 paths, but B = H1 contains P1 and P2',
                 'E1', cid, ': unknown fact')
            need(not has_linkage(B, ('s1', 't1'), ('s2', 't2')), 'E1', cid, ': B has disjoint s1-t1 and s2-t2 paths')
        else:
            need(False, 'E1', cid, ': unknown verdict', v)
        counts[v] += 1
    distinct_rooted(bridges, A, 'E1: the bridges')
    compare(summary, {'candidates': len(cands), 'verdicts': dict(counts), 'remaining': 0,
                      'proper minors that are T_i (copies)': dict(targets),
                      'one-edge minors of G with a minor model of some T_i': dict(ops)}, 'E1')
    return '%d bridges, distinct up to relabelling the interior; %s' % (len(cands), nums(counts))


# ---------------------------------------------------------------- E2

def e2_family():
    """The 48 graphs of input.json 'the graphs', as edge sets, with their kinds."""
    psi = [(w, x) for w in 'pq' for x in 'ace'] + [(w, x) for w in 'rz' for x in 'bdf'] + \
          [('a', 'b'), ('c', 'd'), ('e', 'f')]
    forms_B = {'S': [('u', 'a'), ('u', 'c')] + [(w, x) for w in ("p'", "q'") for x in ('b', 'd', 'u')],
               'T': [('u', 'b'), ('u', 'd')] + [(w, x) for w in ("p'", "q'") for x in ('a', 'c', 'u')]}
    others = []
    for trip in ('abe', 'abf', 'cde', 'cdf', 'aef', 'bef', 'cef', 'def'):
        path = [x for x in trip if {x, {'a': 'b', 'b': 'a', 'c': 'd', 'd': 'c', 'e': 'f', 'f': 'e'}[x]} <= set(trip)]
        for g in path:
            rest = [x for x in trip if x != g]
            others.append(('another three-attachment component',
                           [("u''", g)] + [(w, x) for w in ("p''", "q''") for x in rest + ["u''"]]))
    for s_side, t_side in ((('a', 'e'), ('b', 'f')), (('c', 'e'), ('d', 'f'))):
        for centre, leaves in ((s_side, t_side), (t_side, s_side)):
            others.append(('a four-attachment component without the disk drawing',
                           [("u''", x) for x in centre] + [(w, x) for w in ("p''", "q''") for x in list(leaves) + ["u''"]]))
    for e1 in (('a', 'f'), ('e', 'b')):
        for e2 in (('c', 'f'), ('e', 'd')):
            others.append(('edges s_i t_j on both index pairs {1,3} and {2,3}', [e1, e2]))
    out = []
    for fb in forms_B.values():
        for kind, es in others:
            out.append((kind, frozenset(frozenset(e) for e in psi + fb + es)))
    return out


def check_E2(root):
    d, inp, cands, V, summary = item_files(root, 'E2')
    fam = e2_family()
    need(len(fam) == 48 and len({e for _, e in fam}) == 48, 'E2: the family of the text should have 48 graphs')
    want = {e: k for k, e in fam}
    seen, kinds, counts = set(), Counter(), Counter()
    for c in cands:
        cid = c['id']
        G = graph(c['graph'])
        check_graph6(G, c['graph']['vertices'], c['graph6'], 'E2 ' + cid)
        es = frozenset(edge_set(G))
        need(es in want and es not in seen, 'E2', cid, ': not one of the graphs of the text, or repeated')
        seen.add(es)
        need(c['form']['kind'] == want[es], 'E2', cid, ': the kind', c['form']['kind'], 'is not', want[es])
        kinds[want[es]] += 1
        r = V[cid]
        w = r['witness']
        need(r['verdict'] == 'proper nontoroidal minor' and w['type'] == 'proper minor', 'E2', cid, ': verdict')
        why = w['why the minor is nontoroidal']
        need(why['reason'] == 'it fails all twelve choices of Lemma 4.3', 'E2', cid, ': reason')
        PSI1_AF.check_contains(G)
        PSI1_AF.all_choices_fail(apply_op(G, w['operation'], w['edge']), why['the twelve choices'])
        counts[r['verdict']] += 1
    need(len(seen) == 48, 'E2: not all 48 graphs are recorded')
    compare(summary, {'graphs': 48, 'by kind': dict(kinds), 'verdicts': dict(counts), 'remaining': 0}, 'E2')
    return '48 graphs, exactly those of the text; each has a one-edge minor failing all twelve choices'


# ---------------------------------------------------------------- E3

ALSO_E3 = re.compile(r'the vertex (\S+) has degree three and its neighbours (\S+), (\S+) are adjacent, so a torus '
                     r'embedding of G - (\S+) would extend to G \(proof of Lemma 3\.3\(a\)\); G is nontoroidal because '
                     r'B \+ ac \+ cd \+ bd is nonplanar')


def check_E3(root):
    d, inp, cands, V, summary = item_files(root, 'E3')
    T = catalogue(root)
    A = ['a', 'b', 'c', 'd']
    E = [('a', 'c'), ('c', 'd'), ('b', 'd')]
    X = graph(inp['the rest of G'])
    need(same_graph(X, from_pairs(psi1_edges(PSI1_AF) + [('a', 'f')])), 'E3: the rest of G is not Psi_1 + af')
    counts, ops, why_c = Counter(), Counter(), Counter()
    bridges = []
    for c in cands:
        cid = c['id']
        B = graph(c['bridge'])
        need(set(B) & set(X) == set(A), 'E3', cid, ': the bridge does not meet the rest of G in a, b, c, d')
        bridge_hypotheses(B, A)
        check_kuratowski(plus(B, E), c['B + ac + cd + bd is nonplanar']['edges'])
        bridges.append((cid, B))
        G = union(X, B)
        r = V[cid]
        w = r['witness']
        need(r['verdict'] == 'proper nontoroidal minor' and w['type'] == 'proper minor', 'E3', cid, ': verdict')
        op, e = w['operation'], w['edge']
        M = apply_op(G, op, e)
        why = w['why the minor is nontoroidal']
        reason = why['reason']
        if reason.startswith("it contains (Psi_1 ∪ B') + af"):
            Bp = allowed_op(B, A, op + ' edge', e)
            check_kuratowski(plus(Bp, E), why['Kuratowski subgraph'])
            need(subgraph_of(union(X, Bp), M), 'E3', cid, ": the minor does not contain (Psi_1 + af) and B'")
            how = 'an allowed reduction keeping B + ac + cd + bd nonplanar'
        elif reason == 'it fails all twelve choices of Lemma 4.3':
            PSI1_AF.all_choices_fail(M, why['the twelve choices'])
            how = reason.replace('it fails', 'the minor fails')
        elif reason.startswith('it has a minor model of'):
            need(why['type'] == 'minor model' and reason == 'it has a minor model of %s' % why['of'], 'E3', cid)
            model_of(M, T[why['of']], why['branch sets'])
            how = 'a minor model of %s in the minor' % why['of']
        else:
            need(reason == 'genus at least two, by the genus computations of Appendix A', 'E3', cid, ': reason', reason)
            m = ALSO_E3.fullmatch(why['also'])
            need(m, 'E3', cid, ': the note is not of the form this check reads')
            v, x, y, xy = m.groups()
            need(op == 'delete' and xy == x + y and sorted(e) == sorted([x, y]), 'E3', cid, ': the note names another edge')
            need(len(G[v]) == 3 and {x, y} <= G[v] and y in G[x], 'E3', cid, ': the note is wrong')
            how = 'a genus computation'
        ops[op] += 1
        why_c[how] += 1
        counts[r['verdict']] += 1
    distinct_rooted(bridges, A, 'E3: the bridges')
    compare(summary, {'candidates': len(cands), 'verdicts': dict(counts), 'one-edge minors': dict(ops),
                      'why nontoroidal': dict(why_c), 'remaining': 0}, 'E3')
    return '%d bridges, distinct up to relabelling the interior; one-edge minors %s (%d left to the genus computation)' % (
        len(cands), dict(ops), why_c['a genus computation'])


# ---------------------------------------------------------------- E7, E8

TESTS = {
    'E7': {'B + C': [('sx', 'tx'), ('tx', 'ty'), ('ty', 'sy'), ('sy', 'sx')],
           'B + sy sx + sx tx + tx ty': [('sy', 'sx'), ('sx', 'tx'), ('tx', 'ty')],
           'B + K4[sx,tx,sy,ty]': [('sx', 'tx'), ('sx', 'sy'), ('sx', 'ty'), ('tx', 'sy'), ('tx', 'ty'), ('sy', 'ty')],
           'B + sx tx + sy ty': [('sx', 'tx'), ('sy', 'ty')]},
    'E8': {'B + ac + cd + db': [('a', 'c'), ('c', 'd'), ('d', 'b')],
           'B + ac + bd': [('a', 'c'), ('b', 'd')],
           'B + K4[a,b,c,d]': [('a', 'b'), ('a', 'c'), ('a', 'd'), ('b', 'c'), ('b', 'd'), ('c', 'd')],
           'B + cd': [('c', 'd')],
           'Psi1 + B': [(w, x) for w in 'pq' for x in 'ace'] + [(w, x) for w in 'rz' for x in 'bdf'] +
                       [('c', 'd'), ('e', 'f')]},
}
MAIN = {'E7': 'B + C', 'E8': 'B + ac + cd + db'}               # nonplanar for B, kept by reductions
OTHER = {'E7': 'B + sy sx + sx tx + tx ty', 'E8': 'B + ac + bd'}  # planar for the bridges that pass
SUBGRAPH_OF_G = {'E7': {'B + sx tx + sy ty'}, 'E8': {'B + cd', 'Psi1 + B'}}   # tests that are minors of G
K4 = {'E7': 'B + K4[sx,tx,sy,ty]', 'E8': 'B + K4[a,b,c,d]'}
FACTS = {
    'E7': {'nonplanar': 'the path completion B + s_y s_x + s_x t_x + t_x t_y is nonplanar',
           'paths': 'B has disjoint s_x-t_y and s_y-t_x paths'},
    'E8': {'nonplanar': 'B + ac + bd is nonplanar'},
}
DEG3 = 'an interior vertex of degree three has two adjacent neighbours'
CUT = {'E7': 'B plus all edges among its four attachments has a cut of at most two vertices',
       'E8': 'B plus all edges among a, b, c, d has a cut of at most two vertices'}
FORMS = {   # the bridges given in the text: Proposition 4.29, and B1-B5 of Lemma 4.30
    'E7': {'the form of Proposition 4.29': [("p'", 'sx'), ("p'", 'sy'), ("p'", 'v'), ("q'", 'sx'), ("q'", 'sy'),
                                           ("q'", 'v'), ('v', 'tx'), ('v', 'ty')]},
    'E8': {'B1': [("p'", 'a'), ("p'", 'c'), ("p'", 'v'), ("q'", 'a'), ("q'", 'c'), ("q'", 'v'),
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
SWAP = {'E7': {'sx': 'tx', 'tx': 'sx', 'sy': 'ty', 'ty': 'sy'}, 'E8': {'a': 'b', 'b': 'a', 'c': 'd', 'd': 'c'}}
EXTERIOR = {'edge eb': ('e', 'b'), 'edge af': ('a', 'f')}


def check_passing(item, B, roots, rec, tests):
    """The record that B passes every test of the case."""
    main, other = MAIN[item], OTHER[item]
    k = rec['%s is nonplanar' % main]
    need(k['type'] == 'Kuratowski subgraph', 'bad certificate')
    check_kuratowski(plus(B, tests[main]), k['edges'])
    e = rec['%s is planar' % other]
    check_embedding(plus(B, tests[other]), e, genus=0)
    listed = {}
    for m in rec['every allowed one-step reduction makes %s planar' % main]:
        H = allowed_op(B, roots, m['operation'], m['arguments'])
        need(genus_of_rotation(plus(H, tests[main]), m['planar embedding of the test after it']) == 0,
             'the main test after', m['operation'], m['arguments'], 'is not embedded in the plane')
        key = op_key(m['operation'], m['arguments'], set(roots))
        need(key not in listed, 'an operation is listed twice')
        listed[key] = m
    allowed = {op_key(op, args, set(roots)) for op, args in all_allowed_ops(B, roots)}
    need(set(listed) == allowed, 'the listed operations are not all the allowed one-step operations')
    for x in rec.get('edge deletions keeping %s nonplanar' % main, []):
        u, w = x['deleted edge']
        try:
            allowed_op(B, roots, 'delete edge', [u, w])
            allowed_ok = True
        except AssertionError:
            allowed_ok = False
        need(not allowed_ok, 'an edge deletion keeping the test nonplanar is allowed')
        check_kuratowski(plus(delete_edge(B, (u, w)), tests[main]), x['Kuratowski subgraph'])
    return len(listed)


def check_E78(root, item):
    d, inp, cands, V, summary = item_files(root, item)
    T = catalogue(root)
    roots = inp['attachments']
    tests = {k: [tuple(e) for e in v['added edges']] for k, v in inp['test graphs'].items()}
    need(set(tests) == set(TESTS[item]) and all(edge_set(from_pairs(tests[k])) == edge_set(from_pairs(TESTS[item][k]))
                                                for k in tests), item, ': the test graphs differ from those of the text')
    forms = {k: graph(v) for k, v in inp['bridges given in the text'].items()}
    need(set(forms) == set(FORMS[item]) and all(edge_set(forms[k]) == edge_set(from_pairs(FORMS[item][k]))
                                                for k in forms), item, ': the bridges given in the text differ')
    counts, passing, whole = Counter(), Counter(), Counter()
    found = Counter()
    bridges = []
    n_ops = 0
    for c in cands:
        cid = c['id']
        b = c['bridge']
        B = graph({'vertices': sorted(set(b['vertices']) | set(roots)), 'edges': b['edges']})
        bridge_hypotheses(B, roots)
        bridges.append((cid, B))
        r = V[cid]
        v, w = r['verdict'], r['witness']
        if v == 'excluded by a hypothesis of the case':
            f = w['fact']
            if f == FACTS[item].get('nonplanar'):
                need(w['test graph'] == OTHER[item], item, cid, ': the fact names the test', w['test graph'])
                check_kuratowski(plus(B, tests[w['test graph']]), w['Kuratowski subgraph'])
            elif f == FACTS[item].get('paths'):
                p, q = w['paths']
                check_path(B, p, 'sx', 'ty', roots)
                check_path(B, q, 'sy', 'tx', roots)
                need(not set(p) & set(q), item, cid, ': the paths meet')
            elif f == DEG3:
                need(w['in'] in SUBGRAPH_OF_G[item] and w['in'] != 'Psi1 + B', item, cid, ': test', w['in'])
                H = plus(B, tests[w['in']])
                x, y = w['neighbours']
                u = w['vertex']
                need(u not in roots and len(B[u]) == 3 and {x, y} <= B[u] and y in H[x], item, cid,
                     ': the degree-three vertex and its neighbours are not as stated')
            elif f == CUT[item]:
                need(w['test graph'] == K4[item], item, cid, ': the cut is in', w['test graph'])
                check_cut(plus(B, tests[K4[item]]), w['cut'])
            else:
                need(False, item, cid, ': unknown fact', f)
        elif v == 'K5 minor':
            need(w['type'] == 'K5 branch sets' and w['in'] in SUBGRAPH_OF_G[item], item, cid, ': K5 in', w.get('in'))
            check_k5(plus(B, tests[w['in']]), w['branch sets'])
        elif v == 'proper nontoroidal minor':
            need(w['type'] == 'allowed reduction' and w['failing test after the reduction'] == MAIN[item], item, cid)
            H = allowed_op(B, roots, w['operation'], w['arguments'])
            check_kuratowski(plus(H, tests[MAIN[item]]), w['Kuratowski subgraph'])
        else:
            m = re.fullmatch(r'isomorphic to (.+)', v)
            need(m and m.group(1) in forms, item, cid, ': unknown verdict', v)
            name = m.group(1)
            need(w['type'] == 'isomorphism' and w['onto'] == name, item, cid, ': bad isomorphism')
            mp = w['map']
            check_iso(B, forms[name], mp)
            att = {x: mp[x] for x in roots}
            need(att == {x: x for x in roots} or att == SWAP[item], item, cid,
                 ': the isomorphism does not fix the attachments or interchange the two sides')
            ip = r['in the paper']
            need(ip['form'] == name and ip['A1 and A2 interchanged'] == (att == SWAP[item]) and
                 all(mp[x] == y for x, y in ip['names of the interior vertices'].items()), item, cid, ': in the paper')
            n_ops += check_passing(item, B, roots, r['also: B passes every test of the case'], tests)
            passing[name] += 1
            for g in r.get('whole graphs', []):
                need(item == 'E8', item, cid, ': whole graphs')
                G = graph(g['graph'])
                want = plus(plus(B, tests['Psi1 + B']), [EXTERIOR[g['exterior']]])
                need(same_graph(G, want), item, cid, ': the whole graph is not (Psi_1 + B) +', g['exterior'])
                check_graph6(G, g['graph']['vertices'], g['graph6'], item + ' ' + cid)
                gv, gw = g['verdict'], g['witness']
                if gv == 'proper nontoroidal minor':
                    need(gw['type'] == 'proper minor' and set(gw) == {'type', 'operation', 'edge'}, item, cid)
                    apply_op(G, gw['operation'], gw['edge'])
                else:
                    tn = ISO.fullmatch(gv).group(1)
                    need(gw['type'] == 'isomorphism' and gw['onto'] == tn, item, cid, ': bad isomorphism')
                    check_iso(G, T[tn], {k: vertex_of(T[tn], x) for k, x in gw['map'].items()})
                    found[tn] += 1
                whole[gv] += 1
        counts[v] += 1
    distinct_rooted(bridges, roots, item + ': the bridges')
    want = {'bridges': len(cands), 'verdicts': dict(counts),
            'bridges that pass, by form in the paper (labelled, so a form and its interchange count separately)':
                dict(passing)}
    if item == 'E8':
        want['whole graphs from the bridges that pass'] = dict(whole)
        want['isomorphic to T_i (copies)'] = dict(found)
    compare(summary, want, item)
    extra = '; whole graphs %s' % dict(whole) if whole else ''
    return '%d bridges, distinct up to relabelling the interior; %s; the bridges that pass are %s, with all their ' \
           'allowed one-step operations (%d in all)%s' % (len(cands), nums(counts), dict(passing), n_ops, extra)


# ---------------------------------------------------------------- E12

def e12_forms():
    base = [(w, x) for w in ("p'", "q'") for x in ('a', 'u', 'v')]
    out = []
    star = [('m', g) for g in 'bcd']
    for J in [[('m', 'u'), ('m', 'v')]] + [[('m', x), (g, y)] for x, y in (('u', 'v'), ('v', 'u')) for g in 'bcd']:
        out.append(base + star + J)
    for g2 in 'bcd':
        g1, g3 = sorted(set('bcd') - {g2})
        path = [(g1, 'm1'), ('m1', g2), (g2, 'm2'), ('m2', g3)]
        joins = [[('m1', 'u'), ('m2', 'v')], [('m1', 'v'), ('m2', 'u')]] + \
                [[('m1', x), ('m2', x), (g, y)] for x, y in (('u', 'v'), ('v', 'u')) for g in 'bcd']
        for J in joins:
            out.append(base + path + J)
    return [from_pairs(es) for es in out]


def check_E12(root):
    d, inp, cands, V, summary = item_files(root, 'E12')
    T = catalogue(root)
    A = ['a', 'b', 'c', 'd']
    tests = [('B + ac + bd', [('a', 'c'), ('b', 'd')]), ('B + ad + bc', [('a', 'd'), ('b', 'c')]),
             ('B + ab + ac + ad', [('a', 'b'), ('a', 'c'), ('a', 'd')]),
             ('B + ab + bc + bd', [('a', 'b'), ('b', 'c'), ('b', 'd')])]
    psi2 = from_pairs([(w, x) for w in 'pq' for x in 'abc'] + [(w, x) for w in 'rz' for x in 'abd'])
    counts, found = Counter(), Counter()
    # the bridges are exactly the 31 forms, up to relabelling p', q' and the vertices of M
    fixed = set('abcduv')
    forms = e12_forms()
    need(len(forms) == 31, 'E12: the text gives 31 forms')
    idx = {}
    for i, F in enumerate(forms):
        idx.setdefault(rooted_hash(F, fixed), []).append(i)
    hit = set()
    for c in cands:
        cid = c['id']
        G = graph(c['graph'])
        check_graph6(G, c['graph']['vertices'], c['graph6'], 'E12 ' + cid)
        B = graph(c['bridge B'])
        need(same_graph(G, union(B, psi2)) and set(B) & set(psi2) == set(A), 'E12', cid, ': G is not Psi_2 and B')
        need(G["p'"] == G["q'"] == {'a', 'u', 'v'}, 'E12', cid, ": N(p') or N(q') is not {a, u, v}")
        bridge_hypotheses(B, A)
        check_embedding(B, c['B is planar'], genus=0)
        got = [x['test'] for x in c['the four tests are nonplanar']]
        need(got == [t for t, _ in tests], 'E12', cid, ': the four tests')
        for x, (_, es) in zip(c['the four tests are nonplanar'], tests):
            check_kuratowski(plus(B, es), x['Kuratowski subgraph'])
        j = [i for i in idx.get(rooted_hash(B, fixed), []) if rooted_iso(B, forms[i], fixed)]
        need(len(j) == 1 and j[0] not in hit, 'E12', cid, ': the bridge is not one of the forms, or repeats one')
        hit.add(j[0])
        r = V[cid]
        v, w = r['verdict'], r['witness']
        if ISO.fullmatch(v):
            tn = ISO.fullmatch(v).group(1)
            need(w['type'] == 'isomorphism' and w['onto'] == tn, 'E12', cid)
            check_iso(G, T[tn], {k: vertex_of(T[tn], x) for k, x in w['map'].items()})
            found[tn] += 1
        elif v == 'K5 minor':
            check_k5(G, w['branch sets'])
        elif v == 'cut of at most two vertices':
            check_cut(G, w['vertices'])
        else:
            need(v == 'proper nontoroidal minor' and w['type'] == 'proper minor' and w['operation'] == 'delete', 'E12', cid)
            Bp = delete_edge(B, tuple(w['edge']))
            bridge_hypotheses(Bp, A)
            why = w['why the minor is nontoroidal']
            got = [x['test'] for x in why['tests']]
            need(got == ["B'" + t[1:] for t, _ in tests], 'E12', cid, ': the four tests of the minor')
            for x, (_, es) in zip(why['tests'], tests):
                check_kuratowski(plus(Bp, es), x['Kuratowski subgraph'])
        counts[v] += 1
    need(len(hit) == 31, 'E12: not every form is recorded')
    compare(summary, {'forms': 31, 'verdicts': dict(counts), 'isomorphic to T_i (copies)': dict(found)}, 'E12')
    return '31 bridges, exactly the forms of the text; %s' % nums(counts)


def rooted_hash(B, fixed):
    G = nx.Graph()
    G.add_nodes_from((v, {'r': v if v in fixed else '.'}) for v in B)
    G.add_edges_from(edges(B))
    return nx.weisfeiler_lehman_graph_hash(G, node_attr='r', iterations=3)


def rooted_iso(B, F, fixed):
    def nxg(H):
        G = nx.Graph()
        G.add_nodes_from((v, {'r': v if v in fixed else '.'}) for v in H)
        G.add_edges_from(edges(H))
        return G
    return nx.is_isomorphic(nxg(B), nxg(F), node_match=lambda p, q: p['r'] == q['r'])


# ---------------------------------------------------------------- E4

S4 = ('s1', 's2', 's3')
T4 = ('t1', 't2', 't3')
MATCH4 = [('s1', 't1'), ('s2', 't2'), ('s3', 't3')]
CROSS4 = [('s1', 't2'), ('s1', 't3'), ('t1', 's2'), ('t1', 's3'), ('s2', 't3'), ('t2', 's3')]
QUADS4 = [('s1', 't1', 's2', 't2'), ('s1', 't1', 's3', 't3'), ('s2', 't2', 's3', 't3')]
OPPOSITE4 = [(('s1', 't2'), ('t1', 's2')), (('s1', 't3'), ('t1', 's3')), (('s2', 't3'), ('t2', 's3'))]
BASE4 = [(w, s) for w in 'pq' for s in S4] + [(w, t) for w in 'rz' for t in T4] + MATCH4
GROUPS4 = ['atom_threshold_monotonicity', 'total_at_most_three_on_each_triple_support',
           'independent_degree_three_neighborhood', 'opposite_cross_connections_of_distinct_bridges',
           'four_root_isolation_remaining_clauses', 'global_total_triples_at_most_five_prefix_counter',
           'five_triples_forbid_quad_and_crossedge_atoms']


def e4_form_edges(form, support):
    """The edges of an atom by its description, interior vertices @h, @d1, @d2, @u, @f."""
    if form == 'one vertex':
        return [('@h', x) for x in support]
    m = re.fullmatch(r'T\((\w+)\)', form)
    if m:
        r = m.group(1)
        need(r in support, 'T(r) with r outside its support')
        others = [x for x in support if x != r]
        return [(d, x) for d in ('@d1', '@d2') for x in others + ['@u']] + [('@u', r)]
    m = re.fullmatch(r'S\((\w+);(\w+)\)', form)
    need(m, 'unknown form', form)
    x, y = m.groups()
    (z,) = [w for w in support if w not in (x, y)]
    return [(d, w) for d in ('@d1', '@d2') for w in (x, y, '@u')] + [('@f', w) for w in (x, z, '@u')]


def e4_atoms_expected():
    """(kind, support, form) of every atom the case allows."""
    out = []
    roots = S4 + T4
    match = {frozenset(e) for e in MATCH4}
    for sup in combinations(roots, 3):
        k = len(set(sup) & set(S4))
        if k not in (1, 2):
            continue
        same = sorted(set(sup) & (set(S4) if k == 2 else set(T4)))
        x, y = same
        for form in ['one vertex'] + ['T(%s)' % r for r in sup] + ['S(%s;%s)' % (x, y), 'S(%s;%s)' % (y, x)]:
            es = e4_form_edges(form, list(sup))
            adj = from_pairs(BASE4 + es)
            bad = any(len(adj[v]) == 3 and any(frozenset((a, b)) in match or b in adj[a]
                                               for a, b in combinations(sorted(adj[v]), 2))
                      for v in adj if v.startswith('@'))
            if not bad:
                out.append(('triple', frozenset(sup), form))
    for q in QUADS4:
        out.append(('quad', frozenset(q), 'one vertex'))
    for e in CROSS4:
        out.append(('crossedge', frozenset(e), 'edge'))
    return out


def e4_parse(a):
    s = a['is']
    m = re.fullmatch(r'the edge (\w\d)(\w\d)', s)
    if m:
        return 'crossedge', frozenset(m.groups()), 'edge', [m.groups()]
    m = re.fullmatch(r'(.+) on \{([^}]*)\}', s)
    need(m, 'E4: cannot read the atom', s)
    sup = [x.strip() for x in m.group(2).split(',')]
    form = m.group(1)
    kind = 'quad' if len(sup) == 4 else 'triple'
    return kind, frozenset(sup), form, e4_form_edges(form, sup)


def e4_base_clauses(atoms, var):
    """The base clauses, rebuilt from the constraints of the case, group by
    group, each clause added once (literals sorted by variable)."""
    clauses, seen, groups = [], set(), {}
    kinds = [a['kind'] for a in atoms]

    def add(cl):
        n = tuple(sorted(set(cl), key=lambda x: (abs(x), x)))
        need(not any(-x in n for x in n), 'a clause with both signs of a variable')
        if n not in seen:
            seen.add(n)
            clauses.append(n)
            return 1
        return 0

    g = 0
    for i, a in enumerate(atoms):
        for k in range(1, a['most']):
            g += add([-var[i][k], var[i][k - 1]])
    groups[GROUPS4[0]] = g

    g = 0
    supports = sorted({a['support'] for a in atoms if a['kind'] == 'triple'}, key=sorted)
    for sup in supports:
        block = [i for i, a in enumerate(atoms) if a['kind'] == 'triple' and a['support'] == sup]
        for cs in compositions(4, len(block), 3):
            g += add([-var[i][n - 1] for i, n in zip(block, cs) if n])
    groups[GROUPS4[1]] = g

    g = 0
    for i, a in enumerate(atoms):
        if a['kind'] != 'triple':
            continue
        adj = from_pairs(a['edges'])
        for j, e in enumerate(atoms):
            if e['kind'] == 'crossedge' and any(len(adj[v]) == 3 and e['support'] <= adj[v]
                                                for v in adj if v.startswith('@')):
                g += add([-var[i][0], -var[j][0]])
    groups[GROUPS4[2]] = g

    g = 0
    for i, j in combinations(range(len(atoms)), 2):
        A_, B_ = atoms[i]['support'], atoms[j]['support']
        if any((set(e) <= A_ and set(f) <= B_) or (set(f) <= A_ and set(e) <= B_) for e, f in OPPOSITE4):
            g += add([-var[i][0], -var[j][0]])
    groups[GROUPS4[3]] = g

    g = 0
    for i, q in enumerate(atoms):
        if q['kind'] == 'quad':
            for j, o in enumerate(atoms):
                if j != i and o['support'] <= q['support']:
                    g += add([-var[i][0], -var[j][0]])
    groups[GROUPS4[4]] = g

    xs = [v for i, a in enumerate(atoms) if a['kind'] == 'triple' for v in var[i]]
    top = max(v for vs in var for v in vs)
    s = {}
    nxt = top + 1
    for i in range(1, len(xs) + 1):
        for j in range(1, min(i, 5) + 1):
            s[i, j] = nxt
            nxt += 1
    g = 0
    for i, x in enumerate(xs, 1):
        g += add([-x, s[i, 1]])
        for j in range(1, min(i - 1, 5) + 1):
            g += add([-s[i - 1, j], s[i, j]])
        for j in range(2, min(i, 5) + 1):
            g += add([-x, -s[i - 1, j - 1], s[i, j]])
        if i > 5:
            g += add([-x, -s[i - 1, 5]])
    groups[GROUPS4[5]] = g
    g = 0
    for i, a in enumerate(atoms):
        if a['kind'] != 'triple':
            g += add([-s[len(xs), 5], -var[i][0]])
    groups[GROUPS4[6]] = g
    return clauses, groups, xs, s, nxt - 1


def compositions(total, length, cap):
    if length == 0:
        if total == 0:
            yield ()
        return
    for n in range(min(cap, total) + 1):
        for tail in compositions(total - n, length - 1, cap):
            yield (n,) + tail


def e4_counter_semantics(xs, s, n_other):
    """By dynamic programming over the inputs x_1, ..., x_n in order: the
    counter clauses (as rebuilt in e4_base_clauses) together with the clauses
    "not s[n,5] or not a" for the other atoms a have an assignment of the
    s-variables exactly when at most five inputs are true, and, when exactly
    five are, every such assignment makes s[n,5] true (so no other atom can be
    present). The registers of layer i are s[i,1], ..., s[i,min(i,5)]; every
    counter clause involves x_i and layers i-1 and i only."""
    n = len(xs)
    need(n >= 6, 'counter: too few inputs')
    trans = {}

    def step(R, x, i):
        """The register tuples of layer i reachable from those of layer i-1 in R."""
        key = (R, x, min(i, 6))
        if key not in trans:
            w = min(i, 5)
            out = set()
            for prev in R:
                for cur in product((0, 1), repeat=w):
                    ok = not (x and not cur[0])                                   # x_i -> s[i,1]
                    ok = ok and all(not prev[j - 1] or cur[j - 1] for j in range(1, min(i - 1, 5) + 1))
                    ok = ok and all(not (x and prev[j - 2]) or cur[j - 1] for j in range(2, w + 1))
                    ok = ok and not (i > 5 and x and prev[4])                     # no sixth true input
                    if ok:
                        out.add(cur)
            trans[key] = frozenset(out)
        return trans[key]

    states = {(0, frozenset([()]))}        # (true inputs so far, capped at 6; reachable register tuples)
    for i in range(1, n + 1):
        states = {(min(c + x, 6), step(R, x, i)) for c, R in states for x in (0, 1)}
    for c, R in states:
        if c <= 4:
            need(any(not r[4] for r in R), 'counter: a count of', c, 'is excluded, or forces s[n,5]')
        elif c == 5:
            need(R and all(r[4] for r in R), 'counter: five true inputs do not force s[n,5]')
        else:
            need(not R, 'counter: six true inputs are not excluded')
    return True


def e4_support_semantics(atoms, var, clauses):
    """Within each triple support: the clauses of group 2 with the chain
    clauses allow a count vector (each count 0..3) exactly when the total is
    at most three."""
    cs = set(clauses)
    supports = sorted({a['support'] for a in atoms if a['kind'] == 'triple'}, key=sorted)
    for sup in supports:
        block = [i for i, a in enumerate(atoms) if a['kind'] == 'triple' and a['support'] == sup]
        bvars = {v for i in block for v in var[i]}
        mine = [cl for cl in cs if all(abs(l) in bvars for l in cl)]
        for counts in product(range(4), repeat=len(block)):
            true = {var[i][k] for i, n in zip(block, counts) for k in range(n)}
            sat = all(any((l > 0) == (abs(l) in true) for l in cl) for cl in mine)
            need(sat == (sum(counts) <= 3), 'E4: the support clauses do not express "at most three"')


def check_E4(root):
    d, inp, cands, V, summary = item_files(root, 'E4')
    T = catalogue(root)
    need(edge_set(from_pairs([tuple(e) for e in inp['the base graph']['edges']])) == edge_set(from_pairs(BASE4)),
         'E4: the base graph is not Psi_1')
    # the atoms
    atoms, var = [], []
    expect = Counter(e4_atoms_expected())
    seen = Counter()
    nxt = 1
    for k, a in enumerate(inp['atoms']):
        need(a['id'] == 'A%d' % (k + 1), 'E4: atom ids')
        kind, sup, form, es = e4_parse(a)
        need(sorted(a['attachments']) == sorted(sup), 'E4: attachments of', a['id'])
        need(edge_set(from_pairs([tuple(e) for e in a['edges']])) == edge_set(from_pairs(es)), 'E4: the edges of',
             a['id'], 'are not those of', a['is'])
        most = 3 if kind == 'triple' else 1
        need(a['most copies'] == most and a['variables'] == list(range(nxt, nxt + most)), 'E4: copies or variables of',
             a['id'])
        nxt += most
        seen[(kind, sup, form)] += 1
        atoms.append({'kind': kind, 'support': sup, 'edges': [tuple(e) for e in a['edges']], 'most': most})
        var.append(a['variables'])
    need(seen == expect, 'E4: the atoms are not exactly the types of the case')
    kinds_order = [a['kind'] for a in atoms]
    need(kinds_order == sorted(kinds_order, key=['triple', 'quad', 'crossedge'].index), 'E4: atom order')
    # the base clauses
    base, groups, xs, s, top = e4_base_clauses(atoms, var)
    F = inp['the formula']
    need(groups == F['base clauses'] and len(base) == F['base clause count'], 'E4: the base clause groups', groups)
    need(xs == list(range(1, len(xs) + 1)), 'E4: the triple variables are not the first ones')
    e4_counter_semantics(xs, s, sum(1 for a in atoms if a['kind'] != 'triple'))
    e4_support_semantics(atoms, var, base)
    # the records
    ids = {('A%d' % (i + 1)): i for i in range(len(atoms))}
    rec_clauses = []
    counts, found, ops_twelve, ops_model = Counter(), Counter(), Counter(), Counter()
    n_tor = n_non = 0
    for c in cands:
        cid = c['id']
        n = [0] * len(atoms)
        for k, m in c['atoms'].items():
            need(k in ids and 1 <= m <= atoms[ids[k]]['most'], 'E4', cid, ': atom count', k, m)
            n[ids[k]] = m
        G = from_pairs(BASE4)
        for i, m in enumerate(n):
            for cp in range(m):
                for x, y in atoms[i]['edges']:
                    nm = lambda z: 'A%d#%d@%s' % (i + 1, cp, z[1:]) if z.startswith('@') else z
                    G.setdefault(nm(x), set()).add(nm(y))
                    G.setdefault(nm(y), set()).add(nm(x))
        H = graph(c['graph'])
        need(same_graph(G, H), 'E4', cid, ': the graph is not the base graph with the stated atoms')
        check_graph6(H, c['graph']['vertices'], c['graph6'], 'E4 ' + cid)
        r = V[cid]
        v, w = r['verdict'], r['witness']
        if c['kind'] == 'toroidal, maximal':
            need(v == 'toroidal', 'E4', cid, ': a toroidal record with verdict', v)
            check_embedding(H, w, at_most=1)
            pc = r['also: a passing choice of Lemma 4.3']
            need(pc['type'] == 'passing choice of Lemma 4.3', 'E4', cid)
            PSI1_ST.choice_passes(H, pc)
            rec_clauses.append(tuple(sorted(var[i][n[i]] for i in range(len(n)) if n[i] < atoms[i]['most'])))
            n_tor += 1
        else:
            need(c['kind'] == 'nontoroidal, minimal', 'E4', cid, ': kind', c['kind'])
            nt = r['nontoroidal']
            need(nt['reason'] == 'it fails all twelve choices of Lemma 4.3', 'E4', cid)
            PSI1_ST.all_choices_fail(H, nt['the twelve choices'])
            if ISO.fullmatch(v):
                tn = ISO.fullmatch(v).group(1)
                need(w['type'] == 'isomorphism' and w['onto'] == tn, 'E4', cid)
                check_iso(H, T[tn], {k: vertex_of(T[tn], x) for k, x in w['map'].items()})
                need(not any(n[i] for i, a in enumerate(atoms) if a['kind'] == 'quad'), 'E4', cid,
                     ': a copy of', tn, 'with a four-attachment bridge')
                found[tn] += 1
            else:
                need(v == 'proper nontoroidal minor' and w['type'] == 'proper minor', 'E4', cid, ': verdict', v)
                M = apply_op(H, w['operation'], w['edge'])
                why = w['why the minor is nontoroidal']
                if why['reason'] == 'it fails all twelve choices of Lemma 4.3':
                    PSI1_ST.all_choices_fail(M, why['the twelve choices'])
                    ops_twelve[w['operation']] += 1
                else:
                    need(why['type'] == 'minor model' and why['reason'] == 'it has a minor model of %s' % why['of'],
                         'E4', cid)
                    model_of(M, T[why['of']], why['branch sets'])
                    ops_model[w['operation']] += 1
            rec_clauses.append(tuple(sorted(-var[i][n[i] - 1] for i in range(len(n)) if n[i] > 0)))
            n_non += 1
        counts[v] += 1
    # the formula
    files = F['files']
    cnf_path = os.path.join(d, files['cover.cnf'])
    nv, cnf = read_cnf(cnf_path)
    need(nv == top, 'E4: cover.cnf has', nv, 'variables, the base clauses use', top)
    need(Counter(cnf) == Counter(tuple(sorted(c)) for c in base) + Counter(rec_clauses),
         'E4: cover.cnf is not the base clauses and one clause per record')
    need(len(cnf) == F['clause count'], 'E4: clause count')
    drat_ok = check_drat(cnf_path, os.path.join(d, files['cover.drat']))
    want = {'recorded graphs': len(cands), 'toroidal, maximal': n_tor, 'nontoroidal, minimal': n_non,
            'verdicts': dict(counts), 'isomorphic to T_i (copies)': dict(found),
            'four-attachment bridges in the copies of T48 and T58': 0,
            'one-edge minors failing all twelve choices': dict(ops_twelve),
            'one-edge minors with a minor model of some T_i': dict(ops_model)}
    compare(summary, want, 'E4')
    need(summary['formula']['variables'] == nv and summary['formula']['clauses'] == len(cnf), 'E4: summary formula')
    return '81 atoms and %d base clauses rebuilt; %d recorded graphs (%d toroidal, %d nontoroidal); %s; %s' % (
        len(base), len(cands), n_tor, n_non, nums(counts),
        'DRAT proof verified by drat-trim' if drat_ok else 'DRAT proof not checked (drat-trim not found)')


def check(root, item):
    if item in ('E7', 'E8'):
        return check_E78(root, item)
    return {'E1': check_E1, 'E2': check_E2, 'E3': check_E3, 'E4': check_E4, 'E12': check_E12}[item](root)
