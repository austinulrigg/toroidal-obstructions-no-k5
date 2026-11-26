"""The enumerations E19-E29 (Proposition 5.33 and Section 6 of the paper).

Each folder holds input.json, candidates.json, verdicts.json and summary.json,
any of them possibly gzip-compressed; E29 also holds bases.json.gz. Only these
files and O1 (for the graphs T1, ..., T74) are read. For every item:

 1. Records. Every candidate has a distinct id and exactly one verdict, and
    every verdict belongs to a candidate.
 2. Graphs. Each candidate's graph is rebuilt from its vertex and edge lists
    (no loops, repeated edges or unknown ends), and its graph6 string, where
    given, encodes the same graph with the vertices in the listed order.
 3. Vertex names, as the item's 'vertex names' describe them. The twins of
    Psi have exactly their stated neighbourhoods: p, q on {a,b,c}, and r, z on
    {a,b,d} in E19 (Psi_2) or on {a,d,e} in E20-E29 (Psi_3). Further vertices
    p3, p4, ... and r3, r4, ... have the neighbourhoods of p and r. In E29 a
    vertex v_xyz has neighbourhood {x,y,z}.
 4. How each candidate is formed.
      E19  The graph is Psi_2 with the forms of its combination added, the
           interior vertex i<j> of the k-th form becoming f<k>_<j>, and with
           the edge cd when the combination has it (see 9).
      E29  The graph is its base graph with n_c - 2 vertices p3, p4 on {a,b,c}
           and n_e - 2 vertices r3, r4 on {a,d,e} added, for the (n_c, n_e) of
           its case (see 10).
      the others  The interior vertices of the bridge (u0, u1, ...) induce a
           connected graph whose neighbours outside are exactly the
           attachments of the case (the roots of its specification in
           input.json). Where the bridge is recorded, it is exactly these
           vertices and attachments with every edge at an interior vertex. All
           candidates of one case (and variant) have the same exterior, the
           graph left when the interior is deleted. Where input.json gives the
           number of bridges the construction found for a case, the case has
           that many candidates for each of its variants.
 5. Verdicts, by the kind of witness.
      toroidal                     an embedding whose traced faces give genus
                                   at most one
      K5 minor                     five disjoint connected branch sets with an
                                   edge between every two
      cut of at most two vertices  one or two vertices whose deletion
                                   disconnects the graph
      isomorphic to Ti             a bijection onto Ti of O1, edges onto edges
      excluded by a hypothesis     the stated fact holds (see 6)
      proper nontoroidal minor     the edge deleted or contracted is an edge of
                                   the graph. The nontoroidality of the minor
                                   is left to the genus computations in the
                                   folder genus, and such witnesses are
                                   counted.
 6. Facts of the hypotheses of the case (field 'also excluded by a hypothesis
    of the case', and the discarded bases of E29) are rechecked on the graph.
    A twin pair is two distinct vertices whose neighbourhood is exactly the
    stated triple of three distinct vertices.
      two twin triples are disjoint / share two vertices / meet in zero or two
      vertices: two twin triples with twin pairs, meeting as stated;
      three twin triples share a pair of vertices: three distinct twin triples
      with twin pairs, all three containing two common vertices;
      a choice of twin pairs leaves a component adjacent to all four (five) of
      their attachments: two twin triples meeting in two vertices (in one
      vertex) with twin pairs; the component is disjoint from the triples and
      twins, connected, has no neighbour outside them (so it is a whole
      component once they are deleted), and is adjacent to every vertex of the
      two triples;
      a choice of twin pairs leaves a component with at least four neighbours:
      as before with triples meeting in one vertex; the recorded neighbours are
      exactly the neighbours of the component, at least four of them;
      a vertex of degree three has two adjacent neighbours: as stated.
 7. E26-E28, the planarity of the bridge completions: each completion is the
    recorded bridge with the edges its name lists added, these are two of the
    three pairs of its attachments, the three completions are the three ways
    of choosing two such pairs, and each answer has its witness: a plane
    embedding (genus zero) or a Kuratowski subgraph of that completion.
 8. Summary. Every field of summary.json is recomputed from the records (the
    counts by verdict, the copies of each Ti, the lists of Ti excluded and not
    excluded by a hypothesis of the case, the counts of each fact, the counts
    by case, and the counts of combinations of E19), and a field the program
    does not recompute is an error. Every quoted passage of the paper in
    input.json cites the item, and the counts it states (the numbers of
    candidates, by case where the text counts by case, and for E29 the kept
    bases and the counts by verdict) are the counts of the records. (The
    passage quoted for E19 ends before its counts, so there only the summary
    is compared with the lattice.)
 9. E19, the forms and the lattice of combinations (Proposition 5.33).
    The 76 labeled forms: each has a triple among abc, abd, acd, bcd, interior
    vertices i0, i1, ... inducing a connected graph, attachments exactly its
    triple and no edge between attachments; forms of the same name are one
    shape with x, y, z placed as their field 'x,y,z' says; the eight shapes,
    with every bijection of x, y, z onto each triple, give exactly 76 classes
    up to isomorphisms fixing the attachments, one recorded form in each.
    The domain is every multiset of one to four forms, not all stars, with cd
    absent or present. Every recorded combination lies in the domain and is
    recorded once. Its lattice status agrees with its verdict: 'violation'
    has a fact (rechecked), 'toroidal, no violation' has none, a torus
    embedding never has status 'nontoroidal', and a proper nontoroidal minor
    or an isomorphism onto Ti always has it. The auxiliary combinations are
    all stars, outside the domain, recorded once and not among the candidates,
    and each graph is rebuilt from its forms; a 'nontoroidal' status is
    witnessed by the recorded isomorphism onto the named Ti or left to the
    genus computations; their other recorded outcomes are checked on the
    named graph. Then every combination of the
    domain is visited, smaller combinations first: each one not recorded must
    contain a smaller combination (fewer forms, or the same forms without cd)
    with status 'violation' or 'nontoroidal', as the text's argument requires;
    a recorded combination with status 'toroidal, no violation' must contain
    none; a recorded combination with a torus embedding must contain no
    smaller combination with status 'nontoroidal' (a subgraph of a toroidal
    graph is toroidal). Exactly four recorded combinations are not excluded by
    a fact, a torus embedding, a cut or a proper nontoroidal minor; all four
    are isomorphic to T15, and one of them is the representative described in
    the text (Psi_2, a vertex h on {a,b,c}, and the bridge with edges au, cp',
    cq', dp', dq', up', uq', cd absent).
10. E29, the bases (Construction 6.24 with n_c = n_e = 2). Each base is
    rebuilt (graph, graph6, twins of Psi_3) and agrees with its description
    (its single vertices v_xyz, its cross edges, which are the only edges
    among a, ..., e, and its forms B<k>_..., each connected with attachments
    exactly its triple, at most two forms, triples and cross edges allowed for
    the matching or the path). A discarded base has its fact rechecked; the
    text says the recorded witnesses use degree-three vertices outside
    {a,b,c,d,e} and that both kinds of discard remain valid for all five
    choices of (n_c, n_e): this is checked by adding the further twins of each
    choice and rechecking the fact (for a component, the component of the
    larger graph containing it still has at least four neighbours). The
    counts of bases.json agree with the bases, the kept bases are exactly
    those used by candidates, and each kept base with each of the five
    choices of (n_c, n_e) gives exactly one candidate.
"""
import os
import re
from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement, permutations, product

from witnesses import (CheckFailed, need, load, find, graph, copy, add_edges, n_edges, same_graph,
                       components, from_graph6, check_embedding, check_k5, check_cut, check_iso,
                       check_kuratowski, as_vertex_of, catalogue)


# ---------------------------------------------------------------- the items

PSI5 = {'p': 'abc', 'q': 'abc', 'r': 'abd', 'z': 'abd'}       # Psi_2, Section 5
PSI6 = {'p': 'abc', 'q': 'abc', 'r': 'ade', 'z': 'ade'}       # Psi_3, Section 6

CONF = {'E%d' % i: {'psi': PSI6, 'formed': 'bridge', 'interior': 'u'} for i in range(20, 29)}
CONF['E19'] = {'psi': PSI5, 'formed': 'forms'}
CONF['E29'] = {'psi': PSI6, 'formed': 'bases', 'single vertices': True}
for _k in ('E26', 'E27', 'E28'):
    CONF[_k]['completions'] = True

VERDICT_KEYS = {'id', 'verdict', 'witness', 'also excluded by a hypothesis of the case',
                'planarity of the bridge completions', 'lattice status'}
ISO = re.compile(r'isomorphic to (T\d+)')
E29_CHOICES = [(2, 2), (2, 3), (2, 4), (3, 2), (4, 2)]


class Tally:
    """What the records of one item say, collected while they are checked."""

    def __init__(self):
        self.n = 0
        self.verdicts = Counter()
        self.iso = Counter()                    # Ti -> copies
        self.iso_free = set()                   # Ti with a copy not excluded by a fact
        self.iso_excl = set()                   # Ti with a copy excluded by a fact
        self.facts = Counter()
        self.case_count = Counter()
        self.case_free = defaultdict(set)
        self.variant_count = Counter()
        self.edge_only = 0                      # proper minors named only by the edge
        self.exteriors = {}
        self.completions = 0
        self.planar = []
        self.e19_rows = []
        self.e19_graphs = {}
        self.e29_pairs = set()
        self.extra = {}
        self.notes = []


def fmt(n):
    return '{:,}'.format(n)


def tsort(names):
    return sorted(names, key=lambda s: int(s[1:]))


def keys_are(d, allowed, what):
    need(isinstance(d, dict), what, 'is not an object')
    need(set(d) == set(allowed), what, 'has fields', sorted(d), 'where', sorted(allowed), 'are expected')


def guarded(where, fn, *args):
    """Run one check, naming the record in any failure."""
    try:
        return fn(*args)
    except CheckFailed as e:
        raise CheckFailed('%s: %s' % (where, e)) from None
    except (KeyError, IndexError, TypeError, ValueError, AttributeError) as e:
        raise CheckFailed('%s: malformed record (%s %s)' % (where, type(e).__name__, e)) from None


# ---------------------------------------------------------------- graphs and names

def record_graph(rec):
    """The graph of a record, with its graph6 string (when given) checked to
    encode it, the vertices numbered in the listed order."""
    adj = graph(rec['graph'])
    if 'graph6' in rec:
        order = rec['graph']['vertices']
        g6 = from_graph6(rec['graph6'])
        need(len(g6) == len(order), 'graph6 has', len(g6), 'vertices and the graph', len(order))
        pos = {v: i for i, v in enumerate(order)}
        for v in order:
            need({pos[u] for u in adj[v]} == g6[pos[v]], 'graph6 does not encode the graph (at vertex %s)' % v)
    return adj


def check_names(adj, conf):
    """The twins of Psi and the other named vertices have the neighbourhoods
    the item's vertex names give them."""
    psi = conf['psi']
    for v, T in psi.items():
        need(v in adj, 'there is no vertex', v)
        need(adj[v] == set(T), v, 'has neighbourhood', sorted(adj[v]), 'and not {%s}' % ','.join(T))
    for v in adj:
        if len(v) > 1 and v[0] in 'pr' and v[1:].isdigit() and int(v[1:]) >= 3:
            need(adj[v] == set(psi[v[0]]), v, 'has neighbourhood', sorted(adj[v]),
                 'and not {%s}' % ','.join(psi[v[0]]))
    if conf.get('single vertices'):
        for v in adj:
            if v.startswith('v_'):
                need(re.fullmatch(r'v_[a-e]{3}', v) and adj[v] == set(v[2:]), v, 'has neighbourhood', sorted(adj[v]))


def is_interior(v, conf):
    return v[:1] == 'u' and v[1:].isdigit()


def roots_of(case, conf, inp, item):
    cases = inp.get('cases of the construction', {})
    need(case in cases, 'case', case, 'is not described in input.json')
    return set(cases[case]['specification']['roots'])


def bridge_formation(adj, c, conf, inp, item, t):
    """The candidate is its exterior plus a bridge whose connected interior
    attaches exactly at the roots of its case; the exterior is the same for
    every candidate of the case and variant."""
    interior = {v for v in adj if is_interior(v, conf)}
    need(interior, 'the graph has no interior vertices')
    need(len(components(adj, interior)) == 1, 'the interior of the bridge is not connected')
    att = {y for x in interior for y in adj[x]} - interior
    case = c.get('case')
    roots = roots_of(case, conf, inp, item)
    need(att == roots, 'the bridge attaches at', sorted(att), 'and not at the roots', sorted(roots), 'of its case')
    if 'bridge' in c:
        B = graph(c['bridge'])
        need(set(B) == interior | att, 'the recorded bridge is not the interior with its attachments')
        be = {frozenset((x, y)) for x in B for y in B[x]}
        ge = {frozenset((x, y)) for x in interior for y in adj[x]}
        need(be == ge, 'the edges of the recorded bridge are not the edges at its interior')
    key = (case, c.get('variant'))
    ext = (frozenset(set(adj) - interior),
           frozenset(frozenset((x, y)) for x in adj if x not in interior for y in adj[x] if y not in interior))
    first = t.exteriors.setdefault(key, ext)
    need(first == ext, 'the exterior differs from that of the other candidates of case', key)
    if 'bd added' in c:
        need(c['bd added'] == ('d' in adj['b']), "the field 'bd added' does not agree with the edge bd")
    return interior, att


# ---------------------------------------------------------------- facts

FACT_DISJOINT = 'two twin triples are disjoint'
FACT_SHARE2 = 'two twin triples share two vertices'
FACT_ZERO_TWO = 'two twin triples meet in zero or two vertices'
FACT_THREE = 'three twin triples share a pair of vertices'
FACT_FOUR_ATT = 'a choice of twin pairs leaves a component adjacent to all four of their attachments'
FACT_FIVE_ATT = 'a choice of twin pairs leaves a component adjacent to all five of their attachments'
FACT_FOUR_NB = 'a choice of twin pairs leaves a component with at least four neighbours'
FACT_ADJ = 'a vertex of degree three has two adjacent neighbours'


def twin_triples(adj, f, k):
    """k triples of three distinct vertices, each with a twin pair: two distinct
    vertices whose neighbourhood is exactly the triple."""
    T, P = f['triples'], f['twins']
    need(isinstance(T, list) and isinstance(P, list) and len(T) == k == len(P), 'the fact needs', k,
         'triples with their twin pairs')
    out = []
    for tri, pair in zip(T, P):
        need(len(tri) == 3 == len(set(tri)), 'a twin triple needs three distinct vertices:', tri)
        need(len(pair) == 2 and pair[0] != pair[1], 'a twin pair needs two distinct vertices:', pair)
        for x in pair:
            need(x in adj and adj[x] == set(tri), x, 'does not have neighbourhood', tri)
        out.append(set(tri))
    return out, [list(p) for p in P]


def chosen_component(adj, f):
    """Two twin triples with twin pairs and a component left when the triples
    and the twins are deleted; returns the triples, the component and its
    neighbours."""
    (A, B), (PA, PB) = twin_triples(adj, f, 2)
    core = A | B | set(PA) | set(PB)
    comp = f['component']
    S = set(comp)
    need(S and len(S) == len(comp) and S <= set(adj), 'the component is not a set of vertices of the graph')
    need(not S & core, 'the component meets the chosen twins or their triples')
    need(len(components(adj, S)) == 1, 'the component is not connected')
    nb = {y for x in S for y in adj[x]} - S
    need(nb <= core, 'the component is not a whole component once the twins and their triples are deleted')
    return A, B, S, nb


def check_fact(adj, f, inp):
    """Recheck a fact of a hypothesis of the case on the graph adj."""
    need(isinstance(f, dict) and isinstance(f.get('fact'), str), 'the fact is not recorded')
    fact = f['fact']
    if fact == FACT_ADJ:
        keys_are(f, {'fact', 'vertex', 'neighbours'}, 'the fact')
        v, (x, y) = f['vertex'], f['neighbours']
        need(v in adj and len(adj[v]) == 3, v, 'is not a vertex of degree three')
        need(x != y and x in adj[v] and y in adj[v] and y in adj[x], x, y, 'are not adjacent neighbours of', v)
    elif fact in (FACT_DISJOINT, FACT_SHARE2, FACT_ZERO_TWO):
        keys_are(f, {'fact', 'triples', 'twins'}, 'the fact')
        (A, B), _ = twin_triples(adj, f, 2)
        k = len(A & B)
        want = {FACT_DISJOINT: (0,), FACT_SHARE2: (2,), FACT_ZERO_TWO: (0, 2)}[fact]
        need(k in want, 'the two twin triples meet in', k, 'vertices')
    elif fact == FACT_THREE:
        keys_are(f, {'fact', 'triples', 'twins'}, 'the fact')
        (A, B, C), _ = twin_triples(adj, f, 3)
        need(len({frozenset(A), frozenset(B), frozenset(C)}) == 3, 'the three twin triples are not distinct')
        need(len(A & B & C) >= 2, 'the three twin triples do not share a pair')
    elif fact in (FACT_FOUR_ATT, FACT_FIVE_ATT):
        keys_are(f, {'fact', 'triples', 'twins', 'component'}, 'the fact')
        A, B, S, nb = chosen_component(adj, f)
        need(len(A & B) == (2 if fact == FACT_FOUR_ATT else 1), 'the twin triples meet in', len(A & B), 'vertices')
        need(nb == A | B, 'the component is not adjacent to every attachment', sorted(A | B))
    elif fact == FACT_FOUR_NB:
        keys_are(f, {'fact', 'triples', 'twins', 'component', 'neighbours'}, 'the fact')
        A, B, S, nb = chosen_component(adj, f)
        need(len(A & B) == 1, 'the twin triples meet in', len(A & B), 'vertices, not one')
        rec = f['neighbours']
        need(len(rec) == len(set(rec)) and set(rec) == nb, 'the neighbours of the component are', sorted(nb),
             'not', rec)
        need(len(nb) >= 4, 'the component has only', len(nb), 'neighbours')
    else:
        raise CheckFailed('unknown fact: %s' % fact)


# ---------------------------------------------------------------- verdicts

def check_proper_minor(adj, w, cat, t, cid):
    need(w.get('type') == 'proper minor', 'the witness of a proper nontoroidal minor is of type', w.get('type'))
    keys_are(w, {'type', 'operation', 'edge'}, 'the witness')
    op, e = w['operation'], w['edge']
    need(op in ('delete', 'contract'), 'unknown operation', op)
    need(isinstance(e, list) and len(e) == 2, 'the operation needs one edge')
    u, v = e
    need(u in adj and v in adj[u], 'the edge', u, v, 'of the minor is not in the graph')
    t.edge_only += 1


def check_verdict(adj, verdict, w, cat, inp, t, cid):
    """Check the witness of a verdict by its kind; returns Ti for an isomorphism."""
    need(isinstance(w, dict), 'the witness is not an object')
    if verdict == 'toroidal':
        keys_are(w, {'type', 'genus', 'rotation'}, 'the embedding')
        need(w['type'] == 'embedding', 'the witness of toroidal is of type', w['type'])
        check_embedding(adj, w, at_most=1)
    elif verdict == 'K5 minor':
        keys_are(w, {'type', 'branch sets'}, 'the K5 witness')
        need(w['type'] == 'K5 branch sets', 'the witness of a K5 minor is of type', w['type'])
        check_k5(adj, w['branch sets'])
    elif verdict == 'cut of at most two vertices':
        keys_are(w, {'type', 'vertices'}, 'the cut')
        need(w['type'] == 'cut', 'the witness of a cut is of type', w['type'])
        check_cut(adj, w['vertices'])
    elif ISO.fullmatch(verdict):
        name = ISO.fullmatch(verdict).group(1)
        need(name in cat, 'unknown graph', name)
        keys_are(w, {'type', 'onto', 'map'}, 'the isomorphism')
        need(w['type'] == 'isomorphism' and w['onto'] == name, 'the isomorphism is not onto', name)
        T = cat[name]
        check_iso(adj, T, {k: as_vertex_of(T, x) for k, x in w['map'].items()})
        return name
    elif verdict == 'proper nontoroidal minor':
        check_proper_minor(adj, w, cat, t, cid)
    elif verdict == 'excluded by a hypothesis of the case':
        check_fact(adj, w, inp)
    else:
        raise CheckFailed('unknown verdict: %s' % verdict)
    return None


def check_completions(c, att, tests, t):
    """E26-E28: the three two-edge completions of the bridge, each formed as
    its name says and with a plane embedding or a Kuratowski subgraph."""
    need('bridge' in c, 'completions without a recorded bridge')
    B = graph(c['bridge'])
    need(len(att) == 3, 'the bridge does not have three attachments')
    pairs = {frozenset(p) for p in combinations(sorted(att), 2)}
    need(isinstance(tests, list), 'the completions are not a list')
    seen = set()
    for x in tests:
        keys_are(x, {'completion', 'added edges', 'answer', 'witness'}, 'a completion')
        name = x['completion']
        parts = name.split('+')
        need(len(parts) == 3 and parts[0] in ('U', 'B') and all(len(s) == 2 for s in parts[1:]),
             'cannot read the completion', name)
        named = {frozenset(s) for s in parts[1:]}
        added = {frozenset(e) for e in x['added edges']}
        need(len(x['added edges']) == 2 and named == added, 'the completion', name, 'adds', x['added edges'])
        need(added <= pairs, 'the completion', name, 'adds an edge not joining two attachments')
        need(frozenset(added) not in seen, 'the completion', name, 'appears twice')
        seen.add(frozenset(added))
        H = add_edges(B, x['added edges'])
        wit = x['witness']
        if x['answer'] == 'planar':
            keys_are(wit, {'type', 'genus', 'rotation'}, 'the plane embedding')
            need(wit['type'] == 'embedding', 'a planar completion needs an embedding')
            check_embedding(H, wit, genus=0)
            t.planar.append((c.get('case'), name))
        elif x['answer'] == 'nonplanar':
            keys_are(wit, {'type', 'edges'}, 'the Kuratowski subgraph')
            need(wit['type'] == 'Kuratowski subgraph', 'a nonplanar completion needs a Kuratowski subgraph')
            check_kuratowski(H, wit['edges'])
        else:
            raise CheckFailed('unknown answer %s' % x['answer'])
        t.completions += 1
    need(seen == {frozenset(s) for s in combinations(pairs, 2)}, 'the completions are not the three two-edge completions')


# ---------------------------------------------------------------- one candidate

def check_candidate(c, v, conf, inp, cat, ctx, t, item):
    adj = record_graph(c)
    check_names(adj, conf)
    att = None
    if conf['formed'] == 'bridge':
        _, att = bridge_formation(adj, c, conf, inp, item, t)
    elif conf['formed'] == 'forms':
        e19_formation(adj, c, ctx['forms'])
    else:
        e29_formation(adj, c, ctx['bases'], t)
    extra = set(v) - VERDICT_KEYS
    need(not extra, 'the verdict has fields the checker does not know:', sorted(extra))
    verdict, w = v['verdict'], v['witness']
    name = check_verdict(adj, verdict, w, cat, inp, t, c['id'])
    fact = v.get('also excluded by a hypothesis of the case')
    if fact is not None:
        check_fact(adj, fact, inp)
        t.facts[fact['fact']] += 1
    if 'planarity of the bridge completions' in v:
        need(conf.get('completions'), 'completions are not expected in this item')
        check_completions(c, att, v['planarity of the bridge completions'], t)
    if item == 'E19':
        need('lattice status' in v, 'no lattice status')
        e19_status(v['lattice status'], verdict, fact is not None)
        ms = [ctx['forms']['index'][f] for f in c['config']['forms']]
        # a status is witnessed in the folder by a rechecked fact or an isomorphism onto some Ti
        witnessed = fact is not None or name is not None
        t.e19_rows.append((c['id'], tuple(sorted(ms)), int(c['config']['cd']), v['lattice status'], verdict,
                           fact is not None, witnessed))
        if name is not None:
            t.e19_graphs[c['id']] = adj
    else:
        need('lattice status' not in v, 'a lattice status outside E19')
    t.n += 1
    t.verdicts[verdict] += 1
    case = c.get('case')
    if case is not None:
        t.case_count[case] += 1
        t.variant_count[(case, c.get('variant'))] += 1
    if name is not None:
        t.iso[name] += 1
        if fact is None:
            t.iso_free.add(name)
            if case is not None:
                t.case_free[case].add(name)
        else:
            t.iso_excl.add(name)


# ---------------------------------------------------------------- E19

def canonical(edges, fixed):
    """A canonical form of a graph given by its edges, up to isomorphisms
    fixing the vertices in fixed: the least sorted edge list over the
    renamings of the other vertices to #0, #1, ... that respect a colour
    refinement started from the fixed neighbours and degrees."""
    fixed = set(fixed)
    adj = {}
    for u, v in edges:
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
    inner = [v for v in adj if v not in fixed]
    col = {v: (tuple(sorted(adj[v] & fixed)), len(adj[v])) for v in inner}
    for _ in range(3):
        col = {v: (col[v], tuple(sorted(col[u] for u in adj[v] if u in col))) for v in inner}
    groups = defaultdict(list)
    for v in inner:
        groups[col[v]].append(v)
    groups = [groups[k] for k in sorted(groups)]
    best = None
    for perm in product(*(permutations(g) for g in groups)):
        order = [v for g in perm for v in g]
        nm = {v: '#%d' % i for i, v in enumerate(order)}
        key = tuple(sorted(tuple(sorted((nm.get(a, a), nm.get(b, b)))) for a, b in edges))
        if best is None or key < best:
            best = key
    return best


def e19_forms(inp):
    """The 76 labeled forms (see 9 in the module text)."""
    L = inp['lattice']
    keys_are(L, {'forms', 'domain', 'auxiliary'}, 'the lattice')
    F = L['forms']
    ids = [f['id'] for f in F]
    need(len(F) == 76 and len(set(ids)) == 76, 'the lattice does not have 76 distinct forms')
    shapes = defaultdict(set)
    recorded = set()
    for f in F:
        keys_are(f, {'id', 'triple', 'name', 'x,y,z', 'edges'}, 'form ' + str(f.get('id')))
        T = f['triple']
        need(T in ('abc', 'abd', 'acd', 'bcd'), 'form', f['id'], 'has triple', T)
        need(sorted(f['x,y,z']) == list(T), 'form', f['id'], ': x, y, z are not placed on its triple')
        adj = {}
        for u, v in f['edges']:
            need(u != v and v not in adj.get(u, ()), 'form', f['id'], 'has a loop or a repeated edge')
            adj.setdefault(u, set()).add(v)
            adj.setdefault(v, set()).add(u)
        inner = [v for v in adj if v not in set(T)]
        need(sorted(inner) == sorted('i%d' % j for j in range(len(inner))), 'form', f['id'], 'has vertices', sorted(adj))
        need(set(adj) - set(inner) == set(T), 'form', f['id'], 'does not attach at all of', T)
        need(all(x in inner or y in inner for x in adj for y in adj[x]), 'form', f['id'], 'has an edge between attachments')
        need(len(components(adj, inner)) == 1, 'the interior of form', f['id'], 'is not connected')
        back = dict(zip(f['x,y,z'], 'xyz'))
        shapes[f['name']].add(canonical([[back.get(x, x) for x in e] for e in f['edges']], 'xyz'))
        recorded.add((T, canonical(f['edges'], T)))
    for name, s in shapes.items():
        need(len(s) == 1, 'the forms named', name, 'are not one shape')
    need(len(shapes) == 8, 'the forms have', len(shapes), 'shapes and not eight')
    need(len(recorded) == 76, 'two recorded forms are the same labeled form')
    classes = set()
    for (shape,) in shapes.values():
        for T in ('abc', 'abd', 'acd', 'bcd'):
            for perm in permutations(T):
                m = dict(zip('xyz', perm))
                classes.add((T, canonical([[m.get(x, x) for x in e] for e in shape], T)))
    need(classes == recorded, 'the recorded forms are not the', len(classes), 'labelings of the eight shapes')
    order = sorted(ids)
    name = {f['id']: f['name'] for f in F}
    return {'by id': {f['id']: f for f in F}, 'index': {fid: i for i, fid in enumerate(order)}, 'ids': order,
            'star': [name[fid] == 'star' for fid in order]}


def e19_graph(cfg, forms):
    """Psi_2 with the forms of the combination and cd added."""
    keys_are(cfg, {'forms', 'cd'}, 'the combination')
    need(isinstance(cfg['forms'], list) and 1 <= len(cfg['forms']) <= 4, 'a combination has one to four forms')
    need(cfg['cd'] in (True, False), 'cd is neither absent nor present')
    adj = {}

    def add(u, v):
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
    for t_, T in PSI5.items():
        for x in T:
            add(t_, x)
    for k, fid in enumerate(cfg['forms']):
        need(fid in forms['by id'], 'unknown form', fid)
        for e in forms['by id'][fid]['edges']:
            u, v = ('f%d_%s' % (k, x[1:]) if x[:1] == 'i' else x for x in e)
            add(u, v)
    if cfg['cd']:
        add('c', 'd')
    return adj


def e19_formation(adj, c, forms):
    cfg = c['config']
    need(same_graph(adj, e19_graph(cfg, forms)), 'the graph is not Psi_2 with the forms of its combination')
    need(c.get('case') == 'size %d' % len(cfg['forms']), 'case', c.get('case'), 'for', len(cfg['forms']), 'forms')
    need((c.get('variant') or '') == ('cd' if cfg['cd'] else ''), 'variant', c.get('variant'), 'with cd', cfg['cd'])


def e19_status(status, verdict, has_fact):
    need(status in ('violation', 'nontoroidal', 'toroidal, no violation'), 'unknown lattice status', status)
    if status == 'violation':
        need(has_fact, 'status violation without a fact')
    if status == 'toroidal, no violation':
        need(not has_fact, 'status toroidal, no violation with a fact')
    if verdict == 'toroidal':
        need(status != 'nontoroidal', 'a torus embedding with status nontoroidal')
    if verdict == 'proper nontoroidal minor' or ISO.fullmatch(verdict):
        need(status == 'nontoroidal', 'verdict', verdict, 'with status', status)


def find_isomorphism(G, H, fixed=()):
    """An isomorphism from G onto H fixing the vertices in fixed, found by
    backtracking, or None."""
    if len(G) != len(H) or n_edges(G) != n_edges(H):
        return None
    if sorted(len(s) for s in G.values()) != sorted(len(s) for s in H.values()):
        return None
    fixed = list(fixed)
    if any(v not in G or v not in H or len(G[v]) != len(H[v]) for v in fixed):
        return None
    # fixed vertices first, then always the vertex with most neighbours already placed
    order, seen = list(fixed), set(fixed)
    rest = sorted((v for v in G if v not in seen), key=str)
    while rest:
        v = max(rest, key=lambda x: (len(G[x] & seen), len(G[x])))
        order.append(v)
        seen.add(v)
        rest.remove(v)
    m, used = {}, set()

    def extend(i):
        if i == len(order):
            return True
        v = order[i]
        cands = [v] if v in fixed else [t_ for t_ in H if t_ not in used and len(H[t_]) == len(G[v])]
        for t_ in cands:
            if t_ in used:
                continue
            if all((u in G[v]) == (tu in H[t_]) for u, tu in m.items()):
                m[v] = t_
                used.add(t_)
                if extend(i + 1):
                    return True
                del m[v]
                used.discard(t_)
        return False
    return dict(m) if extend(0) else None


def e19_auxiliary(inp, forms, cat, rec_keys, cand_ids):
    """Check the stars-only combinations recorded for the exclusions and return their keys and flags."""
    out = []
    seen = set()
    for a in inp['lattice']['auxiliary']:
        aid = a.get('id')

        def one():
            cfg = a['config']
            G = e19_graph(cfg, forms)
            ms = tuple(sorted(forms['index'][f] for f in cfg['forms']))
            need(all(forms['star'][x] for x in ms), 'the combination is not all stars')
            need(a['in domain'] is False, 'it is marked as in the domain')
            key = (int(cfg['cd']), ms)
            need(key not in rec_keys and key not in seen, 'the combination is recorded twice')
            need(aid not in cand_ids, 'its id is also a candidate id')
            seen.add(key)
            need(a.get('case') == 'size %d' % len(ms), 'case', a.get('case'))
            need((a.get('variant') or '') == ('cd' if cfg['cd'] else ''), 'variant', a.get('variant'))
            adj = graph(a['graph'])
            need(same_graph(adj, G), 'the graph is not Psi_2 with the forms of its combination')
            status = a['lattice status']
            need(status in ('nontoroidal', 'toroidal, no violation'), 'status', status, 'has no fact to check')
            o = a['outcome']
            k = o.get('outcome')
            iso_ok = False
            if k == 'cut of order at most two':
                check_cut(adj, o['cut'])
            elif k == 'proper nontoroidal minor':
                need(status == 'nontoroidal', 'a proper nontoroidal minor with status', status)
                need(o['operation'] in ('delete', 'contract'), 'unknown operation')
                u, v = o['edge']
                need(u in adj and v in adj[u], 'the edge', u, v, 'is not in the graph')
            elif k == 'toroidal':
                need(status != 'nontoroidal', 'a torus embedding with status nontoroidal')
                rot, g = {}, 0
                for p in o['embedding']:
                    rot.update(p['rotation'])
                    g += p['genus']
                check_embedding(adj, {'type': 'embedding', 'genus': g, 'rotation': rot}, at_most=1)
            elif k == 'catalogue graph':
                need(status == 'nontoroidal', 'an isomorphism onto Ti with status', status)
                need(o['catalogue'] in cat, 'unknown graph', o['catalogue'])
                T = cat[o['catalogue']]
                check_iso(adj, T, {x: as_vertex_of(T, y) for x, y in o['map'].items()})
                iso_ok = True
            else:
                raise CheckFailed('unknown outcome %s' % k)
            return key, status, iso_ok
        out.append(guarded('E19 auxiliary combination %s' % aid, one))
    return out


def e19_lattice(inp, forms, cat, t):
    """Visit every combination of the domain (see 9 in the module text)."""
    BAD, NT, WIT, REC = 1, 2, 4, 8
    R = 77
    star = forms['star']
    nf = len(forms['ids'])

    def code(ms):
        c = 0
        for x in ms:
            c = c * R + x + 1
        return c

    def decode(c):
        out = []
        while c:
            c, d = divmod(c, R)
            out.append(forms['ids'][d - 1])
        return list(reversed(out))
    flags = (bytearray(R ** 4), bytearray(R ** 4))
    rec = ({}, {})
    for cid, ms, cd, status, verdict, has_fact, witnessed in t.e19_rows:
        need(not all(star[x] for x in ms), 'E19 %s: the combination is all stars, outside the domain' % cid)
        c = code(ms)
        need(c not in rec[cd], 'E19: combination %s with cd %s is recorded twice' % (decode(c), bool(cd)))
        rec[cd][c] = (cid, status, verdict)
        if status == 'violation':
            flags[cd][c] |= BAD | WIT | REC
        elif status == 'nontoroidal':
            flags[cd][c] |= BAD | NT | REC | (WIT if witnessed else 0)
    rec_keys = {(r[2], r[1]) for r in t.e19_rows}
    cand_ids = {r[0] for r in t.e19_rows}
    aux = e19_auxiliary(inp, forms, cat, rec_keys, cand_ids)
    aux_nt = aux_nt_wit = 0
    for (cd, ms), status, iso_ok in aux:
        if status == 'nontoroidal':
            flags[cd][code(ms)] |= BAD | NT | (WIT if iso_ok else 0)
            aux_nt += 1
            aux_nt_wit += bool(iso_ok)
    cnt = Counter()

    def visit(c, subs, all_stars):
        for cd in (0, 1):
            arr = flags[cd]
            s = 0
            for sc in subs:
                s |= arr[sc]
            if cd:
                s |= flags[0][c]
            arr[c] |= s
            if all_stars:
                continue
            cnt['domain'] += 1
            r = rec[cd].get(c)
            if r is None:
                need(s & BAD, 'E19: combination %s with cd %s is neither recorded nor contains a smaller combination '
                     'with status violation or nontoroidal' % (decode(c), bool(cd)))
                cnt['excluded'] += 1
                cnt['excluded, witnessed'] += bool(s & WIT)
            else:
                cid, status, verdict = r
                cnt['recorded'] += 1
                cnt['recorded, already excluded'] += bool(s & BAD)
                need(not (status == 'toroidal, no violation' and s & BAD),
                     'E19 %s: status toroidal, no violation, but it contains a smaller combination with status '
                     'violation or nontoroidal' % cid)
                need(not (verdict == 'toroidal' and s & NT),
                     'E19 %s: a torus embedding, but it contains a smaller combination with status nontoroidal' % cid)
    for i in range(nf):
        visit(i + 1, (), star[i])
    for i in range(nf):
        for j in range(i, nf):
            visit((i + 1) * R + j + 1, (i + 1, j + 1), star[i] and star[j])
    for i in range(nf):
        for j in range(i, nf):
            cij = (i + 1) * R + j + 1
            for k in range(j, nf):
                visit(cij * R + k + 1, ((j + 1) * R + k + 1, (i + 1) * R + k + 1, cij),
                      star[i] and star[j] and star[k])
    for i in range(nf):
        for j in range(i, nf):
            cij = (i + 1) * R + j + 1
            for k in range(j, nf):
                cijk = cij * R + k + 1
                cjk = (j + 1) * R + k + 1
                cik = (i + 1) * R + k + 1
                for l in range(k, nf):
                    visit(cijk * R + l + 1, (cjk * R + l + 1, cik * R + l + 1, cij * R + l + 1, cijk),
                          star[i] and star[j] and star[k] and star[l])
    need(cnt['recorded'] == len(t.e19_rows), 'E19: a recorded combination lies outside the domain')
    # the combinations that remain: no fact, no torus embedding, no cut, no proper nontoroidal minor
    remaining = [r for r in t.e19_rows if not r[5]
                 and r[4] not in ('toroidal', 'cut of at most two vertices', 'proper nontoroidal minor')]
    need(len(remaining) == 4 and all(r[4] == 'isomorphic to T15' for r in remaining),
         'E19: the combinations that remain are', [(r[0], r[4]) for r in remaining], 'and not four copies of T15')
    rep = {}
    for a_, b_ in ([(x, y) for x in 'pq' for y in 'abc'] + [(x, y) for x in 'rz' for y in 'abd'] +
                   [('h', y) for y in 'abc'] +
                   [('a', 'u'), ('c', "p'"), ('c', "q'"), ('d', "p'"), ('d', "q'"), ('u', "p'"), ('u', "q'")]):
        rep.setdefault(a_, set()).add(b_)
        rep.setdefault(b_, set()).add(a_)
    hit = [r[0] for r in remaining if find_isomorphism(t.e19_graphs[r[0]], rep, 'abcdpqrz') is not None]
    need(hit, 'E19: the representative described in the text is none of the four remaining combinations')
    t.extra.update({'in the domain': cnt['domain'], 'recorded': cnt['recorded'],
                    'excluded': cnt['excluded'], 'auxiliary': len(aux)})
    rec_nt = sum(1 for r in t.e19_rows if r[3] == 'nontoroidal')
    rec_nt_wit = sum(1 for r in t.e19_rows if r[3] == 'nontoroidal' and r[6])
    t.notes.append('76 forms = labelings of 8 shapes; %s combinations of the domain, %s recorded, %s excluded through '
                   'a smaller combination (%s of them through one whose violation or nontoroidality is witnessed '
                   'in the folder), %s recorded combinations also contain one; four remain, all T15, one the '
                   'representative of the text (%s)'
                   % (fmt(cnt['domain']), fmt(cnt['recorded']), fmt(cnt['excluded']), fmt(cnt['excluded, witnessed']),
                      fmt(cnt['recorded, already excluded']), hit[0]))
    t.notes.append('nontoroidal statuses left to the genus computations: %d of %d recorded, %d of %d auxiliary'
                   % (rec_nt - rec_nt_wit, rec_nt, aux_nt - aux_nt_wit, aux_nt))


# ---------------------------------------------------------------- E29

def with_twins(adj, nc, ne):
    """adj with nc - 2 vertices p3, p4 on {a,b,c} and ne - 2 vertices r3, r4 on {a,d,e}."""
    H = copy(adj)
    for i in range(3, nc + 1):
        need('p%d' % i not in H, 'the base already has p%d' % i)
        H = add_edges(H, [('p%d' % i, x) for x in 'abc'])
    for i in range(3, ne + 1):
        need('r%d' % i not in H, 'the base already has r%d' % i)
        H = add_edges(H, [('r%d' % i, x) for x in 'ade'])
    return H


def still_excluded(adj, f):
    """The discard of a base remains valid in adj (the base with further twins)."""
    if f['fact'] == FACT_ADJ:
        check_fact(adj, f, {})
        return
    (A, B), (PA, PB) = twin_triples(adj, f, 2)
    core = A | B | set(PA) | set(PB)
    comp = set(f['component'])
    big, stack = set(comp), list(comp)
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in core and y not in big:
                big.add(y)
                stack.append(y)
    nb = {y for x in big for y in adj[x]} - big
    need(len(nb) >= 4, 'with further twins the component has only', len(nb), 'neighbours')


def base_description(adj, d):
    """A base agrees with its description."""
    keys_are(d, {'kind', 'forms', 'form triples', 'singletons', 'cross edges'}, 'the description')
    allowed = {'matching': ({'abd', 'ace'}, {'bd', 'ce'}),
               'path': ({'abd', 'acd', 'ace', 'bcd', 'cde'}, {'bd', 'cd', 'ce'})}
    need(d['kind'] in allowed, 'unknown kind', d['kind'])
    trip_ok, cross_ok = allowed[d['kind']]
    need(len(d['forms']) == len(d['form triples']) <= 2, 'more than two forms, or forms without triples')
    need(set(d['form triples']) | set(d['singletons']) <= trip_ok, 'an attachment triple not allowed for the', d['kind'])
    need(set(d['cross edges']) <= cross_ok and len(set(d['cross edges'])) == len(d['cross edges']),
         'cross edges not allowed for the', d['kind'])
    cross = {frozenset((x, y)) for x in adj for y in adj[x] if x in 'abcde' and y in 'abcde'}
    need(cross == {frozenset(s) for s in d['cross edges']}, 'the edges among a, ..., e are not the cross edges')
    single = {v for v in adj if v.startswith('v_')}
    need(len(d['singletons']) == len(set(d['singletons'])) and single == {'v_' + s for s in d['singletons']},
         'the single vertices are not those of the description')
    used = set('abcdepqrz') | single
    for k, T in enumerate(d['form triples']):
        fv = {v for v in adj if v.startswith('B%d_' % k)}
        need(fv and len(components(adj, fv)) == 1, 'form', k, 'is empty or not connected')
        need({y for x in fv for y in adj[x]} - fv == set(T), 'form', k, 'does not attach exactly at', T)
        used |= fv
    need(set(adj) == used, 'vertices outside the description:', sorted(set(adj) - used))


def e29_bases(B, conf, inp, t):
    """The bases of E29 (see 10 in the module text)."""
    keys_are(B, {'item', 'what', 'counts', 'bases'}, 'bases.json')
    need(B['item'] == 'E29', 'bases.json is not of E29')
    ids = [b['id'] for b in B['bases']]
    need(len(set(ids)) == len(ids), 'bases.json: repeated base id')
    kept, counts = {}, Counter()
    for b in B['bases']:
        def one():
            adj = record_graph(b)
            check_names(adj, conf)
            base_description(adj, b['description'])
            fields = {'id', 'description', 'graph', 'verdict'}
            if b['verdict'] == 'kept':
                need(fields <= set(b) <= fields | {'graph6'}, 'a kept base has fields', sorted(b))
                kept[b['id']] = adj
                counts['kept'] += 1
                return
            need(b['verdict'] == 'excluded by a hypothesis of the case', 'unknown verdict', b['verdict'])
            need(fields | {'witness'} <= set(b) <= fields | {'witness', 'graph6'}, 'a discarded base has fields', sorted(b))
            f = b['witness']
            need(f['fact'] in (FACT_ADJ, FACT_FOUR_NB), 'a base discarded by', f['fact'])
            check_fact(adj, f, inp)
            if f['fact'] == FACT_ADJ:
                need(f['vertex'] not in ('a', 'b', 'c', 'd', 'e'), 'the vertex of degree three is one of a, ..., e')
            for nc, ne in E29_CHOICES[1:]:
                still_excluded(with_twins(adj, nc, ne), f)
            counts[f['fact']] += 1
        guarded('E29 base %s' % b.get('id'), one)
    need(B['counts'] == dict(counts), 'bases.json: counts', B['counts'], 'but the bases give', dict(counts))
    need(fmt(len(ids)) in B['what'] and fmt(counts['kept']) in B['what'],
         'bases.json: its description does not state the numbers of bases and of kept bases')
    t.extra['kept'] = counts['kept']
    t.notes.append('%s bases rebuilt and matched to their descriptions: %s kept, %s discarded with the fact rechecked '
                   'and rechecked again for the other four choices of (n_c, n_e)'
                   % (fmt(len(ids)), fmt(counts['kept']), fmt(len(ids) - counts['kept'])))
    return kept


def e29_formation(adj, c, kept, t):
    need(c.get('base') in kept, 'the base', c.get('base'), 'is not a kept base')
    m = re.fullmatch(r'k=(\d+),l=(\d+)', c.get('case') or '')
    need(m, 'case', c.get('case'), 'does not give (n_c, n_e)')
    nc, ne = int(m.group(1)), int(m.group(2))
    need((nc, ne) in E29_CHOICES, '(n_c, n_e) =', (nc, ne), 'is not one of the five choices')
    need(same_graph(adj, with_twins(kept[c['base']], nc, ne)), 'the graph is not its base with', nc - 2, 'and',
         ne - 2, 'further twins')
    need((c['base'], nc, ne) not in t.e29_pairs, 'two candidates from base', c['base'], 'with', (nc, ne))
    t.e29_pairs.add((c['base'], nc, ne))


# ---------------------------------------------------------------- summary and text

def quoted_counts(item, t):
    """The counts the quoted text of the paper states for the item."""
    per_case = sorted(set(t.case_count.values()))
    if item in ('E20', 'E21', 'E22', 'E23', 'E25'):
        return [t.n]
    if item in ('E24', 'E26', 'E27', 'E28'):
        return per_case
    if item == 'E29':
        return [t.extra['kept'], t.n, t.verdicts['toroidal'], t.verdicts['proper nontoroidal minor'], sum(t.iso.values())]
    return []


def check_text(item, inp, t):
    P = inp.get('paper text')
    need(isinstance(P, list) and P, 'input.json quotes no text of the paper')
    cite = re.compile(r'\\supp\{(?:[^}]*[ ,~])?%s[,}]' % item)
    for p in P:
        need(cite.search(p['text']), 'a quoted passage does not cite', item)
    text = '\n'.join(p['text'] for p in P)
    nums = quoted_counts(item, t)
    for n in nums:
        need('$%s$' % fmt(n).replace(',', '{,}') in text, 'the quoted text does not state the count', fmt(n))
    return nums


def check_summary(item, inp, s, t):
    exp = {'item': item, 'where': inp['where'], 'candidates': t.n, 'verdicts': dict(t.verdicts),
           'isomorphic to T_i (copies)': dict(t.iso),
           'graphs T_i not excluded by a hypothesis of the case': tsort(t.iso_free)}
    if t.facts:
        exp['also excluded by a hypothesis of the case'] = dict(t.facts)
    if t.iso_excl:
        exp['graphs T_i excluded by a hypothesis of the case'] = tsort(t.iso_excl)
    if len(t.case_count) > 1:
        exp['by case'] = {k: {'candidates': n, 'graphs T_i not excluded by a hypothesis of the case': tsort(t.case_free[k])}
                          for k, n in t.case_count.items()}
    if item == 'E19':
        exp['combinations'] = None
    for k in s:
        need(k in exp, 'summary.json has a field the checker does not recompute:', k)
    for k, x in exp.items():
        need(k in s, 'summary.json lacks the field', k)
        if k == 'combinations':
            cm = s[k]
            want = {'in the domain': t.extra['in the domain'], 'recorded': t.extra['recorded'],
                    'excluded through a recorded smaller combination with a violation or nontoroidal status': t.extra['excluded'],
                    'recorded stars-only combinations used only for these exclusions': t.extra['auxiliary']}
            need(set(cm) == set(want) | {'checked'} and isinstance(cm['checked'], str),
                 'summary.json: the counts of combinations have fields', sorted(cm))
            for kk, n in want.items():
                need(cm[kk] == n, 'summary.json: combinations,', kk, 'is', cm[kk], 'but the lattice gives', n)
        else:
            need(s[k] == x, 'summary.json:', k, 'is', s[k], 'but the records give', x)


def check_construction_counts(inp, t):
    """Where input.json gives the number of bridges the construction found for a
    case, the case has that many candidates for each of its variants."""
    unused = []
    for k, x in inp.get('cases of the construction', {}).items():
        n = x.get('parts up to rooted isomorphism', x.get('parts up to label-preserving isomorphism'))
        if n is None:
            continue
        if k not in t.case_count:
            unused.append(k)
            continue
        variants = {v for (case, v) in t.variant_count if case == k}
        for v in variants:
            need(t.variant_count[(k, v)] == n, 'case', k, 'variant', v, 'has', t.variant_count[(k, v)],
                 'candidates, but input.json gives', n, 'bridges')
    if unused:
        t.notes.append('input.json also describes cases with no candidates: ' + ', '.join(sorted(unused)))


# ---------------------------------------------------------------- the item

def check(root, item):
    """Check one of E19-E29; returns a one-line report."""
    conf = CONF[item]
    folder = os.path.join(root, 'enumerations', item)
    need(os.path.isdir(folder), 'there is no folder', folder)

    def read(name):
        p = find(os.path.join(folder, name))
        need(p is not None, '%s: there is no %s' % (item, name))
        return load(p)
    inp, summ = read('input.json'), read('summary.json')
    need(inp.get('item') == item and summ.get('item') == item, 'input.json or summary.json is not of', item)
    cat = catalogue(root)
    t = Tally()
    ctx = {}
    if item == 'E19':
        ctx['forms'] = guarded('E19 forms', e19_forms, inp)
    if item == 'E29':
        ctx['bases'] = e29_bases(read('bases.json'), conf, inp, t)
    cf, vf = read('candidates.json'), read('verdicts.json')
    need(cf.get('item') == item and vf.get('item') == item, 'candidates.json or verdicts.json is not of', item)
    cands, verds = cf['candidates'], vf['verdicts']
    ids = [c['id'] for c in cands]
    need(len(set(ids)) == len(ids), item + ': repeated candidate id')
    vby = {}
    for v in verds:
        need(v['id'] not in vby, item + ': two verdicts for', v['id'])
        vby[v['id']] = v
    missing = [i for i in ids if i not in vby]
    need(not missing, item + ': candidates without a verdict:', missing[:5])
    stray = set(vby) - set(ids)
    need(not stray, item + ': verdicts for unknown candidates:', sorted(stray)[:5])
    del verds
    for c in cands:
        guarded('%s %s' % (item, c['id']), check_candidate, c, vby[c['id']], conf, inp, cat, ctx, t, item)
    del cands, vby
    if item == 'E19':
        e19_lattice(inp, ctx['forms'], cat, t)
    if item == 'E29':
        want = {(b, nc, ne) for b in ctx['bases'] for nc, ne in E29_CHOICES}
        need(t.e29_pairs == want, 'E29: the candidates are not the kept bases with the five choices of (n_c, n_e)')
    if conf['formed'] == 'bridge':
        check_construction_counts(inp, t)
    guarded(item + ' summary.json', check_summary, item, inp, summ, t)
    nums = guarded(item + ' input.json', check_text, item, inp, t)
    return report(item, t, nums)


def report(item, t, nums):
    v = t.verdicts
    seg = []
    if v.get('toroidal'):
        seg.append('%s toroidal' % fmt(v['toroidal']))
    if v.get('K5 minor'):
        seg.append('%s K5 minors' % fmt(v['K5 minor']))
    if v.get('cut of at most two vertices'):
        seg.append('%s cuts' % fmt(v['cut of at most two vertices']))
    if t.iso:
        seg.append('%d isomorphisms (%s)' % (sum(t.iso.values()), ', '.join(tsort(t.iso))))
    pm = v.get('proper nontoroidal minor', 0)
    if pm:
        seg.append('%s proper minors, named by the edge, left to the genus computations' % fmt(pm))
    out = '%s candidates, one verdict each [%s]' % (fmt(t.n), '; '.join(seg))
    if t.facts:
        out += '; %s facts rechecked' % fmt(sum(t.facts.values()))
    if t.exteriors:
        out += '; bridges attach at the roots, one exterior per case and variant (%d)' % len(t.exteriors)
    if t.completions:
        out += '; %d completions (planar: %s)' % (t.completions, ', '.join('%s %s' % x for x in t.planar) or 'none')
    for n in t.notes:
        out += '; ' + n
    out += '; summary agrees'
    if nums:
        out += ', and the quoted text states %s' % ', '.join(fmt(n) for n in nums)
    return out
