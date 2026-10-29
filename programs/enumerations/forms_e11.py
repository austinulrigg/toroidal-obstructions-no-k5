"""The forms used by E11 (Proposition 5.14, the cut X = {a, b, w}, Section 5.4.1).

X = {a, b, w} is a three-cut of J; C and D are the components of B - X
containing c and d (D_c and D_d in the paper; the program writes P for B). For C, N_C is its set of neighbors in X (w is always
one: each component of C - c meets c and w, and C = {c} forces cw), and
H_C = C + {a, b, w} + ab + aw + bw. By Lemma 2.5
the tests depend on C only through N_C and the planarity of the four graphs
H_C + E, E a subset of {ac, bc}; H_C itself is planar, as the component I
on abw can replace the triangle. So C may be replaced by any graph with the
same signature (N_C, four planarity values) in which every component of
C - c meets c and w, and an obstruction uses a C with no proper minor of
the same signature.

c_forms() generates a set of C containing all such minimal ones:
  - no subdivision (all four planar): at most two single-vertex components
    on c, w and some of a, b, with or without cw;
  - a K3,3 subdivision of H_C + E0 for E0 = {ac}, {bc} or {ac, bc}, which
    must use E0 (the smaller completions are planar): Construction 3.17
    with roots a, b, c, w, added edges ab, aw, bw and E0, required E0, and
    the original pair cw taken as an edge and as a subdivided segment;
  - H_C + ac and H_C + bc both nonplanar: Lemma 2.7
    makes (H_C + ac + bc)/ab nonplanar; Construction 3.17 on the quotient
    with roots r, c, w, added rc, rw, required rc, original cw, and every
    lift of each incidence at r to a, b or both.
In each, the interior pieces are divided into components in every way and
joined by trees (Lemma 5.8); a component without
an edge to w, or to c, gets one in every way; and a or b without an
interior neighbor gets nothing, one edge to an interior vertex, or a further
single-vertex component on c, w and it (or one on c, w, a, b), as Lemma
5.8 allows two further components for neighbors at a, b.
Then every form with a one-step proper minor of the same signature that
still has every component of C - c meeting c and w is removed (a superset
of the minimal forms remains), and the rest are taken up to isomorphism
fixing a, b, c, w.

u_forms() gives the graphs U of Proposition 5.7: the
candidate shapes with every labeling, kept when U + ab is planar and U + aw
or U + bw is planar (sixteen graphs).

r_forms() gives the components R containing both c and d when there are two
single vertices on abw (case (c)).
"""

import itertools

import networkx as nx

import construction as c2
import psi1_domain as sd

ROOTS = ("a", "b", "c", "w")


def planar(G, extra=()):
    H = G.copy()
    H.add_edges_from(extra)
    return nx.check_planarity(H)[0]


def inner(H, roots=ROOTS):
    return [v for v in H if v not in roots]


def n_c(C):
    out = {x for x in "abw" if any(v not in ROOTS for v in C[x])}
    if C.has_edge("c", "w"):
        out.add("w")
    return "".join(sorted(out))


def signature(C):
    H = C.copy()
    H.add_edges_from([("a", "b"), ("a", "w"), ("b", "w")])
    return (n_c(C), planar(H), planar(H, [("a", "c")]), planar(H, [("b", "c")]),
            planar(H, [("a", "c"), ("b", "c")]))


def valid(C):
    for x, y in C.edges():
        if x in ROOTS and y in ROOTS and {x, y} != {"c", "w"}:
            return False
    comps = c2.interior_components(C, ROOTS)
    for comp in comps:
        if not any(C.has_edge(v, "c") for v in comp) or not any(C.has_edge(v, "w") for v in comp):
            return False
    if not (comps or C.has_edge("c", "w")):
        return False
    H = C.copy()
    H.add_edges_from([("a", "b"), ("a", "w"), ("b", "w")])
    return planar(H)          # H_C is planar: the component(s) on abw replace the triangle


def allowed_pair(x, y):
    return {x, y} == {"c", "w"}


def fix_components(H, comps, need=("w", "c")):
    """Every way of giving each component an edge to each vertex of need
    that it misses."""
    tasks = []
    for comp in comps:
        for x in need:
            if not any(H.has_edge(v, x) for v in comp):
                tasks.append([(x, v) for v in comp])
    for pick in itertools.product(*tasks):
        M = H.copy()
        M.add_edges_from(pick)
        yield M


def ab_options(H):
    """a or b without an interior neighbor: nothing, an edge to an interior
    vertex, or a further vertex on c, w and it; both: also one further
    vertex on c, w, a, b."""
    ins = inner(H)
    lonely = [x for x in "ab" if not any(v in ins for v in H[x])]
    per = []
    for x in lonely:
        per.append([None] + [("edge", x, v) for v in ins] + [("vertex", (x,))])
    seen = []
    for pick in itertools.product(*per):
        seen.append(pick)
    if len(lonely) == 2:
        seen.append((("vertex", ("a", "b")),))
    for pick in seen:
        M = H.copy()
        k = 0
        for opt in pick:
            if opt is None:
                continue
            if opt[0] == "edge":
                M.add_edge(opt[1], opt[2])
            else:
                f = "f%d" % k
                k += 1
                M.add_edges_from([(f, "c"), (f, "w")] + [(f, x) for x in opt[1]])
        yield M


def no_subdivision_forms():
    out = []
    sets = [(), ("a",), ("b",), ("a", "b")]
    for cw in (False, True):
        for k in range(3):
            for S in itertools.combinations_with_replacement(range(4), k):
                C = nx.Graph()
                C.add_nodes_from(ROOTS)
                if cw:
                    C.add_edge("c", "w")
                for j, s in enumerate(S):
                    C.add_edges_from([("f%d" % j, "c"), ("f%d" % j, "w")] + [("f%d" % j, x) for x in sets[s]])
                out.append((C, {"origin": "no subdivision"}))
    return out


def subdivision_forms(E0, stats):
    added = [("a", "b"), ("a", "w"), ("b", "w")] + list(E0)
    for K, part, rec in c2.retained(ROOTS, added=added, required=E0, original=[("c", "w")]):
        stats["retained"] = stats.get("retained", 0) + 1
        for H, comps in c2.grouped(part, ROOTS):
            for M in fix_components(H, comps):
                for C in ab_options(M):
                    stats["forms"] = stats.get("forms", 0) + 1
                    yield C, {"origin": "subdivision of H_C + %s" % "+".join("".join(e) for e in E0)}


def lifts_ab(Q):
    """Lift each incidence at r of the quotient part Q to a, b or both."""
    ins = [v for v in Q if v not in ("r", "c", "w")]
    at_r = sorted(v for v in Q["r"] if v in ins)
    base = Q.copy()
    base.remove_node("r")
    base.add_nodes_from(ROOTS)
    for choice in itertools.product(("a", "b", "ab"), repeat=len(at_r)):
        M = base.copy()
        for x, v in zip(choice, at_r):
            M.add_edges_from((y, v) for y in x)
        yield M


def quotient_forms(stats):
    QR = ("r", "c", "w")
    for K, part, rec in c2.retained(QR, added=[("r", "c"), ("r", "w")], required=[("r", "c")],
                                    original=[("c", "w")]):
        stats["retained"] = stats.get("retained", 0) + 1
        for H, comps in c2.grouped(part, QR):
            for M in fix_components(H, comps):
                for L in lifts_ab(M):
                    for C in ab_options(L):
                        stats["forms"] = stats.get("forms", 0) + 1
                        yield C, {"origin": "quotient (H_C + ac + bc)/ab"}


def c_forms():
    stats = {}
    raw = no_subdivision_forms()
    for E0 in ([("a", "c")], [("b", "c")], [("a", "c"), ("b", "c")]):
        st = {}
        raw += list(subdivision_forms(E0, st))
        stats["subdivision " + "+".join("".join(e) for e in E0)] = st
    st = {}
    raw += list(quotient_forms(st))
    stats["quotient"] = st
    stats["generated"] = len(raw)
    ok = [(C, i) for C, i in raw if valid(C)]
    stats["valid"] = len(ok)
    ok = c2.dedupe(ok, ROOTS)
    stats["valid up to isomorphism"] = len(ok)
    kept = c2.local_minimal(ok, ROOTS, signature, valid, allowed_pair)
    stats["without a one-step minor of the same signature"] = len(kept)
    return kept, stats


# ----------------------------------------------------------------------------
# U: Proposition 5.7

def u_forms():
    A = ("a", "b", "w")
    shapes = []
    shapes.append(("one vertex", sd.form_edges("single", "a", "b", "w")))
    for t in A:
        o = [x for x in A if x != t]
        shapes.append(("T(%s)" % t, sd.form_edges("Tx", t, o[0], o[1])))
    for x, y in itertools.permutations(A, 2):
        z = [v for v in A if v not in (x, y)][0]
        shapes.append(("S(%s;%s)" % (x, y), sd.form_edges("Sxy", x, y, z)))
    for x in A:
        y, z = [v for v in A if v != x]
        shapes.append(("Q(%s;%s,%s)" % (x, y, z), sd.form_edges("Q", y, z, x)))
    for x, y, z in itertools.permutations(A):
        shapes.append(("four-cycle, %s pendant at the pair of %s" % (x, z), sd.form_edges("FC1", x, y, z)))
    k23 = [("@m%d" % i, "@n%d" % j) for i in range(2) for j in range(3)] + \
        [("@n0", "a"), ("@n1", "b"), ("@n2", "w")]
    shapes.append(("pendant K2,3", k23))
    out = []
    for name, E in shapes:
        U = nx.Graph()
        U.add_nodes_from(A)
        U.add_edges_from((x.lstrip("@") if x.startswith("@") else x, y.lstrip("@") if y.startswith("@") else y)
                         for x, y in E)
        if planar(U, [("a", "b")]) and (planar(U, [("a", "w")]) or planar(U, [("b", "w")])):
            out.append((name, U))
    return out


# ----------------------------------------------------------------------------
# R containing c and d, two single vertices on abw (case (c))

R_ROOTS = ("a", "c", "d", "w")


def r_signature(R):
    nb = "".join(sorted(x for x in "aw" if any(v not in R_ROOTS for v in R[x]) or
                        (x == "w" and (R.has_edge("c", "w") or R.has_edge("d", "w")))))
    return (nb, planar(R, [("a", "c"), ("a", "d"), ("a", "w")]))


def r_valid(R):
    for x, y in R.edges():
        if x in R_ROOTS and y in R_ROOTS and {x, y} not in ({"c", "w"}, {"d", "w"}):
            return False
    comps = c2.interior_components(R, R_ROOTS)
    for comp in comps:
        if not any(R.has_edge(v, "w") for v in comp):
            return False
        if not any(R.has_edge(v, "c") or R.has_edge(v, "d") for v in comp):
            return False
    S = R.subgraph([v for v in R if v not in ("a", "w")])
    return nx.has_path(S, "c", "d") if "c" in S and "d" in S else False


def r_allowed_pair(x, y):
    return {x, y} in ({"c", "w"}, {"d", "w"})


def r_forms():
    stats = {"retained": 0, "forms": 0}
    raw = []
    for K, part, rec in c2.retained(R_ROOTS, added=[("a", "c"), ("a", "d"), ("a", "w")],
                                    original=[("c", "w"), ("d", "w")]):
        if not (K.has_edge("a", "c") or K.has_edge("a", "d")):
            continue
        stats["retained"] += 1
        for H, comps in c2.grouped(part, R_ROOTS):
            tasks = []
            for comp in comps:
                if not any(H.has_edge(v, "w") for v in comp):
                    tasks.append([("w", v) for v in comp])
                if not any(H.has_edge(v, "c") or H.has_edge(v, "d") for v in comp):
                    tasks.append([(x, v) for x in "cd" for v in comp])
            for pick in itertools.product(*tasks):
                M = H.copy()
                M.add_edges_from(pick)
                ins = inner(M, R_ROOTS)
                S = M.subgraph([v for v in M if v not in ("a", "w")])
                join = [None] if nx.has_path(S, "c", "d") else [("vertex", ("c", "d"))]
                lonely_a = not any(v in ins for v in M["a"])
                aopt = [None] + ([("edge", v) for v in ins] + [("vertex", ("a", "c")), ("vertex", ("a", "d"))]
                                 if lonely_a else [])
                for j, ao in itertools.product(join, aopt):
                    R = M.copy()
                    if j is not None:
                        R.add_edges_from([("g0", "c"), ("g0", "d"), ("g0", "w")])
                    if ao is not None:
                        if ao[0] == "edge":
                            R.add_edge("a", ao[1])
                        else:
                            R.add_edges_from([("g1", "a"), ("g1", "w"), ("g1", ao[1][1])])
                    stats["forms"] += 1
                    raw.append((R, {"origin": "subdivision of R + {a,w} + ac + ad + aw"}))
    ok = [(R, i) for R, i in raw if r_valid(R)]
    stats["valid"] = len(ok)
    ok = c2.dedupe(ok, R_ROOTS)
    stats["valid up to isomorphism"] = len(ok)
    kept = c2.local_minimal(ok, R_ROOTS, r_signature, r_valid, r_allowed_pair)
    stats["without a one-step minor of the same signature"] = len(kept)
    return kept, stats
