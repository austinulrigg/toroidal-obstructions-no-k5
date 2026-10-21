"""The forms used by E14 (Proposition 5.13, the cut X = {a, u, v} with two
other components, Section 5.4.2).

X = {a, u, v}, p and q on X, Psi2 = twins on abc and abd. The other
components of P - X are R, containing x and y, and D, containing r, with
{x, y, r} = {b, c, d}; up to c <-> d the partitions are
  P1: R contains b, c and D contains d;
  P2: R contains c, d and D contains b.
Neighborhoods N(R), N(D) are subsets of X with one or two vertices,
meeting {u, v}, and together covering u, v.

The four tests and Lemma 2.5. With p, q adjacent to
all of X, a test is nonplanar iff a third bridge is adjacent to all of X,
or a bridge B with two attachments has B + (their edge) nonplanar, or a
bridge with one attachment is nonplanar. A test edge between R and D joins
them into one bridge, which then meets a (every such test also adds an edge
at a, or the joined bridge already meets it) and u, v, so that test fails.
  P1: ac+bd and ab+bc+bd join R, D and fail. ad+bc gives R + bc and D + ad;
      ab+ac+ad gives R + ab + ac and D + ad.
  P2: ac+bd, ad+bc and ab+bc+bd fail. ab+ac+ad gives R + ac + ad and D + ab.
Hence the tests depend on R only through
  P1: N(R), planarity of R + bc + K[N(R)], and of R + ab + ac + K[N(R)+a]
      (when N(R) + a has two vertices; with three it is a third full bridge);
  P2: N(R), planarity of R + ac + ad + K[N(R)+a];
and on D only through N(D) and planarity of D + ar + K[N(D)+a].
These are the signatures below. Replacing R or D by a proper minor with the
same signature, which keeps every component meeting u or v and the
component connected, gives a proper nontoroidal minor of G. So in an
obstruction R and D have no such minor, and the generated sets need only
contain those.

R generators (the text's constructions, with mirrors u <-> v):
  G0  no subdivision: a vertex t on x, y and some of a, u, v; possibly a
      second single vertex; any original edges xu, xv, yu, yv (the seven
      path forms of the text are among these);
  G1  P1, N(R) = uv, R + uv + bc nonplanar: Construction 3.16 with roots
      b, c, u, v, added bc, uv, required bc, original pairs bu, bv, cu, cv;
  G2  P1, the tests of equation (2): roots a, b, c, u, added bc, au,
      required bc, original bu, cu, the interior with u connected;
  G3  P2, N(R) in {u, au}, R + au + ac + ad nonplanar: roots a, c, d, u,
      added au, ac, ad, the subdivision using ac or ad, original cu, du.
D generators: the nine forms on (s, t, r) for st in {au, av, uv}, two single
vertices on rst, or empty interior, with original rs, rt other than ar.
"""

import itertools

import networkx as nx

import construction as c2
import psi1_domain as sd

PARTITIONS = {"P1": (("b", "c"), "d"), "P2": (("c", "d"), "b")}
X = ("a", "u", "v")


def planar(G, extra=()):
    H = G.copy()
    H.add_edges_from(extra)
    return nx.check_planarity(H)[0]


def nbhd(H, members):
    out = set()
    for m in members:
        if m in H:
            out |= {z for z in H[m] if z in X}
    for v in H:
        if v not in X and v not in "abcd" and v in H:
            out |= {z for z in H[v] if z in X}
    return "".join(sorted(out))


def k_on(S):
    return list(itertools.combinations(sorted(S), 2))


def r_signature(R, part):
    (x, y), r = PARTITIONS[part]
    N = nbhd(R, (x, y))
    if part == "P1":
        s1 = planar(R, [("b", "c")] + k_on(N))
        Na = set(N) | {"a"}
        s2 = planar(R, [("a", "b"), ("a", "c")] + k_on(Na)) if len(Na) <= 2 else None
        return (N, s1, s2)
    Na = set(N) | {"a"}
    s = planar(R, [("a", "c"), ("a", "d")] + k_on(Na)) if len(Na) <= 2 else None
    return (N, s)


def d_signature(D, part):
    (x, y), r = PARTITIONS[part]
    N = nbhd(D, (r,))
    Na = set(N) | {"a"}
    return (N, planar(D, [("a", r)] + k_on(Na)) if len(Na) <= 2 else None)


def interior(H, members):
    return [v for v in H if v not in X and v not in "abcd"]


def component_ok(H, members, allowed):
    """Edges between named vertices only in allowed; the component (members
    and interior) connected; every interior component meets u or v; the
    neighborhood has one or two vertices and meets u, v."""
    named = set(X) | set(members)
    for e in H.edges():
        if e[0] in named and e[1] in named and frozenset(e) not in allowed:
            return False
        if any(z in "abcd" and z not in members for z in e if z not in X):
            return False
    ins = interior(H, members)
    S = H.subgraph(list(members) + ins)
    if S.number_of_nodes() == 0 or not nx.is_connected(S):
        return False
    for comp in nx.connected_components(H.subgraph(ins)):
        if not any(z in ("u", "v") for v in comp for z in H[v]):
            return False
    N = nbhd(H, members)
    return 1 <= len(N) <= 2 and bool(set(N) & {"u", "v"})


def r_allowed(part):
    (x, y), r = PARTITIONS[part]
    return {frozenset(p) for p in itertools.product((x, y), ("u", "v"))}


def d_allowed(part):
    (x, y), r = PARTITIONS[part]
    return {frozenset((r, "u")), frozenset((r, "v"))}


def mirror(H):
    return nx.relabel_nodes(H, {"u": "v", "v": "u"})


def lonely_edges(H, who, ins):
    """Every way of giving each vertex of who without an interior neighbor
    one edge to the interior."""
    lonely = [z for z in who if z in H and not any(v in ins for v in H[z])] + \
             [z for z in who if z not in H]
    for nbrs in itertools.product(ins, repeat=len(lonely)):
        M = H.copy()
        M.add_edges_from(zip(lonely, nbrs))
        yield M


# ----------------------------------------------------------------------------

def g0(part):
    (x, y), r = PARTITIONS[part]
    subsets = [S for k in range(4) for S in itertools.combinations(X, k)]
    second = [None] + [(A, B) for k in (1, 2) for A in itertools.combinations((x, y), k)
                       for j in (1, 2, 3) for B in itertools.combinations(X, j)]
    origs = [E for k in range(5) for E in itertools.combinations(list(itertools.product((x, y), ("u", "v"))), k)]
    for S in subsets:
        for sec in second:
            for E in origs:
                H = nx.Graph()
                H.add_nodes_from((x, y))
                H.add_edges_from([("t", x), ("t", y)] + [("t", z) for z in S])
                if sec is not None:
                    H.add_edges_from([("t2", z) for z in sec[0] + sec[1]])
                H.add_edges_from(E)
                yield H, {"origin": "G0 no subdivision"}


def subdivision_forms(roots, added, required, original, members, need_uv, connect_with_u, tag, use_one_of=None):
    for K, part, rec in c2.retained(roots, added=added, required=required, original=original):
        if use_one_of and not any(K.has_edge(*e) for e in use_one_of):
            continue
        if connect_with_u:
            groupings = [(H, comps) for H, comps in c2.grouped(part, tuple(z for z in roots if z != "u"))
                         if len(comps) == 1]
        else:
            groupings = list(c2.grouped(part, roots))
        for H, comps in groupings:
            # every interior component must meet u or v
            tasks = []
            for comp in comps:
                comp = [v for v in comp if v not in X]
                if comp and not any(z in ("u", "v") for v in comp for z in H[v]):
                    tasks.append([(z, v) for z in ("u", "v") for v in comp])
            for pick in itertools.product(*tasks):
                M = H.copy()
                M.add_edges_from(pick)
                ins = interior(M, members)
                for L in lonely_edges(M, members, ins):
                    S = L.subgraph(list(members) + ins)
                    joins = [None] if nx.has_path(S, members[0], members[1]) else \
                        [("g0", members[0]), ("g0", members[1])]
                    for j in ([None] if joins == [None] else [z for z in ("u", "v")]):
                        L2 = L.copy()
                        if j is not None:
                            L2.add_edges_from([("g0", members[0]), ("g0", members[1]), ("g0", j)])
                        ins2 = interior(L2, members)
                        lacking = [z for z in need_uv if not any(w in L2[z] for w in list(ins2) + list(members))] \
                            if need_uv else []
                        opts = [[None] + [("edge", z, w) for w in ins2] + [("vertex", z, m) for m in members]
                                for z in lacking]
                        for pk in itertools.product(*opts):
                            L3 = L2.copy()
                            for k, o in enumerate(pk):
                                if o is None:
                                    continue
                                if o[0] == "edge":
                                    L3.add_edge(o[1], o[2])
                                else:
                                    L3.add_edges_from([("h%d" % k, o[1]), ("h%d" % k, o[2])])
                            yield L3, {"origin": tag}


def r_forms(part):
    (x, y), r = PARTITIONS[part]
    raw = list(g0(part))
    if part == "P1":
        raw += list(subdivision_forms(("b", "c", "u", "v"), [("b", "c"), ("u", "v")], [("b", "c")],
                                      [("b", "u"), ("b", "v"), ("c", "u"), ("c", "v")], ("b", "c"),
                                      ("u", "v"), False, "G1 subdivision of R + uv + bc"))
        raw += list(subdivision_forms(("a", "b", "c", "u"), [("b", "c"), ("a", "u")], [("b", "c")],
                                      [("b", "u"), ("c", "u")], ("b", "c"), None, True,
                                      "G2 hard tests, subdivision of H + bc"))
    else:
        raw += list(subdivision_forms(("a", "c", "d", "u"), [("a", "u"), ("a", "c"), ("a", "d")], [],
                                      [("c", "u"), ("d", "u")], ("c", "d"), None, False,
                                      "G3 subdivision of R + au + ac + ad",
                                      use_one_of=[("a", "c"), ("a", "d")]))
    raw += [(mirror(H), dict(i, mirrored=True)) for H, i in raw if i["origin"] != "G0 no subdivision"]
    allowed = r_allowed(part)
    ok = []
    for H, i in raw:
        H = H.copy()
        H.remove_nodes_from([z for z in X if z in H and H.degree(z) == 0])
        if component_ok(H, (x, y), allowed):
            ok.append((H, i))
    roots = (x, y, "a", "u", "v")
    norm = []
    for H, i in ok:
        H = H.copy()
        H.add_nodes_from(roots)
        norm.append((H, i))
    stats = {"generated": len(raw), "valid": len(norm)}
    norm = c2.dedupe(norm, roots)
    stats["valid up to isomorphism"] = len(norm)
    kept = c2.local_minimal(norm, roots, lambda H: r_signature(H, part),
                            lambda H: component_ok(H, (x, y), allowed),
                            lambda p, q: frozenset((p, q)) in allowed)
    stats["without a one-step minor of the same signature"] = len(kept)
    return kept, stats


def d_forms(part):
    (x, y), r = PARTITIONS[part]
    raw = []
    for s, t in (("a", "u"), ("a", "v"), ("u", "v")):
        bodies = [(f, sd.form_edges(f, s, t, r)) for f in
                  ("single", "Tx", "Ty", "Tz", "Sxy", "Syx", "Q", "FC1", "FC2")]
        bodies.append(("two single vertices", [(z, w) for z in ("@k1", "@k2") for w in (s, t, r)]))
        edges = [(r, z) for z in (s, t) if z != "a"]
        for name, body in bodies:
            for k in range(len(edges) + 1):
                for E in itertools.combinations(edges, k):
                    H = nx.Graph()
                    H.add_edges_from((p.replace("@", "k_"), q.replace("@", "k_")) for p, q in body)
                    H.add_edges_from(E)
                    raw.append((H, {"origin": "%s on %s%s%s" % (name, r, s, t)}))
    for E in ([(r, "u")], [(r, "v")], [(r, "u"), (r, "v")]):
        H = nx.Graph()
        H.add_edges_from(E)
        raw.append((H, {"origin": "empty interior"}))
    allowed = d_allowed(part)
    roots = (r, "a", "u", "v")
    norm = []
    for H, i in raw:
        if component_ok(H, (r,), allowed):
            H = H.copy()
            H.add_nodes_from(roots)
            norm.append((H, i))
    stats = {"generated": len(raw), "valid": len(norm)}
    norm = c2.dedupe(norm, roots)
    kept = c2.local_minimal(norm, roots, lambda H: d_signature(H, part),
                            lambda H: component_ok(H, (r,), allowed),
                            lambda p, q: frozenset((p, q)) in allowed)
    stats["without a one-step minor of the same signature"] = len(kept)
    return kept, stats
