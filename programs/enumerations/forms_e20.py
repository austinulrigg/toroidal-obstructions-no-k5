"""The eight forms of Section 5.6 (Lemma 5.26) and the structural conditions
of Sections 5.5 and 5.6 used to discard configurations in E20.

A form is a bridge U on a triple {x, y, z} with U + xy + xz planar, built
from Proposition 3.10 (prop:three-root-blocks): the blocks of U - x are edges and
four-cycles in a sequence from y to z; in each four-cycle the two vertices
meeting the neighboring blocks (or equal to y, z) are opposite and the
other two have degree three and are adjacent to x; paths of edge blocks
have length 0, 1 or 2, a middle vertex of a length-two path being adjacent
to x. The text's eight possibilities:
  star         one vertex adjacent to x, y, z
  one cycle    end lengths (0,1), (1,0)   = T(z), T(y)
               end lengths (0,2), (2,0)   = S(x;y), S(x;z)
               end lengths (1,1)          = Q(x;y,z)
  two cycles   end lengths (0,0), joined by a path of length 1 or 2.
"""
import itertools

import networkx as nx


def cycle_chain(x, y, z, ends, joins):
    """Blocks of U - x: end path from y of length ends[0], cycles joined by
    paths of lengths joins, end path to z of length ends[1]. Interior
    vertices are named i0, i1, ... ."""
    G = nx.Graph()
    count = [0]

    def new():
        v = "i%d" % count[0]
        count[0] += 1
        return v

    def path(a, length):
        """Path of the given length from a; returns its far end."""
        if length == 0:
            return a
        if length == 1:
            b = new()
            G.add_edge(a, b)
            return b
        m, b = new(), new()
        G.add_edges_from([(a, m), (m, b), (m, x)])
        return b

    def cycle(s1):
        p, q, s2 = new(), new(), new()
        G.add_edges_from([(s1, p), (p, s2), (s2, q), (q, s1), (p, x), (q, x)])
        return s2

    cur = path(y, ends[0])
    for k in range(len(joins) + 1):
        cur = cycle(cur)
        if k < len(joins):
            cur = path(cur, joins[k])
    # end path to z: identify the far end with z
    last = ends[1]
    if last == 0:
        G = nx.relabel_nodes(G, {cur: z})
    elif last == 1:
        G.add_edge(cur, z)
    else:
        m = new()
        G.add_edges_from([(cur, m), (m, z), (m, x)])
    return G


def forms_xyz(x, y, z):
    """The eight forms for the attachment x and ends y, z."""
    star = nx.Graph()
    star.add_edges_from([("i0", x), ("i0", y), ("i0", z)])
    out = [("star", star)]
    for name, ends in (("T(z)", (0, 1)), ("T(y)", (1, 0)), ("S(x;y)", (0, 2)),
                       ("S(x;z)", (2, 0)), ("Q(x;y,z)", (1, 1))):
        out.append((name, cycle_chain(x, y, z, ends, [])))
    for j in (1, 2):
        out.append(("two cycles, join %d" % j, cycle_chain(x, y, z, (0, 0), [j])))
    return out


TRIPLES = ("abc", "abd", "acd", "bcd")
TWINS = {"x1": "abc", "x2": "abc", "y1": "abd", "y2": "abd"}


def form_library():
    """The 76 labeled forms: every form for every choice of x and order of
    y, z on each triple, up to label-preserving isomorphism. Returns a list
    of dicts with id, triple, name, (x,y,z) and edges (interior i0, i1, ...)."""
    import enumeration as el
    lib = []
    for T in TRIPLES:
        seen = []
        for perm in itertools.permutations(T):
            for name, U in forms_xyz(*perm):
                if any(el.rooted_isomorphic(U, V, tuple(T)) for V in seen):
                    continue
                seen.append(U)
                lib.append({"id": "F%02d" % len(lib), "triple": T, "name": name,
                            "x,y,z": list(perm),
                            "edges": sorted(sorted(e) for e in U.edges())})
    return lib


def config_graph(lib, config, cd):
    """Psi2 (x1, x2 on abc; y1, y2 on abd), the forms in config (a tuple of
    indices into lib, repetition allowed), and cd if requested. The interior
    vertex i of the k-th form is named f<k>_<i>."""
    G = nx.Graph()
    G.add_nodes_from("abcd")
    for t, nb in TWINS.items():
        G.add_edges_from((t, v) for v in nb)
    for k, idx in enumerate(config):
        for u, v in lib[idx]["edges"]:
            G.add_edge(u if u in "abcd" else "f%d_%s" % (k, u[1:]),
                       v if v in "abcd" else "f%d_%s" % (k, v[1:]))
    if cd:
        G.add_edge("c", "d")
    return G


def twin_triples(G):
    """Neighborhoods shared by at least two vertices of degree three."""
    by = {}
    for v in G:
        if G.degree(v) == 3:
            by.setdefault(frozenset(G[v]), []).append(v)
    return {N: sorted(vs, key=str) for N, vs in by.items() if len(vs) >= 2}


def structural_violation(G):
    """The structural conditions of Sections 5.5 and 5.6 that persist when
    forms are added (see generate_E20.py). Returns a witness or None."""
    tt = twin_triples(G)
    trip = sorted(tt, key=lambda s: sorted(map(str, s)))
    for A, B in itertools.combinations(trip, 2):
        if not (A & B):
            return {"condition": "disjoint twin triples",
                    "triples": [sorted(A), sorted(B)], "twins": [tt[A][:2], tt[B][:2]]}
    for A, B, C in itertools.combinations(trip, 3):
        if len(A & B & C) >= 2:
            return {"condition": "three twin triples share a pair",
                    "triples": [sorted(A), sorted(B), sorted(C)],
                    "twins": [tt[A][:2], tt[B][:2], tt[C][:2]]}
    for A, B in itertools.combinations(trip, 2):
        if len(A & B) == 2:
            att = A | B
            psi = set(att) | set(tt[A][:2]) | set(tt[B][:2])
            H = G.subgraph([v for v in G if v not in psi])
            for comp in nx.connected_components(H):
                nb = {u for v in comp for u in G[v]} & att
                if nb == att:
                    return {"condition": "component adjacent to all four attachments",
                            "triples": [sorted(A), sorted(B)], "twins": [tt[A][:2], tt[B][:2]],
                            "component": sorted(comp, key=str)}
    return None
