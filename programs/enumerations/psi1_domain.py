"""The family of bridge assemblies on Psi_1 of Proposition 4.24 (E5).

Names: a = s1, b = t1, c = s2, d = t2, e = s3, f = t3; p, q are adjacent
to S = {a, c, e} and r, z to T = {b, d, f}. The bridges are those on
S u T in G - {p, q, r, z} other than the paths (Section 4.4), so they
attach only to S u T.

Atoms and the facts of the paper that give them:
  matching edges ab, cd, ef: the path P_m is the edge s_m t_m, or it lies
      in a bridge (Proposition 4.12), and then s_m t_m is not an edge
      (equation (1)).
  cross edges ad, af, bc, be, cf, de, never both s_i t_j and s_j t_i
      (Lemma 4.7).
  three-attachment bridges on a support {x, y, z} with x, y on one side and
      z on the other (Lemma 4.15; supports inside S or T give T5 by
      Lemma 4.5, and are not in the family). A bridge with a planar
      triangle completion is a single vertex; otherwise, as it has a planar
      relevant completion (Corollary 4.23), it is one of the forms of
      Proposition 3.14: T(x), T(y), T(z), S(x;y), S(y;x), and, for the
      bridge containing a path P_m when the support is {s_m, t_m, e} (roles
      g, e, h with g the end of P_m on the side of e), also Q(h;g,e) and the
      two four-cycle graphs of Lemma 3.13 with g, e, h as x, y, z.
  four-attachment bridges, single vertices (the hypothesis of
      Proposition 4.24): on {s_i, t_i, s_j, t_j} (types (a), (b)), and on
      {s_m, t_m, s_x, t_y} with {m, x, y} = {1, 2, 3} only as the bridge
      containing P_m (type (c); outside the paths it is excluded by
      Lemma 4.6). Other four-sets give T4 or T5 (Lemmas 4.4, 4.5).

Clauses (each holds in every obstruction of Proposition 4.24):
  (a) Lemma 3.4(a): a degree-three vertex has independent neighbours, so an
      atom whose degree-three interior vertex sees both ends of an edge of
      S u T excludes that edge.
  (b) at most three bridges outside the paths on one attachment set
      (Lemma 4.16), plus the bridge containing P_m when the set contains
      s_m, t_m and s_m t_m is not an edge.
  (c) a form allowed only for the bridge containing P_m (Q, four-cycle, type
      (c)) excludes s_m t_m, and there is at most one bridge containing P_m.
  (d) Psi_1 exists: s_m t_m is an edge or some bridge attaches to both s_m
      and t_m.

form_edges() gives the forms of Proposition 3.14 and Lemma 3.13; it is also
used by forms_e11.py, forms_e14.py and generate_E13.py.
"""

from itertools import combinations

import networkx as nx
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

S, T = ("a", "c", "e"), ("b", "d", "f")
MATCH = {1: ("a", "b"), 2: ("c", "d"), 3: ("e", "f")}
PARTNER = {"a": "b", "b": "a", "c": "d", "d": "c", "e": "f", "f": "e"}
INDEX = {"a": 1, "b": 1, "c": 2, "d": 2, "e": 3, "f": 3}
CROSS = [("a", "d"), ("a", "f"), ("b", "c"), ("b", "e"), ("c", "f"), ("d", "e")]
OPPOSITE = [(("a", "d"), ("b", "c")), (("a", "f"), ("b", "e")), (("c", "f"), ("d", "e"))]


def base_graph():
    G = nx.Graph()
    G.add_edges_from((v, s) for v in "pq" for s in S)
    G.add_edges_from((v, t) for v in "rz" for t in T)
    return G


def form_edges(form, x, y, z):
    """Interior vertices start with '@'."""
    if form == "single":
        return [("@h", x), ("@h", y), ("@h", z)]
    if form.startswith("T"):          # T(w): u adjacent to w, d1, d2 to the other two and u
        w = {"Tx": x, "Ty": y, "Tz": z}[form]
        o = [v for v in (x, y, z) if v != w]
        return [(d, v) for d in ("@d1", "@d2") for v in o + ["@u"]] + [("@u", w)]
    if form in ("Sxy", "Syx"):        # S(x;y): d's on x, y, u; f on x, z, u
        p, q = (x, y) if form == "Sxy" else (y, x)
        return [(d, v) for d in ("@d1", "@d2") for v in (p, q, "@u")] + [("@f", v) for v in (p, z, "@u")]
    if form == "Q":                   # Q(z;x,y): d's on z, u, v; u on x; v on y
        return [(d, v) for d in ("@d1", "@d2") for v in (z, "@u", "@v")] + [("@u", x), ("@v", y)]
    if form in ("FC1", "FC2"):        # cycle a1 b1 a2 b2, z a1, z a2, then y b1, y b2, x a2 or x b1, x b2, y a2
        E = [("@a1", "@b1"), ("@b1", "@a2"), ("@a2", "@b2"), ("@b2", "@a1"), (z, "@a1"), (z, "@a2")]
        if form == "FC1":
            E += [(y, "@b1"), (y, "@b2"), (x, "@a2")]
        else:
            E += [(x, "@b1"), (x, "@b2"), (y, "@a2")]
        return E
    raise ValueError(form)


def degree3_conflicts(edges, attachments):
    """Pairs of attachments that are both neighbors of a degree-three
    interior vertex, in sorted order."""
    adj = {}
    for u, v in edges:
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
    out = set()
    for v, nb in adj.items():
        if v.startswith("@") and len(nb) == 3:
            for a, b in combinations(sorted(nb), 2):
                if not a.startswith("@") and not b.startswith("@"):
                    out.add(tuple(sorted((a, b))))
    return sorted(out)


def internal_violation(edges):
    adj = {}
    for u, v in edges:
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
    return any(v.startswith("@") and len(nb) == 3 and any(b in adj.get(a, ()) for a, b in combinations(nb, 2))
               for v, nb in adj.items())


def atoms(forms=("single", "Tx", "Ty", "Tz", "Sxy", "Syx", "Q", "FC1", "FC2")):
    out = []
    for m, (s, t) in MATCH.items():
        out.append({"id": "M%d" % m, "kind": "matching", "edges": [(s, t)], "attachments": [s, t], "m": m})
    for e in CROSS:
        out.append({"id": "X" + "".join(e), "kind": "cross", "edges": [e], "attachments": list(e)})
    for side, other in ((S, T), (T, S)):
        for x, y in combinations(side, 2):
            for z in other:
                support = tuple(sorted((x, y, z)))
                m = INDEX[z] if PARTNER[z] in (x, y) else None
                if m is not None:
                    g = PARTNER[z]                  # the end of P_m on the side of x, y
                    e_ = y if g == x else x
                    roles = (g, e_, z)              # g, e, h as x, y, z
                else:
                    roles = (x, y, z)
                for form in forms:
                    if form in ("Q", "FC1", "FC2") and m is None:
                        continue
                    E = form_edges(form, *roles)
                    if internal_violation(E):
                        continue
                    out.append({"id": "B%s_%s" % ("".join(support), form), "kind": "triple",
                                "support": support, "form": form, "edges": E,
                                "attachments": list(support), "path m": m,
                                "path only": form in ("Q", "FC1", "FC2"), "roles": list(roles)})
    for i, j in combinations((1, 2, 3), 2):
        sup = tuple(sorted(MATCH[i] + MATCH[j]))
        out.append({"id": "Q%s" % "".join(sup), "kind": "quad", "support": sup, "form": "single",
                    "edges": [("@h", v) for v in sup], "attachments": list(sup), "path m": None,
                    "path only": False})
    for m in (1, 2, 3):
        x, y = [k for k in (1, 2, 3) if k != m]
        for sx, ty in ((MATCH[x][0], MATCH[y][1]), (MATCH[y][0], MATCH[x][1])):
            sup = tuple(sorted(MATCH[m] + (sx, ty)))
            out.append({"id": "C%s" % "".join(sup), "kind": "quad-c", "support": sup, "form": "single",
                        "edges": [("@h", v) for v in sup], "attachments": list(sup), "path m": m,
                        "path only": True})
    return out


def maxima(atom_list):
    mx = []
    for a in atom_list:
        if a["kind"] in ("matching", "cross"):
            mx.append(1)
        elif a.get("path only"):
            mx.append(1)
        elif a["kind"] == "triple":
            mx.append(4 if a["path m"] is not None else 3)
        elif a["kind"] == "quad":
            mx.append(3)
        else:
            mx.append(1)
    return mx


def clauses(atom_list, mx, var, extra=None):
    """The clauses (a)-(d) above, plus extra(atom_list, var) if given.
    Cardinality constraints use auxiliary variables numbered after the
    threshold variables; returns (clauses, top variable)."""
    nthr = sum(mx)
    pool = IDPool(start_from=nthr + 1)
    idx = {a["id"]: i for i, a in enumerate(atom_list)}
    M = {m: var(idx["M%d" % m], 1) for m in (1, 2, 3)}
    edge_var = {tuple(sorted(MATCH[m])): M[m] for m in (1, 2, 3)}
    for e in CROSS:
        edge_var[tuple(sorted(e))] = var(idx["X" + "".join(e)], 1)
    out = []
    # (a) Lemma 3.4(a)
    for i, a in enumerate(atom_list):
        if a["kind"] in ("matching", "cross"):
            continue
        for pair in degree3_conflicts(a["edges"], a["attachments"]):
            if pair in edge_var:
                out.append([-var(i, 1), -edge_var[pair]])
    # Lemma 4.7: not both s_i t_j and s_j t_i
    for e1, e2 in OPPOSITE:
        out.append([-edge_var[tuple(sorted(e1))], -edge_var[tuple(sorted(e2))]])
    # (b) bridges on one triple support
    supports = {}
    for i, a in enumerate(atom_list):
        if a["kind"] == "triple":
            supports.setdefault(a["support"], []).append(i)
    for sup, ids in supports.items():
        lits = [var(i, k) for i in ids for k in range(1, mx[i] + 1)]
        m = atom_list[ids[0]]["path m"]
        for bound, guard in ((4, None), (3, M[m] if m is not None else None)):
            if len(lits) <= bound:
                continue
            enc = CardEnc.atmost(lits=lits, bound=bound, vpool=pool, encoding=EncType.seqcounter)
            for cl in enc.clauses:
                out.append(cl + ([-guard] if guard is not None else []))
    # (c) forms only for the bridge containing P_m
    for m in (1, 2, 3):
        only = [i for i, a in enumerate(atom_list) if a.get("path only") and a["path m"] == m]
        for i in only:
            out.append([-var(i, 1), -M[m]])
        for i, j in combinations(only, 2):
            out.append([-var(i, 1), -var(j, 1)])
    # (d) Psi1 exists
    for m in (1, 2, 3):
        s, t = MATCH[m]
        out.append([M[m]] + [var(i, 1) for i, a in enumerate(atom_list)
                             if a["kind"] not in ("matching", "cross") and s in a["attachments"]
                             and t in a["attachments"]])
    if extra:
        out += extra(atom_list, var, idx)
    return out, pool.top


def label_symmetries():
    """The twelve relabelings of Psi1 preserving the domain: a permutation of
    the indices 1, 2, 3 applied to S and T together, possibly followed by
    interchanging S and T (with p, q and r, z)."""
    from itertools import permutations
    out = []
    for perm in permutations(range(3)):
        for swap in (False, True):
            m = {}
            for i in range(3):
                j = perm[i]
                if not swap:
                    m[S[i]], m[T[i]] = S[j], T[j]
                else:
                    m[S[i]], m[T[i]] = T[j], S[j]
            m.update({"p": "p", "q": "q", "r": "r", "z": "z"} if not swap else {"p": "r", "q": "z", "r": "p", "z": "q"})
            out.append(m)
    return out


def atom_permutations(atom_list, maps, maxima):
    """For each label map, the permutation of atom indices it induces (the
    image of every atom's edge set, interior names kept, is an atom's edge
    set), checked to preserve the maxima."""
    key = lambda E: frozenset(frozenset(e) for e in E)
    idx = {key(a["edges"]): i for i, a in enumerate(atom_list)}
    perms = []
    for m in maps:
        f = lambda v: v if v.startswith("@") else m[v]
        perm = [idx[key([(f(x), f(y)) for x, y in a["edges"]])] for a in atom_list]
        assert sorted(perm) == list(range(len(atom_list)))
        assert all(maxima[i] == maxima[perm[i]] for i in range(len(atom_list)))
        perms.append(perm)
    return perms
