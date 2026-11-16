"""The forms of the component containing two of the attachments, for E10
(Proposition 5.16, the component D_1 containing c and d) and E13
(Proposition 5.24, the component D_1 containing x and y), with the
completion step of Construction 3.18 carried out in full: after the
placement of the attachments on a K3,3 subdivision and the grouping of its
interior pieces into components (construction.retained and
construction.grouped), add at most two further single-vertex components with
every admissible neighbourhood, and every inclusion-minimal set of attaching
edges from the retained components (and, in E13, of original edges between
the attachments) that meets the connectivity and attachment requirements of
the case; the whole configuration is required to be minimal
(form_requirements.completions). Forms with a one-step proper minor of the same
signature that meets the requirements are then removed, as before.

The rest (placements, groupings, signatures, validity, the no-subdivision
forms of E13, the forms of D_c, D_d and D_2, the assembly and the
classification) is that of forms_e10.py and forms_e13.py."""
import itertools

import networkx as nx

import construction as c2
import forms_e10 as f17
import forms_e13 as f19
from form_requirements import Requirement, completions


def subsets_with(base, optional):
    out = []
    for k in range(len(optional) + 1):
        for S in itertools.combinations(sorted(optional), k):
            out.append(frozenset(base) | frozenset(S))
    return out


# ---------------------------------------------------------------------- E10

def e10_r_forms():
    """D_1 on the attachments a, c, d, w: added edges ac, ad, aw, the
    subdivision using ac or ad, original pairs cw, dw as an edge or a once-
    subdivided segment. Requirements: every interior component adjacent to w
    and to c or d; c and d joined through the interior (D_1 - {a, w}
    connected); N(D_1) = w (no interior neighbour at a) or aw (one at a)."""
    roots = f17.R_ROOTS
    stats = {"retained subdivisions": 0, "groupings": 0, "completed forms": 0}
    singles_w = [frozenset(s) for s in ({"w", "c"}, {"w", "d"}, {"w", "c", "d"})]
    singles_aw = singles_w + [s | {"a"} for s in singles_w]
    req_w = Requirement(per_component=[{"w"}, {"c", "d"}], members=("c", "d"), forbid={"a"})
    req_aw = Requirement(per_component=[{"w"}, {"c", "d"}], members=("c", "d"),
                         global_adjacent=[("a", False, False)])
    raw = []
    for K, part, rec in c2.retained(roots, added=[("a", "c"), ("a", "d"), ("a", "w")],
                                    original=[("c", "w"), ("d", "w")]):
        if not (K.has_edge("a", "c") or K.has_edge("a", "d")):
            continue
        stats["retained subdivisions"] += 1
        for H, comps in c2.grouped(part, roots):
            stats["groupings"] += 1
            for req, singles, targets in ((req_w, singles_w, ("c", "d", "w")),
                                          (req_aw, singles_aw, ("a", "c", "d", "w"))):
                for R, info in completions(H, comps, roots, req, singles, 2, targets):
                    stats["completed forms"] += 1
                    raw.append((R, {"origin": "subdivision of R + {a,w} + ac + ad + aw", **info}))
    ok = [(R, i) for R, i in raw if f17.r_valid(R)]
    stats["valid"] = len(ok)
    ok = c2.dedupe(ok, roots)
    stats["valid up to isomorphism"] = len(ok)
    kept = c2.local_minimal(ok, roots, f17.r_signature, f17.r_valid, f17.r_allowed_pair)
    stats["without a one-step minor of the same signature"] = len(kept)
    return kept, stats


# ---------------------------------------------------------------------- E13

def _component_lists(H, inner_roots):
    return c2.interior_components(H, inner_roots)


def e13_subdivision_forms(part):
    """G1, G2 (partition P1) and G3 (partition P2), each with the full
    completion step."""
    out = []
    stats = {}
    if part == "P1":
        # G1: N(D_1) = uv and B_1 + uv + bc nonplanar.
        roots = ("b", "c", "u", "v")
        req = Requirement(per_component=[{"u", "v"}, {"b", "c"}], members=("b", "c"),
                          global_adjacent=[("u", True, False), ("v", True, False)])
        singles = sorted({frozenset(x) | frozenset(y) for x in subsets_with((), ("b", "c")) if x
                          for y in subsets_with((), ("u", "v")) if y}, key=sorted)
        st = {"retained subdivisions": 0, "groupings": 0, "completed forms": 0}
        for K, prt, rec in c2.retained(roots, added=[("b", "c"), ("u", "v")], required=[("b", "c")],
                                       original=[("b", "u"), ("b", "v"), ("c", "u"), ("c", "v")]):
            st["retained subdivisions"] += 1
            for H, comps in c2.grouped(prt, roots):
                st["groupings"] += 1
                for R, info in completions(H, comps, roots, req, singles, 2, roots,
                                           orig_options=[("b", "u"), ("b", "v"), ("c", "u"), ("c", "v")]):
                    st["completed forms"] += 1
                    out.append((R, {"origin": "G1 subdivision of R + uv + bc", **info}))
        stats["G1"] = st
        # G2: the hard tests; the interior with u joined into one component.
        roots = ("a", "b", "c", "u")
        req = Requirement(per_component=[{"u"}], members=("b", "c"))
        singles = [frozenset(s) for s in ({"b", "u"}, {"c", "u"}, {"b", "c", "u"})]
        st = {"retained subdivisions": 0, "groupings": 0, "completed forms": 0}
        for K, prt, rec in c2.retained(roots, added=[("b", "c"), ("a", "u")], required=[("b", "c")],
                                       original=[("b", "u"), ("c", "u")]):
            st["retained subdivisions"] += 1
            for H, comps in c2.grouped(prt, tuple(z for z in roots if z != "u")):
                if len(comps) != 1:
                    continue
                st["groupings"] += 1
                inner = _component_lists(H, roots)
                for R, info in completions(H, inner, roots, req, singles, 2, ("b", "c")):
                    st["completed forms"] += 1
                    out.append((R, {"origin": "G2 hard tests, subdivision of H + bc", **info}))
        stats["G2"] = st
    else:
        # G3: N(D_1) in {u, au} and B_1 + au + ac + ad nonplanar.
        roots = ("a", "c", "d", "u")
        req_u = Requirement(per_component=[{"u"}], members=("c", "d"), forbid={"a"})
        req_au = Requirement(per_component=[{"u"}], members=("c", "d"), global_adjacent=[("a", False, True)])
        singles = [frozenset(s) for s in ({"c", "u"}, {"d", "u"}, {"c", "d", "u"})]
        st = {"retained subdivisions": 0, "groupings": 0, "completed forms": 0}
        for K, prt, rec in c2.retained(roots, added=[("a", "u"), ("a", "c"), ("a", "d")],
                                       original=[("c", "u"), ("d", "u")]):
            if not (K.has_edge("a", "c") or K.has_edge("a", "d")):
                continue
            st["retained subdivisions"] += 1
            for H, comps in c2.grouped(prt, roots):
                st["groupings"] += 1
                for req, targets in ((req_u, ("c", "d", "u")), (req_au, ("a", "c", "d", "u"))):
                    for R, info in completions(H, comps, roots, req, singles, 2, targets):
                        st["completed forms"] += 1
                        out.append((R, {"origin": "G3 subdivision of R + au + ac + ad", **info}))
        stats["G3"] = st
    return out, stats


def e13_r_forms(part):
    (x, y), r = f19.PARTITIONS[part]
    raw = list(f19.g0(part))
    sub, sub_stats = e13_subdivision_forms(part)
    raw += sub
    raw += [(f19.mirror(H), dict(i, mirrored=True)) for H, i in sub]
    allowed = f19.r_allowed(part)
    ok = []
    for H, i in raw:
        H = H.copy()
        H.remove_nodes_from([z for z in f19.X if z in H and H.degree(z) == 0])
        if f19.component_ok(H, (x, y), allowed):
            ok.append((H, i))
    roots = (x, y, "a", "u", "v")
    norm = []
    for H, i in ok:
        H = H.copy()
        H.add_nodes_from(roots)
        norm.append((H, i))
    stats = {"generated": len(raw), "valid": len(norm), "construction": sub_stats}
    norm = c2.dedupe(norm, roots)
    stats["valid up to isomorphism"] = len(norm)
    kept = c2.local_minimal(norm, roots, lambda H: f19.r_signature(H, part),
                            lambda H: f19.component_ok(H, (x, y), allowed),
                            lambda p, q: frozenset((p, q)) in allowed)
    stats["without a one-step minor of the same signature"] = len(kept)
    return kept, stats
