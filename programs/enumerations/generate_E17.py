"""E17: a cut {a, c, w} with b and d in one component (Lemma 5.20).

The program's names: x1, x2 are the twins on abc, y1, y2 those on abd, and
p, q the twins on {a, c, w} (p', q' in the paper). Here b and d lie in one
component of B - X, X = {a, c, w}. Deleting a and the six vertices of
degree three on abc, abd, acw leaves the paper's H_1 on the attachments
b, c, d, w, with interior I, subject to equation (2) of the paper:
H_1[{b, d} u I] connected, H_1[{w} u I] connected, and c has a neighbour
in I. Only bw and dw can join two attachments. G is nontoroidal exactly
when H_1 + bd + cw is nonplanar, and H_1 + cw is planar.

Construction 3.17 with attachments b, c, d, w and added edges bd, cw; the
retained K3,3 subdivision must use bd, and original edges bw, dw are
permitted on it (both as a direct edge and as a segment with one interior
vertex). As the proof states, by Lemma 5.8:

  * the pieces of the subdivision minus {b, c, d, w} are split into
    interior components in every way (every set partition), and the
    pieces of each component are joined by every tree, with every choice
    of ends;
  * at most two further components are added, each a single vertex s;
    by equation (2), N(s) contains w and one of b, d, and every such set
    is allowed: wb, wd, wbd, wbc, wdc, wbcd;
  * then every inclusion-minimal set of attaching edges that, with the
    further components, satisfies equation (2) is added, and the whole
    configuration is required to be minimal.

Every resulting graph G is classified (enumeration.classify), and a
candidate whose B is nonplanar is also recorded with a Kuratowski subgraph
of B, as Section 5.4 assumes B planar. The proof finds no obstruction here.

Run as: python generate_E17.py FOLDER
"""
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import networkx as nx

import enumeration as el
import psi2_runs as ext

ROOTS = ("b", "c", "d", "w")
TWINS = {"x1": "abc", "x2": "abc", "y1": "abd", "y2": "abd"}
SINGLE_OPTIONS = [frozenset(s) for s in ("wb", "wd", "wbd", "wbc", "wdc", "wbcd")]


def requirement(atts):
    """(r2:cross-connectivity) at the level of interior components, given
    the set of attachments each component is adjacent to."""
    if not atts or not all("w" in a for a in atts) or not any("c" in a for a in atts):
        return False
    return bd_connected(atts)


def bd_connected(atts):
    H = nx.Graph()
    H.add_nodes_from(["b", "d"])
    for i, a in enumerate(atts):
        H.add_node(i)
        for r in a & {"b", "d"}:
            H.add_edge(i, r)
    return nx.is_connected(H)


def additions(groups_att, singles):
    n = len(groups_att)
    allatts = list(groups_att) + list(singles)
    wadd = [(g, "w") for g in range(n) if "w" not in groups_att[g]]
    if any("c" in a for a in allatts):
        cadds = [[]]
    else:
        cadds = [[(g, "c")] for g in range(n)]
    bd_cands = [(g, r) for g in range(n) for r in "bd"]
    bdsets = ext.minimal_group_additions(groups_att, singles, bd_connected, bd_cands)
    for c in cadds:
        for bd in bdsets:
            A = wadd + c + bd
            if ext.globally_minimal(groups_att, A, singles, requirement):
                yield A


def generate():
    spec = el.Spec(ROOTS, c_edges=[("b", "d"), ("c", "w")], required=[("b", "d")])
    reps, raw = el.placement_classes(4)
    stats = Counter()
    stats["root placements"] = raw
    stats["root placement classes"] = len(reps)
    buckets, parts = {}, []
    for pl0, orbit in reps:
        for K, part, used_orig in ext.build_retained_variants(pl0, spec, [("b", "w"), ("d", "w")]):
            stats["valid retained subdivisions"] += 1
            pieces = el.interior_pieces(part, ROOTS)
            for partition in ext.set_partitions(range(len(pieces))):
                groups = [[pieces[i] for i in grp] for grp in partition]
                for forest in ext.group_forests(part, ROOTS, groups):
                    base = part.copy()
                    base.add_edges_from(forest)
                    comps = ext.components_attach(base, ROOTS)
                    gatt = [a for _, a in comps]
                    gverts = [c for c, _ in comps]
                    for k in range(3):
                        for singles in itertools.combinations_with_replacement(SINGLE_OPTIONS, k):
                            for A in additions(gatt, singles):
                                for verts in itertools.product(*[gverts[g] for g, r in A]):
                                    stats["recipes"] += 1
                                    P = base.copy()
                                    P.add_edges_from((v, r) for v, (g, r) in zip(verts, A))
                                    for j, s in enumerate(singles):
                                        P.add_edges_from(("s%d" % j, r) for r in s)
                                    key = el.rooted_key(P, ROOTS)
                                    if any(el.rooted_isomorphic(P, parts[i][0], ROOTS)
                                           for i in buckets.get(key, [])):
                                        continue
                                    buckets.setdefault(key, []).append(len(parts))
                                    parts.append((el.relabel_interior(P, ROOTS), {
                                        "placement": [list(pl0[0]), [list(p) for p in pl0[1]]],
                                        "K'": sorted(sorted(map(str, e)) for e in K.edges()),
                                        "original edges on K'": sorted(sorted(e) for e in used_orig),
                                        "components": [[list(map(str, p)) for p in g] for g in groups],
                                        "further components": [sorted(s) for s in singles],
                                        "added attaching edges": [[str(v), r] for v, (g, r) in zip(verts, A)]}))
    stats["parts up to rooted isomorphism"] = len(parts)
    return parts, dict(stats), spec


def assemble(P):
    G = nx.Graph()
    for t, nb in TWINS.items():
        G.add_edges_from((t, v) for v in nb)
    for v in ("p", "q"):
        G.add_edges_from((v, x) for x in ("a", "c", "w"))
    G.add_edges_from(P.edges())
    return G


STRUCT = {"P nonplanar": {"kind": "nonplanar subgraph", "remove vertices": sorted(TWINS),
                          "remove edges within": ["a", "b", "c", "d"],
                          "meaning": "P = G minus the four twins, with no edges among a,b,c,d; "
                                     "Section 5.4 assumes P planar"}}

if __name__ == "__main__":
    out = Path(sys.argv[1])
    t0 = time.time()
    parts, st, spec = generate()
    print("generated", st, round(time.time() - t0, 1), flush=True)
    pay = dict(spec.payload(), **{
        "permitted original edges": [["b", "w"], ["d", "w"]],
        "further components": "at most two single vertices, N(s) one of " +
                              ", ".join("".join(sorted(s)) for s in SINGLE_OPTIONS),
        "requirements": "K[{b,d} u I] connected, K[{w} u I] connected, c has a neighbor in I",
        "interior components": "every partition of the pieces of K' - R, each joined by every tree"})
    chk = ext.nonplanar_subgraph_check("P nonplanar", lambda G: [v for v in G if v not in TWINS],
                                       ["a", "b", "c", "d"])
    s = ext.run_whole_item("E17", [("R_bd", parts, st, pay, assemble)], [], out, structural=[chk],
                           structural_payload=STRUCT,
                           note="x1,x2 on abc; y1,y2 on abd; p,q on acw; K on b,c,d,w")
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes",
                                         "whole-graph outcomes before structural exclusions",
                                         "catalogue graphs among structurally excluded candidates",
                                         "survivors", "agrees with the text", "seconds")}, indent=1))
