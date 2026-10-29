"""E16: a cut {a, c, w} with b and d in different components, D meeting a (Proposition 5.21).

The program's names: x1, x2 are the twins on abc, y1, y2 those on abd, and
p, q the twins on {a, c, w} (p', q' in the paper). B is the bridge of the
proof, on a, b, c, d; it has the three-cut X = {a, c, w}, and its other
components are the paper's D_b, containing b, with neighbourhood {c, w}
in X, and D_d, containing d, with neighbourhood {a, w}. The paper's H_b and
H_d are the graphs on {b, c, w} and {d, a, w}.

By the proof, B + ad + bc is planar exactly when H_b + cw + bc and
H_d + aw + ad are both planar, so one of them is nonplanar. Both choices
are enumerated:

  altB: H_b + cw + bc nonplanar. Construction 3.17 with attachments
        b, c, w, added edges cw, bc, the edge bc required; the other side
        contracted to one vertex zD with N(zD) = {a, d, w}.
  altD: H_d + aw + ad nonplanar. Construction 3.17 with attachments
        a, d, w, added edges aw, ad, the edge ad required; the other side
        contracted to one vertex zB with N(zB) = {b, c, w}.

An original edge bw (altB) or dw (altD) may lie on the retained
subdivision; the construction takes both a direct edge and a segment with
one interior vertex there (psi2_runs.generate_parts_ext). Every resulting
graph G is classified (enumeration.classify), and a candidate whose B is
nonplanar is also recorded with a Kuratowski subgraph of B, as Section 5.4
assumes B planar. The proof finds no obstruction here.

Run as: python generate_E16.py FOLDER
"""
import json
import sys
from pathlib import Path

import networkx as nx

import enumeration as el
import psi2_runs as ext

TWINS = {"x1": "abc", "x2": "abc", "y1": "abd", "y2": "abd"}


def base():
    G = nx.Graph()
    for t, nb in TWINS.items():
        G.add_edges_from((t, v) for v in nb)
    for v in ("p", "q"):
        G.add_edges_from((v, x) for x in ("a", "c", "w"))
    return G


def assemble_alt_B(P):
    G = base()
    G.add_edges_from(P.edges())
    G.add_edges_from([("zD", "a"), ("zD", "d"), ("zD", "w")])
    return G


def assemble_alt_D(P):
    G = base()
    G.add_edges_from(P.edges())
    G.add_edges_from([("zB", "b"), ("zB", "c"), ("zB", "w")])
    return G


def p_vertices(G):
    return [v for v in G if v not in TWINS]


STRUCT = {"P nonplanar": {"kind": "nonplanar subgraph", "remove vertices": sorted(TWINS),
                          "remove edges within": ["a", "b", "c", "d"],
                          "meaning": "P = G minus the four twins, with no edges among a,b,c,d; "
                                     "Section 5.4 assumes P planar"}}

if __name__ == "__main__":
    out = Path(sys.argv[1])
    specB = el.Spec(("b", "c", "w"), c_edges=[("c", "w"), ("b", "c")], required=[("b", "c")])
    specD = el.Spec(("a", "d", "w"), c_edges=[("a", "w"), ("a", "d")], required=[("a", "d")])
    partsB, stB = ext.generate_parts_ext(specB, originals=[("b", "w")])
    partsD, stD = ext.generate_parts_ext(specD, originals=[("d", "w")])
    payB = dict(specB.payload(), **{"permitted original edges": [["b", "w"]],
                                    "contracted side": "zD adjacent to a, d, w"})
    payD = dict(specD.payload(), **{"permitted original edges": [["d", "w"]],
                                    "contracted side": "zB adjacent to b, c, w"})
    chk = ext.nonplanar_subgraph_check("P nonplanar", p_vertices, ["a", "b", "c", "d"])
    s = ext.run_whole_item("E16", [("altB", partsB, stB, payB, assemble_alt_B), ("altD", partsD, stD, payD, assemble_alt_D)],
        [], out, structural=[chk], structural_payload=STRUCT,
        note="x1,x2 on abc; y1,y2 on abd; p,q on acw; B on bcw, D on adw")
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes",
                                         "whole-graph outcomes before structural exclusions",
                                         "catalogue graphs among structurally excluded candidates",
                                         "survivors", "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["recipes"], v["parts up to rooted isomorphism"], v["outcomes"])
