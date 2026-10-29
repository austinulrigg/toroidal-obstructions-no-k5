"""E32: a bridge on {a, b, d} and another bridge on abc, ade, acd or cde (Lemma 6.17).

The program's names: p1, p2 are the twins on abc and q1, q2 those on ade.
As in the proof, the other bridge is contracted to a vertex v, and the
exterior is Psi_3 + ce + v with N(v) = abc, ade or acd, or Psi_3 + v with
N(v) = cde. The faces of these exteriors are recorded in R10; G is
toroidal only if B + ab + ad is planar, so B is a part of Construction 3.17
with attachments a, b, d and added edges ab, ad. Every graph is
classified, first by looking for a proper nontoroidal minor (multi_genus)
and then by enumeration.classify, and each catalogue graph found is
checked against the conditions G satisfies here, as in E31. The proof
finds a proper nontoroidal minor for every graph except two in the cde
case, both with B + ad + bd planar.

Run as: python generate_E32.py FOLDER
"""
import json
import sys
from pathlib import Path

import enumeration as el

import psi3_runs as ext

R = ("a", "b", "d")
D = [("a", "b"), ("a", "d")]
TESTS = [("U+ab+ad", [("a", "b"), ("a", "d")]), ("U+ab+bd", [("a", "b"), ("b", "d")]),
         ("U+ad+bd", [("a", "d"), ("b", "d")])]
EXTERIORS = [("v_on_abc", "abc", True), ("v_on_ade", "ade", True), ("v_on_acd", "acd", True),
             ("v_on_cde", "cde", False)]


def exterior(vt, ce):
    G = ext.psi3()
    if ce:
        G.add_edge("c", "e")
    ext.add_vertex(G, "v", vt)
    return G


def case(name, vt, ce):
    spec = el.Spec(R, c_edges=D, name=name)

    def assemble(P):
        G = exterior(vt, ce)
        G.add_edges_from(P.edges())
        return G
    return (name, spec, assemble, {"min_interior": 1})


if __name__ == "__main__":
    out = sys.argv[1]
    s = ext.run_item_minor_first("E32", [case(*c) for c in EXTERIORS], [], out,
                                 note="exterior Psi3 + ce + v (v on abc, ade, acd) or Psi3 + v (v on cde) "
                                      "+ part U of Construction 3.17 with R = abd, D = ab, ad")
    found = ext.annotate(out, {c[0]: (R, TESTS) for c in EXTERIORS}, want_row=3)
    s["catalogue graphs and their structural exclusions"] = found
    s["agrees with the text"] = s["complete"] and all(
        (x["structural exclusion"] or x["planar hypothesis tests"]) for xs in found.values() for x in xs)
    json.dump(s, open(Path(out) / "summary.json", "w"), indent=1)
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "catalogue graphs and their structural exclusions",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["recipes"], v["parts up to rooted isomorphism"], v["outcomes"])
