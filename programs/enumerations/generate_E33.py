"""E33: a bridge on {a, b, d} with B + ad nonplanar (Lemma 6.17).

The program's names: p1, p2 are the twins on abc and q1, q2 those on ade.
As in the proof, G = Psi_3 + ce + t + B with N(t) = {b, c, d}, or
G = Psi_3 + ce + cd + B, and by equation (5) G is toroidal only if
B + ab + ad or B + ad + bd is planar. When B + ad is nonplanar, B is a part
of Construction 3.17 with attachments a, b, d and the added edge ad. Every
graph is classified, first by looking for a proper nontoroidal minor
(multi_genus) and then by enumeration.classify, and each catalogue graph
found is checked against the conditions G satisfies here, as in E31. The
proof finds T52, T57, T58, T59 in the cd case, the first three with two
twin triples meeting twice or not at all and the last with a choice of
twin pairs leaving a component with four neighbours.

Run as: python generate_E33.py FOLDER
"""
import json
import sys
from pathlib import Path

import enumeration as el

import psi3_runs as ext

R = ("a", "b", "d")
D = [("a", "d")]
TESTS = [("U+ab+ad", [("a", "b"), ("a", "d")]), ("U+ab+bd", [("a", "b"), ("b", "d")]),
         ("U+ad+bd", [("a", "d"), ("b", "d")])]
EXTERIORS = [("t_on_bcd", True), ("edge_cd", False)]


def case(name, with_t):
    spec = el.Spec(R, c_edges=D, name=name)

    def assemble(P):
        G = ext.psi3()
        G.add_edge("c", "e")
        if with_t:
            ext.add_vertex(G, "t", "bcd")
        else:
            G.add_edge("c", "d")
        G.add_edges_from(P.edges())
        return G
    return (name, spec, assemble, {"min_interior": 1})


if __name__ == "__main__":
    out = sys.argv[1]
    s = ext.run_item_minor_first("E33", [case(*c) for c in EXTERIORS], ["T59"], out,
                                 note="Psi3 + ce + t (N(t) = bcd) or Psi3 + ce + cd, "
                                      "+ part U of Construction 3.17 with R = abd, D = ad")
    found = ext.annotate(out, {c[0]: (R, TESTS) for c in EXTERIORS}, want_row=3)
    s["catalogue graphs and their structural exclusions"] = found
    s["agrees with the text"] = s["complete"] and set(found) == {"T59"} and all(
        x["id"].startswith("E33-edge_cd") and x["row of Theorem 6.1"] == 2 for x in found["T59"])
    s["every catalogue graph is from an earlier case"] = all(
        x["structural exclusion"] in ("earlier row", "twin triples meet other than once")
        for xs in found.values() for x in xs)
    s["difference from the text"] = sorted(set(found) - {"T59"}, key=lambda x: int(x[1:]))
    json.dump(s, open(Path(out) / "summary.json", "w"), indent=1)
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "catalogue graphs and their structural exclusions",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["recipes"], v["parts up to rooted isomorphism"], v["outcomes"])
