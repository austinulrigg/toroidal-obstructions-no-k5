"""E24: three vertices with neighbourhood {a, d, e} (Lemma 6.15).

The program's names: p1, p2 are the twins on abc and q1, q2, q3 the
vertices on ade. A component W has neighbourhood {a, c, d, e} with
bridge B. As in the proof, the rest of G is two vertices on abc, three on
ade, the edge bd and possibly cd, and G is toroidal exactly when B + E is
planar for the single edge set E = {ac, ad, de, ae} without cd, or
E = {ac, cd, de, ae} with cd. B is a part of Construction 3.18 with
attachments a, c, d, e and these added edges, and G = Psi_3 + q3 + bd
(+ cd) + B. Every graph is classified, first by looking for a proper
nontoroidal minor (multi_genus) and then by enumeration.classify
(psi3_runs.run_item_minor_first). The proof finds no obstruction here.

Run as: python generate_E24.py FOLDER
"""
import json
import sys

import enumeration as el

import psi3_runs as ext

R = ("a", "c", "d", "e")
CASES = [("no_cd", [("a", "c"), ("a", "d"), ("d", "e"), ("a", "e")], [("b", "d")]),
         ("with_cd", [("a", "c"), ("c", "d"), ("d", "e"), ("a", "e")], [("b", "d"), ("c", "d")])]


def case(name, D, cross):
    spec = el.Spec(R, c_edges=D, name=name)

    def assemble(P):
        G = ext.psi3()
        ext.add_vertex(G, "q3", "ade")
        G.add_edges_from(cross)
        G.add_edges_from(P.edges())
        return G
    return (name, spec, assemble, {"min_interior": 1})


if __name__ == "__main__":
    out = sys.argv[1]
    s = ext.run_item_minor_first("E24", [case(*c) for c in CASES], [], out,
                                 note="Psi3 + q3 on ade + bd (+ cd) + part P of Construction 3.18 with "
                                      "R = acde and the single test of the table's last column")
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["recipes"], v["parts up to rooted isomorphism"], v["outcomes"])
