"""E26: a component on {b, c, d, e} with other bridges on both abd and ace (Proposition 6.12).

The program's names: p1, p2 are the twins on abc and q1, q2 those on ade.
As in the proof, G = Psi_3 + w + v + C with N(w) = {b, c, d, e}, v the
contracted component with neighbourhood {a, c, e}, and C a bridge on
{a, b, d} with nonplanar triangle completion: a part of Construction 3.16
with attachments a, b, d and the added triangle. The symmetric case, C on
{a, c, e} and v on {a, b, d}, is run as well. Every graph is classified
(enumeration.classify). The proof finds no obstruction here.

Run as: python generate_E26.py FOLDER
"""
import json
import sys

import enumeration as el

import psi3_runs as ext


def case(name, C_triple, v_triple):
    spec = el.Spec(C_triple, c_edges=ext.triangle(C_triple), name=name)

    def assemble(P):
        G = ext.psi3()
        ext.add_vertex(G, "w", "bcde")
        ext.add_vertex(G, "v", v_triple)
        G.add_edges_from(P.edges())
        return G
    return (name, spec, assemble, {"min_interior": 1})


if __name__ == "__main__":
    out = sys.argv[1]
    s = el.run_item("E26", [case("C_on_abd", ("a", "b", "d"), "ace"),
                     case("C_on_ace", ("a", "c", "e"), "abd")],
                    [], out,
                    note="Psi3 (p1,p2 on abc; q1,q2 on ade) + w on bcde + v on the opposite triple "
                         "+ part C of Construction 3.16 with R = the triple and D = its triangle")
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["parts up to rooted isomorphism"], v["outcomes"])
