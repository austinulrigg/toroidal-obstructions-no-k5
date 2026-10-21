"""E28: another bridge with all three two-edge completions nonplanar (Lemma 6.13).

The program's names: p1, p2 are the twins on abc and q1, q2 those on ade.
A component W has neighbourhood {a, c, d, e}, contracted to a vertex w.
For another bridge C on each triple, the proof gives an exterior and a
completion of C that must be planar if G is toroidal:

    C on   G                    completion that must be planar
    abd    Psi_3 + w + C        C + ab + ad
    ade    Psi_3 + w + C + bd   C + ad + ae
    acd    Psi_3 + w + C + bd   C + ac + ad
    bcd    Psi_3 + w + C        C + bc + bd + cd

(C on abc is excluded in the proof by the two-cut {a, c}.) C is a part of
Construction 3.16 with the triple as attachments and the added edges of
the table. Every graph is classified (enumeration.classify); the data
folder keeps the case C on abd, which the proof cites. The proof finds no
obstruction here.

Run as: python generate_E28.py FOLDER
"""
import json
import sys

import enumeration as el

import psi3_runs as ext

CASES = [("C_on_abd", ("a", "b", "d"), [("a", "b"), ("a", "d")], False),
         ("C_on_ade", ("a", "d", "e"), [("a", "d"), ("a", "e")], True),
         ("C_on_acd", ("a", "c", "d"), [("a", "c"), ("a", "d")], True),
         ("C_on_bcd", ("b", "c", "d"), [("b", "c"), ("b", "d"), ("c", "d")], False)]


def case(name, R, D, bd):
    spec = el.Spec(R, c_edges=D, name=name)

    def assemble(P):
        G = ext.psi3()
        ext.add_vertex(G, "w", "acde")
        if bd:
            G.add_edge("b", "d")
        G.add_edges_from(P.edges())
        return G
    return (name, spec, assemble, {"min_interior": 1})


if __name__ == "__main__":
    out = sys.argv[1]
    s = el.run_item("E28", [case(*c) for c in CASES], [], out,
                    note="Psi3 + w on acde (+ bd for C on ade, acd) + part C of Construction 3.16 "
                         "with the added edges of the table of the proof")
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["parts up to rooted isomorphism"], v["outcomes"])
