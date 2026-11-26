"""E23: the second bridge C in an earlier proof of Lemma 6.14.

Here p1, p2 are the twins on {a, b, c}, and q1, q2 are the twins on {a, d, e}.
Contract the component on {a, c, d, e} to w. For each attachment triple of C,
Construction 3.18 uses the added edges shown below.

    Triple    Graph                 Completion
    abd       Psi_3 + w + C         C + ab + ad
    ade       Psi_3 + w + C + bd    C + ad + ae
    acd       Psi_3 + w + C + bd    C + ac + ad
    bcd       Psi_3 + w + C         C + bc + bd + cd

The generator classifies the graphs for all four triples. The published
data in E23 contain the abd case. The current proof uses R10.

    python generate_E23.py FOLDER"""
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
    s = el.run_item("E23", [case(*c) for c in CASES], [], out,
                    note="Psi3 + w on acde (+ bd for C on ade, acd) + part C of Construction 3.18 "
                         "with the added edges of the table of the proof")
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["parts up to rooted isomorphism"], v["outcomes"])
