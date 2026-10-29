"""E30: the remaining component with neighbourhood {a, c, d, e} (Proposition 6.15).

The program's names: p1, p2 are the twins on abc and q1, q2 those on ade.
As in the proof, G = Psi_3 + B + bd + Q with Q empty or {cd}, and G is
toroidal exactly when one of B + ac + Q + de and B + ac + Q + ad + ae is
planar; both fail. Put U = B + ac + Q. Four cases:

  U planar (with or without cd): Construction 3.17 with attachments
    a, c, d, e, added edges ac, de (+ cd), the edge de required
    (Lemmas 3.15 and 2.6);
  U nonplanar (with or without cd): Construction 3.17 with added edges
    ac (+ cd), nothing required.

Every graph is classified, first by looking for a proper nontoroidal minor
(multi_genus) and then by enumeration.classify, and each catalogue graph
found is placed in its row of Theorem 6.1. The proof finds T14, T37, T41,
T59.

Run as: python generate_E30.py FOLDER
"""
import json
import sys

import enumeration as el

import psi3_runs as ext

R = ("a", "c", "d", "e")
AC, DE, CD = ("a", "c"), ("d", "e"), ("c", "d")
CASES = [("U_planar_cd", [AC, DE, CD], [DE], [CD]),
         ("U_planar_no_cd", [AC, DE], [DE], []),
         ("U_nonplanar_cd", [AC, CD], [], [CD]),
         ("U_nonplanar_no_cd", [AC], [], [])]


def case(name, D, req, Q):
    spec = el.Spec(R, c_edges=D, required=req, name=name)

    def assemble(P):
        G = ext.psi3()
        G.add_edge("b", "d")
        G.add_edges_from(Q)
        G.add_edges_from(P.edges())
        return G
    return (name, spec, assemble, {"min_interior": 1})


if __name__ == "__main__":
    out = sys.argv[1]
    s = ext.run_item_minor_first("E30", [case(*c) for c in CASES], ["T14", "T37", "T41", "T59"], out,
                                 note="Psi3 + bd (+ cd) + part P of Construction 3.17 with R = acde; "
                                      "D = ac, de (+cd) with de required when U = P + ac + Q is planar, "
                                      "D = ac (+cd) when U is nonplanar")
    # every catalogue graph found must satisfy the conditions of this case (row 2 of Theorem 6.1)
    found = ext.annotate(out, {c[0]: (R, None) for c in CASES}, want_row=2)
    s["catalogue graphs and their rows of Theorem 6.1"] = found
    json.dump(s, open(out + "/summary.json", "w"), indent=1)
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "agrees with the text", "seconds",
                                         "catalogue graphs and their rows of Theorem 6.1")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["recipes"], v["parts up to rooted isomorphism"], v["outcomes"])
