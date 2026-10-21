"""E9: a bridge with four attachments and nonplanar B or B + ab + ac + bc (Proposition 5.9).

The names: Psi_2 is a, b, c, d with n_c twins on abc and n_d twins on abd
(the program calls the twins p0, p1, ... and r0, r1, ...), and B is the
bridge with attachments a, b, c, d. Two cases of the proof:

  A: B nonplanar. Construction 3.16 with attachments a, b, c, d and no
     added edges. Each part is completed with (n_c, n_d) = (2,2), (3,2)
     and (3,3); Lemma 5.6 excludes (3,3), which is included anyway.
  B: (n_c, n_d) = (3,2) and B + ab + ac + bc nonplanar. Construction 3.16
     with added edges ab, ac, bc, completed with (3,2).

Every resulting graph is classified (enumeration.classify). The proof finds
T17, T18, T19, T22, T24, T26, T30, T34, T39, T48, T52, T57, T66, T74 in
case A and T11 in case B.

Run as: python generate_E9.py FOLDER
"""
import json
import sys

import networkx as nx
import enumeration as el

R = ("a", "b", "c", "d")


def psi2(k, l):
    G = nx.Graph()
    for i in range(k):
        G.add_edges_from(("p%d" % i, x) for x in "abc")
    for i in range(l):
        G.add_edges_from(("r%d" % i, x) for x in "abd")
    return G


def case(name, c_edges, kls):
    spec = el.Spec(R, c_edges=c_edges)

    def assemble(P):
        out = []
        for k, l in kls:
            G = psi2(k, l)
            G.add_edges_from(P.edges())
            out.append(G)
        return out
    return (name, spec, assemble, {"min_interior": 1})


if __name__ == "__main__":
    out = sys.argv[1]
    s = el.run_item("E9", [case("A", [], [(2, 2), (3, 2), (3, 3)]),
                     case("B", [("a", "b"), ("a", "c"), ("b", "c")], [(3, 2)])],
                    ["T11", "T17", "T18", "T19", "T22", "T24", "T26", "T30", "T34", "T39", "T48", "T52",
                     "T57", "T66", "T74"], out,
                    note="assembled graphs: case A with (k,l) = (2,2), (3,2), (3,3) in that order; case B with (3,2)")
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found", "catalogue graphs by case",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["parts up to rooted isomorphism"], v["outcomes"])
