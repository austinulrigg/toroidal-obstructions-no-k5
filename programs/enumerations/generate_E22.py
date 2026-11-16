"""E22: a component on {b, c, d, e} with one other bridge, on abd, and no cross edges (Proposition 6.13).

The program's names: p1, p2 are the twins on abc and q1, q2 those on ade.
As in the proof, the component on b, c, d, e is a single vertex w, and the
other bridge C lies on {a, b, d} with C + ab + ad nonplanar. The case
one_bridge is Construction 3.18 with attachments a, b, d and added edges
ab, ad, and G = Psi_3 + w + C.

Every graph is classified (enumeration.classify), and a candidate
violating a condition of Section 6 is also recorded with a certificate.
The proof finds no obstruction here.

Run as: python generate_E22.py FOLDER
"""
import json
import sys
from pathlib import Path

import networkx as nx

import enumeration as el
import psi2_runs as ext

PSI3 = {"p1": "abc", "p2": "abc", "q1": "ade", "q2": "ade"}


def make_assemble(extra):
    def assemble(P):
        G = nx.Graph()
        for t, nb in PSI3.items():
            G.add_edges_from((t, v) for v in nb)
        G.add_edges_from(("w", v) for v in "bcde")
        G.add_edges_from(P.edges())
        G.add_edges_from(extra)
        return G
    return assemble


if __name__ == "__main__":
    out = Path(sys.argv[1])
    spec1 = el.Spec(("a", "b", "d"), c_edges=[("a", "b"), ("a", "d")])
    parts1, st1 = el.generate_parts(spec1)
    s = ext.run_whole_item("E22", [("one_bridge", parts1, st1, dict(spec1.payload(), exterior="Psi3 + w, N(w)=bcde"), make_assemble([]))],
        [], out, structural=[lambda G: ext.section6_violation(G, forbid_five=True)],
                           structural_payload=ext.SECTION6_STRUCT,
        note="p1,p2 on abc; q1,q2 on ade; w on bcde; C on abd")
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "survivors",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["recipes"], v["parts up to rooted isomorphism"], v["outcomes"])
