"""E21: a component with five neighbours and another bridge C on abd (Proposition 6.6).

The program's names: p1, p2 are the twins on abc and q1, q2 those on ade.
As in the proof, after relabelling the only cross pair is bd, the interior
of the bridge with five attachments is a single vertex w with
N(w) = {a, b, c, d, e}, and the other bridge C lies on {a, b, d} with
nonplanar triangle completion. C is a part of Construction 3.16 with
attachments a, b, d and the added triangle. G = Psi_3 + w + C, with bd
absent or present. Every graph is classified (enumeration.classify); a
candidate with two twin triples disjoint or sharing two vertices is also
recorded with a certificate, as Section 6 assumes every two twin triples
share exactly one vertex. The proof finds no obstruction here.

Run as: python generate_E21.py FOLDER
"""
import json
import sys
from pathlib import Path

import networkx as nx

import enumeration as el
import psi2_runs as ext

PSI3 = {"p1": "abc", "p2": "abc", "q1": "ade", "q2": "ade"}


def assemble(P):
    out = []
    for bd in (False, True):
        G = nx.Graph()
        for t, nb in PSI3.items():
            G.add_edges_from((t, v) for v in nb)
        G.add_edges_from(("w", v) for v in "abcde")
        G.add_edges_from(P.edges())
        if bd:
            G.add_edge("b", "d")
        out.append(("bd" if bd else "nobd", G))
    return out


if __name__ == "__main__":
    out = Path(sys.argv[1])
    spec = el.Spec(("a", "b", "d"), c_edges=[("a", "b"), ("a", "d"), ("b", "d")])
    parts, st = el.generate_parts(spec)
    s = ext.run_whole_item("E21", [("C_on_abd", parts, st, dict(spec.payload(), W="one vertex w, N(w)=abcde",
                                                         bd="absent or present"), assemble)],
                           [], out, structural=[lambda G: ext.section6_violation(G, forbid_five=False)],
                           structural_payload=ext.SECTION6_STRUCT,
                           note="p1,p2 on abc; q1,q2 on ade; w on abcde; C on abd")
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "survivors",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["recipes"], v["parts up to rooted isomorphism"], v["outcomes"])
