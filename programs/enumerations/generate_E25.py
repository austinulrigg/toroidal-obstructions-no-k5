"""E25: the cross edge bd with one other bridge C on abd (Proposition 6.12).

The program's names: p1, p2 are the twins on abc, q1, q2 those on ade, and
h1, h2 the twins with neighbourhood {c, e, w} (p', q' in the paper). As in
the proof, the bridge on b, c, d, e is reduced to w, h1, h2 with
N(w) = {b, d, h1, h2}, N(h1) = N(h2) = {c, e, w}, and the other bridge C
lies on {a, b, d} with nonplanar triangle completion: a part of
Construction 3.17 with attachments a, b, d and the added triangle.
G = Psi_3 + bd + w + h1 + h2 + C. Every graph is classified
(enumeration.classify). A candidate is also recorded with a certificate
when C + ab + ad is nonplanar (it must be planar when a cross edge is
present) or when it violates a condition of Section 6. The proof finds no
obstruction here.

Run as: python generate_E25.py FOLDER
"""
import json
import sys
from pathlib import Path

import networkx as nx

import enumeration as el
import psi2_runs as ext

PSI3 = {"p1": "abc", "p2": "abc", "q1": "ade", "q2": "ade"}
EXTERIOR = ["p1", "p2", "q1", "q2", "c", "e", "w", "h1", "h2"]


def assemble(P):
    G = nx.Graph()
    for t, nb in PSI3.items():
        G.add_edges_from((t, v) for v in nb)
    G.add_edge("b", "d")
    G.add_edges_from([("w", "b"), ("w", "d"), ("w", "h1"), ("w", "h2")])
    for h in ("h1", "h2"):
        G.add_edges_from((h, x) for x in "ce")
    G.add_edges_from(P.edges())
    return G


def b_check(G):
    H = G.subgraph([v for v in G if v not in EXTERIOR]).copy()
    H.remove_edges_from([("a", "b"), ("a", "d"), ("b", "d")])
    H.add_edges_from([("a", "b"), ("a", "d")])
    H.remove_edges_from([("b", "d")])
    ok, cert = el.planarity_certificate(H)
    if ok:
        return None
    return {"condition": "C + ab + ad nonplanar", "kind": "nonplanar subgraph", "kuratowski": cert}


STRUCT = {"C + ab + ad nonplanar": {"kind": "nonplanar subgraph", "remove vertices": EXTERIOR,
                                    "remove edges within": ["a", "b", "d"],
                                    "add edges": [["a", "b"], ["a", "d"]],
                                    "meaning": "B = C + ab + ad must be planar when a cross edge is present"}}

if __name__ == "__main__":
    out = Path(sys.argv[1])
    spec = el.Spec(("a", "b", "d"), c_edges=[("a", "b"), ("a", "d"), ("b", "d")])
    parts, st = el.generate_parts(spec)
    s = ext.run_whole_item("E25", [("C_on_abd", parts, st, dict(spec.payload(),
                                                         exterior="Psi3 + bd + w + h1 + h2, N(w)=b,d,h1,h2, "
                                                                  "N(h1)=N(h2)=c,e,w"), assemble)],
                           [], out, structural=[b_check, lambda G: ext.section6_violation(G, forbid_five=True)],
                           structural_payload=dict(STRUCT, **ext.SECTION6_STRUCT),
                           note="p1,p2 on abc; q1,q2 on ade; h1,h2 on cew; C on abd")
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "survivors",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["recipes"], v["parts up to rooted isomorphism"], v["outcomes"])
