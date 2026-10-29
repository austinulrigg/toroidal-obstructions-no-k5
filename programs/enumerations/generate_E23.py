"""E23: a component with neighbourhood {b, c, d, e} and no other bridges (Lemma 6.10).

The program's names: p1, p2 are the twins on abc and q1, q2 those on ade.
The bridge B of the component has attachments b, c, d, e; there are no
other bridges, and, as in the proof, no cross edge. The graph H_2 = B + bc
+ de is nonplanar, so B is a part of Construction 3.17 with attachments
b, c, d, e and added edges bc, de, and G = Psi_3 + B. Every graph is
classified: the graph and all its one-edge minors are first tested with
multi_genus (psi2_runs.fast_outcomes), the rest by enumeration.classify. A
candidate violating a condition of Section 6 (two twin triples disjoint or
sharing two vertices, or a choice of twin pairs leaving a component
adjacent to all five of their attachments) is also recorded with a
certificate. The proof finds T54.

Run as: python generate_E23.py FOLDER
"""
import json
import sys
from pathlib import Path

import networkx as nx

import enumeration as el
import psi2_runs as ext

PSI3 = {"p1": "abc", "p2": "abc", "q1": "ade", "q2": "ade"}
CROSS = [None]


def assemble(P):
    out = []
    for x in CROSS:
        G = nx.Graph()
        for t, nb in PSI3.items():
            G.add_edges_from((t, v) for v in nb)
        G.add_edges_from(P.edges())
        if x:
            G.add_edge(*x)
        out.append(("cross_" + "".join(x) if x else "nocross", G))
    return out


if __name__ == "__main__":
    out = Path(sys.argv[1])
    spec = el.Spec(("b", "c", "d", "e"), c_edges=[("b", "c"), ("d", "e")])
    parts, st = el.generate_parts(spec)
    s = ext.run_whole_item("E23", [("W_bcde", parts, st, dict(spec.payload(),
                                                      completions="no cross edge"),
                             assemble)],
                           ["T54"], out, structural=[lambda G: ext.section6_violation(G, forbid_five=True)],
                           structural_payload=ext.SECTION6_STRUCT, fast=True, slim=True,
                           note="p1,p2 on abc; q1,q2 on ade; W on bcde; no other bridges")
    d = json.load(open(out / "candidates.json"))
    from collections import Counter
    by = Counter((r["variant"], r["outcome"]["outcome"] + (" " + r["outcome"].get("catalogue", "")
                                                          if r["outcome"]["outcome"] == "catalogue graph" else ""))
                 for r in d["candidates"])
    s["outcomes by completion"] = {"%s: %s" % k: v for k, v in sorted(by.items())}
    json.dump(s, open(out / "summary.json", "w"), indent=1)
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "outcomes by completion", "survivors",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["recipes"], v["parts up to rooted isomorphism"])
