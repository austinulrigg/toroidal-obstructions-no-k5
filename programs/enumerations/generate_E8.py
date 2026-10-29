"""E8: bridges of type (b) with at least two interior vertices (Lemma 4.29 and Proposition 4.30).

The names: a = s1, b = t1, c = s2, d = t2, e = s3, f = t3, with twins p, q
on {a, c, e} and r, z on {b, d, f}. B contains P1 and has attachments
a, b, c, d; P_x = cd and P_w = ef are edges. By the proof of Lemma 4.29,
B is a part of Construction 3.17 with these attachments and added edges
ac, cd, db, the retained K3,3 subdivision uses cd, and B has at least two
interior vertices (the program names the bridge U). Each part is tested
against the conditions of the proof, in this order: B + ac + cd + db is
nonplanar; B + ac + bd is planar; Psi_1 + B has no K5 minor (it is a minor
of G); B + K4[a,b,c,d] has no cut of at most two vertices; no interior
vertex of degree three is adjacent to both c and d (Lemma 3.4(a), as cd is
an edge of G); and every one-step allowed operation makes B + ac + cd + db
planar. The proof finds nine labeled bridges, five up to interchanging A1
and A2.

For each bridge that passes, the two graphs of Lemma 4.29(i),
G = (Psi_1 + B) + eb and G = (Psi_1 + B) + af, are classified
(enumeration.classify). Proposition 4.30 uses that they give only T66.

Run as: python generate_E8.py FOLDER
"""
import json
import sys
import time
from collections import Counter
from itertools import combinations
from pathlib import Path

import networkx as nx
import enumeration as el

R = ("a", "b", "c", "d")
D = [("a", "c"), ("c", "d"), ("d", "b")]
MATCH = [("a", "c"), ("b", "d")]
PSI1_EDGES = ([(v, x) for v in "pq" for x in "ace"] + [(v, x) for v in "rz" for x in "bdf"]
              + [("c", "d"), ("e", "f")])


def plus(B, edges):
    H = B.copy()
    H.add_edges_from(edges)
    return H


def classify_bridge(U):
    ok, kur = el.planarity_certificate(plus(U, D))
    if ok:
        return {"outcome": "U + ac + cd + db planar", "rotation": kur}
    ok, cert = el.planarity_certificate(plus(U, MATCH))
    if not ok:
        return {"outcome": "U + ac + bd nonplanar", "kuratowski": cert}
    # Psi1 + U is a subgraph of G in every exterior (R_x = cd, R_w = ef)
    sets = el.gr.find_k5_minor(plus(U, PSI1_EDGES), trials=3000, seed=1)
    if sets:
        return {"outcome": "Psi1 + U has a K5 minor", "branch_sets": [[str(v) for v in b] for b in sets]}
    match_rot = cert
    cut = el.small_cut(plus(U, combinations(R, 2)))
    if cut is not None:
        return {"outcome": "U + K4[abcd] has a cut of order at most two", "cut": [str(v) for v in cut]}
    for v in sorted((v for v in U if v not in R), key=str):
        if U.degree(v) == 3 and "c" in U[v] and "d" in U[v]:
            return {"outcome": "degree-three vertex with adjacent neighbors (Lemma 3.4(a))",
                    "vertex": str(v), "adjacent neighbors": ["c", "d"]}
    planar_minors = []
    for (op, args), H in el.allowed_operations(U, R):
        ok, cert = el.planarity_certificate(plus(H, D))
        if not ok:
            return {"outcome": "allowed minor keeps U + ac + cd + db nonplanar", "operation": op,
                    "arguments": args, "kuratowski": cert}
        planar_minors.append({"operation": op, "arguments": args, "rotation": cert})
    return {"outcome": "survivor", "B + C kuratowski": kur, "path completion rotation": match_rot,
            "allowed minors": planar_minors}


def psi1():
    G = nx.Graph()
    G.add_edges_from((v, s) for v in "pq" for s in "ace")
    G.add_edges_from((v, t) for v in "rz" for t in "bdf")
    G.add_edges_from([("c", "d"), ("e", "f")])
    return G


# the two graphs of Lemma 4.29(i)
EXTERIORS = {"edge eb": [("e", "b")], "edge af": [("a", "f")]}


if __name__ == "__main__":
    t0 = time.time()
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    spec = el.Spec(R, c_edges=D, required=[("c", "d")])
    parts, st = el.generate_parts(spec, min_interior=2)
    records = []
    for k, (U, recipe) in enumerate(parts):
        res = classify_bridge(U)
        rec = {"id": "E8-%04d" % k, "part": el.graph_record(U), "recipe": recipe, "outcome": res}
        if res["outcome"] == "survivor":
            rec["whole graphs"] = []
            for name, edges in EXTERIORS.items():
                G = psi1()
                G.add_edges_from(U.edges())
                G.add_edges_from(edges)
                rec["whole graphs"].append({"exterior": name, "graph": el.graph_record(G),
                                            "outcome": el.classify(G, seed=k)})
        records.append(rec)
    counts = Counter(r["outcome"]["outcome"] for r in records)
    surv = [r for r in records if r["outcome"]["outcome"] == "survivor"]
    whole = Counter(w["outcome"]["outcome"] for r in surv for w in r["whole graphs"])
    found = sorted({w["outcome"]["catalogue"] for r in surv for w in r["whole graphs"]
                    if w["outcome"]["outcome"] == "catalogue graph"})
    summary = {"item": "E8", "where": "Lemma 4.29 and Proposition 4.30",
               "specification": spec.payload(), "generation": st, "candidates": len(records),
               "bridge outcomes": dict(counts), "surviving bridges": len(surv),
               "whole-graph outcomes for the surviving bridges": dict(whole),
               "catalogue graphs found": found, "targets in the text": ["T66"],
               "survivors": [{"id": r["id"], "edges": r["part"]["edges"],
                              "exteriors": {w["exterior"]: w["outcome"]["outcome"] +
                                            (" " + w["outcome"].get("catalogue", "") if w["outcome"]["outcome"] == "catalogue graph" else "")
                                            for w in r["whole graphs"]}} for r in surv],
               "seconds": round(time.time() - t0, 1)}
    json.dump({"item": "E8", "roots": list(R),
               "test graphs": {"U + ac + cd + db": D, "U + ac + bd": MATCH,
                               "U + K4[abcd]": [list(e) for e in combinations(R, 2)],
                               "U + cd": [["c", "d"]],
                               "Psi1 + U": [list(e) for e in PSI1_EDGES]},
               "outcome tests": {
                   "U + ac + cd + db planar": ["planar", "U + ac + cd + db"],
                   "U + ac + bd nonplanar": ["nonplanar", "U + ac + bd"],
                   "Psi1 + U has a K5 minor": ["K5 minor", "Psi1 + U"],
                   "U + K4[abcd] has a cut of order at most two": ["cut", "U + K4[abcd]"],
                   "degree-three vertex with adjacent neighbors (Lemma 3.4(a))": ["degree three", "U + cd"],
                   "allowed minor keeps U + ac + cd + db nonplanar": ["allowed minor nonplanar", "U + ac + cd + db"],
                   "survivor": ["survivor", "U + ac + cd + db", "U + ac + bd"]},
               "candidates": records}, open(out / "candidates.json", "w"))
    json.dump(summary, open(out / "summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))
