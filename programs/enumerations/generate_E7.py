"""E7: bridges of type (a) with at least two interior vertices (Proposition 4.28).

B has attachments sx, tx, sy, ty (s_x, t_x, s_y, t_y in the paper), and C is
the cycle sx tx ty sy sx. By the proof, B is a part of Construction 3.17
with these attachments and added edges E(C), and the retained K3,3
subdivision uses sy ty. Each part is tested against the conditions of the
proof, in this order: B + C is nonplanar; the path completion
B + sy sx + sx tx + tx ty is planar; B + K4[sx,tx,sy,ty] has no cut of at
most two vertices; no interior vertex of degree three has two adjacent
neighbours in B + sx tx + sy ty (Lemma 3.4(a), as P_x and P_y are edges);
B has no disjoint sx-ty and sy-tx paths; and every one-step allowed
operation of Section 3 makes B + C planar. A part failing a condition is
recorded with a certificate. The proof finds one bridge, up to
interchanging A1 and A2: twins p', q' on {sx, sy, v}, with v also adjacent
to tx and ty.

Run as: python generate_E7.py FOLDER
"""
import json
import sys
from itertools import combinations
from pathlib import Path

import enumeration as el

R = ("sx", "tx", "sy", "ty")
CYCLE = [("sx", "tx"), ("tx", "ty"), ("ty", "sy"), ("sy", "sx")]
PATH = [("sy", "sx"), ("sx", "tx"), ("tx", "ty")]


def plus(B, edges):
    H = B.copy()
    H.add_edges_from(edges)
    return H


def classify_bridge(B):
    ok, cert = el.planarity_certificate(plus(B, CYCLE))
    if ok:
        return {"outcome": "B + C planar", "rotation": cert}
    kur = cert
    ok, cert = el.planarity_certificate(plus(B, PATH))
    if not ok:
        return {"outcome": "path completion nonplanar", "kuratowski": cert}
    path_rot = cert
    K4 = plus(B, combinations(R, 2))
    cut = el.small_cut(K4)
    if cut is not None:
        return {"outcome": "B + K4[R] has a cut of order at most two", "cut": [str(v) for v in cut]}
    # Lemma 3.4(a): a vertex of degree three has independent neighbors. In G
    # the paths P_x = sx tx and P_y = sy ty are edges.
    Gpart = plus(B, [("sx", "tx"), ("sy", "ty")])
    for v in sorted((v for v in B if v not in R), key=str):
        if B.degree(v) == 3:
            nb = sorted(B[v], key=str)
            bad = [[str(a), str(b)] for a, b in combinations(nb, 2) if Gpart.has_edge(a, b)]
            if bad:
                return {"outcome": "degree-three vertex with adjacent neighbors (Lemma 3.4(a))",
                        "vertex": str(v), "adjacent neighbors": bad[0]}
    dp = el.disjoint_paths(B, R, [("sx", "ty"), ("sy", "tx")])
    if dp:
        return {"outcome": "disjoint sx-ty and sy-tx paths", "paths": [[str(v) for v in p] for p in dp]}
    planar_minors = []
    for (op, args), H in el.allowed_operations(B, R):
        ok, cert = el.planarity_certificate(plus(H, CYCLE))
        if not ok:
            return {"outcome": "allowed minor keeps B + C nonplanar", "operation": op, "arguments": args,
                    "kuratowski": cert}
        planar_minors.append({"operation": op, "arguments": args, "rotation": cert})
    return {"outcome": "survivor", "B + C kuratowski": kur, "path completion rotation": path_rot,
            "allowed minors": planar_minors}


if __name__ == "__main__":
    import time
    from collections import Counter
    t0 = time.time()
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    spec = el.Spec(R, c_edges=CYCLE, required=[("sy", "ty")])
    parts, st = el.generate_parts(spec, min_interior=2)
    records = []
    for k, (B, recipe) in enumerate(parts):
        res = classify_bridge(B)
        rec = {"id": "E7-%04d" % k, "part": el.graph_record(B), "recipe": recipe, "outcome": res}
        if res["outcome"] == "survivor":
            # the edges whose deletion keeps B + C nonplanar
            excl = []
            for e in sorted(B.edges(), key=el.edge_key):
                H = B.copy(); H.remove_edge(*e)
                ok, cert = el.planarity_certificate(plus(H, CYCLE))
                if not ok:
                    excl.append({"deleted edge": list(el.edge_key(e)), "kuratowski": cert})
            rec["edge deletions keeping B + C nonplanar"] = excl
        records.append(rec)
    counts = Counter(r["outcome"]["outcome"] for r in records)
    surv = [r for r in records if r["outcome"]["outcome"] == "survivor"]
    summary = {"item": "E7", "where": "Proposition 4.28",
               "specification": spec.payload(), "generation": st, "candidates": len(records),
               "outcomes": dict(counts),
               "survivors": [{"id": r["id"], "edges": r["part"]["edges"],
                              "deletions keeping B + C nonplanar": [x["deleted edge"] for x in
                                                                    r["edge deletions keeping B + C nonplanar"]]}
                             for r in surv],
               "seconds": round(time.time() - t0, 1)}
    json.dump({"item": "E7", "roots": list(R),
               "test graphs": {"B + C": CYCLE, "B + path": PATH,
                               "B + K4[R]": [list(e) for e in combinations(R, 2)],
                               "B + Rx + Ry": [["sx", "tx"], ["sy", "ty"]]},
               "outcome tests": {
                   "B + C planar": ["planar", "B + C"],
                   "path completion nonplanar": ["nonplanar", "B + path"],
                   "B + K4[R] has a cut of order at most two": ["cut", "B + K4[R]"],
                   "degree-three vertex with adjacent neighbors (Lemma 3.4(a))": ["degree three", "B + Rx + Ry"],
                   "disjoint sx-ty and sy-tx paths": ["disjoint paths", [["sx", "ty"], ["sy", "tx"]]],
                   "allowed minor keeps B + C nonplanar": ["allowed minor nonplanar", "B + C"],
                   "survivor": ["survivor", "B + C", "B + path"]},
               "candidates": records}, open(out / "candidates.json", "w"))
    json.dump(summary, open(out / "summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))
