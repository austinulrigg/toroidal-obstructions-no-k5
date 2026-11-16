"""E13: the cut X = {a, u, v} with two other components (Proposition 5.24).

The program's names: t0, t1 are the twins on abc and t2, t3 those on abd
(p, q, r, z in the paper), and p, q have neighbourhood X = {a, u, v} (p', q'
in the paper). G = Psi_2 + p, q + D_1 + D_2, where D_1 contains two of b, c,
d and D_2 the third (the program writes R for D_1 and D for D_2). Up to
interchanging c and d the partitions are P1 (D_1 contains b, c) and P2 (D_1
contains c, d). The forms of D_1 come from forms_one_component.e13_r_forms
(the forms with no K3,3 subdivision, and the constructions of the proof with
the full completion step, with their mirror images under u <-> v), and those
of D_2 from forms_e13.d_forms; forms with a one-step proper minor of the
same signature are removed first. Every pair whose neighbourhoods together
contain u and v is assembled; the graphs are compared up to isomorphism,
each class is classified (enumeration.classify), and each is placed in the
rows of the table of Lemma 5.19:
  row 1  P1, N(D_1) in {u, au}, N(D_2) = uv
  row 2  P1, N(D_1) in {u, au}, N(D_2) in {v, av}   (with its condition)
  row 3  P1, N(D_1) = uv, N(D_2) = uv
  row 4  P1, N(D_1) = uv, N(D_2) in {v, av}
  row 5  P2, one neighbourhood uv
  row 6  P2, N(D_1) in {u, au}, N(D_2) in {v, av}
each up to u <-> v. The proof finds no obstruction here (1,229 candidates).

The forms are kept in forms.pkl in the output folder, and the
classification in progress.jsonl, so that a second start continues.
Run as: python generate_E13.py FOLDER
"""
import json
import pickle
import sys
from collections import Counter
from pathlib import Path

import networkx as nx

import enumeration as el
import forms_e13 as f
import forms_one_component as fx


def _lookup(part, NR, ND):
    if part == "P1":
        if NR in ("u", "au") and ND == "uv":
            return "1"
        if NR in ("u", "au") and ND in ("v", "av"):
            return "2"
        if NR == "uv" and ND == "uv":
            return "3"
        if NR == "uv" and ND in ("v", "av"):
            return "4"
    else:
        if "uv" in (NR, ND):
            return "5"
        if NR in ("u", "au") and ND in ("v", "av"):
            return "6"
    return None


def row(part, NR, ND):
    sw = lambda N: "".join(sorted({"u": "v", "v": "u", "a": "a"}[z] for z in N))
    return _lookup(part, NR, ND) or _lookup(part, sw(NR), sw(ND)) or "outside the table"


def prefixed(H, keep, tag):
    return nx.relabel_nodes(H, {v: v if v in keep else "%s%s" % (tag, v) for v in H})


def base():
    G = nx.Graph()
    G.add_edges_from((t, x) for t in ("t0", "t1") for x in "abc")
    G.add_edges_from((t, x) for t in ("t2", "t3") for x in "abd")
    G.add_edges_from((p, x) for p in "pq" for x in "auv")
    return G


def build(forms):
    graphs = []
    named = set("abcduv")
    for part, (R, rst, D, dst) in forms.items():
        for i, (RH, rinfo) in enumerate(R):
            NR = f.nbhd(RH, f.PARTITIONS[part][0])
            for j, (DH, dinfo) in enumerate(D):
                ND = f.nbhd(DH, (f.PARTITIONS[part][1],))
                if not {"u", "v"} <= set(NR) | set(ND):
                    continue
                G = base()
                G.add_edges_from(prefixed(RH, named, "R").edges())
                G.add_edges_from(prefixed(DH, named, "D").edges())
                graphs.append((G, {"partition": part, "R form": i, "D form": j, "N(R)": NR, "N(D)": ND,
                                   "row": row(part, NR, ND), "R origin": rinfo["origin"],
                                   "D origin": dinfo["origin"]}))
    return graphs


def dedupe(graphs):
    buckets, out = {}, []
    for G, info in graphs:
        key = (G.number_of_nodes(), G.number_of_edges(), nx.weisfeiler_lehman_graph_hash(G, iterations=4))
        hit = None
        for j in buckets.get(key, []):
            if nx.is_isomorphic(G, out[j][0]):
                hit = j
                break
        if hit is None:
            buckets.setdefault(key, []).append(len(out))
            out.append((G, {"first": info, "also": []}))
        else:
            out[hit][1]["also"].append(info)
    return out


if __name__ == "__main__":
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    cache = out / "forms.pkl"
    if cache.exists():
        with open(cache, "rb") as fh:
            forms = pickle.load(fh)
    else:
        forms = {}
        for part in ("P1", "P2"):
            R, rst = fx.e13_r_forms(part)
            D, dst = f.d_forms(part)
            forms[part] = (R, rst, D, dst)
            print(part, "D_1 forms", len(R), "D_2 forms", len(D), flush=True)
        with open(cache, "wb") as fh:
            pickle.dump(forms, fh)
    record = {}
    for part, (R, rst, D, dst) in forms.items():
        record[part] = {
            "R forms": [{"index": i, "signature": list(f.r_signature(H, part)), "origin": inf["origin"],
                         "mirrored": inf.get("mirrored", False), "graph": el.graph_record(H)}
                        for i, (H, inf) in enumerate(R)],
            "D forms": [{"index": i, "signature": list(f.d_signature(H, part)), "origin": inf["origin"],
                         "graph": el.graph_record(H)} for i, (H, inf) in enumerate(D)],
            "R statistics": rst, "D statistics": dst}
    with open(out / "forms.json", "w") as fh:
        json.dump(record, fh, indent=1)
    graphs = build(forms)
    print("assembled", len(graphs), Counter(i["row"] for _, i in graphs), flush=True)
    uniq = dedupe(graphs)
    print("up to isomorphism", len(uniq), flush=True)
    s = el.run_item("E13", [("assemblies", None, lambda G: G,
                             {"parts": uniq, "stats": {"assembled": len(graphs), "up to isomorphism": len(uniq)},
                              "payload": "see forms.json"})],
                    [], str(out), note="G = Psi2 + p, q on auv + D_1 + D_2")
    data = json.load(open(out / "candidates.json"))
    rows, found = Counter(), {}
    for r in data["candidates"]:
        rs = {r["recipe"]["first"]["row"]} | {x["row"] for x in r["recipe"]["also"]}
        for x in rs:
            rows[x] += 1
            if r["outcome"]["outcome"] == "catalogue graph":
                found.setdefault(x, set()).add(r["outcome"]["catalogue"])
    s["candidates by table row"] = dict(rows)
    s["catalogue graphs by table row"] = {k: sorted(v, key=lambda x: int(x[1:])) for k, v in found.items()}
    with open(out / "summary.json", "w") as fh:
        json.dump(s, fh, indent=1)
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "candidates by table row", "catalogue graphs by table row", "seconds")},
                     indent=1))
