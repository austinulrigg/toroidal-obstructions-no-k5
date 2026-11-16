"""E10: the cut X = {a, b, w} (Proposition 5.16).

The program's names: t0, t1 are the twins on abc and t2, t3 those on abd
(p, q, r, z in the paper). As in the proof, G is Psi_2 together with
  (a) one component D_0 whose bridge B_0 is one of the eleven forms of the
      proof, and the components D_c, D_d, with every subset of the edges
      aw, bw;
  (b) two single vertices k1, k2 on abw, and D_c, D_d;
  (c) two single vertices k1, k2 on abw, and the component D_1 containing
      both c and d, with neighbourhood w or aw in X (bw is the mirror
      image).
The forms of D_c (and of D_d, with d for c) come from forms_e10.c_forms and
those of D_1 from forms_one_component.e10_r_forms; forms with a one-step
proper minor of the same signature are removed first. Every assembled graph
is compared with the others up to isomorphism, and each class is classified
(enumeration.classify). Each candidate also records whether B (G minus the
twins, with no edges among a, b, c, d) is planar, as Section 5.4 retains
the graphs with B planar. The proof finds 33 forms of D_c and 54 of D_1,
and T56, T72, T73 among the graphs with B planar.

The forms are kept in forms.pkl in the output folder, and the
classification in progress.jsonl, so that a second start continues.
Run as: python generate_E10.py FOLDER
"""
import json
import pickle
import sys
from collections import Counter
from pathlib import Path

import networkx as nx

import enumeration as el
import forms_e10 as f
import forms_one_component as fx

ELEVEN = ['one vertex', 'T(a)', 'T(b)', 'T(w)', 'S(a;b)', 'S(a;w)', 'S(b;a)', 'S(b;w)',
          'Q(a;b,w)', 'Q(b;a,w)', 'Q(w;a,b)']
TWINS = ("t0", "t1", "t2", "t3")


def prefixed(H, roots, tag, ren=None):
    ren = ren or {}
    m = {v: ren.get(v, v) if v in roots else "%s%s" % (tag, v) for v in H}
    return nx.relabel_nodes(H, m)


def psi2():
    G = nx.Graph()
    G.add_edges_from((t, x) for t in ("t0", "t1") for x in "abc")
    G.add_edges_from((t, x) for t in ("t2", "t3") for x in "abd")
    return G


def build(cforms, rforms):
    U = [(n, H) for n, H in f.u_forms() if n in ELEVEN]
    assert [n for n, _ in U] == ELEVEN
    C = [(i, prefixed(H, f.ROOTS, "C")) for i, (H, _) in enumerate(cforms)]
    D = [(i, prefixed(H, f.ROOTS, "D", {"c": "d"})) for i, (H, _) in enumerate(cforms)]
    graphs = []
    for un, UG in U:
        Ug = prefixed(UG, ("a", "b", "w"), "U")
        for ci, CG in C:
            for di, DG in D:
                for E in ([], [("a", "w")], [("b", "w")], [("a", "w"), ("b", "w")]):
                    G = psi2()
                    G.add_edges_from(Ug.edges())
                    G.add_edges_from(CG.edges())
                    G.add_edges_from(DG.edges())
                    G.add_edges_from(E)
                    graphs.append((G, {"case": "a", "U": un, "C form": ci, "D form": di,
                                       "edges in X": ["".join(e) for e in E]}))
    for ci, CG in C:
        for di, DG in D:
            G = psi2()
            G.add_edges_from((k, x) for k in ("k1", "k2") for x in "abw")
            G.add_edges_from(CG.edges())
            G.add_edges_from(DG.edges())
            graphs.append((G, {"case": "b", "C form": ci, "D form": di}))
    for ri, (RG, _) in enumerate(rforms):
        G = psi2()
        G.add_edges_from((k, x) for k in ("k1", "k2") for x in "abw")
        G.add_edges_from(prefixed(RG, f.R_ROOTS, "R").edges())
        graphs.append((G, {"case": "c", "R form": ri, "N(R)": f.r_signature(RG)[0]}))
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


def b_planar(G):
    """Whether B = G minus the twins, with no edges among a, b, c, d, is planar."""
    X = set("abcd")
    B = nx.Graph()
    B.add_nodes_from(v for v in G if v not in TWINS)
    B.add_edges_from(e for e in G.edges() if not set(e) & set(TWINS) and not set(e) <= X)
    return nx.check_planarity(B)[0]


def forms(out):
    cache = out / "forms.pkl"
    if cache.exists():
        with open(cache, "rb") as fh:
            return pickle.load(fh)
    cforms, cst = f.c_forms()
    rforms, rst = fx.e10_r_forms()
    with open(cache, "wb") as fh:
        pickle.dump((cforms, cst, rforms, rst), fh)
    return cforms, cst, rforms, rst


if __name__ == "__main__":
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    cforms, cst, rforms, rst = forms(out)
    print("D_c forms", len(cforms), "D_1 forms", len(rforms), flush=True)
    graphs = build(cforms, rforms)
    print("assembled", len(graphs), Counter(i["case"] for _, i in graphs), flush=True)
    uniq = dedupe(graphs)
    print("up to isomorphism", len(uniq), flush=True)
    forms_record = {"C forms": [{"index": i, "signature": list(f.signature(H)), "origin": inf["origin"],
                                 "graph": el.graph_record(H)} for i, (H, inf) in enumerate(cforms)],
                    "R forms": [{"index": i, "signature": list(f.r_signature(H)), "origin": inf["origin"],
                                 "graph": el.graph_record(H)} for i, (H, inf) in enumerate(rforms)],
                    "U forms": [{"name": n, "graph": el.graph_record(H)} for n, H in f.u_forms() if n in ELEVEN],
                    "C form statistics": cst, "R form statistics": rst}
    with open(out / "forms.json", "w") as fh:
        json.dump(forms_record, fh, indent=1)
    s = el.run_item("E10", [("assemblies", None, lambda G: G,
                             {"parts": uniq, "stats": {"assembled": len(graphs), "up to isomorphism": len(uniq)},
                              "payload": "see forms.json"})],
                    ["T56", "T72", "T73"], str(out),
                    note="G = Psi2 + (B_0, or k1, k2 on abw) + D_c + D_d (+ subsets of aw, bw in case (a)), "
                         "or Psi2 + k1, k2 + D_1")
    data = json.load(open(out / "candidates.json"))
    for r in data["candidates"]:
        G = nx.Graph()
        G.add_edges_from(r["graph"]["edges"])
        r["P planar"] = b_planar(G)
    with open(out / "candidates.json", "w") as fh:
        json.dump(data, fh)
    # the text concerns the graphs with B planar (a class counts when some member has B planar;
    # package.py records such a member)
    planar_found = set()
    by_key = {json.dumps(i, sort_keys=True): G for G, i in graphs}
    for r in data["candidates"]:
        if r["outcome"]["outcome"] == "catalogue graph":
            infos = [r["recipe"]["first"]] + r["recipe"]["also"]
            if any(b_planar(by_key[json.dumps(i, sort_keys=True)]) for i in infos):
                planar_found.add(r["outcome"]["catalogue"])
    s["catalogue graphs with B planar"] = sorted(planar_found, key=lambda x: int(x[1:]))
    s["agrees with the text"] = s["complete"] and s["catalogue graphs with B planar"] == ["T56", "T72", "T73"]
    with open(out / "summary.json", "w") as fh:
        json.dump(s, fh, indent=1)
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "catalogue graphs with B planar", "agrees with the text", "seconds")},
                     indent=1))
