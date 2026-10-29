"""E27: other bridges whose triangle completion is nonplanar (Lemma 6.13).

The program's names: p1, p2 are the twins on abc and q1, q2 those on ade.
A component W has neighbourhood {a, c, d, e}, contracted to a vertex w.
Let U be another bridge, on a triple t, with nonplanar triangle completion
and some two-edge completion U + xy + xz planar. By Lemma 3.8, U has as an
allowed minor one of the forms T, S, Q of Section 3, whose twins must meet
abc and ade in exactly one vertex each. This program lists the pairs
(triple, form) allowed by that condition, over the five triples abc, ade,
abd, acd, bcd, every ordering of the triple and every case of Lemma 3.8;
the proof uses the triples ade, acd and abd, and the data folder keeps
those.
For each pair it builds the minor of G of the proof:
H = Psi_3 + w + F, with the edge bd added when t is abc, ade or acd (an
outside b-d path avoiding U, shortened). Each H is classified with
enumeration.classify, without the test for a cut of at most two vertices
(H is a minor of G, not G itself), and each nontoroidal H gets a sequence
of edge deletions and contractions down to a graph isomorphic to some T_i.
When the recorded minor of T_i does not survive the operation that
classify names, the first step of the sequence is named instead.
The proof finds nine forms on these triples, eight of them with a proper
nontoroidal minor; the ninth, Q(a;b,d), gives a toroidal graph.

Run as: python generate_E27.py FOLDER
"""
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import networkx as nx
import enumeration as el

import psi3_runs as ext

TRIPLES = ["abc", "ade", "abd", "acd", "bcd"]
ABC, ADE = set("abc"), set("ade")


def forms_on(t):
    """Every (name, graph, twin triple) from Lemma 3.8 on the triple t."""
    out = {}
    for x, y, z in itertools.permutations(t):
        cands = [("T(%s)" % z, ext.form_T(z, x, y), {x, y}),          # case (ii), D0 = {z}
                 ("S(%s;%s)" % (x, y), ext.form_S(x, y, z), {x, y}),   # case (ii)
                 ("T(%s)" % y, ext.form_T(y, x, z), {x, z}),          # case (iii), D0 = {y}
                 ("S(%s;%s)" % (x, z), ext.form_S(x, z, y), {x, z}),   # case (iii)
                 ("Q(%s;%s,%s)" % (x, *sorted((y, z))), ext.form_Q(x, *sorted((y, z))), {x})]  # case (i)
        for name, F, Ytop in cands:
            out.setdefault(name, (F, Ytop))
    return out


def survives(H, op, edge, steps, name, mapping):
    """Whether the minor of T_name given by steps (with mapping from its
    vertices to those of T_name) is still a minor of H after the named
    deletion or contraction of edge: the branch sets the steps define,
    changed by the operation, form a model of T_name in the result."""
    cur = H.copy()
    blocks = {v: {v} for v in H}
    for o, (u, v) in steps:
        if o == "delete":
            cur.remove_edge(u, v)
        else:
            cur = el.gr.contract_edge(cur, (u, v))
            blocks[u] |= blocks.pop(v)
    T = dict(el.catalog())[name]
    model = {mapping[x]: blocks[x] for x in cur if cur.degree(x) > 0}
    u, v = edge
    if op == "delete":
        H2 = el.gr.delete_edge(H, (u, v))
    else:
        H2 = el.gr.contract_edge(H, (u, v))
        model = {t: {u if y == v else y for y in B} for t, B in model.items()}
    used = [y for B in model.values() for y in B]
    if len(used) != len(set(used)) or set(model) != set(T):
        return False
    if not all(nx.is_connected(H2.subgraph(B)) for B in model.values()):
        return False
    return all(any(H2.has_edge(x, y) for x in model[s] for y in model[t]) for s, t in T.edges())


def admissible(Ytop):
    """Y = Ytop plus interior vertices meets abc and ade exactly once each."""
    return len(Ytop & ABC) == 1 and len(Ytop & ADE) == 1


if __name__ == "__main__":
    t0 = time.time()
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    records, rejected = [], []
    for t in TRIPLES:
        for name, (F, Ytop) in sorted(forms_on(t).items()):
            if not admissible(Ytop):
                rejected.append([t, name, "".join(sorted(Ytop))])
                continue
            H = ext.psi3()
            ext.add_vertex(H, "w", "acde")
            H.add_edges_from(F.edges())
            bd = not ({"b", "d"} <= set(t))
            if bd:
                H.add_edge("b", "d")
            res = el.classify(H, seed=len(records), need_three_connected=False)
            rec = {"id": "E27-%s-%s" % (t, name), "case": t, "form": name,
                   "bd added": bd, "part": el.graph_record(F), "graph": el.graph_record(H),
                   "outcome": res}
            if res["outcome"] != "toroidal" and res["outcome"] != "K5 minor":
                steps, k, m = ext.obstruction_minor(H, seed=len(records))
                rec["obstruction minor"] = {"steps": steps, "catalogue": k, "map": m}
                if res["outcome"] == "proper nontoroidal minor" and not survives(
                        H, res["operation"], res["edge"], steps, k, m):
                    # name the first step of the sequence, which the minor survives
                    res["operation"], res["edge"] = steps[0][0], list(steps[0][1])
            records.append(rec)
    counts = Counter(r["outcome"]["outcome"] for r in records)
    # the proof uses the triples ade, acd and abd: nine forms, eight excluded
    used = [r for r in records if r["case"] in ("ade", "acd", "abd")]
    not_excluded = [r["form"] + " on " + r["case"] for r in used
                    if r["outcome"]["outcome"] not in ("K5 minor", "proper nontoroidal minor")]
    target = ["Q(a;b,d) on abd"]
    summary = {"item": "E27", "where": "Lemma 6.13",
               "forms admitted (twin triple meets abc and ade once each)": len(records),
               "forms rejected by the twin-triple condition": rejected,
               "outcomes": dict(counts),
               "table": [{"triple": r["case"], "form": r["form"], "bd added": r["bd added"],
                          "outcome": r["outcome"]["outcome"],
                          "obstruction minor": (r.get("obstruction minor") or {}).get("catalogue")}
                         for r in records],
               "not excluded": not_excluded, "exceptions in the text": target,
               "agrees with the text": sorted(not_excluded) == sorted(target),
               "seconds": round(time.time() - t0, 1)}
    json.dump({"item": "E27", "candidates": records}, open(out / "candidates.json", "w"))
    json.dump(summary, open(out / "summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))
