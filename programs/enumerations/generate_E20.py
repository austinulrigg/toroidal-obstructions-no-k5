"""E20: combinations of the bridge forms on three attachments (Proposition 5.20).

The program's names: x1, x2 are the twins on abc and y1, y2 those on abd.
By Lemma 5.19 every bridge on a triple is one of eight forms on (x; y, z):
a single vertex, T(z), T(y), S(x;y), S(x;z), Q(x;y,z), and two four-cycles
joined by a path of length one or two (forms_e20.py). Placed with every
bijection of (x, y, z) onto a triple of {a, b, c, d}, they give 76 labeled
forms, 19 on each triple. A combination is a multiset of one to four
labeled forms, not all single vertices, with cd absent or present.

As the proof states, a combination cannot give G if it has two disjoint
twin triples, three distinct twin triples sharing a pair, or a bridge with
four attachments for some choice of twin pairs, and these conditions
remain when another form or cd is added; a combination containing a
smaller nontoroidal combination has a proper nontoroidal minor. So the
combinations are processed by size: a combination is built and classified
only when every combination with one form fewer, or without cd, is toroidal
and satisfies the three conditions. Its status (nontoroidal, violating a
condition, or neither) is recorded for the later sizes; the combinations
of single vertices alone are classified as well, as smaller combinations,
and recorded under 'lattice' with the forms. Every classified combination
is classified as a whole graph (enumeration.classify). The proof finds T15.

Run as: python generate_E20.py FOLDER
"""
import itertools
import json
import sys
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import networkx as nx

import common
import enumeration as el
import forms_e20 as f
import multigenus as mg
import psi2_runs as ext


def immediate_subs(c, cd):
    out = [(c[:i] + c[i + 1:], cd) for i in range(len(c)) if i == 0 or c[i] != c[i - 1]]
    if cd:
        out.append((c, False))
    return out


def _classify(job):
    idx, edges, nodes, seed = job
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    return idx, el.classify(G, seed=seed)


if __name__ == "__main__":
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    lib = f.form_library()
    N = len(lib)
    stars = {i for i, F in enumerate(lib) if F["name"] == "star"}
    status = {((), False): "toroidal, no violation", ((), True): "toroidal, no violation"}
    records = []
    lattice = Counter()
    for size in (1, 2, 3, 4):
        for cd in (False, True):
            batch = []
            for c in itertools.combinations_with_replacement(range(N), size):
                in_domain = not set(c) <= stars
                if all(status.get(s) == "toroidal, no violation" for s in immediate_subs(c, cd)):
                    batch.append((c, cd))
                elif in_domain:
                    lattice["excluded through a sub-configuration"] += 1
            Gs = [f.config_graph(lib, c, d) for c, d in batch]
            over = mg.nontoroidal(Gs) if Gs else []
            for (c, d), G, o in zip(batch, Gs, over):
                viol = f.structural_violation(G)
                st = "nontoroidal" if o else ("violation" if viol else "toroidal, no violation")
                status[(c, d)] = st
                in_domain = not set(c) <= stars
                lattice["classified" if in_domain else "auxiliary (stars only)"] += 1
                rec = {"id": "E20-%s-%s" % ("-".join(lib[i]["id"] for i in c), "cd" if d else "nocd"),
                       "case": "size %d" % size, "variant": "cd" if d else "",
                       "config": {"forms": [lib[i]["id"] for i in c], "cd": d},
                       "in domain": in_domain, "lattice status": st,
                       "part": {"vertices": [], "edges": []}, "recipe": {},
                       "graph": el.graph_record(G)}
                if viol:
                    rec["structural exclusion"] = dict(viol, kind="twin structure")
                records.append(rec)
    jobs = [(i, [list(e) for e in nx.Graph([tuple(e) for e in r["graph"]["edges"]]).edges()],
             r["graph"]["vertices"], i) for i, r in enumerate(records)]
    with Pool(common.workers()) as pool:
        for idx, res in pool.imap_unordered(_classify, jobs, chunksize=4):
            records[idx]["outcome"] = res
    # auxiliary star-only records are kept for the lattice, but are not candidates of this item
    cand = [r for r in records if r["in domain"]]
    aux = [r for r in records if not r["in domain"]]
    total = 2 * sum(len(list(itertools.combinations_with_replacement(range(N), k))) for k in (1, 2, 3, 4)) \
        - 2 * sum(len(list(itertools.combinations_with_replacement(range(len(stars)), k))) for k in (1, 2, 3, 4))
    STRUCT = {"disjoint twin triples": {"kind": "twin structure"},
              "three twin triples share a pair": {"kind": "twin structure"},
              "component adjacent to all four attachments": {"kind": "twin structure"}}
    summary = ext.finish_whole_item(
        "E20", cand,
        {"forms": {"specification": "76 labeled forms, forms_e20.form_library()", "labeled forms": N}},
        ["T15"], out, STRUCT, "x1,x2 on abc; y1,y2 on abd; forms named f<k>_<i>", t0,
        extra_summary={"domain configurations": total, "lattice": dict(lattice),
                       "lattice statuses of classified configurations":
                           dict(Counter(r["lattice status"] for r in cand))})
    # rewrite candidates.json with the lattice data needed by the checker
    data = json.load(open(out / "candidates.json"))
    data["lattice"] = {"forms": lib, "domain": "multisets of 1-4 forms, not all stars, cd absent or present",
                       "auxiliary": aux}
    json.dump(data, open(out / "candidates.json", "w"))
    json.dump(lib, open(out / "forms.json", "w"), indent=1)
    print(json.dumps({k: summary[k] for k in ("domain configurations", "lattice", "candidates",
                                               "outcomes", "lattice statuses of classified configurations",
                                               "survivors", "agrees with the text", "seconds")}, indent=1))
