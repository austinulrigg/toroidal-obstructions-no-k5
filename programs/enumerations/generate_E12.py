"""E12: the cut X = {a, u, v} with three other components (Proposition 5.24).

The program's names: t0, t1 are the twins on abc and t2, t3 those on abd
(p, q, r, z in the paper), and p, q have neighbourhood X = {a, u, v} (p', q'
in the paper). G = Psi_2 + p, q + R_b + R_c + R_d, where R_g is the
component of P - X containing g, for g in {b, c, d}. Either R_g is g alone,
adjacent to u, v or both, or R_g has neighbourhood st in {au, av, uv} in X
and R_g - g is one of the nine forms on (s, t, g) as (x, y, z) of
Proposition 3.15 and Lemma 3.14 (one vertex, T(s), T(t), T(g), S(s;t),
S(t;s), Q(g;s,t) and the two four-cycle forms), or two single vertices on
{g, s, t}; the edges gs, gt of the original graph may be present, but never
ag.

Atoms: for each g, each pair st and each form, the form on {g, s, t} (at
most two copies of the one-vertex form, which give the two single
vertices, and one copy of the others); and the edges gu, gv. Clauses: for
each g at most one form, a form on {g, s, t} excludes the edge from g to the
third vertex of X, and some atom is present for each g. The family is
covered as in covering.py (one process, the SAT solver CaDiCaL 1.5.3); the
final formula is written to cover.cnf and a DRAT proof of its
unsatisfiability, found by Glucose 4, to cover.drat.

Each maximal toroidal graph gets a torus embedding (PAGE,
graphs.page_embedding), and each minimal nontoroidal graph that satisfies
the clauses is classified (enumeration.classify; the random search for a K5
minor in the k-th one, counting from 0 over all nontoroidal graphs, uses the
seed k). The proof finds T51 and T72.

The covering is kept in cover.pkl in the output folder and the
classification in progress.jsonl, so that a second start continues.
Run as: python generate_E12.py FOLDER
"""
from collections import Counter
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
import json
import pickle
import sys
import time

import networkx as nx
from pysat.formula import IDPool

import common
import covering as cv
import enumeration as el
import graphs as gr
import psi1_domain as sd

FORMS9 = ("single", "Tx", "Ty", "Tz", "Sxy", "Syx", "Q", "FC1", "FC2")


def family():
    base = nx.Graph()
    base.add_edges_from((t, x) for t in ("t0", "t1") for x in "abc")
    base.add_edges_from((t, x) for t in ("t2", "t3") for x in "abd")
    base.add_edges_from((p, x) for p in "pq" for x in "auv")
    atoms = []
    for r in "bcd":
        for s, t in (("a", "u"), ("a", "v"), ("u", "v")):
            for form in FORMS9:
                atoms.append({"id": "%s_%s%s_%s" % (r, s, t, form), "kind": "form", "r": r,
                              "support": [s, t], "form": form, "edges": sd.form_edges(form, s, t, r),
                              "attachments": [s, t, r]})
        for s in "uv":
            atoms.append({"id": "%s_edge_%s" % (r, s), "kind": "edge", "r": r, "end": s,
                          "edges": [(r, s)], "attachments": [r, s]})
    mx = [2 if a.get("form") == "single" else 1 for a in atoms]

    def clauses(var):
        pool = IDPool(start_from=sum(mx) + 1)
        out = []
        for r in "bcd":
            forms = [i for i, a in enumerate(atoms) if a["r"] == r and a["kind"] == "form"]
            edges = {a["end"]: i for i, a in enumerate(atoms) if a["r"] == r and a["kind"] == "edge"}
            for i, j in combinations(forms, 2):
                # one support and one form (two single vertices are two copies of one atom)
                out.append([-var(i, 1), -var(j, 1)])
            for i in forms:
                s, t = atoms[i]["support"]
                for end, k in edges.items():
                    if end not in (s, t):
                        out.append([-var(i, 1), -var(k, 1)])
            out.append([var(i, 1) for i in forms] + [var(k, 1) for k in edges.values()])
        return out, pool.top

    return atoms, mx, base, clauses


def two_single(atoms, n):
    return any(c == 2 for i, c in enumerate(n) if atoms[i].get("form") == "single")


def label(atoms, n):
    return {atoms[i]["id"]: c for i, c in enumerate(n) if c}


def main(out):
    t0 = time.time()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    log = lambda m: print("%s (%.0fs)" % (m, time.time() - t0), flush=True)
    atoms, mx, base, clauses = family()
    dom = cv.Domain(atoms, base, mx, clauses, "E12")
    log("%d atoms, %d variables, %d base clauses" % (len(atoms), dom.nvars, len(dom.base_clauses)))
    cache = out / "cover.pkl"
    if cache.exists():
        with open(cache, "rb") as f:
            P, F, cnf = pickle.load(f)
    else:
        P, F, cnf = cv.cover(dom, log=log)
        with open(cache, "wb") as f:
            pickle.dump((P, F, cnf), f)
    log("covering: %d maximal toroidal, %d minimal nontoroidal" % (len(P), len(F)))
    cv.write_dimacs(out / "cover.cnf", dom.nvars, cnf)
    verified = cv.fresh_proof(out / "cover.cnf", out / "cover.drat")
    log("drat-trim: %s" % {True: "verified", False: "NOT VERIFIED", None: "not found, proof not checked"}[verified])
    if verified is False:
        raise SystemExit("the DRAT proof was not verified")
    passing = []
    for k, n in enumerate(P):
        g, parts = gr.page_embedding(dom.graph(n))
        assert g <= 1
        passing.append({"id": "P%05d" % k, "vector": label(atoms, n), "embedding": el.emb_record(parts)})
    log("torus embeddings of %d graphs" % len(passing))
    failing, jobs = [], []
    for k, n in enumerate(F):
        rec = {"id": "F%05d" % k, "vector": label(atoms, n), "in domain": cv.in_domain(dom, n)}
        if rec["in domain"]:
            G = dom.graph(n)
            rec["graph"] = el.graph_record(G)
            rec["two single vertices"] = two_single(atoms, n)
            jobs.append((k, [list(e) for e in G.edges()], k))
        failing.append(rec)
    progress = out / "progress.jsonl"
    done = {}
    if progress.exists():
        for line in progress.read_text().splitlines():
            try:
                r = json.loads(line)
            except ValueError:
                continue        # an interrupted last line
            done[r["index"]] = r
    todo = []
    for job in jobs:
        r = done.get(job[0])
        if r is not None and r["edges"] == sorted(sorted(e) for e in job[1]):
            failing[job[0]]["outcome"] = r["outcome"]
        else:
            todo.append(job)
    log("classification: %d done before, %d to do" % (len(jobs) - len(todo), len(todo)))
    with progress.open("a") as pf, Pool(common.workers()) as pool:
        for n, (k, res) in enumerate(pool.imap_unordered(el._classify_job, todo, chunksize=4), 1):
            failing[k]["outcome"] = res
            pf.write(json.dumps({"index": k, "edges": sorted(sorted(e) for e in failing[k]["graph"]["edges"]),
                                 "outcome": res}) + "\n")
            pf.flush()
            if n % 500 == 0:
                log("classified %d of %d" % (n + len(jobs) - len(todo), len(jobs)))
    outcomes = Counter(r["outcome"]["outcome"] for r in failing if "outcome" in r)
    log("classified: %s" % dict(outcomes))
    found = sorted({r["outcome"]["catalogue"] for r in failing
                    if r.get("outcome", {}).get("outcome") == "catalogue graph"}, key=lambda x: int(x[1:]))
    with open(out / "domain.json", "w") as f:
        json.dump({"atoms": atoms, "maxima": mx, "base edges": sorted(sorted(e) for e in base.edges())}, f)
    with open(out / "passing.json", "w") as f:
        json.dump({"passing": passing}, f)
    with open(out / "failing.json", "w") as f:
        json.dump({"failing": failing}, f)
    summary = {"item": "E12", "where": "Proposition 5.24", "atoms": len(atoms), "variables": dom.nvars,
               "base clauses": len(dom.base_clauses), "maximal toroidal": len(P), "minimal nontoroidal": len(F),
               "minimal nontoroidal satisfying the clauses": len(jobs), "outcomes": dict(outcomes),
               "catalogue graphs found": found,
               "with two single vertices": sorted({r["outcome"]["catalogue"] for r in failing
                                                   if r.get("two single vertices")
                                                   and r["outcome"]["outcome"] == "catalogue graph"}),
               "drat-trim": {True: "verified", None: "not run"}[verified],
               "agrees with the text": found == ["T51", "T72"] and not outcomes.get("unresolved")
               and not outcomes.get("new obstruction"),
               "seconds": round(time.time() - t0, 1)}
    with open(out / "summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
