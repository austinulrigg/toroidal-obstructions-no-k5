"""E5: assemblies of bridges on Psi_1 when some path is not an edge and
every four-attachment bridge has one interior vertex (Proposition 4.25).

The family is that of psi1_domain.py: 162 atoms (the nine edges s_i t_j, the
144 labeled three-attachment forms and the nine single vertices with four
attachments) and the clauses of the proof. It is covered as in covering.py
by maximal toroidal graphs and minimal nontoroidal graphs, with worker
processes (parallel_cover); every vector found is recorded with its images
under the twelve symmetries of Psi_1. The final formula is written to
cover.cnf, with a DRAT proof of its unsatisfiability found by Glucose 4 in
cover.drat. Then
  - each maximal toroidal graph gets a torus embedding (embeddings.py,
    multi_genus); an image of another one gets the image of its embedding;
  - each minimal nontoroidal graph that satisfies the clauses is classified:
    a cut of at most two vertices (enumeration.small_cut), else
    isomorphism with some T_i (enumeration.identify), else
    enumeration.classify with the seed k for the k-th nontoroidal vector
    (counting from 0); an image of another such graph gets the image of
    its outcome.
The proof finds T32, T43, T61 in part (a), T38, T45, T48, T66, T67 in
part (b), none in part (c), and T48, T58 among the graphs with all three
edges s_m t_m present.

Which vectors the covering finds depends on the order in which the worker
processes finish, so two runs record different, equally valid, families of
graphs, and the genus computations that the proper nontoroidal minors of a
new run need are not those of the data folder. With the option
--use-covering-of FOLDER the covering is not run: the vectors and their
symmetry images are read from the records candidates.json(.gz) of FOLDER
(the folder enumerations/E5 of the data), and everything after the covering
is done as above. This reproduces that folder.

The covering state is saved in cover_checkpoint.pkl, and the embeddings
and classifications in embedding_jobs.jsonl and classification_jobs.jsonl,
so that a second start continues.
Run as: python generate_E5.py FOLDER [--workers N] [--cover-seconds S]
                                     [--use-covering-of FOLDER]
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import gzip
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import sys
import time

import common
import covering as cv
import enumeration as el
import psi1_domain as sd

TEXT = {"(a)": ["T32", "T43", "T61"], "(b)": ["T38", "T45", "T48", "T66", "T67"], "(c)": [],
        "all three edges s_m t_m present": ["T48", "T58"]}


def atomic_json(path, obj):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as f:
        json.dump(obj, f, separators=(",", ":"))
        f.write("\n")
    os.replace(tmp, path)


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def label(atoms, n):
    return {atoms[i]["id"]: c for i, c in enumerate(n) if c}


def kind_of(atoms, n):
    forms = {atoms[i].get("form") for i, c in enumerate(n) if c and atoms[i]["kind"] in ("triple", "quad", "quad-c")}
    paths_edges = all(n[i] for i, a in enumerate(atoms) if a["kind"] == "matching")
    return {"all paths edges": paths_edges, "forms": sorted(f for f in forms if f),
            "singletons only": forms <= {"single"},
            "uses T, S or Q": bool(forms & {"Tx", "Ty", "Tz", "Sxy", "Syx", "Q"}),
            "uses a four-cycle form": bool(forms & {"FC1", "FC2"})}


def part(shape):
    if shape["all paths edges"]:
        return "all three edges s_m t_m present"
    if shape["singletons only"]:
        return "(a)"
    tsq, fc = shape["uses T, S or Q"], shape["uses a four-cycle form"]
    return "(b) and (c)" if tsq and fc else ("(b)" if tsq else "(c)")


def do_job(job):
    mode, idx, edges = job
    import networkx as nx
    import embeddings
    G = nx.Graph()
    G.add_edges_from(edges)
    if mode == "embedding":
        g, parts = embeddings.toroidal_embedding(G)
        if g > 1:
            raise RuntimeError("a toroidal graph got an embedding of genus larger than one")
        return idx, el.emb_record(parts)
    cut = el.small_cut(G)
    if cut is not None:
        return idx, {"outcome": "cut of order at most two", "cut": [str(v) for v in cut]}
    name, mapping = el.identify(G)
    if name is not None:
        return idx, {"outcome": "catalogue graph", "catalogue": name, "map": mapping}
    return idx, el.classify(G, seed=idx)


def cached_jobs(jobs, cache_path, workers, log, phase):
    """Runs the jobs with worker processes, appending each result to
    cache_path as it finishes; results already there are reused."""
    results = {}
    if cache_path.exists():
        lines = cache_path.read_text().splitlines()
        for pos, line in enumerate(lines):
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                if pos != len(lines) - 1:
                    raise
                break       # only an interrupted last line is dropped
            results[rec["key"]] = rec["result"]
        with cache_path.open("w") as f:
            for key, value in results.items():
                f.write(json.dumps({"key": key, "result": value}, separators=(",", ":")) + "\n")
    answers, todo = {}, []
    for job in jobs:
        key = digest(job)
        if key in results:
            answers[job[1]] = results[key]
        else:
            todo.append((key, job))
    log("%s: %d done before, %d to do" % (phase, len(answers), len(todo)))
    if not todo:
        return answers
    with cache_path.open("a") as f, ProcessPoolExecutor(
            max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        remaining = iter(todo)
        pending = {}

        def fill():
            while len(pending) < workers * 2:
                try:
                    key, job = next(remaining)
                except StopIteration:
                    break
                pending[pool.submit(do_job, job)] = key
        fill()
        last = time.monotonic()
        while pending:
            done, _ = wait(pending, timeout=15, return_when=FIRST_COMPLETED)
            for fut in done:
                key = pending.pop(fut)
                idx, result = fut.result()
                answers[idx] = result
                f.write(json.dumps({"key": key, "result": result}, separators=(",", ":")) + "\n")
                f.flush()
            fill()
            if time.monotonic() - last >= 60 or not pending:
                log("%s: %d of %d done" % (phase, len(answers), len(jobs)))
                last = time.monotonic()
    return answers


def covering_of(folder, atoms, perms, mx):
    """The vectors P, F and their symmetry origins recorded in the
    candidates of a data folder of E5."""
    path = Path(folder) / "candidates.json"
    if path.exists():
        with path.open() as f:
            cands = json.load(f)["candidates"]
    else:
        with gzip.open(str(path) + ".gz", "rt") as f:
            cands = json.load(f)["candidates"]
    P, F, origins = [], [], {}
    for c in cands:
        kind = "P" if c["kind"] == "toroidal, maximal" else "F"
        lst = P if kind == "P" else F
        assert c["id"] == "%s%05d" % (kind, len(lst)), c["id"]
        n = [0] * len(atoms)
        for a, k in c["atoms"].items():
            n[int(a[1:]) - 1] = k
        assert all(0 <= x <= m for x, m in zip(n, mx))
        n = tuple(n)
        t = len(lst)
        if c.get("image of"):
            rep = int(c["image of"][1:])
            k = next(k for k, p in enumerate(perms) if cv.image(lst[rep], p) == n)
        else:
            rep, k = t, 0
        origins[(kind, t)] = (rep, k)
        lst.append(n)
    return P, F, origins


def run(args):
    start = time.monotonic()
    out = Path(args.folder).resolve()
    out.mkdir(parents=True, exist_ok=True)

    def log(msg):
        print("%s (%.0fs)" % (msg, time.monotonic() - start), flush=True)
    atoms = sd.atoms()
    mx = sd.maxima(atoms)
    dom = cv.Domain(atoms, sd.base_graph(), mx, lambda v: sd.clauses(atoms, mx, v), "E5")
    maps = sd.label_symmetries()
    perms = sd.atom_permutations(atoms, maps, mx)
    atomic_json(out / "domain.json", {"atoms": atoms, "maxima": mx,
                                      "base edges": sorted(sorted(e) for e in dom.base.edges())})
    log("%d atoms, %d variables, %d base clauses; %d workers, %d symmetries" %
        (len(atoms), dom.nvars, len(dom.base_clauses), args.workers, len(perms)))
    if args.use_covering_of:
        P, F, origins = covering_of(args.use_covering_of, atoms, perms, mx)
        cnf = list(dom.base_clauses) + [dom.above_clause(v) for v in P] + [dom.below_clause(v) for v in F]
        # The recorded formula lists the same clauses; its order is kept, so
        # that the DRAT proof found for it is that of the folder.
        recorded = Path(args.use_covering_of) / "cover.cnf"
        nv, rec = cv.read_dimacs(recorded if recorded.exists() else str(recorded) + ".gz")
        key = lambda cl: Counter(tuple(sorted(c)) for c in cl)
        if nv != dom.nvars or key(rec) != key(cnf):
            raise SystemExit("the recorded cover.cnf is not the formula of the recorded covering")
        cnf = rec
        log("covering read from %s" % args.use_covering_of)
    else:
        result = cv.parallel_cover(dom, out, workers=args.workers, symmetries=perms, log=log,
                                   max_seconds=args.cover_seconds)
        if not result["complete"]:
            log("the covering is not complete; its state is saved. Run the same command again to continue.")
            return 75
        P, F, origins, cnf = (result[k] for k in ("P", "F", "origins", "cnf"))
    cv.write_dimacs(out / "cover.cnf", dom.nvars, cnf)
    log("covering: %d maximal toroidal, %d minimal nontoroidal; proving that cover.cnf is unsatisfiable"
        % (len(P), len(F)))
    verified = cv.fresh_proof(out / "cover.cnf", out / "cover.drat")
    log("drat-trim: %s" % {True: "verified", False: "NOT VERIFIED", None: "not found, proof not checked"}[verified])
    if verified is False:
        raise SystemExit("the DRAT proof was not verified")

    def relabeler(k):
        root_map, pi = maps[k], perms[k]
        ids = {atoms[i]["id"]: atoms[pi[i]]["id"] for i in range(len(atoms))}

        def f(v):
            if "#" not in v:
                return root_map[v]
            head, tail = v.split("#", 1)
            return ids[head] + "#" + tail
        return f

    jobs = [("embedding", t, sorted(sorted(e) for e in dom.graph(n).edges()))
            for t, n in enumerate(P) if origins[("P", t)][0] == t]
    embedded = cached_jobs(jobs, out / "embedding_jobs.jsonl", args.workers, log, "torus embeddings")
    passing = []
    for t, n in enumerate(P):
        rep, k = origins[("P", t)]
        if rep == t:
            emb = embedded[t]
        else:
            f = relabeler(k)
            emb = [{"vertices": [f(v) for v in p["vertices"]], "genus": p["genus"],
                    "rotation": {f(v): [f(x) for x in nb] for v, nb in p["rotation"].items()}}
                   for p in passing[rep]["embedding"]]
        passing.append({"id": "P%05d" % t, "vector": label(atoms, n), "embedding": emb,
                        "image of": None if rep == t else "P%05d" % rep})
    atomic_json(out / "passing.json", {"passing": passing})
    failing, jobs = [], []
    for t, n in enumerate(F):
        rep, k = origins[("F", t)]
        rec = {"id": "F%05d" % t, "vector": label(atoms, n), "in domain": cv.in_domain(dom, n),
               "image of": None if rep == t else "F%05d" % rep}
        if rec["in domain"]:
            G = dom.graph(n)
            rec["graph"] = el.graph_record(G)
            rec["shape"] = kind_of(atoms, n)
            if rep == t or not failing[rep]["in domain"]:
                jobs.append(("classify", t, rec["graph"]["edges"]))
        failing.append(rec)
    outcomes = cached_jobs(jobs, out / "classification_jobs.jsonl", args.workers, log, "classifications")
    for t, rec in enumerate(failing):
        if not rec["in domain"]:
            continue
        if t in outcomes:
            rec["outcome"] = outcomes[t]
            continue
        rep, k = origins[("F", t)]
        o, f = dict(failing[rep]["outcome"]), relabeler(k)
        if o["outcome"] == "K5 minor":
            o["branch_sets"] = [[f(v) for v in b] for b in o["branch_sets"]]
        elif o["outcome"] == "cut of order at most two":
            o["cut"] = [f(v) for v in o["cut"]]
        elif o["outcome"] == "catalogue graph":
            o["map"] = {f(v): w for v, w in o["map"].items()}
        elif o["outcome"] == "proper nontoroidal minor":
            o["edge"] = sorted(f(v) for v in o["edge"])
        else:
            o = el.classify(dom.graph(F[t]), seed=t)
        rec["outcome"] = o
    atomic_json(out / "failing.json", {"failing": failing})
    counts = Counter(r["outcome"]["outcome"] for r in failing if r["in domain"])
    found = {}
    for r in failing:
        if r.get("outcome", {}).get("outcome") == "catalogue graph":
            found.setdefault(part(r["shape"]), set()).add(r["outcome"]["catalogue"])
    found = {p: sorted(v, key=lambda x: int(x[1:])) for p, v in sorted(found.items())}
    by_part = {p: [] for p in TEXT}
    for p, names in found.items():
        for q in ("(a)", "(b)", "(c)") if p == "(b) and (c)" else (p,):
            by_part[q] = sorted(set(by_part[q]) | set(names), key=lambda x: int(x[1:]))
    summary = {"item": "E5", "where": "Proposition 4.25", "atoms": len(atoms), "variables": dom.nvars,
               "base clauses": len(dom.base_clauses), "maximal toroidal": len(P), "minimal nontoroidal": len(F),
               "minimal nontoroidal satisfying the clauses": sum(r["in domain"] for r in failing),
               "outcomes": dict(counts), "catalogue graphs by part": found,
               "drat-trim": {True: "verified", None: "not run"}[verified],
               "agrees with the text": by_part == TEXT and not counts.get("unresolved")
               and not counts.get("new obstruction"),
               "seconds": round(time.monotonic() - start, 1)}
    with (out / "summary.json").open("w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("folder")
    parser.add_argument("--workers", type=int, default=common.workers())
    parser.add_argument("--cover-seconds", type=float)
    parser.add_argument("--use-covering-of")
    args = parser.parse_args()
    raise SystemExit(run(args))
