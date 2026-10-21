"""Helpers for the enumerations of Section 6 (E21-E34), added to the
shared enumeration.py.

Vertex labels (Section 6): the attachments are a, b, c, d, e; the twins
p1, p2 have neighborhood abc and q1, q2 have neighborhood ade. Psi3 is these
nine vertices with the twelve edges at p1, p2, q1, q2.
"""

import itertools
import json
from pathlib import Path

import networkx as nx

import common
import enumeration as el

ATT = ("a", "b", "c", "d", "e")


def psi3():
    """Psi_3: twins p1, p2 on abc and q1, q2 on ade."""
    G = nx.Graph()
    G.add_edges_from((p, x) for p in ("p1", "p2") for x in "abc")
    G.add_edges_from((q, x) for q in ("q1", "q2") for x in "ade")
    return G


def add_vertex(G, name, nbrs):
    G.add_edges_from((name, x) for x in nbrs)
    return G


def with_edges(G, edges):
    H = G.copy()
    H.add_edges_from(edges)
    return H


def triangle(t):
    return [(t[0], t[1]), (t[0], t[2]), (t[1], t[2])]


# ----------------------------------------------------------------------------
# The bridge forms T, S, Q of Lemma 3.7 (core:three-root-cut).
# Interior vertices get the prefix pre so that several forms can coexist.

def form_T(z, x, y, pre="U"):
    """T(z): d1, d2 joined to x, y, u, and the edge uz."""
    d1, d2, u = pre + "d1", pre + "d2", pre + "u"
    B = nx.Graph()
    B.add_edges_from((d, w) for d in (d1, d2) for w in (x, y, u))
    B.add_edge(u, z)
    return B


def form_S(x, y, z, pre="U"):
    """S(x;y): d1, d2 joined to x, y, u; f joined to x, z, u."""
    d1, d2, u, f = pre + "d1", pre + "d2", pre + "u", pre + "f"
    B = nx.Graph()
    B.add_edges_from((d, w) for d in (d1, d2) for w in (x, y, u))
    B.add_edges_from((f, w) for w in (x, z, u))
    return B


def form_Q(x, y, z, pre="U"):
    """Q(x;y,z): four-cycle u d1 v d2 and edges x d1, x d2, u y, v z."""
    d1, d2, u, v = pre + "d1", pre + "d2", pre + "u", pre + "v"
    B = nx.Graph()
    B.add_edges_from([(u, d1), (d1, v), (v, d2), (d2, u), (x, d1), (x, d2), (u, y), (v, z)])
    return B


# ----------------------------------------------------------------------------
# Certificates for minors of G that need not equal G

def obstruction_minor(H, seed=0):
    """For a graph H of genus at least two, a sequence of edge deletions and
    contractions down to a graph all of whose one-step minors are toroidal,
    together with its identification in T1..T74 (after deleting isolated
    vertices). Uses multi_genus for every genus question. Returns
    (steps, name, map) or (steps, None, None)."""
    import random
    rng = random.Random(seed)
    cur = H.copy()
    steps = []
    assert el.genus_above_one([cur])[0]
    while True:
        es = sorted(cur.edges(), key=el.edge_key)
        rng.shuffle(es)
        cands = [("delete", e, el.gr.delete_edge(cur, e)) for e in es] + \
                [("contract", e, el.gr.contract_edge(cur, e)) for e in es]
        flags = el.genus_above_one([g for _, _, g in cands])
        nxt = None
        for (op, e, g), f in zip(cands, flags):
            if f:
                nxt = (op, e, g)
                break
        if nxt is None:
            break
        op, e, g = nxt
        steps.append([op, [str(x) for x in e]])
        cur = g
    core = cur.copy()
    core.remove_nodes_from([v for v in list(core) if core.degree(v) == 0])
    name, m = el.identify(core)
    return steps, name, m


def save_json(obj, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=1)


# ----------------------------------------------------------------------------
# Minor models of a fixed graph

def catalogue_graph(name):
    for nm, T in el.catalog():
        if nm == name:
            return nx.relabel_nodes(T, {v: str(v) for v in T})
    raise KeyError(name)


def _forests(H, j):
    es = sorted(H.edges(), key=el.edge_key)
    for F in itertools.combinations(es, j):
        parent = {}

        def find(v):
            while parent.get(v, v) != v:
                v = parent[v]
            return v
        ok = True
        for u, v in F:
            a, b = find(u), find(v)
            if a == b:
                ok = False
                break
            parent[a] = b
        if ok:
            yield F


def find_minor_model(H, T):
    """Branch sets of a minor model of T in H, or None if there is none.
    Exhaustive: every model arises from deleting a set S of vertices and
    contracting a forest F of H - S with |S| + |F| = |H| - |T|, after which
    T is a spanning subgraph of the result. Returns {T vertex: [H vertices]}."""
    n, t = H.number_of_nodes(), T.number_of_nodes()
    if n < t or H.number_of_edges() < T.number_of_edges():
        return None
    r = n - t
    tdeg = sorted((d for _, d in T.degree()), reverse=True)
    tm = T.number_of_edges()
    for i in range(r, -1, -1):          # deletions first: cheapest to test
        for S in itertools.combinations(sorted(H, key=str), i):
            H1 = H.copy()
            H1.remove_nodes_from(S)
            for F in _forests(H1, r - i):
                parent = {v: v for v in H1}

                def find(v):
                    while parent[v] != v:
                        v = parent[v]
                    return v
                for u, v in F:
                    parent[find(v)] = find(u)
                Q = nx.Graph()
                Q.add_nodes_from({find(v) for v in H1})
                Q.add_edges_from((find(u), find(v)) for u, v in H1.edges() if find(u) != find(v))
                if Q.number_of_edges() < tm:
                    continue
                qdeg = sorted((d for _, d in Q.degree()), reverse=True)
                if any(a < b for a, b in zip(qdeg, tdeg)):
                    continue
                GM = nx.algorithms.isomorphism.GraphMatcher(Q, T)
                for m in GM.subgraph_monomorphisms_iter():
                    bags = {}
                    for v in H1:
                        bags.setdefault(find(v), []).append(str(v))
                    return {str(tv): sorted(bags[qv]) for qv, tv in m.items()}
    return None


def check_minor_model(H, T, bags):
    """Check a minor model: disjoint connected bags, an edge between the bags
    of adjacent vertices of T."""
    used = [v for b in bags.values() for v in b]
    assert len(used) == len(set(used)) and set(used) <= set(map(str, H))
    Hs = nx.relabel_nodes(H, {v: str(v) for v in H})
    for b in bags.values():
        assert b and nx.is_connected(Hs.subgraph(b))
    for x, y in T.edges():
        assert any(Hs.has_edge(u, v) for u in bags[str(x)] for v in bags[str(y)])
    return True


# ----------------------------------------------------------------------------
# A faster whole-graph classification for large domains

def classify_minor_first(G, seed=0, need_three_connected=True):
    """The outcomes of enumeration.classify, tried in a different order: first
    a proper minor G - e or G / e of genus at least two (one multi_genus
    batch), and only if there is none, enumeration.classify itself. Every
    outcome carries the same certificate as in enumeration.classify, and is a
    valid exclusion for a candidate that must equal G."""
    es = sorted(G.edges(), key=el.edge_key)
    minors = [("delete", e, el.gr.delete_edge(G, e)) for e in es] + \
             [("contract", e, el.gr.contract_edge(G, e)) for e in es]
    flags = el.genus_above_one([H for _, _, H in minors])
    for (op, e, H), f in zip(minors, flags):
        if f:
            return {"outcome": "proper nontoroidal minor", "operation": op,
                    "edge": [str(x) for x in el.edge_key(e)]}
    return el.classify(G, seed=seed, need_three_connected=need_three_connected)


def _job_minor_first(job):
    idx, edges, seed = job
    G = nx.Graph()
    G.add_edges_from(edges)
    return idx, classify_minor_first(G, seed=seed)


def run_item_minor_first(item, cases, targets, outdir, note=""):
    """enumeration.run_item with classify_minor_first in place of classify;
    the files written have the same format."""
    import time
    from collections import Counter
    from multiprocessing import Pool
    t0 = time.time()
    workers = common.workers()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    records, jobs, case_stats = [], [], {}
    for name, spec, assemble, extra in cases:
        if spec is None:
            parts, st = extra["parts"], extra.get("stats", {})
        else:
            parts, st = el.generate_parts(spec, **extra)
        case_stats[name] = {"specification": spec.payload() if spec is not None else extra.get("payload"), **st}
        for k, (P, recipe) in enumerate(parts):
            Gs = assemble(P)
            if not isinstance(Gs, list):
                Gs = [Gs]
            for j, G in enumerate(Gs):
                rid = "%s-%s-%04d%s" % (item, name, k, "" if len(Gs) == 1 else "-%d" % j)
                records.append({"id": rid, "case": name, "part": el.graph_record(P), "recipe": recipe,
                                "graph": el.graph_record(G)})
                jobs.append((len(records) - 1, [list(e) for e in G.edges()], len(records)))
    with Pool(workers) as pool:
        for idx, res in pool.imap_unordered(_job_minor_first, jobs, chunksize=8):
            records[idx]["outcome"] = res
    counts = Counter(r["outcome"]["outcome"] for r in records)
    by_case = {}
    for r in records:
        by_case.setdefault(r["case"], Counter())[r["outcome"]["outcome"]] += 1
    found = sorted({r["outcome"]["catalogue"] for r in records
                    if r["outcome"]["outcome"] == "catalogue graph"}, key=lambda x: int(x[1:]))
    complete = counts.get("unresolved", 0) == 0 and counts.get("new obstruction", 0) == 0
    summary = {"item": item, "note": note,
               "classification order": "proper nontoroidal minor first, then enumeration.classify",
               "cases": {k: dict(v, outcomes=dict(by_case.get(k, {}))) for k, v in case_stats.items()},
               "candidates": len(records), "outcomes": dict(counts),
               "catalogue graphs found": found, "targets in the text": targets,
               "agrees with the text": complete and found == sorted(targets, key=lambda x: int(x[1:])),
               "complete": complete, "seconds": round(time.time() - t0, 1),
               "tools": {"networkx": nx.__version__,
                         "genus lower bounds": "multi_genus (package graph-genus)",
                         "embeddings": "PAGE (package graph-genus)"}}
    with open(outdir / "candidates.json", "w") as f:
        json.dump({"item": item, "candidates": records}, f)
    with open(outdir / "summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    return summary


# ----------------------------------------------------------------------------
# Structural conditions of Section 6 for a whole graph G

def twin_classes(G):
    """{neighborhood triple (sorted tuple): [degree-three vertices with it]},
    for triples with at least two such vertices."""
    by = {}
    for v in G:
        if G.degree(v) == 3:
            by.setdefault(tuple(sorted(map(str, G[v]))), []).append(str(v))
    return {t: sorted(vs) for t, vs in by.items() if len(vs) >= 2}


def adjacent_neighbors(G):
    """A degree-three vertex with two adjacent neighbors (violating Lemma
    3.3(a), lem:independent-neighborhood), as (v, x, y), or None."""
    for v in sorted(G, key=str):
        if G.degree(v) == 3:
            nb = sorted(G[v], key=str)
            for x, y in itertools.combinations(nb, 2):
                if G.has_edge(x, y):
                    return [str(v), str(x), str(y)]
    return None


def section6_row(G):
    """For a graph G all of whose twin triples meet pairwise in exactly one
    vertex: over every choice of two twin pairs on distinct twin triples,
    the earliest row of Theorem 6.1 (thm:one-common) reached, with a witness.
    Rows: 0 a component with five neighbors (Section 6.2), 1 four neighbors
    not including the common vertex (6.3), 2 four including it (6.4),
    3 every component has at most three neighbors (6.5-6.6).
    Returns (row, witness) or (None, reason) if the twin-triple condition
    fails."""
    tw = twin_classes(G)
    trip = sorted(tw)
    if len(trip) < 2:
        return None, {"reason": "fewer than two twin triples", "twin triples": trip}
    for A, B in itertools.combinations(trip, 2):
        k = len(set(A) & set(B))
        if k != 1:
            return None, {"reason": "two twin triples meet in %d vertices" % k, "triples": [A, B],
                          "twins": [tw[A], tw[B]]}
    best = (3, None)
    for A, B in itertools.permutations(trip, 2):
        if A > B:
            continue
        common = (set(A) & set(B)).pop()
        att = set(A) | set(B)
        for TA in itertools.combinations(tw[A], 2):
            for TB in itertools.combinations(tw[B], 2):
                core = att | set(TA) | set(TB)
                H = G.subgraph([v for v in G if str(v) not in core])
                for comp in nx.connected_components(H):
                    nb = {str(u) for v in comp for u in G[v]} & att
                    if len(nb) == 5:
                        row = 0
                    elif len(nb) == 4:
                        row = 1 if common not in nb else 2
                    else:
                        continue
                    if row < best[0]:
                        best = (row, {"triples": [A, B], "twins": [list(TA), list(TB)],
                                      "component": sorted(map(str, comp)), "neighbors": sorted(nb)})
    return best


def structural_certificate(G, want_row=3):
    """A checkable reason why the whole graph G cannot be an obstruction of
    the case at hand, or None: a degree-three vertex with adjacent
    neighbors (Lemma 3.3(a)); two twin triples meeting other than once (an
    earlier case, Sections 4-5); or, for want_row = 3 (Sections 6.5-6.6) or
    2 (Section 6.4), a choice of twin pairs giving an earlier row of
    Theorem 6.1."""
    adj = adjacent_neighbors(G)
    if adj:
        return {"kind": "adjacent neighbors", "vertex": adj[0], "neighbors": adj[1:]}
    row, wit = section6_row(G)
    if row is None:
        if "twins" in wit:
            return {"kind": "twin triples meet other than once", "triples": [list(t) for t in wit["triples"]],
                    "twins": wit["twins"]}
        return None
    if row < want_row:
        return {"kind": "earlier row", "row": row, "triples": [list(t) for t in wit["triples"]],
                "twins": wit["twins"], "component": wit["component"], "neighbors": wit["neighbors"]}
    return None


def planarity_tests(U, tests):
    """[{name, added edges, answer, certificate}] for each named edge set added
    to U: a plane rotation system or the edges of a Kuratowski subgraph."""
    out = []
    for name, edges in tests:
        ok, cert = el.planarity_certificate(with_edges(U, edges))
        out.append({"name": name, "added edges": [list(e) for e in edges],
                    "answer": "planar" if ok else "nonplanar", "certificate": cert})
    return out


def annotate(outdir, case_info, want_row, all_records=False):
    """Post-process candidates.json of a whole-graph item. case_info maps a
    case name to (roots, hypothesis), hypothesis a list of (name, edges) or
    None. For each record that is a catalogue graph (or every record, with
    all_records): attach a structural exclusion if there is one, and, with a
    hypothesis, the planarity of U + edges (U = the record's part, roots
    kept) for each listed edge set. Returns a summary of the catalogue
    graphs found."""
    outdir = Path(outdir)
    data = json.load(open(outdir / "candidates.json"))
    found = {}
    for r in data["candidates"]:
        o = r["outcome"]
        if not (all_records or o["outcome"] == "catalogue graph"):
            continue
        roots, hypothesis = case_info[r["case"]]
        G = nx.Graph()
        G.add_edges_from(r["graph"]["edges"])
        cert = structural_certificate(G, want_row)
        if cert:
            r["structural exclusion"] = cert
        if hypothesis:
            U = nx.Graph()
            U.add_nodes_from(r["part"]["vertices"])
            U.add_edges_from(r["part"]["edges"])
            U.add_nodes_from(roots)
            r["roots"] = list(roots)
            r["bridge tests"] = planarity_tests(U, hypothesis)
        if o["outcome"] == "catalogue graph":
            row, _ = section6_row(G)
            planar_tests = [t["name"] for t in r.get("bridge tests", []) if t["answer"] == "planar"]
            found.setdefault(o["catalogue"], []).append(
                {"id": r["id"], "row of Theorem 6.1": row,
                 "structural exclusion": cert["kind"] if cert else None,
                 "planar hypothesis tests": planar_tests})
    json.dump(data, open(outdir / "candidates.json", "w"))
    return found
