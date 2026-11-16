"""Shared tools for the enumerations built with Construction 3.18.

Construction 3.18 (common:marked-construction) lists candidates for a part P
of G with labeled attachments R, as in Lemma 3.17 (core:retained-subdivision):
a K3,3 subdivision K' of P + C, where C is a planar graph meeting P in R
(edges within R, and outside vertices Y whose neighbors are all in R), with
every segment shortened as in Lemma 3.17(ii), then the interior pieces joined
by forests and each attachment with no interior neighbor given one edge
(Lemma 3.17(iii)).

This module generates such parts, compares them up to label-preserving
isomorphism, and classifies whole graphs G built from them. Classification
answers whether G is an obstruction and records a certificate:

  K5 minor          five branch sets (checked by is_k5_model)
  cut               a vertex set of size at most two separating G
  catalogue graph   an isomorphism onto one of T1, ..., T74
  toroidal          a rotation system of genus at most one (from PAGE)
  proper minor      an edge e with G - e or G / e of genus at least two
                    (multi_genus)
  new obstruction   none of these: G has genus at least two, every G - e and
                    G / e has genus at most one, and a Wagner decomposition
                    shows G has no K5 minor. This would contradict the paper.
"""

import hashlib
import itertools
import json
from pathlib import Path

import networkx as nx

import common
import graphs as gr

BRANCH_EDGES = tuple((a, b) for a in range(3) for b in range(3, 6))


# ----------------------------------------------------------------------------
# K3,3 placements

def k33_automorphisms():
    for left in itertools.permutations(range(3)):
        for right in itertools.permutations(range(3)):
            for swap in (False, True):
                if not swap:
                    yield left + tuple(v + 3 for v in right)
                else:
                    yield tuple(v + 3 for v in left) + right


AUTS = tuple(k33_automorphisms())


def raw_placements(m):
    """Every placement of m marked items on K3,3: each at a branch vertex, at
    a position inside one of the nine branch paths (in order), or off.
    A placement is (branches, paths): branches[v] is the item at branch
    vertex v or -1; paths[k] is the tuple of items inside branch path k, read
    from its left end (0, 1 or 2) to its right end (3, 4 or 5)."""
    branches = [-1] * 6
    paths = [[] for _ in range(9)]

    def extend(item):
        if item == m:
            yield (tuple(branches), tuple(tuple(p) for p in paths))
            return
        yield from extend(item + 1)
        for v in range(6):
            if branches[v] == -1:
                branches[v] = item
                yield from extend(item + 1)
                branches[v] = -1
        for path in paths:
            for pos in range(len(path) + 1):
                path.insert(pos, item)
                yield from extend(item + 1)
                path.pop(pos)

    yield from extend(0)


def transform(pl, perm):
    branches, paths = pl
    nb = [-1] * 6
    np_ = [()] * 9
    for src, dst in enumerate(perm):
        nb[dst] = branches[src]
    for k, (a, b) in enumerate(BRANCH_EDGES):
        x, y = perm[a], perm[b]
        seq = paths[k]
        if x >= 3:
            x, y = y, x
            seq = seq[::-1]
        np_[x * 3 + y - 3] = seq
    return (tuple(nb), tuple(np_))


_class_cache = {}


def placement_classes(m):
    """One placement from each orbit of the automorphism group of K3,3 acting
    on the raw placements of m items, with the orbit sizes."""
    if m in _class_cache:
        return _class_cache[m]
    seen, reps = set(), []
    total = 0
    for pl in raw_placements(m):
        total += 1
        if pl in seen:
            continue
        orbit = {transform(pl, a) for a in AUTS}
        seen |= orbit
        reps.append((min(orbit), len(orbit)))
    assert total == len(seen)
    reps.sort()
    _class_cache[m] = (reps, total)
    return reps, total


def extend_placement(pl, base, extra):
    """All placements of items base, ..., base + extra - 1 added to pl, in
    every free branch vertex or path position (no 'off', no symmetry)."""
    def rec(pl, item):
        if item == base + extra:
            yield pl
            return
        branches, paths = pl
        for v in range(6):
            if branches[v] == -1:
                nb = list(branches)
                nb[v] = item
                yield from rec((tuple(nb), paths), item + 1)
        for k in range(9):
            for pos in range(len(paths[k]) + 1):
                np_ = list(paths)
                np_[k] = paths[k][:pos] + (item,) + paths[k][pos:]
                yield from rec((branches, tuple(np_)), item + 1)
    yield from rec(pl, base)


# ----------------------------------------------------------------------------
# Building the retained subdivision and the part

class Spec:
    """Parameters of one application of Construction 3.18.

    roots      labels of R
    outside    {label: tuple of neighbors in R} for the vertices of C - R
    c_edges    edges of C within R
    required   edges that K' must use (pairs of labels)
    forbid_root_edges  original edges within R are never in P (P is a bridge)
    """

    def __init__(self, roots, outside=None, c_edges=(), required=(), name=""):
        self.roots = tuple(roots)
        self.outside = dict(outside or {})
        self.c_edges = {frozenset(e) for e in c_edges}
        self.required = {frozenset(e) for e in required}
        self.name = name
        self.marked = self.roots + tuple(self.outside)
        for y, nb in self.outside.items():
            assert set(nb) <= set(self.roots)

    def payload(self):
        return {"roots": list(self.roots),
                "outside": {k: list(v) for k, v in self.outside.items()},
                "C edges within R": sorted(sorted(e) for e in self.c_edges),
                "required edges of K'": sorted(sorted(e) for e in self.required)}

    def c_adjacent(self, x, y):
        if x in self.outside:
            return y in self.outside[x]
        if y in self.outside:
            return x in self.outside[y]
        return frozenset((x, y)) in self.c_edges


def build_retained(pl, spec):
    """K' for the placement pl of spec.marked, or None if a segment at an
    outside vertex is not an edge of C or a required edge is not used.
    Returns (Kp, part) where part = K' - C - Y with all roots kept."""
    branches, paths = pl
    names = spec.marked
    label = []
    count = 0
    for v in range(6):
        if branches[v] >= 0:
            label.append(names[branches[v]])
        else:
            label.append("i%d" % count)
            count += 1
    K = nx.Graph()
    used_c = set()
    for k, (u, v) in enumerate(BRANCH_EDGES):
        chain = [label[u]] + [names[i] for i in paths[k]] + [label[v]]
        for x, y in zip(chain, chain[1:]):
            xm, ym = x in names, y in names
            if (x in spec.outside) or (y in spec.outside):
                if not (xm and ym) or not spec.c_adjacent(x, y):
                    return None
                K.add_edge(x, y)
                used_c.add(frozenset((x, y)))
            elif x in spec.roots and y in spec.roots:
                if spec.c_adjacent(x, y):
                    K.add_edge(x, y)
                    used_c.add(frozenset((x, y)))
                else:
                    mid = "i%d" % count
                    count += 1
                    K.add_edge(x, mid)
                    K.add_edge(mid, y)
            else:
                K.add_edge(x, y)
    if not spec.required <= {frozenset(e) for e in K.edges()}:
        return None
    part = K.copy()
    part.remove_edges_from([tuple(e) for e in used_c])
    part.remove_nodes_from([y for y in spec.outside if y in part])
    part.add_nodes_from(spec.roots)
    return K, part


def interior_pieces(part, roots):
    inner = [v for v in part if v not in roots]
    return [sorted(c) for c in sorted(nx.connected_components(part.subgraph(inner)), key=lambda c: sorted(c))]


def trees_on(k):
    if k == 1:
        yield ()
    elif k == 2:
        yield ((0, 1),)
    else:
        for seq in itertools.product(range(k), repeat=k - 2):
            T = nx.from_prufer_sequence(list(seq))
            yield tuple(sorted(tuple(sorted(e)) for e in T.edges()))


def completions(part, roots):
    """Lemma 3.17(iii): join the interior pieces by a forest of interior
    edges, with every choice of ends, and give each root with no interior
    neighbor one edge to an interior vertex, in every way."""
    pieces = interior_pieces(part, roots)
    inner = sorted(v for v in part if v not in roots)
    lonely = [r for r in roots if not any(u not in roots for u in part[r])]
    for tree in trees_on(len(pieces)):
        ends = [list(itertools.product(pieces[a], pieces[b])) for a, b in tree]
        for forest in itertools.product(*ends):
            for nbrs in itertools.product(inner, repeat=len(lonely)):
                P = part.copy()
                P.add_edges_from(forest)
                P.add_edges_from(zip(lonely, nbrs))
                yield P


def rooted_key(P, roots):
    """A label-preserving invariant: Weisfeiler-Lehman hash with each root
    colored by its own label and interior vertices uncolored."""
    H = P.copy()
    for v in H:
        H.nodes[v]["c"] = v if v in roots else "*"
    return nx.weisfeiler_lehman_graph_hash(H, node_attr="c", iterations=4)


def rooted_isomorphic(P, Q, roots):
    if P.number_of_nodes() != Q.number_of_nodes() or P.number_of_edges() != Q.number_of_edges():
        return False
    nm = lambda a, b: a["c"] == b["c"]
    A, B = P.copy(), Q.copy()
    for H in (A, B):
        for v in H:
            H.nodes[v]["c"] = v if v in roots else "*"
    return nx.is_isomorphic(A, B, node_match=nm)


def relabel_interior(P, roots):
    """Interior vertices renamed u0, u1, ... in sorted order."""
    inner = sorted((v for v in P if v not in roots), key=str)
    return nx.relabel_nodes(P, {v: "u%d" % i for i, v in enumerate(inner)})


def generate_parts(spec, min_interior=1, max_interior=None, root_classes=None):
    """All parts P from Construction 3.18 for spec, up to label-preserving
    isomorphism. Roots are placed up to the symmetry of K3,3, outside
    vertices in every position. Returns (parts, stats); each part is
    (P, first recipe)."""
    r, y = len(spec.roots), len(spec.outside)
    reps, raw = placement_classes(r)
    buckets, parts = {}, []
    stats = {"root placements": raw, "root placement classes": len(reps),
             "full placements": 0, "valid retained subdivisions": 0, "recipes": 0}
    for pl0, orbit in reps:
        for pl in extend_placement(pl0, r, y):
            stats["full placements"] += 1
            # outside vertices must lie on K'
            res = build_retained(pl, spec)
            if res is None:
                continue
            K, part = res
            stats["valid retained subdivisions"] += 1
            for P in completions(part, spec.roots):
                stats["recipes"] += 1
                n_in = P.number_of_nodes() - r
                if n_in < min_interior or (max_interior and n_in > max_interior):
                    continue
                key = rooted_key(P, spec.roots)
                found = False
                for j in buckets.get(key, []):
                    if rooted_isomorphic(P, parts[j][0], spec.roots):
                        found = True
                        break
                if not found:
                    buckets.setdefault(key, []).append(len(parts))
                    parts.append((relabel_interior(P, spec.roots),
                                  {"placement": [list(pl[0]), [list(p) for p in pl[1]]],
                                   "K'": sorted(sorted(map(str, e)) for e in K.edges())}))
    stats["parts up to rooted isomorphism"] = len(parts)
    return parts, stats


# ----------------------------------------------------------------------------
# Classifying whole graphs

_catalog = None


def catalog():
    global _catalog
    if _catalog is None:
        _catalog = [(r["name"], gr.from_graph6(r["graph6"])) for r in common.catalogue_records()]
    return _catalog


def identify(G):
    for name, T in catalog():
        if T.number_of_nodes() == G.number_of_nodes() and T.number_of_edges() == G.number_of_edges():
            m = gr.isomorphism(G, T)
            if m is not None:
                return name, {str(k): v for k, v in m.items()}
    return None, None


def emb_record(parts):
    return [{"vertices": [str(v) for v in p["vertices"]], "genus": p["genus"],
             "rotation": {str(k): [str(x) for x in v] for k, v in p["rotation"].items()}}
            for p in parts]


def small_cut(G):
    if not nx.is_connected(G):
        return []
    for v in G:
        H = G.copy(); H.remove_node(v)
        if not nx.is_connected(H):
            return [v]
    for u, v in itertools.combinations(sorted(G, key=str), 2):
        H = G.copy(); H.remove_nodes_from([u, v])
        if not nx.is_connected(H):
            return [u, v]
    return None


def genus_above_one(graphs):
    """For each graph, whether its genus is at least two (multi_genus, with
    the genus added over connected components)."""
    import multigenus as mg
    jobs, owners = [], []
    for i, G in enumerate(graphs):
        for comp in nx.connected_components(G):
            H = G.subgraph(comp)
            if H.number_of_edges() >= 9:
                jobs.append(H)
                owners.append(i)
    over1 = mg.nontoroidal(jobs) if jobs else []
    over0 = mg.nonplanar(jobs) if jobs else []
    result = [False] * len(graphs)
    nonplanar_count = [0] * len(graphs)
    for i, a, b in zip(owners, over1, over0):
        if a:
            result[i] = True
        if b:
            nonplanar_count[i] += 1
    return [r or c >= 2 for r, c in zip(result, nonplanar_count)]


def edge_key(e):
    return tuple(sorted(map(str, e)))


def classify(G, seed=0, k5_trials=2000, need_three_connected=True):
    """Classify the graph G; see the module docstring. Genus lower bounds
    come from multi_genus, embeddings from PAGE."""
    out = {}
    sets = gr.find_k5_minor(G, trials=k5_trials, seed=seed)
    if sets:
        out.update(outcome="K5 minor", branch_sets=[[str(v) for v in s] for s in sets])
        return out
    if need_three_connected:
        cut = small_cut(G)
        if cut is not None:
            out.update(outcome="cut of order at most two", cut=[str(v) for v in cut])
            return out
    name, m = identify(G)
    if name is not None:
        out.update(outcome="catalogue graph", catalogue=name, map=m)
        return out
    if not genus_above_one([G])[0]:
        g, parts = gr.page_embedding(G)
        assert g <= 1
        out.update(outcome="toroidal", embedding=emb_record(parts))
        return out
    es = sorted(G.edges(), key=edge_key)
    minors = [("delete", e, gr.delete_edge(G, e)) for e in es] + \
             [("contract", e, gr.contract_edge(G, e)) for e in es]
    flags = genus_above_one([H for _, _, H in minors])
    for (op, e, H), f in zip(minors, flags):
        if f:
            out.update(outcome="proper nontoroidal minor", operation=op, edge=[str(x) for x in edge_key(e)])
            return out
    w = gr.wagner_decomposition(G)
    if w is None:
        sets = gr.find_k5_minor(G, trials=0, deletion_trials=300000, seed=seed)
        if sets:
            out.update(outcome="K5 minor", branch_sets=[[str(v) for v in s] for s in sets])
            return out
        out.update(outcome="unresolved")
        return out
    out.update(outcome="new obstruction", wagner=w)
    return out


def graph_record(G):
    return {"vertices": sorted(map(str, G)), "edges": sorted(sorted(map(str, e)) for e in G.edges())}


def sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


# ----------------------------------------------------------------------------
# Running an item

def _classify_job(job):
    """job = (index, edges, seed) or (index, edges, seed, vertices): the
    graph is built from the vertices in the given order (when given) and
    then the edges; the random search for a K5 minor depends on that order."""
    idx, edges, seed = job[:3]
    G = nx.Graph()
    if len(job) > 3:
        G.add_nodes_from(job[3])
    G.add_edges_from(edges)
    return idx, classify(G, seed=seed)


def run_item(item, cases, targets, outdir, note="", seed_of=None, nodes_first=None):
    """cases: list of (name, spec, assemble, extra) where assemble(P) returns
    the whole graph G for the part P (or a list of graphs). Generates every
    part, classifies every assembled graph, writes candidates.json and
    summary.json in outdir, and returns the summary. targets: the graphs T_i
    the text says the item leaves. The random search for a K5 minor in the
    n-th graph (counting from 1) uses the seed n, or seed_of(case, k, n) for
    the k-th part of the case when seed_of is given. When nodes_first(case)
    is true, the graphs of the case are handed to the classification with
    their vertices in the order of the assembled graph."""
    import time
    from collections import Counter
    from multiprocessing import Pool
    workers = common.workers()
    t0 = time.time()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    records, jobs = [], []
    case_stats = {}
    for name, spec, assemble, extra in cases:
        if spec is None:
            parts, st = extra["parts"], extra.get("stats", {})
        else:
            parts, st = generate_parts(spec, **extra)
        case_stats[name] = {"specification": spec.payload() if spec is not None else extra.get("payload"),
                            **st}
        for k, (P, recipe) in enumerate(parts):
            Gs = assemble(P)
            if not isinstance(Gs, list):
                Gs = [Gs]
            for j, G in enumerate(Gs):
                rid = "%s-%s-%04d%s" % (item, name, k, "" if len(Gs) == 1 else "-%d" % j)
                records.append({"id": rid, "case": name, "part": graph_record(P), "recipe": recipe,
                                "graph": graph_record(G)})
                seed = len(records) if seed_of is None else seed_of(name, k, len(records))
                job = (len(records) - 1, [list(e) for e in G.edges()], seed)
                if nodes_first is not None and nodes_first(name):
                    job += ([str(v) for v in G.nodes()],)
                jobs.append(job)
    # Each verdict is appended to progress.jsonl as soon as it is found, so an
    # interrupted run continues where it stopped when it is started again.
    progress = outdir / "progress.jsonl"
    keys = {job[0]: hashlib.sha256(json.dumps([sorted(sorted(map(str, e)) for e in job[1])] + list(job[3:]) if
                                              len(job) > 3 else sorted(sorted(map(str, e)) for e in job[1]))
                                   .encode()).hexdigest() for job in jobs}
    done = {}
    if progress.exists():
        with open(progress) as pf:
            for line in pf:
                try:
                    r = json.loads(line)
                except ValueError:
                    continue
                done[r["index"]] = r
    todo = []
    for job in jobs:
        r = done.get(job[0])
        if r is not None and r["key"] == keys[job[0]] and r["seed"] == job[2]:
            records[job[0]]["outcome"] = r["outcome"]
        else:
            todo.append(job)
    if done:
        print("resuming: %d of %d verdicts already found" % (len(jobs) - len(todo), len(jobs)), flush=True)
    with open(progress, "a") as pf, Pool(workers) as pool:
        for n, (idx, res) in enumerate(pool.imap_unordered(_classify_job, todo, chunksize=8), 1):
            records[idx]["outcome"] = res
            pf.write(json.dumps({"index": idx, "key": keys[idx], "seed": jobs[idx][2], "outcome": res}) + "\n")
            pf.flush()
            if n % 500 == 0:
                print("classified %d of %d" % (n + len(jobs) - len(todo), len(jobs)), flush=True)
    counts = Counter(r["outcome"]["outcome"] for r in records)
    by_case = {}
    for r in records:
        by_case.setdefault(r["case"], Counter())[r["outcome"]["outcome"]] += 1
    found = sorted({r["outcome"]["catalogue"] for r in records
                    if r["outcome"]["outcome"] == "catalogue graph"}, key=lambda x: int(x[1:]))
    found_by_case = {}
    for r in records:
        if r["outcome"]["outcome"] == "catalogue graph":
            found_by_case.setdefault(r["case"], set()).add(r["outcome"]["catalogue"])
    found_by_case = {k: sorted(v, key=lambda x: int(x[1:])) for k, v in found_by_case.items()}
    complete = counts.get("unresolved", 0) == 0 and counts.get("new obstruction", 0) == 0
    summary = {"item": item, "note": note,
               "cases": {k: dict(v, outcomes=dict(by_case.get(k, {}))) for k, v in case_stats.items()},
               "candidates": len(records), "outcomes": dict(counts),
               "catalogue graphs found": found, "catalogue graphs by case": found_by_case,
               "targets in the text": targets,
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
# Bridge-level tests: planarity with certificates, allowed operations

def planarity_certificate(G):
    """(True, rotation) for a planar graph, with a plane rotation system;
    (False, edges of a Kuratowski subgraph) otherwise."""
    ok, cert = nx.check_planarity(G, counterexample=True)
    if ok:
        return True, {str(v): [str(u) for u in cert.neighbors_cw_order(v)] for v in cert}
    return False, sorted(sorted(map(str, e)) for e in cert.edges())


def allowed_operations(B, roots):
    """Every one-step allowed operation on the bridge B with attachments
    roots (Section 3.1): delete an interior vertex or an edge, or contract
    an edge with an interior end and delete the new edges between
    attachments. The result must have a nonempty connected interior and an
    interior neighbor at every attachment. Yields (description, B')."""
    roots = set(roots)

    def valid(H):
        inner = [v for v in H if v not in roots]
        if not inner or not nx.is_connected(H.subgraph(inner)):
            return False
        return all(any(u not in roots for u in H[r]) for r in roots if r in H) and roots <= set(H)

    for v in sorted((v for v in B if v not in roots), key=str):
        H = B.copy(); H.remove_node(v)
        if valid(H):
            yield ("delete vertex", [str(v)]), H
    for e in sorted(B.edges(), key=edge_key):
        H = B.copy(); H.remove_edge(*e)
        if valid(H):
            yield ("delete edge", list(edge_key(e))), H
    for e in sorted(B.edges(), key=edge_key):
        u, v = e
        if u in roots and v in roots:
            continue
        if v in roots:
            u, v = v, u          # contract the interior vertex v into u
        H = nx.contracted_nodes(B, u, v, self_loops=False, copy=True)
        if u in roots:
            H.remove_edges_from([(u, w) for w in list(H[u]) if w in roots])
        H = nx.Graph(H)
        if valid(H):
            yield ("contract edge", [str(u), str(v)]), H


def disjoint_paths(B, roots, pairs):
    """Vertex-disjoint paths joining each pair, with internal vertices
    outside roots, or None. Brute force over simple paths (small B)."""
    inner = [v for v in B if v not in roots]

    def paths(s, t, banned):
        H = B.subgraph([v for v in inner if v not in banned] + [s, t])
        if s in H and t in H:
            for p in nx.all_simple_paths(H, s, t):
                if all(x not in roots for x in p[1:-1]):
                    yield p

    (s1, t1), (s2, t2) = pairs
    for p in paths(s1, t1, set()):
        for q in paths(s2, t2, set(p)):
            return [p, q]
    return None
