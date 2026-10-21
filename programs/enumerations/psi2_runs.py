"""Extensions of enumeration.py for the enumerations of Section 5 that
classify whole graphs (E16-E25).

1. Permitted original edges between attachments (Construction 3.16(i)):
   when two attachments x, y are consecutive on a branch path of K', are
   not joined by an added edge, and xy is a permitted original edge, the
   segment is taken both as the original edge xy (kept in the part) and
   with one interior vertex. generate_parts_ext does this and otherwise
   follows enumeration.generate_parts exactly (roots up to the 72 symmetries
   of K3,3, outside vertices in every position, segments shortened as in
   Lemma 3.15(ii), interior pieces joined by every forest, one edge at each
   attachment with no interior neighbor).

2. Structural exclusions recorded beside a whole-graph outcome. Every
   candidate gets an outcome of enumeration.classify, and a candidate that
   violates a structural condition of G at that point of the text also gets
   a field "structural exclusion" with a certificate.

3. run_whole_item: like enumeration.run_item, with the structural checks.
"""

import itertools
import json
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import networkx as nx

import common
import enumeration as el


# ----------------------------------------------------------------------------
# Construction 3.16 with permitted original edges between attachments

def build_retained_variants(pl, spec, originals=()):
    """Like enumeration.build_retained, but a segment between consecutive
    attachments x, y with xy in originals (and no added edge xy) is taken
    both as the original edge and with one interior vertex. Returns a list
    of (K', part, original edges used by K')."""
    originals = {frozenset(e) for e in originals}
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
    chains = []
    for k, (u, v) in enumerate(el.BRANCH_EDGES):
        chains.append([label[u]] + [names[i] for i in paths[k]] + [label[v]])
    # validity of segments at outside vertices does not depend on the choices
    choice_segments = []
    for chain in chains:
        for x, y in zip(chain, chain[1:]):
            xm, ym = x in names, y in names
            if (x in spec.outside) or (y in spec.outside):
                if not (xm and ym) or not spec.c_adjacent(x, y):
                    return []
            elif x in spec.roots and y in spec.roots and not spec.c_adjacent(x, y) \
                    and frozenset((x, y)) in originals:
                choice_segments.append(frozenset((x, y)))
    results = []
    for choice in itertools.product((False, True), repeat=len(choice_segments)):
        direct = {s for s, c in zip(choice_segments, choice) if c}
        K = nx.Graph()
        used_c, used_orig = set(), set()
        cnt = count
        for chain in chains:
            for x, y in zip(chain, chain[1:]):
                if (x in spec.outside) or (y in spec.outside):
                    K.add_edge(x, y)
                    used_c.add(frozenset((x, y)))
                elif x in spec.roots and y in spec.roots:
                    if spec.c_adjacent(x, y):
                        K.add_edge(x, y)
                        used_c.add(frozenset((x, y)))
                    elif frozenset((x, y)) in direct:
                        K.add_edge(x, y)
                        used_orig.add(frozenset((x, y)))
                    else:
                        mid = "i%d" % cnt
                        cnt += 1
                        K.add_edge(x, mid)
                        K.add_edge(mid, y)
                else:
                    K.add_edge(x, y)
        if not spec.required <= {frozenset(e) for e in K.edges()}:
            continue
        part = K.copy()
        part.remove_edges_from([tuple(e) for e in used_c])
        part.remove_nodes_from([y for y in spec.outside if y in part])
        part.add_nodes_from(spec.roots)
        results.append((K, part, used_orig))
    return results


def generate_parts_ext(spec, originals=(), min_interior=1, max_interior=None):
    """enumeration.generate_parts with permitted original edges between
    attachments (see build_retained_variants). Each part is (P, recipe);
    the original edges on K' are edges of P."""
    r, y = len(spec.roots), len(spec.outside)
    reps, raw = el.placement_classes(r)
    buckets, parts = {}, []
    stats = {"root placements": raw, "root placement classes": len(reps),
             "full placements": 0, "valid retained subdivisions": 0, "recipes": 0}
    for pl0, orbit in reps:
        for pl in el.extend_placement(pl0, r, y):
            stats["full placements"] += 1
            for K, part, used_orig in build_retained_variants(pl, spec, originals):
                stats["valid retained subdivisions"] += 1
                for P in el.completions(part, spec.roots):
                    stats["recipes"] += 1
                    n_in = P.number_of_nodes() - r
                    if n_in < min_interior or (max_interior and n_in > max_interior):
                        continue
                    key = el.rooted_key(P, spec.roots)
                    found = False
                    for j in buckets.get(key, []):
                        if el.rooted_isomorphic(P, parts[j][0], spec.roots):
                            found = True
                            break
                    if not found:
                        buckets.setdefault(key, []).append(len(parts))
                        parts.append((el.relabel_interior(P, spec.roots),
                                      {"placement": [list(pl[0]), [list(p) for p in pl[1]]],
                                       "K'": sorted(sorted(map(str, e)) for e in K.edges()),
                                       "original edges on K'": sorted(sorted(e) for e in used_orig)}))
    stats["parts up to rooted isomorphism"] = len(parts)
    return parts, stats


# ----------------------------------------------------------------------------
# Structural checks on whole graphs

def nonplanar_subgraph_check(name, keep_vertices, drop_edges_within=()):
    """A check: the subgraph of G on keep_vertices(G), minus the edges with
    both ends in drop_edges_within, is nonplanar. Returns a function giving
    a certificate (Kuratowski edges) or None."""
    drop = set(drop_edges_within)

    def check(G):
        vs = keep_vertices(G)
        H = G.subgraph(vs).copy()
        H.remove_edges_from([(u, v) for u, v in list(H.edges()) if u in drop and v in drop])
        ok, cert = el.planarity_certificate(H)
        if ok:
            return None
        return {"condition": name, "kind": "nonplanar subgraph",
                "kuratowski": cert}
    return check


def _job(job):
    idx, edges, nodes, seed = job
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    return idx, el.classify(G, seed=seed)


def run_whole_item(item, cases, targets, outdir, structural=(),
                   structural_payload=None, note="", extra_summary=None, fast=False,
                   slim=False):
    """cases: list of (name, parts, stats, payload, assemble) where parts is
    a list of (P, recipe) and assemble(P) returns a whole graph or a list of
    (label, graph). Every graph is classified (enumeration.classify), and the
    structural checks (list of functions G -> certificate or None) are
    applied; a candidate is a survivor if it is a catalogue graph, a new
    obstruction or unresolved, and passes every structural check."""
    t0 = time.time()
    workers = common.workers()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    records, jobs, case_stats = [], [], {}
    for name, parts, st, payload, assemble in cases:
        case_stats[name] = {"specification": payload, **st}
        for k, (P, recipe) in enumerate(parts):
            Gs = assemble(P)
            if not isinstance(Gs, list):
                Gs = [("", Gs)]
            for lab, G in Gs:
                rid = "%s-%s-%04d%s" % (item, name, k, ("-" + lab) if lab else "")
                rec = {"id": rid, "case": name, "variant": lab, "part": el.graph_record(P),
                       "recipe": recipe, "graph": el.graph_record(G)}
                if slim:
                    # the part is the graph minus the fixed exterior; keep only the placement
                    del rec["part"]
                    rec["recipe"] = {"placement": recipe.get("placement")}
                for chk in structural:
                    cert = chk(G)
                    if cert is not None:
                        rec["structural exclusion"] = cert
                        break
                records.append(rec)
                if fast:
                    jobs.append(G)
                else:
                    jobs.append((len(records) - 1, [list(e) for e in G.edges()], list(G.nodes()), len(records)))
    if fast:
        for rec, res in zip(records, fast_outcomes(jobs, workers=workers)):
            rec["outcome"] = res
    else:
        with Pool(workers) as pool:
            for idx, res in pool.imap_unordered(_job, jobs, chunksize=8):
                records[idx]["outcome"] = res
    return finish_whole_item(item, records, case_stats, targets, outdir,
                             structural_payload, note, t0, extra_summary)


SURVIVOR_KINDS = ("catalogue graph", "new obstruction", "unresolved")


def finish_whole_item(item, records, case_stats, targets, outdir,
                      structural_payload, note, t0, extra_summary=None):
    def final(r):
        if "structural exclusion" in r:
            return "excluded: " + r["structural exclusion"]["condition"]
        return r["outcome"]["outcome"]
    counts = Counter(final(r) for r in records)
    raw = Counter(r["outcome"]["outcome"] for r in records)
    by_case = {}
    for r in records:
        by_case.setdefault(r["case"], Counter())[final(r)] += 1
    surv = [r for r in records if "structural exclusion" not in r
            and r["outcome"]["outcome"] in SURVIVOR_KINDS]
    found = sorted({r["outcome"]["catalogue"] for r in surv
                    if r["outcome"]["outcome"] == "catalogue graph"}, key=lambda x: int(x[1:]))
    open_ = [r["id"] for r in records if r["outcome"]["outcome"] in ("unresolved", "new obstruction")
             and "structural exclusion" not in r]
    excluded_catalogue = Counter(r["outcome"]["catalogue"] for r in records
                                 if "structural exclusion" in r and r["outcome"]["outcome"] == "catalogue graph")
    complete = not open_
    tg = sorted(targets, key=lambda x: int(x[1:]))
    summary = {"item": item, "note": note,
               "cases": {k: dict(v, outcomes=dict(by_case.get(k, {}))) for k, v in case_stats.items()},
               "structural conditions": structural_payload,
               "candidates": len(records), "outcomes": dict(counts),
               "whole-graph outcomes before structural exclusions": dict(raw),
               "catalogue graphs among structurally excluded candidates": dict(excluded_catalogue),
               "survivors": [{"id": r["id"], "catalogue": r["outcome"].get("catalogue"),
                              "outcome": r["outcome"]["outcome"]} for r in surv],
               "catalogue graphs found": found, "targets in the text": targets,
               "agrees with the text": complete and found == tg,
               "complete": complete, "open candidates": open_,
               "seconds": round(time.time() - t0, 1),
               "tools": {"networkx": nx.__version__,
                         "genus lower bounds": "multi_genus (package graph-genus)",
                         "embeddings": "PAGE (package graph-genus)"}}
    if extra_summary:
        summary.update(extra_summary)
    with open(outdir / "candidates.json", "w") as f:
        json.dump({"item": item, "structural conditions": structural_payload,
                   "candidates": records}, f)
    with open(outdir / "summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    return summary


# ----------------------------------------------------------------------------
# Several interior components, further singleton components, and minimal
# sets of attaching edges meeting connectivity requirements
# (Construction 3.16(ii) with h > 0, Lemma r2:retained-components)

def set_partitions(items):
    items = list(items)
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for p in set_partitions(rest):
        yield [[first]] + p
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]


def group_forests(part, roots, groups):
    """For a partition of the interior pieces into groups (future
    components), every choice of a tree joining the pieces of each group,
    with every choice of ends. Yields the list of forest edges."""
    per_group = []
    for grp in groups:
        opts = []
        for tree in el.trees_on(len(grp)):
            ends = [list(itertools.product(grp[a], grp[b])) for a, b in tree]
            for chosen in itertools.product(*ends):
                opts.append(list(chosen))
        per_group.append(opts)
    for combo in itertools.product(*per_group):
        yield [e for es in combo for e in es]


def components_attach(H, roots):
    """Interior components of H with the set of roots each is adjacent to.

    The components are listed in a fixed order: when the interior vertices
    are fewer than half of the vertices of H, by their largest vertex name,
    from the last name to the first; otherwise in the order in which the
    vertices of H first meet them. (networkx lists them, in the first case,
    in the order of a Python set of the interior vertices, which can change
    from one Python process to another. Any fixed order gives a valid
    enumeration; this one gives the data folder.)"""
    inner = [v for v in H if v not in roots]
    comps = list(nx.connected_components(H.subgraph(inner)))
    if 2 * len(inner) < len(H):
        comps.sort(key=lambda c: max(map(str, c)), reverse=True)
    else:
        pos = {v: i for i, v in enumerate(inner)}
        comps.sort(key=lambda c: min(pos[v] for v in c))
    out = []
    for comp in comps:
        att = {u for v in comp for u in H[v] if u in roots}
        out.append((sorted(comp, key=str), att))
    return out


def minimal_group_additions(groups_att, single_atts, requirement, candidates):
    """groups_att: list of root sets (one per existing interior component);
    single_atts: root sets of the further singleton components;
    requirement(list of root sets) -> bool, monotone in adding roots;
    candidates: list of (group index, root) pairs that may be added.
    Returns every inclusion-minimal set of candidate pairs whose addition
    satisfies the requirement (with the singletons present)."""
    found = []
    cands = [(g, r) for g, r in candidates if r not in groups_att[g]]
    for k in range(len(cands) + 1):
        for S in itertools.combinations(range(len(cands)), k):
            Sset = set(S)
            if any(f <= Sset for f in found):
                continue
            atts = [set(a) for a in groups_att]
            for i in S:
                g, r = cands[i]
                atts[g].add(r)
            if requirement(atts + [set(s) for s in single_atts]):
                found.append(Sset)
    return [[cands[i] for i in sorted(f)] for f in found]


def globally_minimal(groups_att, added, single_atts, requirement):
    """Check that no proper subset of the additions (added pairs, singleton
    edges, whole singletons) satisfies the requirement. By monotonicity
    for a fixed set of singletons, single removals suffice (see E17)."""
    def atts_with(drop_pair=None, drop_single=None, drop_single_root=None):
        atts = [set(a) for a in groups_att]
        for p in added:
            if p != drop_pair:
                atts[p[0]].add(p[1])
        ss = []
        for i, s in enumerate(single_atts):
            if i == drop_single:
                continue
            s = set(s)
            if drop_single_root and drop_single_root[0] == i:
                s.discard(drop_single_root[1])
            ss.append(s)
        return atts + ss
    for p in added:
        if requirement(atts_with(drop_pair=p)):
            return False
    for i, s in enumerate(single_atts):
        if requirement(atts_with(drop_single=i)):
            return False
        for r in s:
            if len(s) > 1 and requirement(atts_with(drop_single_root=(i, r))):
                return False
    return True


# ----------------------------------------------------------------------------
# Twin structure (Sections 5 and 6)

def twin_triples(G):
    """Neighborhoods shared by at least two vertices of degree three."""
    by = {}
    for v in G:
        if G.degree(v) == 3:
            by.setdefault(frozenset(G[v]), []).append(v)
    return {N: sorted(vs, key=str) for N, vs in by.items() if len(vs) >= 2}


def section6_violation(G, forbid_five=False):
    """Hypotheses of Section 6 on the twin structure of G: any two distinct
    twin triples share exactly one vertex (Theorem thm:one-common); with
    forbid_five, also no two twin pairs on triples sharing one vertex leave
    a component of G - V(Psi3') adjacent to all five of their vertices (the
    hypothesis of Sections 6.3-6.5, after Section 6.2). Returns a witness
    or None."""
    tt = twin_triples(G)
    trip = sorted(tt, key=lambda s: sorted(map(str, s)))
    for A, B in itertools.combinations(trip, 2):
        k = len(A & B)
        if k == 0:
            return {"condition": "disjoint twin triples", "kind": "twin structure",
                    "triples": [sorted(A), sorted(B)], "twins": [tt[A][:2], tt[B][:2]]}
        if k == 2:
            return {"condition": "two twin triples share two vertices", "kind": "twin structure",
                    "triples": [sorted(A), sorted(B)], "twins": [tt[A][:2], tt[B][:2]]}
    if forbid_five:
        for A, B in itertools.combinations(trip, 2):
            att = A | B
            psi = set(att) | set(tt[A][:2]) | set(tt[B][:2])
            H = G.subgraph([v for v in G if v not in psi])
            for comp in nx.connected_components(H):
                if {u for v in comp for u in G[v]} & att == att:
                    return {"condition": "component adjacent to all five attachments",
                            "kind": "twin structure", "triples": [sorted(A), sorted(B)],
                            "twins": [tt[A][:2], tt[B][:2]], "component": sorted(comp, key=str)}
    return None


SECTION6_STRUCT = {
    "disjoint twin triples": {"kind": "twin structure",
                              "meaning": "Section 6 assumes any two twin triples share exactly one vertex"},
    "two twin triples share two vertices": {"kind": "twin structure",
                                            "meaning": "such graphs belong to Section 5"},
    "component adjacent to all five attachments": {"kind": "twin structure",
                                                   "meaning": "such graphs belong to Section 6.2"}}


# ----------------------------------------------------------------------------
# Faster classification for large items: batch genus tests first

def fast_outcomes(graphs, workers=None, seed0=0):
    """Outcomes with the same certificates as enumeration.classify. First
    multi_genus decides which graphs have genus at least two; for those,
    every edge deletion and contraction is tested in one batch, and the
    first one of genus at least two gives the outcome "proper nontoroidal
    minor". Every other graph goes through enumeration.classify."""
    workers = workers or common.workers()
    over = el.genus_above_one(graphs)
    out = [None] * len(graphs)
    chunk = 400
    idx_nt = [i for i, o in enumerate(over) if o]
    for s in range(0, len(idx_nt), chunk):
        block = idx_nt[s:s + chunk]
        minors, owner = [], []
        for i in block:
            G = graphs[i]
            es = sorted(G.edges(), key=el.edge_key)
            for e in es:
                minors.append(("delete", e, el.gr.delete_edge(G, e)))
                owner.append(i)
            for e in es:
                minors.append(("contract", e, el.gr.contract_edge(G, e)))
                owner.append(i)
        flags = el.genus_above_one([H for _, _, H in minors])
        for (op, e, H), i, fl in zip(minors, owner, flags):
            if fl and out[i] is None:
                out[i] = {"outcome": "proper nontoroidal minor", "operation": op,
                          "edge": [str(x) for x in el.edge_key(e)]}
    rest = [i for i in range(len(graphs)) if out[i] is None]
    jobs = [(i, [list(e) for e in graphs[i].edges()], list(graphs[i].nodes()), seed0 + i) for i in rest]
    with Pool(workers) as pool:
        for i, res in pool.imap_unordered(_job, jobs, chunksize=4):
            out[i] = res
    return out
