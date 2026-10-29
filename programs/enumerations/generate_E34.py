"""E34: the combinations of Construction 6.21 (Proposition 6.22).

The program's names: p1, p2, p3, ... are the vertices on abc and
q1, q2, ... those on ade. As in Sections 6.5 to 6.7, whichever two twin
pairs define Psi_3, every component of G - V(Psi_3) has three neighbours,
and the cross pairs are the matching {bd, ce} or the path {bd, cd, ce}
(equation (4)), with new triples abd, ace, or abd, acd, ace, bcd, cde.
Construction 6.21 as implemented here:

(1) At most two bridges with nonplanar triangle completion, on new
    triples. For a center x other than a, the forms are the chains of
    Proposition 3.10 whose twins meet abc and ade once each: one four-cycle
    with end lengths (0,1) or (0,2), where the end of length zero, y,
    makes {x, y} a cross pair; or two four-cycles with ends of length zero
    joined by a path of length one or two, with {x, y} and {x, z} both
    cross pairs. The middle vertex of a segment of length two is adjacent
    to x. For x = a the only form is that of Lemma 6.20. Every two twin
    neighbourhoods meet in exactly one vertex, and a bridge with x = a lies
    on a triple a y z whose pair yz is joined by no other edge or
    component, with no other edge or component joining y or z to the
    opposite pair.
(2) Optionally one vertex on each new triple and any subset of the cross
    edges of the chosen matching or path. The cross pairs realized (edges,
    and pairs inside the neighbourhood of a component with three
    neighbours) must be exactly the chosen set.
(3) k vertices on abc and l on ade with (k, l) = (2,2), (2,3), (2,4), (3,2)
    or (4,2).

The bases (k = l = 2) are compared up to isomorphism fixing a, b, c, d, e.
A base is discarded if a vertex of degree three has two adjacent
neighbours, if two twin triples meet other than once, or if some choice of
two twin pairs leaves a component with four or five neighbours; these
persist when vertices are added on abc and ade. The bases are written to
bases.json.gz. Each remaining base with each (k, l) is a candidate,
classified first by looking for a proper nontoroidal minor (multi_genus)
and then by enumeration.classify, and each catalogue graph found is placed
in its row of Theorem 6.1. The proof finds T7, T33, T35, T49, T70.

Run as: python generate_E34.py FOLDER
"""
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import networkx as nx

import common
import enumeration as el

import psi3_runs as ext

SIDE = {"b": 0, "c": 0, "d": 1, "e": 1}
AVAIL = {"matching": [("b", "d"), ("c", "e")], "path": [("b", "d"), ("c", "d"), ("c", "e")]}
NEW = {"matching": ["abd", "ace"], "path": ["abd", "acd", "ace", "bcd", "cde"]}
MULT = [(2, 2), (2, 3), (2, 4), (3, 2), (4, 2)]


def cross(x, y):
    return x in SIDE and y in SIDE and SIDE[x] != SIDE[y]


def pairs_in(t):
    return {tuple(sorted((x, y))) for x, y in itertools.combinations(t, 2) if cross(x, y)}


def forms_on(t):
    """(name, edges with interior names starting with '#', twin triples,
    x == a) for every form on the triple t."""
    out = []
    for x in t:
        if x == "a":
            y, z = [v for v in t if v != "a"]
            E = [(y, "#s"), ("#s", "#p"), ("#p", "#t"), ("#t", "#q"), ("#q", "#s"), ("#t", z),
                 ("a", "#p"), ("a", "#q")]
            out.append(("A(%s)" % t, E, [{"a", "#s", "#t"}], True))
            continue
        others = [v for v in t if v != x]
        for y in others:
            if not cross(x, y):
                continue
            z = [v for v in others if v != y][0]
            cyc = [(y, "#p"), ("#p", "#t"), ("#t", "#q"), ("#q", y), (x, "#p"), (x, "#q")]
            out.append(("C01(x=%s,y=%s,z=%s)" % (x, y, z), cyc + [("#t", z)], [{x, y, "#t"}], False))
            out.append(("C02(x=%s,y=%s,z=%s)" % (x, y, z), cyc + [("#t", "#m"), ("#m", z), ("#m", x)],
                        [{x, y, "#t"}], False))
        if all(cross(x, v) for v in others):
            y, z = sorted(others)
            c1 = [(y, "#p"), ("#p", "#t"), ("#t", "#q"), ("#q", y), (x, "#p"), (x, "#q")]
            c2 = [("#s", "#p2"), ("#p2", z), (z, "#q2"), ("#q2", "#s"), (x, "#p2"), (x, "#q2")]
            out.append(("CC1(x=%s)" % x, c1 + [("#t", "#s")] + c2, [{x, y, "#t"}, {x, "#s", z}], False))
            out.append(("CC2(x=%s)" % x, c1 + [("#t", "#m"), ("#m", "#s"), ("#m", x)] + c2,
                        [{x, y, "#t"}, {x, "#s", z}], False))
    return out


def rename(edges, triples, pre):
    f = lambda v: pre + v[1:] if v.startswith("#") else v
    return [(f(u), f(v)) for u, v in edges], [{f(v) for v in T} for T in triples]


def meets_once(triples):
    return all(len(A & B) == 1 for A, B in itertools.combinations(triples, 2))


def bases(kind):
    """Every labeled base (k = l = 2) of Construction 6.21 for the matching
    or the path, with the restrictions applied before the structural ones."""
    avail = {tuple(sorted(p)) for p in AVAIL[kind]}
    forms = [(t,) + f for t in NEW[kind] for f in forms_on(t)]
    stats = Counter()
    for nb in range(3):
        for combo in itertools.combinations_with_replacement(range(len(forms)), nb):
            fs = [forms[i] for i in combo]
            fedges, ftrip, ftri = [], [{"a", "b", "c"}, {"a", "d", "e"}], []
            for j, (t, name, E, T, isa) in enumerate(fs):
                E2, T2 = rename(E, T, "B%d_" % j)
                fedges += E2
                ftrip += T2
                ftri.append(t)
            if not meets_once(ftrip):
                stats["form twins meet other than once"] += 1
                continue
            stats["form choices"] += 1
            for singles in itertools.product((0, 1), repeat=len(NEW[kind])):
                strip = [t for t, s in zip(NEW[kind], singles) if s]
                for ce in itertools.product((0, 1), repeat=len(AVAIL[kind])):
                    cedges = [p for p, s in zip(AVAIL[kind], ce) if s]
                    realized = {tuple(sorted(p)) for p in cedges}
                    for t in ftri + strip:
                        realized |= pairs_in(t)
                    if realized != avail:
                        stats["realized cross pairs differ"] += 1
                        continue
                    bad = False
                    for j, (t, name, E, T, isa) in enumerate(fs):
                        if not isa:
                            continue
                        y, z = [v for v in t if v != "a"]
                        other = {tuple(sorted(p)) for p in cedges}
                        for jj, t2 in enumerate(ftri):
                            if jj != j:
                                other |= pairs_in(t2)
                        for t2 in strip:
                            other |= pairs_in(t2)
                        if any(y in p or z in p for p in other) or t in strip or \
                                any(jj != j and t2 == t for jj, t2 in enumerate(ftri)):
                            bad = True
                    if bad:
                        stats["x = a bridge not isolated"] += 1
                        continue
                    G = ext.psi3()
                    G.add_edges_from(fedges)
                    for t in strip:
                        ext.add_vertex(G, "v_" + t, t)
                    G.add_edges_from(cedges)
                    stats["constructions"] += 1
                    yield {"kind": kind, "forms": [f[1] for f in fs], "form triples": ftri,
                           "singletons": strip, "cross edges": ["".join(p) for p in cedges]}, G
    bases.stats[kind] = dict(stats)


bases.stats = {}


def with_mult(G, k, l):
    H = G.copy()
    for i in range(3, k + 1):
        ext.add_vertex(H, "p%d" % i, "abc")
    for i in range(3, l + 1):
        ext.add_vertex(H, "q%d" % i, "ade")
    return H


if __name__ == "__main__":
    t0 = time.time()
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    roots = list("abcde")
    reps, buckets = [], {}
    for kind in ("matching", "path"):
        for desc, G in bases(kind):
            key = el.rooted_key(G, roots)
            dup = None
            for j in buckets.get(key, []):
                if el.rooted_isomorphic(G, reps[j][1], roots):
                    dup = j
                    break
            if dup is None:
                buckets.setdefault(key, []).append(len(reps))
                reps.append((desc, G, [desc]))
            else:
                reps[dup][2].append(desc)
    base_records, kept = [], []
    for i, (desc, G, all_desc) in enumerate(reps):
        cert = ext.structural_certificate(G, want_row=3)
        rec = {"id": "E34-base-%04d" % i, "description": desc, "presentations": len(all_desc),
               "graph": el.graph_record(G)}
        if cert:
            rec["structural exclusion"] = cert
        else:
            kept.append((i, G))
        base_records.append(rec)
    from multiprocessing import Pool
    records, jobs = [], []
    for i, G in kept:
        for k, l in MULT:
            H = with_mult(G, k, l)
            records.append({"id": "E34-%04d-k%dl%d" % (i, k, l), "case": "k=%d,l=%d" % (k, l),
                            "base": "E34-base-%04d" % i, "graph": el.graph_record(H)})
            jobs.append((len(records) - 1, [list(e) for e in H.edges()], len(records)))
    with Pool(common.workers()) as pool:
        for idx, res in pool.imap_unordered(ext._job_minor_first, jobs, chunksize=4):
            records[idx]["outcome"] = res
    found = {}
    for r in records:
        if r["outcome"]["outcome"] == "catalogue graph":
            H = nx.Graph()
            H.add_edges_from(r["graph"]["edges"])
            row, _ = ext.section6_row(H)
            cert = ext.structural_certificate(H, want_row=3)
            if cert:
                r["structural exclusion"] = cert
            found.setdefault(r["outcome"]["catalogue"], []).append(
                {"id": r["id"], "row of Theorem 6.1": row, "structural exclusion": cert["kind"] if cert else None})
    counts = Counter(r["outcome"]["outcome"] for r in records)
    struct = Counter(r["structural exclusion"]["kind"] for r in base_records if "structural exclusion" in r)
    names = sorted(found, key=lambda x: int(x[1:]))
    complete = counts.get("unresolved", 0) == 0 and counts.get("new obstruction", 0) == 0
    summary = {"item": "E34", "where": "Proposition 6.22",
               "generation": bases.stats, "labeled constructions (k=l=2)": sum(
                   v.get("constructions", 0) for v in bases.stats.values()),
               "bases up to isomorphism fixing a..e": len(reps),
               "bases excluded by structural conditions": dict(struct),
               "bases kept": len(kept), "multiplicities": MULT,
               "candidates": len(records), "outcomes": dict(counts),
               "catalogue graphs found": names, "catalogue graphs and their rows of Theorem 6.1": found,
               "targets in the text": ["T7", "T33", "T35", "T49", "T70"],
               "agrees with the text": complete and names == ["T7", "T33", "T35", "T49", "T70"]
               and all(x["row of Theorem 6.1"] == 3 for xs in found.values() for x in xs),
               "complete": complete, "seconds": round(time.time() - t0, 1)}
    json.dump({"item": "E34", "candidates": records}, open(out / "candidates.json", "w"))
    import gzip
    with gzip.open(out / "bases.json.gz", "wt") as f:
        json.dump({"item": "E34", "candidates": base_records}, f, separators=(",", ":"))
    json.dump(summary, open(out / "summary.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "catalogue graphs and their rows of Theorem 6.1"},
                     indent=1))
