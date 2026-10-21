"""Construction 3.16 with original edges between attachments and
disconnected interiors (used by E11 and E14).

retained(roots, added, required, original)
    Construction 3.16(i): every placement of the roots on a K3,3
    subdivision (branch vertex, inside a branch path in order, or off), up
    to the symmetry of K3,3; segments shortened as in Lemma 3.15(ii). A
    segment between two roots is the added edge when its ends are joined by
    one, is taken both as the original edge and with one interior vertex
    when its ends form a permitted original pair, and otherwise has one
    interior vertex. Yields (K', part), part = K' minus the added edges,
    with every root present.

grouped(part, roots)
    The interior pieces (components of part - roots) are divided among the
    components of the interior in every way (every set partition), and the
    pieces of each component are joined by a spanning tree with every
    choice of ends (Lemma 5.8: the vertices of a
    component off the subdivision are contracted onto a forest joining its
    retained pieces). Yields (graph, list of components).

local_minimal(forms, roots, signature, valid, allowed_pair)
    Keeps the forms with no one-step proper minor (deleting an edge, or
    contracting an edge with an interior end without creating a forbidden
    edge between roots) that is valid and has the same signature. Every form
    with no proper minor of the same signature is kept, so the result
    contains every minimal form.
"""

import itertools

import networkx as nx

import enumeration as el


def retained(roots, added=(), required=(), original=()):
    added = {frozenset(e) for e in added}
    required = {frozenset(e) for e in required}
    original = {frozenset(e) for e in original}
    reps, raw = el.placement_classes(len(roots))
    for pl, orbit in reps:
        branches, paths = pl
        label, count = [], 0
        for v in range(6):
            if branches[v] >= 0:
                label.append(roots[branches[v]])
            else:
                label.append("i%d" % count)
                count += 1
        pairs = []
        for k, (u, v) in enumerate(el.BRANCH_EDGES):
            chain = [label[u]] + [roots[i] for i in paths[k]] + [label[v]]
            pairs += list(zip(chain, chain[1:]))
        choices = []
        for x, y in pairs:
            if x in roots and y in roots:
                key = frozenset((x, y))
                if key in added:
                    choices.append(("added",))
                elif key in original:
                    choices.append(("original", "subdivided"))
                else:
                    choices.append(("subdivided",))
            else:
                choices.append(("edge",))
        for pick in itertools.product(*choices):
            K = nx.Graph()
            K.add_nodes_from(label)
            used, mid = set(), count
            for (x, y), c in zip(pairs, pick):
                if c == "subdivided":
                    m = "i%d" % mid
                    mid += 1
                    K.add_edge(x, m)
                    K.add_edge(m, y)
                else:
                    K.add_edge(x, y)
                    if c == "added":
                        used.add(frozenset((x, y)))
            if not required <= used:
                continue
            part = K.copy()
            part.remove_edges_from([tuple(e) for e in used])
            part.add_nodes_from(roots)
            yield K, part, {"placement": [list(pl[0]), [list(p) for p in pl[1]]], "segments": list(pick)}


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for p in set_partitions(rest):
        yield [[first]] + p
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]


def grouped(part, roots):
    inner = [v for v in part if v not in roots]
    pieces = [sorted(c) for c in nx.connected_components(part.subgraph(inner))]
    pieces.sort()
    for groups in set_partitions(list(range(len(pieces)))):
        options = []
        for g in groups:
            trees = []
            for tree in el.trees_on(len(g)):
                ends = [list(itertools.product(pieces[g[a]], pieces[g[b]])) for a, b in tree]
                for edges in itertools.product(*ends):
                    trees.append(edges)
            options.append(trees)
        for pick in itertools.product(*options):
            H = part.copy()
            for edges in pick:
                H.add_edges_from(edges)
            comps = [sorted(v for i in g for v in pieces[i]) for g in groups]
            yield H, comps


# The order of the interior components when the interior vertices are fewer
# than half of the vertices (see interior_components).
COMPONENT_RANK = {"i1": 0, "i0": 1, "i2": 2}


def interior_components(H, roots):
    """The components of H minus the roots (sorted vertex lists), in a fixed
    order: when the interior vertices are fewer than half of the vertices of
    H, by the rank of their vertices in COMPONENT_RANK (i1, i0, i2, then the
    other names in order); otherwise in the order in which the vertices of H
    first meet them. (networkx lists them, in the first case, in the order of
    a Python set of the interior vertices, which can change from one Python
    process to another. Any fixed order gives a valid enumeration; this one
    gives the data folder.)"""
    inner = [v for v in H if v not in roots]
    comps = list(nx.connected_components(H.subgraph(inner)))
    if 2 * len(inner) < len(H):
        rank = lambda v: (COMPONENT_RANK.get(v, len(COMPONENT_RANK)), str(v))
        comps.sort(key=lambda c: min(rank(v) for v in c))
    else:
        pos = {v: i for i, v in enumerate(inner)}
        comps.sort(key=lambda c: min(pos[v] for v in c))
    return [sorted(c) for c in comps]


def one_step_minors(H, roots, allowed_pair):
    """Every graph from H by deleting an edge, or contracting an edge with
    an interior end into its other end, when no new edge joins two roots
    unless allowed_pair says it may."""
    rs = set(roots)
    for e in list(H.edges()):
        M = H.copy()
        M.remove_edge(*e)
        yield ("delete", e), M
        for x, y in (e, e[::-1]):
            if x in rs:
                continue                      # contract the interior end x into y
            if y not in rs and x > y:
                continue                      # interior-interior: once
            new = [z for z in H[x] if z != y]
            if y in rs and any(z in rs and not allowed_pair(y, z) for z in new):
                continue
            M = nx.contracted_nodes(H, y, x, self_loops=False)
            M = nx.Graph(M)
            yield ("contract", (x, y)), M


def local_minimal(forms, roots, signature, valid, allowed_pair):
    kept = []
    for H, info in forms:
        sig = signature(H)
        smaller = None
        for op, M in one_step_minors(H, roots, allowed_pair):
            if valid(M) and signature(M) == sig:
                smaller = op
                break
        if smaller is None:
            kept.append((H, info))
    return kept


def dedupe(forms, roots):
    buckets, out = {}, []
    for H, info in forms:
        key = el.rooted_key(H, roots)
        if any(el.rooted_isomorphic(H, out[j][0], roots) for j in buckets.get(key, [])):
            continue
        buckets.setdefault(key, []).append(len(out))
        out.append((H, info))
    return out
