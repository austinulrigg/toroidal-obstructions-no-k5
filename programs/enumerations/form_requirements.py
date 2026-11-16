"""Construction 3.18, the completion step: further single-vertex components
and every inclusion-minimal set of attaching edges (and, where permitted,
original edges between attachments) that meets the requirements of the case,
with the whole configuration minimal (Lemma 5.11).

A configuration consists of the retained interior components (fixed vertex
sets, each with the set of attachments it is already adjacent to), further
single vertices (each given by its set of attachments), and the chosen
additions. Additions are
    ("edge", g, r)   an attaching edge from some vertex of retained component g
                     to attachment r (every choice of the vertex is made later);
    ("orig", m, r)   an original edge between attachments m and r.
A requirement is a monotone predicate on the configuration. Every minimal set
of additions is found by branching: at each step one violated condition names
a set of additions of which every satisfying superset must contain one.
"""
import itertools

import networkx as nx


class Config:
    def __init__(self, comp_atts, singles, origs=()):
        self.comp_atts = [set(a) for a in comp_atts]     # retained components
        self.singles = [set(s) for s in singles]         # further single vertices
        self.origs = set(origs)                          # original edges (m, r)

    def with_additions(self, adds):
        c = Config(self.comp_atts, self.singles, self.origs)
        for a in adds:
            if a[0] == "edge":
                c.comp_atts[a[1]].add(a[2])
            else:
                c.origs.add((a[1], a[2]))
        return c

    def all_atts(self):
        return self.comp_atts + self.singles


def member_connectivity(cfg, members):
    """The graph on the members and all interior components (retained and
    single), with an edge from a component to each member it is adjacent to.
    Returns (connected?, the node set of the part containing members[0])."""
    G = nx.Graph()
    G.add_nodes_from(("m", x) for x in members)
    for i, a in enumerate(cfg.all_atts()):
        G.add_node(("c", i))
        for x in members:
            if x in a:
                G.add_edge(("c", i), ("m", x))
    part = nx.node_connected_component(G, ("m", members[0]))
    return len(part) == G.number_of_nodes(), part


class Requirement:
    """conditions:
        per_component: list of root sets; every interior component (retained or
            single) must be adjacent to at least one root of each set;
        members: the attachments that, with the interior, must form one
            connected component (or None);
        global_adjacent: list of (root, allow_orig, retained_only): some
            interior vertex (of a retained component when retained_only) must
            be adjacent to root, or, when allow_orig, some member must be
            joined to root by an original edge;
        forbid: roots no interior vertex may be adjacent to."""

    def __init__(self, per_component=(), members=None, global_adjacent=(), forbid=()):
        self.per_component = [set(s) for s in per_component]
        self.members = tuple(members) if members else None
        self.global_adjacent = list(global_adjacent)
        self.forbid = set(forbid)

    def violated(self, cfg, n_retained):
        """None if satisfied, else the list of additions one of which every
        satisfying superset must contain ([] if no addition can help)."""
        for i, a in enumerate(cfg.all_atts()):
            if a & self.forbid:
                return []
            for S in self.per_component:
                if not (a & S):
                    if i < n_retained:
                        return [("edge", i, r) for r in sorted(S)]
                    return []
        for root, allow_orig, retained_only in self.global_adjacent:
            atts = cfg.comp_atts if retained_only else cfg.all_atts()
            ok = any(root in a for a in atts) or (allow_orig and any(r == root for _, r in cfg.origs))
            if not ok:
                adds = [("edge", i, root) for i in range(n_retained)]
                if allow_orig:
                    adds += [("orig", m, root) for m in self.members]
                return adds
        if self.members:
            ok, part = member_connectivity(cfg, self.members)
            if not ok:
                adds = []
                for i in range(n_retained):
                    for x in self.members:
                        if (("c", i) in part) != (("m", x) in part) and x not in cfg.comp_atts[i]:
                            adds.append(("edge", i, x))
                return adds
        return None


def minimal_addition_sets(base, req, n_retained, allowed):
    """Every inclusion-minimal set of additions (from allowed, a set) for
    which req holds."""
    found = set()

    def rec(chosen):
        cfg = base.with_additions(chosen)
        v = req.violated(cfg, n_retained)
        if v is None:
            found.add(frozenset(chosen))
            return
        for a in v:
            if a in allowed and a not in chosen:
                rec(chosen | {a})
    rec(frozenset())
    out = [s for s in found if not any(t < s for t in found)]
    return sorted((sorted(s) for s in out), key=lambda s: (len(s), s))


def globally_minimal(base, adds, req, n_retained):
    """No addition, further single vertex, or edge of a further single vertex
    can be dropped with req still holding."""
    full = base.with_additions(adds)
    for a in adds:
        rest = [b for b in adds if b != a]
        if req.violated(base.with_additions(rest), n_retained) is None:
            return False
    for i, s in enumerate(full.singles):
        c = Config(full.comp_atts, full.singles[:i] + full.singles[i + 1:], full.origs)
        if req.violated(c, n_retained) is None:
            return False
        for r in s:
            if len(s) > 1:
                s2 = set(s) - {r}
                c = Config(full.comp_atts, full.singles[:i] + [s2] + full.singles[i + 1:], full.origs)
                if req.violated(c, n_retained) is None:
                    return False
    return True


def completions(H, comps, roots, req, single_options, max_singles, edge_targets, orig_options=()):
    """For a part H (a networkx graph containing the attachments roots) whose
    interior components are comps (vertex lists), yield every completed graph:
    further single vertices s0, s1, ... with neighbourhoods from
    single_options (up to max_singles of them), and every minimal set of
    additions, attaching edges going to edge_targets from every choice of
    vertex in the component. orig_options lists permitted original edges
    (m, r) between attachments that may be added."""
    comp_atts = [{r for v in comp for r in H[v] if r in roots} for comp in comps]
    n = len(comps)
    allowed = {("edge", i, r) for i in range(n) for r in edge_targets if r not in comp_atts[i]}
    allowed |= {("orig", m, r) for m, r in orig_options if not H.has_edge(m, r)}
    origs0 = {(m, r) for m, r in orig_options if H.has_edge(m, r)}
    seen = set()
    for k in range(max_singles + 1):
        for singles in itertools.combinations_with_replacement(sorted(single_options, key=sorted), k):
            base = Config(comp_atts, singles, origs0)
            for adds in minimal_addition_sets(base, req, n, allowed):
                if not globally_minimal(base, adds, req, n):
                    continue
                edge_adds = [a for a in adds if a[0] == "edge"]
                orig_adds = [a for a in adds if a[0] == "orig"]
                for verts in itertools.product(*[comps[a[1]] for a in edge_adds]):
                    M = H.copy()
                    M.add_edges_from((v, a[2]) for v, a in zip(verts, edge_adds))
                    M.add_edges_from((a[1], a[2]) for a in orig_adds)
                    for j, s in enumerate(singles):
                        M.add_edges_from(("s%d" % j, r) for r in sorted(s))
                    key = (frozenset(frozenset(map(str, e)) for e in M.edges()))
                    if key in seen:
                        continue
                    seen.add(key)
                    yield M, {"further single vertices": [sorted(s) for s in singles],
                              "added attaching edges": [[str(v), a[2]] for v, a in zip(verts, edge_adds)],
                              "added original edges": [[a[1], a[2]] for a in orig_adds]}
