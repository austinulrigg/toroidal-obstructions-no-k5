"""Searches for minors that return a checked witness or nothing (E1, E2, E3).

find_k5_minor tries random contractions with a fixed seed, an exhaustive
search over small components and, if asked, a SAT formulation with lazy
connectivity cuts (PySAT). known_minor looks for a model of a graph of the
catalogue after a bounded number of contractions and edge deletions. None
always means that no witness was returned, never that none exists. Every
returned model is checked (validate_model) before it is returned.
"""

from collections import Counter, deque
from copy import deepcopy
import itertools as it
import random

import networkx as nx



class InvalidCertificate(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise InvalidCertificate(message)


def checked_graph(vertices, edges):
    """Adjacency sets of a simple graph with string vertex names."""
    g = {str(v): set() for v in vertices}
    for u, v in edges:
        u, v = str(u), str(v)
        require(u != v and u in g and v in g, "bad edge")
        g[u].add(v)
        g[v].add(u)
    return g


def _connected(g, vs):
    vs = set(vs)
    if not vs:
        return False
    start = next(iter(vs))
    seen, stack = {start}, [start]
    while stack:
        x = stack.pop()
        for y in g[x]:
            if y in vs and y not in seen:
                seen.add(y)
                stack.append(y)
    return seen == vs


def _edges(g):
    return {tuple(sorted((u, v))) for u in g for v in g[u]}


def minor_model(g, h, raw_sets):
    """Raise unless raw_sets (vertex of h -> branch set in g) is a minor model."""
    require(isinstance(raw_sets, dict), "minor model must map target vertices to branch sets")
    sets = {str(v): list(map(str, group)) for v, group in raw_sets.items()}
    require(set(sets) == set(h), "minor model target vertex mismatch")
    all_vertices = [v for group in sets.values() for v in group]
    require(len(all_vertices) == len(set(all_vertices)), "minor branch sets overlap or repeat")
    require(set(all_vertices) <= set(g), "minor branch set uses unknown vertex")
    for v, group in sets.items():
        require(_connected(g, group), "empty or disconnected branch set for " + v)
    for u, v in _edges(h):
        require(any(g[x] & set(sets[v]) for x in sets[u]), "missing target adjacency")


def isomorphism(g, h, raw_mapping):
    """Raise unless raw_mapping is an isomorphism from g onto h."""
    mapping = {str(v): str(w) for v, w in raw_mapping.items()}
    require(set(mapping) == set(g), "isomorphism source mismatch")
    require(len(set(mapping.values())) == len(h) and set(mapping.values()) == set(h),
            "isomorphism is not a bijection")
    image_edges = {tuple(sorted((mapping[u], mapping[v]))) for u, v in _edges(g)}
    require(image_edges == _edges(h), "isomorphism does not preserve exactly the edges")


def _graph(G):
    """Copy a simple undirected graph with stable, JSON-safe string labels."""
    require(not G.is_directed() and not G.is_multigraph(), "need a simple undirected graph")
    labels = {v: str(v) for v in G}
    require(len(set(labels.values())) == len(labels), "vertex labels collide as strings")
    require(nx.number_of_selfloops(G) == 0, "self-loops are unsupported")
    H = nx.Graph()
    H.add_nodes_from(sorted(labels.values()))
    H.add_edges_from(sorted(tuple(sorted((labels[u], labels[v]))) for u, v in G.edges()))
    return H


def _adj(G):
    return checked_graph(G.nodes(), G.edges())


def _k5():
    return nx.complete_graph([str(i) for i in range(5)])


def validate_model(G, H, model):
    """Raise on an invalid model; otherwise return a canonical branch mapping."""
    G, H = _graph(G), _graph(H)
    if isinstance(model, list):
        require(set(H) == set(map(str, range(len(model)))), "list model needs indexed target")
        model = {str(i): bag for i, bag in enumerate(model)}
    bags = {str(v): sorted(map(str, bag)) for v, bag in model.items()}
    minor_model(_adj(G), _adj(H), bags)
    return bags


def graph_invariant(G):
    """Cheap isomorphism invariant, used only to reject impossible matches."""
    return (len(G), G.number_of_edges(), tuple(sorted(d for _, d in G.degree())),
            tuple(sorted(len(c) for c in nx.connected_components(G))),
            tuple(sorted(nx.triangles(G).values())))


def first_isomorphism(G, H):
    """The isomorphism from G onto H whose list of images of the vertices of G,
    taken in the order of their names, comes first; None if there is none. A
    search may find the isomorphisms in an order that changes from run to
    run; this choice does not depend on it."""
    order = sorted(G)
    best = None
    for iso in nx.isomorphism.GraphMatcher(G, H).isomorphisms_iter():
        key = [iso[v] for v in order]
        if best is None or key < best:
            best = key
    return None if best is None else dict(zip(order, best))


def canonical_labels(H, model):
    """The model obtained from model (vertex of H -> branch set) by the
    automorphism of H that makes the list of labels of the model's vertices,
    taken in the order of their names, come first."""
    owner = {x: h for h, bag in model.items() for x in bag}
    order = sorted(owner)
    best = None
    for aut in nx.isomorphism.GraphMatcher(H, H).isomorphisms_iter():
        key = [int(aut[owner[x]]) for x in order]
        if best is None or key < best[0]:
            best = (key, aut)
    aut = best[1]
    return {aut[h]: bag for h, bag in model.items()}


def transport_model(G, H, model, mapping=None, target=None):
    """Transport a checked model from G to isomorphic H; return None if no iso.

    A list of five bags defaults to target K5. Other targets must be supplied.
    An explicit mapping is checked as an exact source-to-destination iso.
    """
    G, H = _graph(G), _graph(H)
    if target is None:
        require(isinstance(model, list) and len(model) == 5, "target graph required")
        target = _k5()
    target = _graph(target)
    source = validate_model(G, target, model)
    if mapping is None:
        if graph_invariant(G) != graph_invariant(H):
            return None
        mapping = first_isomorphism(G, H)
        if mapping is None:
            return None
    mapping = {str(a): str(b) for a, b in mapping.items()}
    isomorphism(_adj(G), _adj(H), mapping)
    result = validate_model(H, target, {v: [mapping[x] for x in bag]
                                        for v, bag in source.items()})
    if isinstance(model, list):
        return [result[str(i)] for i in range(len(result))]
    return result


class PositiveModelCache:
    """Cache only positive certificates, grouped by global graph invariants."""

    def __init__(self):
        self.buckets = {}

    def add(self, G, target, model, metadata=None):
        G, target = _graph(G), _graph(target)
        bags = validate_model(G, target, model)
        self.buckets.setdefault(graph_invariant(G), []).append(
            (G, target, bags, deepcopy(metadata or {})))

    def find(self, G):
        G = _graph(G)
        for source, target, bags, metadata in self.buckets.get(graph_invariant(G), []):
            mapped = transport_model(source, G, bags, target=target)
            if mapped is not None:
                return {"branch_sets": mapped, "metadata": deepcopy(metadata)}
        return None


def _connected_masks(G):
    vertices = sorted(G)
    index = {v: i for i, v in enumerate(vertices)}
    adjacency = [sum(1 << index[w] for w in G[v]) for v in vertices]
    cache = {0: False}

    def connected(mask):
        if mask not in cache:
            reached = 0
            todo = mask & -mask
            while todo:
                reached |= todo
                frontier = 0
                pending = todo
                while pending:
                    bit = pending & -pending
                    pending -= bit
                    frontier |= adjacency[bit.bit_length() - 1]
                todo = frontier & mask & ~reached
            cache[mask] = reached == mask
        return cache[mask]

    return vertices, adjacency, connected


def exact_k5_minor(G, max_order=10, stats=None):
    """Exhaust partitions on each small connected component, returning a witness.

    In a connected graph any K5 model extends to a partition into five
    connected bags. Restricted-growth partitions avoid color permutations.
    Components larger than max_order are skipped, not certified negative.
    """
    G = _graph(G)
    stats = stats if stats is not None else {}
    stats.update(method="exact_partitions", partitions_tested=0, skipped_large_components=0)
    for component in sorted(nx.connected_components(G), key=lambda c: tuple(sorted(c))):
        if len(component) < 5:
            continue
        if len(component) > max_order:
            stats["skipped_large_components"] += 1
            continue
        Q = G.subgraph(sorted(component))
        vertices, adj, connected = _connected_masks(Q)
        n, bags = len(Q), [0] * 5

        def visit(i, count):
            if n - i < 5 - count:
                return None
            if i == n:
                if count != 5:
                    return None
                stats["partitions_tested"] += 1
                if not all(connected(bag) for bag in bags):
                    return None
                neighborhoods = []
                for bag in bags:
                    neighborhood, pending = 0, bag
                    while pending:
                        bit = pending & -pending
                        pending -= bit
                        neighborhood |= adj[bit.bit_length() - 1]
                    neighborhoods.append(neighborhood)
                if all(neighborhoods[a] & bags[b] for a, b in it.combinations(range(5), 2)):
                    return [[v for j, v in enumerate(vertices) if bag >> j & 1] for bag in bags]
                return None
            for color in range(min(count + 1, 5)):
                bags[color] |= 1 << i
                found = visit(i + 1, max(count, color + 1))
                bags[color] ^= 1 << i
                if found is not None:
                    return found
            return None

        found = visit(0, 0)
        if found is not None:
            checked = validate_model(G, _k5(), found)
            return [checked[str(i)] for i in range(5)]
    return None


def random_k5_minor(G, trials=256, seed=0, stats=None):
    """Uniform current-edge contractions, in stable order, with a local RNG."""
    G = _graph(G)
    require(isinstance(seed, int), "seed must be an integer")
    require(isinstance(trials, int) and trials >= 0, "invalid trial budget")
    stats = stats if stats is not None else {}
    stats.update(method="random_contractions", seed=seed, trials_requested=trials, trials_run=0)
    rng = random.Random(seed)
    for component in sorted(nx.connected_components(G), key=lambda c: tuple(sorted(c))):
        if len(component) < 5:
            continue
        names = sorted(component)
        index = {v: i for i, v in enumerate(names)}
        initial = {index[v]: set(index[w] for w in G[v]) for v in names}
        # Choosing from sorted edges makes seeds independent of insertion order.
        for _ in range(trials):
            stats["trials_run"] += 1
            adj = {v: set(ws) for v, ws in initial.items()}
            bags = {v: [names[v]] for v in adj}
            while len(adj) > 5:
                choices = [(u, v) for u in sorted(adj) for v in sorted(adj[u]) if u < v]
                if not choices:
                    break
                u, v = choices[rng.randrange(len(choices))]
                neighbors = (adj[u] | adj[v]) - {u, v}
                for w in neighbors:
                    adj[w].discard(v)
                    adj[w].add(u)
                adj[u] = neighbors
                del adj[v]
                bags[u].extend(bags.pop(v))
            if len(adj) == 5 and all(len(ws) == 4 for ws in adj.values()):
                found = [sorted(bags[v]) for v in sorted(bags)]
                checked = validate_model(G, _k5(), found)
                return [checked[str(i)] for i in range(5)]
    return None


def sat_k5_minor(G, conflict_budget=20000, max_rounds=200, solver_name="g3", stats=None):
    """Optional exact SAT formulation, bounded search, positive output only.

    Every vertex of one connected component is assigned to one of five bags.
    Edge witnesses enforce pairwise bag adjacency. Connectivity is enforced
    lazily: a selected vertex on each side of a disconnected bag requires a
    selected vertex in the component's external neighborhood. Restricted
    growth orders bag minima and removes the 5! color symmetry. The full
    formulation is exact, but no solver-negative result is exported as proof.
    conflict_budget is per SAT call; max_rounds bounds connectivity refinement.
    """
    from pysat.solvers import Solver

    require(conflict_budget > 0 and max_rounds > 0, "SAT budgets must be positive")
    G = _graph(G)
    stats = stats if stats is not None else {}
    stats.update(method="sat_connectivity_cuts", solver=solver_name,
                 conflict_budget_per_call=conflict_budget, max_rounds=max_rounds,
                 solver_calls=0, connectivity_cuts=0, component_results=[])
    for component in sorted(nx.connected_components(G), key=lambda c: tuple(sorted(c))):
        if len(component) < 5:
            continue
        Q = G.subgraph(sorted(component))
        vertices, n = sorted(Q), len(Q)
        index = {v: i for i, v in enumerate(vertices)}
        x = lambda i, color: 1 + 5 * i + color
        top, clauses = 5 * n, []
        for i in range(n):
            clauses.append([x(i, c) for c in range(5)])
            clauses.extend([-x(i, a), -x(i, b)] for a, b in it.combinations(range(5), 2))
            for c in range(1, 5):
                clauses.append([-x(i, c)] + [x(j, c - 1) for j in range(i)])
        for c in range(5):
            clauses.append([x(i, c) for i in range(n)])
        for a, b in it.combinations(range(5), 2):
            witnesses = []
            for u, v in sorted(tuple(sorted(e)) for e in Q.edges()):
                for p, q in ((index[u], index[v]), (index[v], index[u])):
                    top += 1
                    witnesses.append(top)
                    clauses.extend(([-top, x(p, a)], [-top, x(q, b)]))
            clauses.append(witnesses)
        result = "round_budget_exhausted"
        with Solver(name=solver_name, bootstrap_with=clauses) as solver:
            for _ in range(max_rounds):
                solver.conf_budget(conflict_budget)
                answer = solver.solve_limited()
                stats["solver_calls"] += 1
                if answer is not True:
                    result = "solver_no_model" if answer is False else "conflict_budget_exhausted"
                    break
                assignment = set(solver.get_model())
                bags = [[v for i, v in enumerate(vertices) if x(i, c) in assignment]
                        for c in range(5)]
                cuts = []
                for c, bag in enumerate(bags):
                    pieces = sorted(nx.connected_components(Q.subgraph(bag)),
                                    key=lambda p: tuple(sorted(p)))
                    if len(pieces) <= 1:
                        continue
                    for piece in pieces:
                        first = min(piece)
                        other = min(set(bag) - piece)
                        boundary = set().union(*(set(Q[v]) for v in piece)) - piece
                        cuts.append([-x(index[first], c), -x(index[other], c)] +
                                    [x(index[v], c) for v in sorted(boundary)])
                if not cuts:
                    checked = validate_model(G, _k5(), bags)
                    stats["component_results"].append("witness_found")
                    return [checked[str(i)] for i in range(5)]
                for clause in cuts:
                    solver.add_clause(clause)
                stats["connectivity_cuts"] += len(cuts)
        stats["component_results"].append(result)
    return None


def find_k5_minor(G, trials=256, seed=0, exact_threshold=10, use_sat=False,
                  sat_conflicts=20000, sat_rounds=200, stats=None):
    """Try cheap reproducible contractions, small exact search, then optional SAT."""
    stats = stats if stats is not None else {}
    stats.update(status="no_witness_returned", stages=[])
    stage = {}
    found = random_k5_minor(G, trials=trials, seed=seed, stats=stage)
    stats["stages"].append(stage)
    if found is None and exact_threshold >= 5:
        stage = {}
        found = exact_k5_minor(G, max_order=exact_threshold, stats=stage)
        stats["stages"].append(stage)
    if found is None and use_sat:
        stage = {}
        found = sat_k5_minor(G, conflict_budget=sat_conflicts, max_rounds=sat_rounds, stats=stage)
        stats["stages"].append(stage)
    if found is not None:
        stats["status"] = "witness_found"
    return found


def _spanning_subgraph(Q, H, max_edge_deletions):
    if len(Q) != len(H):
        return None
    excess = Q.number_of_edges() - H.number_of_edges()
    if not 0 <= excess <= max_edge_deletions:
        return None
    qdegrees = sorted((d for _, d in Q.degree()), reverse=True)
    hdegrees = sorted((d for _, d in H.degree()), reverse=True)
    if any(a < b for a, b in zip(qdegrees, hdegrees)):
        return None
    target_invariant = graph_invariant(H)
    ordered_edges = sorted(tuple(sorted(e)) for e in Q.edges())
    for deleted in it.combinations(ordered_edges, excess):
        R = Q.copy()
        R.remove_edges_from(deleted)
        if sorted((d for _, d in R.degree()), reverse=True) != hdegrees:
            continue
        if graph_invariant(R) != target_invariant:
            continue
        matcher = nx.isomorphism.GraphMatcher(R, H)
        if matcher.is_isomorphic():
            return matcher.mapping, [list(e) for e in deleted]
    return None


def known_minor(G, catalog, max_contractions=2, max_edge_deletions=2,
                max_states=5000, stats=None):
    """Bounded contractions plus same-order edge deletions against known targets.

    ``catalog`` is an iterable of (name, NetworkX graph). The source must be
    strictly larger than the target in order or size. No catalogue membership
    or obstruction claim is established here, only the checked minor model.
    """
    require(max_contractions >= 0 and max_edge_deletions >= 0 and max_states > 0,
            "invalid known-minor budget")
    G = _graph(G)
    targets = {}
    for name, raw in sorted(catalog, key=lambda row: str(row[0])):
        H = _graph(raw)
        if (0 <= len(G) - len(H) <= max_contractions and
                (len(G) > len(H) or G.number_of_edges() > H.number_of_edges())):
            targets.setdefault(len(H), []).append((str(name), H))
    stats = stats if stats is not None else {}
    stats.update(method="bounded_contractions_spanning_deletions", states_checked=0,
                 max_states=max_states, max_contractions=max_contractions,
                 max_edge_deletions=max_edge_deletions, status="no_witness_returned")
    if not targets:
        return None
    queue = deque([(G, {v: [v] for v in G}, [])])
    seen = {tuple((v,) for v in sorted(G))}
    min_order = min(targets)
    while queue and stats["states_checked"] < max_states:
        Q, bags, trace = queue.popleft()
        stats["states_checked"] += 1
        for name, H in targets.get(len(Q), []):
            found = _spanning_subgraph(Q, H, max_edge_deletions)
            if found is not None:
                mapping, deleted = found
                model = canonical_labels(H, {mapping[v]: bags[v] for v in Q})
                checked = validate_model(G, H, model)
                stats["status"] = "witness_found"
                return {"type": "proper_obstruction_minor", "target": name,
                        "branch_sets": checked, "contraction_trace": trace,
                        "quotient_deleted_edges": deleted}
        if len(trace) >= max_contractions or len(Q) <= min_order:
            continue
        for u, v in sorted(tuple(sorted(e)) for e in Q.edges()):
            merged = {w: list(bag) for w, bag in bags.items() if w != v}
            merged[u] = sorted(bags[u] + bags[v])
            key = tuple(sorted(tuple(bag) for bag in merged.values()))
            if key in seen:
                continue
            seen.add(key)
            R = nx.contracted_nodes(Q, u, v, self_loops=False)
            R.graph.clear()
            for w in R:
                R.nodes[w].clear()
            queue.append((R, merged, trace + [[u, v]]))
    if queue:
        stats["status"] = "state_budget_exhausted"
    return None
