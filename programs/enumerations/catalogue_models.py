"""A search for a model of a graph of the catalogue in an assembly of E4.

The search contracts a forest of the graph (forests at the fixed Psi_1
first, then every set of edges of the given size) and looks for a spanning
monomorphism of a graph T_i into the quotient. A model found is checked
(minor_search.validate_model) before it is returned; None means only that
this bounded search returned no witness.
"""
import itertools as it
import time

import networkx as nx

from minor_search import validate_model


def embedding(source, target):
    """Return a target-to-source spanning monomorphism of adjacency mappings."""
    if len(source) != len(target):
        return None
    sv, tv = sorted(source), sorted(target)
    n = len(sv)
    si = {v: i for i, v in enumerate(sv)}
    ti = {v: i for i, v in enumerate(tv)}
    sa = [sum(1 << si[u] for u in source[v]) for v in sv]
    ta = [sum(1 << ti[u] for u in target[v]) for v in tv]
    domains = [sum(1 << j for j in range(n)
                   if sa[j].bit_count() >= ta[i].bit_count()) for i in range(n)]
    earlier, later = [[] for _ in range(n)], [[] for _ in range(n)]
    for i in range(n):
        for j in range(i):
            if (ta[i] & ~(1 << j)) == (ta[j] & ~(1 << i)):
                earlier[i].append(j)
                later[j].append(i)
    mapping = [-1] * n

    def visit(remaining_domains, todo, used):
        if not todo:
            return mapping.copy()
        available = {}
        for i in todo:
            mask = remaining_domains[i] & ~used
            for j in earlier[i]:
                if mapping[j] >= 0:
                    mask &= ~((1 << (mapping[j] + 1)) - 1)
            for j in later[i]:
                if mapping[j] >= 0:
                    mask &= (1 << mapping[j]) - 1
            if not mask:
                return None
            available[i] = mask
        i = min(todo, key=lambda j: (available[j].bit_count(), -ta[j].bit_count(), j))
        others = todo - {i}
        mask = available[i]
        while mask:
            bit = mask & -mask
            mask ^= bit
            j = bit.bit_length() - 1
            revised = list(remaining_domains)
            for k in others:
                if ta[i] & (1 << k):
                    revised[k] &= sa[j]
            mapping[i] = j
            result = visit(revised, others, used | bit)
            if result is not None:
                return result
            mapping[i] = -1
        return None

    result = visit(domains, set(range(n)), 0)
    return None if result is None else {tv[k]: sv[result[k]] for k in range(n)}


def _normalized(graph):
    if graph.is_directed() or graph.is_multigraph() or nx.number_of_selfloops(graph):
        raise ValueError('expected a simple undirected graph')
    labels = {v: str(v) for v in graph}
    if len(set(labels.values())) != len(labels):
        raise ValueError('vertex labels collide after conversion to strings')
    return nx.relabel_nodes(graph, labels, copy=True)


def _priority_forests(source, count):
    """Yield useful fixed-core forests first; all choices are checked below."""
    if count == 3:
        for roots, twins in (("ace", "pq"), ("bdf", "rz")):
            for a, b in it.combinations(roots, 2):
                yield ((a, twins[0]), (a, twins[1]), (b, twins[0]))
    if count == 5:
        for a, b in (("a", "b"), ("c", "d"), ("e", "f")):
            yield ((a, b), (a, "p"), (a, "q"), (b, "r"), (b, "z"))
    if count > 0:
        for center in "abcdef":
            if center in source:
                for neighbors in it.combinations(sorted(source[center]), count):
                    yield tuple((center, neighbor) for neighbor in neighbors)


def search(G, catalog, max_contractions=5, stats=None):
    """Return a checked proper catalogue model, or no witness from this search.

    ``catalog`` is an iterable of (target_name, NetworkX_graph). The default
    depth order is 3,4,5,2,1,0. The source and catalogue are not mutated.
    """
    source = _normalized(G)
    if type(max_contractions) is not int or max_contractions < 0:
        raise ValueError('max_contractions must be a nonnegative integer')
    stats = {} if stats is None else stats
    stats.update(status='no_witness_returned', max_contractions=max_contractions,
                 forests_tried=0, distinct_quotients=0, monomorphism_tests=0)
    started = time.monotonic()
    vertices = sorted(source)
    original_edges = sorted(tuple(sorted(e)) for e in source.edges())
    edge_set = set(original_edges)
    targets = {}
    for name, raw_target in catalog:
        target = _normalized(raw_target)
        count = len(source) - len(target)
        if not (0 <= count <= max_contractions):
            continue
        if count == 0 and target.number_of_edges() >= source.number_of_edges():
            continue
        adjacency = {v: set(target[v]) for v in target}
        degrees = sorted(map(len, adjacency.values()), reverse=True)
        targets.setdefault(count, []).append((name, target, adjacency, degrees))
    order = [d for d in (3, 4, 5, 2, 1, 0) if d in targets]
    order.extend(sorted(set(targets) - set(order)))
    stats['contraction_depth_order'] = order
    for count in order:
        seen = set()
        choices = it.chain(_priority_forests(source, count), it.combinations(original_edges, count))
        for chosen in choices:
            selected = tuple(sorted(tuple(sorted(e)) for e in chosen))
            if len(set(selected)) != count or any(e not in edge_set for e in selected):
                continue
            stats['forests_tried'] += 1
            parent = {v: v for v in vertices}

            def root(v):
                while parent[v] != v:
                    parent[v] = parent[parent[v]]
                    v = parent[v]
                return v

            for u, v in selected:
                u, v = root(u), root(v)
                if u == v:
                    break
                parent[max(u, v)] = min(u, v)
            else:
                bags = {}
                for v in vertices:
                    bags.setdefault(root(v), []).append(v)
                key = tuple(sorted(tuple(bag) for bag in bags.values()))
                if key in seen:
                    continue
                seen.add(key)
                stats['distinct_quotients'] += 1
                quotient = {v: set() for v in bags}
                for u, v in original_edges:
                    u, v = root(u), root(v)
                    if u != v:
                        quotient[u].add(v)
                        quotient[v].add(u)
                degrees = sorted(map(len, quotient.values()), reverse=True)
                for name, target, target_adjacency, target_degrees in targets[count]:
                    if any(a < b for a, b in zip(degrees, target_degrees)):
                        continue
                    stats['monomorphism_tests'] += 1
                    mapping = embedding(quotient, target_adjacency)
                    if mapping is None:
                        continue
                    model = validate_model(source, target,
                        {v: bags[u] for v, u in mapping.items()})
                    stats.update(status='witness_found', target=name,
                                 seconds=time.monotonic() - started)
                    return {'type': 'proper_obstruction_minor', 'target': name,
                            'branch_sets': model,
                            'contracted_original_forest': [list(e) for e in selected],
                            'independently_validated': True}
    stats['seconds'] = time.monotonic() - started
    return None
