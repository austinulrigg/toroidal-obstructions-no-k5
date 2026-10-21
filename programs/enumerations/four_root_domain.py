"""Construction 3.16 for four attachments with no outside vertices, no
required edges and one connected interior (E1 and E3).

The four attachments are placed on a K3,3 (at branch vertices, inside branch
paths, or off it) up to the 72 symmetries of K3,3. A segment between
consecutive attachments that are not joined by an added edge gets one
interior vertex; the added edges are removed; the interior pieces are joined
by every tree with every choice of ends, and each attachment without an
interior neighbour gets one in every way. The resulting parts are compared
up to isomorphism fixing each attachment, and the first recipe of each is
kept. domain(roots, added_edges, prefix) returns the parts in the order in
which they are found, with the counts of the construction.
"""
import hashlib
import itertools as it
import json
import math
from dataclasses import dataclass

import networkx as nx

BRANCH_EDGES = tuple((a, b) for a in range(3) for b in range(3, 6))


@dataclass(frozen=True, order=True)
class Placement:
    branches: tuple
    paths: tuple

    def payload(self):
        return {"branches": list(self.branches), "paths": [list(p) for p in self.paths]}

    @property
    def id(self):
        data = json.dumps(self.payload(), separators=(",", ":"), sort_keys=True)
        return "placement-" + hashlib.sha256(data.encode()).hexdigest()[:20]


def raw_placements():
    branches = [-1] * 6
    paths = [[] for _ in range(9)]

    def extend(root):
        if root == 4:
            yield Placement(tuple(branches), tuple(tuple(p) for p in paths))
            return
        yield from extend(root + 1)  # this attachment is off the subdivision
        for vertex in range(6):
            if branches[vertex] == -1:
                branches[vertex] = root
                yield from extend(root + 1)
                branches[vertex] = -1
        for path in paths:
            for position in range(len(path) + 1):
                path.insert(position, root)
                yield from extend(root + 1)
                path.pop(position)

    yield from extend(0)


def automorphisms():
    for left in it.permutations(range(3)):
        for right in it.permutations(range(3)):
            for swap in (False, True):
                if not swap:
                    yield left + tuple(v + 3 for v in right)
                else:
                    yield tuple(v + 3 for v in left) + right


AUTOMORPHISMS = tuple(automorphisms())


def transform(placement, permutation):
    branches = [-1] * 6
    paths = [()] * 9
    for old, new in enumerate(permutation):
        branches[new] = placement.branches[old]
    for old_index, (a, b) in enumerate(BRANCH_EDGES):
        x, y = permutation[a], permutation[b]
        sequence = placement.paths[old_index]
        if x >= 3:
            x, y = y, x
            sequence = sequence[::-1]
        paths[x * 3 + y - 3] = sequence
    return Placement(tuple(branches), tuple(paths))


def placement_classes():
    covered, classes, raw_count = set(), [], 0
    for placement in raw_placements():
        raw_count += 1
        if placement in covered:
            continue
        orbit = {transform(placement, a) for a in AUTOMORPHISMS}
        covered.update(orbit)
        classes.append(min(orbit))
    assert raw_count == len(covered)
    classes.sort()
    expected_raw = sum(math.comb(4, b) * math.perm(6, b) * math.comb(4 - b, i) *
                       math.factorial(i) * math.comb(i + 8, 8) for b in range(5) for i in range(5 - b))
    assert raw_count == expected_raw
    return classes, raw_count


def edge_list(graph):
    return [list(e) for e in sorted(tuple(sorted(e)) for e in graph.edges())]


def color_graph(graph, roots):
    for vertex in graph:
        graph.nodes[vertex]["root_color"] = vertex if vertex in roots else "interior"
    return graph


def build_subdivision(placement, roots, added_edges):
    """The K3,3 subdivision F of the placement and the part P = F minus the
    added edges, with all four attachments kept."""
    F = nx.Graph()
    branch_labels = []
    interior_count = 0
    for root_index in placement.branches:
        if root_index >= 0:
            branch_labels.append(roots[root_index])
        else:
            branch_labels.append(f"i{interior_count}")
            interior_count += 1
    F.add_nodes_from(branch_labels)
    D = {frozenset(e) for e in added_edges}
    for index, (a, b) in enumerate(BRANCH_EDGES):
        sequence = [branch_labels[a], *(roots[r] for r in placement.paths[index]), branch_labels[b]]
        for x, y in it.pairwise(sequence):
            if x in roots and y in roots and frozenset((x, y)) not in D:
                middle = f"i{interior_count}"
                interior_count += 1
                F.add_edges_from(((x, middle), (middle, y)))
            else:
                F.add_edge(x, y)
    P = F.copy()
    P.remove_edges_from(added_edges)
    P.add_nodes_from(roots)
    assert 1 <= interior_count <= 9
    assert not any(P.has_edge(a, b) for a, b in it.combinations(roots, 2))
    return color_graph(F, roots), color_graph(P, roots)


def interior_components(P, roots):
    interior = sorted(set(P) - set(roots))
    return [tuple(sorted(c)) for c in sorted(nx.connected_components(P.subgraph(interior)),
                                             key=lambda c: tuple(sorted(c)))]


def component_trees(k):
    if k == 1:
        yield ()
    elif k == 2:
        yield ((0, 1),)
    else:
        for sequence in it.product(range(k), repeat=k - 2):
            tree = nx.from_prufer_sequence(list(sequence))
            yield tuple(sorted(tuple(sorted(e)) for e in tree.edges()))


def iter_recipes(placement, roots, added_edges):
    _, P = build_subdivision(placement, roots, added_edges)
    components = interior_components(P, roots)
    interior = tuple(sorted(set(P) - set(roots)))
    missing = tuple(r for r in roots if P.degree(r) == 0)
    recipe_index = 0
    for tree in component_trees(len(components)):
        endpoint_choices = [tuple(it.product(components[a], components[b])) for a, b in tree]
        for forest in it.product(*endpoint_choices):
            for neighbors in it.product(interior, repeat=len(missing)):
                attachments = tuple(zip(missing, neighbors))
                graph = P.copy()
                graph.add_edges_from(forest)
                graph.add_edges_from(attachments)
                yield graph, {"placement_id": placement.id, "placement": placement.payload(),
                              "recipe_index": recipe_index,
                              "forest_edges": [list(e) for e in forest],
                              "restored_attachment_edges": [list(e) for e in attachments]}
                recipe_index += 1


def graph_bucket(graph):
    digest = nx.weisfeiler_lehman_graph_hash(graph, node_attr="root_color", iterations=4)
    return f"{len(graph)}:{graph.number_of_edges()}:{digest}"


def rooted_isomorphic(left, right):
    return nx.is_isomorphic(left, right, node_match=lambda x, y: x["root_color"] == y["root_color"])


def domain(roots, added_edges, prefix):
    """The parts, in order, with the counts of the construction. The id of a
    part is prefix followed by the first 24 hexadecimal digits of the SHA-256
    of its sorted edge list."""
    placements, raw_count = placement_classes()
    buckets, records, recipes = {}, [], 0
    for placement in placements:
        for graph, recipe in iter_recipes(placement, roots, added_edges):
            recipes += 1
            bucket = graph_bucket(graph)
            if any(rooted_isomorphic(graph, old) for old in buckets.get(bucket, [])):
                continue
            buckets.setdefault(bucket, []).append(graph)
            edges = edge_list(graph)
            identifier = prefix + hashlib.sha256(json.dumps(edges, separators=(",", ":")).encode()).hexdigest()[:24]
            records.append({"id": identifier, "vertices": sorted(graph), "edges": edges,
                            "roots": list(roots), "interior_vertices": sorted(set(graph) - set(roots)),
                            "bucket": bucket, "first_recipe": recipe})
    stats = {"raw_placements": raw_count, "placement_symmetry_classes": len(placements),
             "expanded_recipes": recipes, "rooted_candidates": len(records)}
    return records, stats


def subdivision_of(record, roots, added_edges):
    """The K3,3 subdivision a part came from (a subgraph of the part plus the
    added edges), as a sorted edge list."""
    placement = record['first_recipe']['placement']
    F, _ = build_subdivision(Placement(tuple(placement['branches']), tuple(map(tuple, placement['paths']))),
                             roots, added_edges)
    return edge_list(F)
