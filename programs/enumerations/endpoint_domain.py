"""The domain of E4 (Theorem 4.22): assemblies of bridges on Psi_1 when
every path P_i is an edge.

The names: a = s1, b = t1, c = s2, d = t2, e = s3, f = t3, with twins p, q
on S = {a, c, e} and r, z on T = {b, d, f}; the paths are the edges ab, cd,
ef. The atoms are the 72 labeled bridge forms on the triples with two
attachments on one side and one on the other (one interior vertex, T(x),
S(x;y), as in the proof), the single vertices on the three sets
{s_i, t_i, s_j, t_j}, and the six edges s_i t_j. A configuration gives the
number of copies of each atom, at most three on one triple and at most five
three-attachment bridges in all (none of the other atoms with five), with
the conditions of the proof (Lemmas 3.4(a), 4.7, 4.16 and 4.20). base_cnf() gives
these conditions as clauses: the first 225 variables say "at least k copies
of atom i", and 1,070 auxiliary variables form a counter for the bound of
five. Every legal count vector has the canonical counter extension of
canonical_assignment(), and every satisfying assignment decodes to a legal
vector.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
import itertools as it
import json
from typing import Iterable, Mapping, Sequence


ROOTS = tuple("abcdef")
S = tuple("ace")
T = tuple("bdf")
CORE_VERTICES = tuple("abcdefpqrz")
MATCHING = (("a", "b"), ("c", "d"), ("e", "f"))
CROSS_EDGES = (("a", "d"), ("a", "f"), ("b", "c"),
               ("b", "e"), ("c", "f"), ("d", "e"))
QUAD_SUPPORTS = (tuple("abcd"), tuple("abef"), tuple("cdef"))
OPPOSITE_PAIRS = ((("a", "d"), ("b", "c")),
                  (("a", "f"), ("b", "e")),
                  (("c", "f"), ("d", "e")))


def edge(u: str, v: str) -> tuple[str, str]:
    return tuple(sorted((u, v)))


CORE_EDGES = tuple(sorted(
    [edge(v, s) for v in ("p", "q") for s in S]
    + [edge(v, t) for v in ("r", "z") for t in T]
    + list(MATCHING)))


@dataclass(frozen=True)
class Atom:
    id: str
    kind: str
    support: tuple[str, ...]
    form: str
    vertices: tuple[str, ...]
    edges: tuple[tuple[str, str], ...]
    max_count: int

    @property
    def interior(self) -> tuple[str, ...]:
        return tuple(v for v in self.vertices if v.startswith("@"))

    def record(self) -> dict:
        return {"id": self.id, "kind": self.kind,
                "support": list(self.support), "form": self.form,
                "vertices": list(self.vertices),
                "edges": [list(e) for e in self.edges],
                "interior": list(self.interior), "max_count": self.max_count}


def _triple(support: tuple[str, ...], form: str,
            distinguished: tuple[str, ...] = ()) -> Atom:
    """Construct the source forms with @-prefixed local interior labels."""
    if form == "singleton":
        raw = [("@h", v) for v in support]
        suffix = "single"
        name = "singleton"
    elif form == "T":
        r, = distinguished
        others = [v for v in support if v != r]
        raw = [(d, v) for d in ("@d1", "@d2") for v in others + ["@u"]]
        raw += [("@u", r)]
        suffix = "T-" + r
        name = f"T({r})"
    elif form == "S":
        x, y = distinguished
        z, = [v for v in support if v not in (x, y)]
        raw = [(d, v) for d in ("@d1", "@d2") for v in (x, y, "@u")]
        raw += [("@f", v) for v in (x, z, "@u")]
        suffix = "S-" + x + "-" + y
        name = f"S({x};{y})"
    else:
        raise ValueError("Unknown triple form")
    edges = tuple(sorted(edge(*e) for e in raw))
    vertices = tuple(sorted(set(support).union(v for e in edges for v in e)))
    return Atom("B-" + "".join(support) + "-" + suffix,
                "triple", support, name, vertices, edges, 3)


def _adjacency(edges: Iterable[tuple[str, str]]) -> dict[str, set[str]]:
    adjacency: dict[str, set[str]] = {}
    for u, v in edges:
        adjacency.setdefault(u, set()).add(v)
        adjacency.setdefault(v, set()).add(u)
    return adjacency


def _fixed_neighbor_violation(atom: Atom) -> bool:
    adjacency = _adjacency(CORE_EDGES + atom.edges)
    return any(len(adjacency[v]) == 3 and any(
        b in adjacency[a] for a, b in it.combinations(sorted(adjacency[v]), 2))
        for v in atom.interior)


@lru_cache(maxsize=1)
def atoms() -> tuple[Atom, ...]:
    """Return 72 triple types, then 3 quad types, then 6 crossedges."""
    triples = []
    for support in it.combinations(ROOTS, 3):
        if len(set(support) & set(S)) not in (1, 2):
            continue
        same = sorted(set(support) & (set(S) if len(set(support) & set(S)) == 2
                                     else set(T)))
        x, y = same
        candidates = [_triple(support, "singleton")]
        candidates += [_triple(support, "T", (r,)) for r in support]
        candidates += [_triple(support, "S", (x, y)),
                       _triple(support, "S", (y, x))]
        triples += [a for a in candidates if not _fixed_neighbor_violation(a)]
    triples.sort(key=lambda a: a.id)
    quads = [Atom("Q-" + "".join(support), "quad", support, "singleton",
                  ("@h",) + support,
                  tuple(sorted(edge("@h", v) for v in support)), 1)
             for support in QUAD_SUPPORTS]
    extras = [Atom("E-" + "".join(e), "crossedge", e, "edge", e, (e,), 1)
              for e in CROSS_EDGES]
    assert len(triples) == 72
    assert Counter(Counter(a.support for a in triples).values()) == {3: 12, 6: 6}
    return tuple(triples + quads + extras)


@lru_cache(maxsize=1)
def threshold_variables() -> dict[str, tuple[int, ...]]:
    result = {}
    next_variable = 1
    for atom in atoms():
        result[atom.id] = tuple(range(next_variable, next_variable + atom.max_count))
        next_variable += atom.max_count
    assert next_variable == 226
    return result


ATOMS = tuple(dict(a.record(), index=i,
                   threshold_variables=list(threshold_variables()[a.id]))
              for i, a in enumerate(atoms()))
MAX_COUNTS = tuple(a.max_count for a in atoms())


def variable(i: int, k: int) -> int:
    """Threshold literal for zero-based atom index and one-based level."""
    if not 0 <= i < len(ATOMS) or not 1 <= k <= MAX_COUNTS[i]:
        raise ValueError("Atom index or count threshold out of range")
    return threshold_variables()[ATOMS[i]["id"]][k - 1]


def normalize_counts(counts: Mapping[str, int] | Sequence[int]) -> dict[str, int]:
    if not isinstance(counts, Mapping):
        if len(counts) != len(ATOMS):
            raise ValueError("A count sequence must have 81 entries")
        counts = dict(zip((a.id for a in atoms()), counts))
    known = {a.id: a for a in atoms()}
    if set(counts) - set(known):
        raise ValueError("Unknown atom ids: " + str(sorted(set(counts) - set(known))))
    result = {}
    for atom in atoms():
        n = counts.get(atom.id, 0)
        if not isinstance(n, int) or not 0 <= n <= atom.max_count:
            raise ValueError(f"Invalid count for {atom.id}: {n!r}")
        result[atom.id] = int(n)
    return result


def configuration_from_model(model: Iterable[int]) -> dict[str, int]:
    """Decode the first 225 variables; reject missing/nonmonotone bits."""
    values = {}
    for literal in model:
        variable = abs(literal)
        value = literal > 0
        if variable in values and values[variable] != value:
            raise ValueError("Conflicting model literals")
        values[variable] = value
    result = {}
    for atom, variables in threshold_variables().items():
        if any(v not in values for v in variables):
            raise ValueError("Missing atom threshold variable")
        bits = [values[v] for v in variables]
        if bits != sorted(bits, reverse=True):
            raise ValueError("Nonmonotone atom threshold variables")
        result[atom] = sum(bits)
    return result


def counts_from_model(model: Iterable[int]) -> tuple[int, ...]:
    decoded = configuration_from_model(model)
    return tuple(decoded[a.id] for a in atoms())


def graph_record(counts: Mapping[str, int]) -> dict:
    """Deterministic labels, not a claim of unrooted isomorphism canonization.

    Local interiors are disjoint.  No construction filter is applied here;
    callers may construct a count vector and then request its violations.
    """
    counts = normalize_counts(counts)
    vertices = set(CORE_VERTICES)
    edges = set(CORE_EDGES)
    instances = []
    for index, atom in enumerate(atoms()):
        for copy in range(counts[atom.id]):
            mapping = {v: (f"v{index:02d}_{copy}_{v[1:]}" if v.startswith("@") else v)
                       for v in atom.vertices}
            vertices.update(mapping.values())
            edges.update(edge(mapping[u], mapping[v]) for u, v in atom.edges)
            instances.append({"atom_id": atom.id, "copy": copy,
                              "interior": [mapping[v] for v in atom.interior],
                              "support": list(atom.support)})
    return {"vertices": sorted(vertices), "edges": [list(e) for e in sorted(edges)],
            "roots": list(ROOTS), "atom_counts": {k: v for k, v in counts.items() if v},
            "instances": instances}


def build_graph(counts: Mapping[str, int]):
    """NetworkX convenience wrapper; all domain/CNF operations are stdlib."""
    import networkx as nx
    record = graph_record(counts)
    graph = nx.Graph()
    graph.add_nodes_from(record["vertices"])
    graph.add_edges_from(record["edges"])
    return graph


graph_from_counts = build_graph


def _opposite_conflict(a: Atom, b: Atom) -> bool:
    A, B = set(a.support), set(b.support)
    return any((set(e) <= A and set(f) <= B) or (set(f) <= A and set(e) <= B)
               for e, f in OPPOSITE_PAIRS)


def _neighbor_conflict(atom: Atom, crossedge: Atom) -> bool:
    adjacency = _adjacency(atom.edges)
    return any(len(adjacency[v]) == 3 and set(crossedge.support) <= adjacency[v]
               for v in atom.interior)


def validate_configuration(counts: Mapping[str, int]) -> list[dict]:
    """Semantic validation independent of clause evaluation, not classification."""
    counts = normalize_counts(counts)
    present = [a for a in atoms() if counts[a.id]]
    violations = []
    totals = Counter()
    for a in present:
        if a.kind == "triple":
            totals[a.support] += counts[a.id]
    for support, count in sorted(totals.items()):
        if count > 3:
            violations.append({"constraint": "support_multiplicity", "support": support,
                               "count": count})
    m = sum(totals.values())
    auxiliaries = [a.id for a in present if a.kind != "triple"]
    if m > 5:
        violations.append({"constraint": "total_triples", "count": m})
    if m >= 5 and auxiliaries:
        violations.append({"constraint": "five_triples_no_auxiliary_atoms",
                           "count": m, "atoms": auxiliaries})
    for a, b in it.combinations(present, 2):
        if _opposite_conflict(a, b):
            violations.append({"constraint": "opposite_distinct_bridge_paths",
                               "atoms": [a.id, b.id]})
    for a in present:
        if a.kind == "quad":
            for b in present:
                if b.id != a.id and set(b.support) <= set(a.support):
                    violations.append({"constraint": "four_root_isolation",
                                       "atoms": [a.id, b.id]})
    record = graph_record(counts)
    adjacency = _adjacency(tuple(map(tuple, record["edges"])))
    for v, neighbors in sorted(adjacency.items()):
        if len(neighbors) == 3:
            for u, w in it.combinations(sorted(neighbors), 2):
                if w in adjacency[u]:
                    violations.append({"constraint": "independent_degree_three_neighborhood",
                                       "vertex": v, "adjacent_neighbors": [u, w]})
    return violations


def legal(counts: Mapping[str, int] | Sequence[int]) -> bool:
    try:
        return not validate_configuration(counts)
    except ValueError:
        return False


def _compositions(total: int, length: int, cap: int = 3):
    if length == 0:
        if total == 0:
            yield ()
        return
    for n in range(min(cap, total) + 1):
        for tail in _compositions(total - n, length - 1, cap):
            yield (n,) + tail


@dataclass
class DomainCNF:
    clauses: list[list[int]]
    variables: dict[int, dict]
    groups: list[dict]
    counter: dict[tuple[int, int], int]
    triple_literals: list[int]

    @property
    def nv(self) -> int:
        return max(self.variables)

    def dimacs(self) -> str:
        return f"p cnf {self.nv} {len(self.clauses)}\n" + "".join(
            " ".join(map(str, c)) + " 0\n" for c in self.clauses)


def domain_cnf() -> DomainCNF:
    """Rebuild every clause deterministically, with no target-based condition."""
    atom_list = atoms()
    thresholds = threshold_variables()
    clauses = []
    seen = set()
    variables = {}
    groups = []
    counter = {}

    def add(clause):
        normalized = tuple(sorted(set(clause), key=lambda x: (abs(x), x)))
        assert not any(-x in normalized for x in normalized)
        if normalized not in seen:
            seen.add(normalized)
            clauses.append(list(normalized))

    def finish_group(name, start):
        groups.append({"name": name, "first_clause_zero_based": start,
                       "clause_count": len(clauses) - start})

    for atom in atom_list:
        for k, variable in enumerate(thresholds[atom.id], 1):
            variables[variable] = {"kind": "atom_threshold", "atom_id": atom.id,
                                   "at_least": k}
    start = len(clauses)
    for atom in atom_list:
        for low, high in zip(thresholds[atom.id], thresholds[atom.id][1:]):
            add([-high, low])
    finish_group("atom_threshold_monotonicity", start)

    start = len(clauses)
    supports = sorted({a.support for a in atom_list if a.kind == "triple"})
    for support in supports:
        block = [a for a in atom_list if a.kind == "triple" and a.support == support]
        for counts in _compositions(4, len(block)):
            add([-thresholds[a.id][n - 1] for a, n in zip(block, counts) if n])
    finish_group("total_at_most_three_on_each_triple_support", start)

    start = len(clauses)
    for atom in atom_list:
        if atom.kind == "triple":
            for extra in atom_list:
                if extra.kind == "crossedge" and _neighbor_conflict(atom, extra):
                    add([-thresholds[atom.id][0], -thresholds[extra.id][0]])
    finish_group("independent_degree_three_neighborhood", start)

    start = len(clauses)
    for a, b in it.combinations(atom_list, 2):
        if _opposite_conflict(a, b):
            add([-thresholds[a.id][0], -thresholds[b.id][0]])
    finish_group("opposite_cross_connections_of_distinct_bridges", start)

    start = len(clauses)
    for quad in atom_list:
        if quad.kind == "quad":
            for other in atom_list:
                if other.id != quad.id and set(other.support) <= set(quad.support):
                    add([-thresholds[quad.id][0], -thresholds[other.id][0]])
    finish_group("four_root_isolation_remaining_clauses", start)

    triple_literals = [v for a in atom_list if a.kind == "triple"
                       for v in thresholds[a.id]]
    assert triple_literals == list(range(1, 217))
    next_variable = 226
    for i in range(1, len(triple_literals) + 1):
        for j in range(1, min(i, 5) + 1):
            counter[i, j] = next_variable
            variables[next_variable] = {"kind": "prefix_count_at_least", "prefix": i,
                                       "at_least": j}
            next_variable += 1
    assert len(counter) == 1070 and next_variable == 1296
    start = len(clauses)
    for i, x in enumerate(triple_literals, 1):
        add([-x, counter[i, 1]])
        for j in range(1, min(i - 1, 5) + 1):
            add([-counter[i - 1, j], counter[i, j]])
        for j in range(2, min(i, 5) + 1):
            add([-x, -counter[i - 1, j - 1], counter[i, j]])
        if i > 5:
            add([-x, -counter[i - 1, 5]])
    finish_group("global_total_triples_at_most_five_prefix_counter", start)
    start = len(clauses)
    for atom in atom_list:
        if atom.kind != "triple":
            add([-counter[len(triple_literals), 5], -thresholds[atom.id][0]])
    finish_group("five_triples_forbid_quad_and_crossedge_atoms", start)
    return DomainCNF(clauses, variables, groups, counter, triple_literals)


def base_cnf() -> tuple[list[list[int]], int]:
    cnf = domain_cnf()
    return cnf.clauses, cnf.nv


def canonical_assignment(counts: Mapping[str, int], cnf: DomainCNF | None = None) -> dict[int, bool]:
    """Exact prefix counts give a canonical auxiliary extension, even if illegal."""
    counts = normalize_counts(counts)
    cnf = domain_cnf() if cnf is None else cnf
    values = {v: counts[atom] >= k for atom, literals in threshold_variables().items()
              for k, v in enumerate(literals, 1)}
    total = 0
    for i, x in enumerate(cnf.triple_literals, 1):
        total += values[x]
        for j in range(1, min(i, 5) + 1):
            values[cnf.counter[i, j]] = total >= j
    return values


def satisfies(clauses: Iterable[Iterable[int]], values: Mapping[int, bool]) -> bool:
    return all(any(values[abs(lit)] == (lit > 0) for lit in clause) for clause in clauses)
