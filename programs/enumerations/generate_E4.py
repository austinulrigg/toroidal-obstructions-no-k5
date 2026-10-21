"""E4: assemblies of bridges on Psi_1 when every path P_i is an edge (Theorem 4.23).

The domain is that of endpoint_domain.py: count vectors of 81 atoms (bridge
forms on triples, single vertices on the sets {s_i, t_i, s_j, t_j}, and the
edges s_i t_j) satisfying the conditions of the proof, given as clauses.

Covering. A SAT solver (PySAT, Glucose 3) proposes a count vector that no
recorded graph covers yet. Its graph is tested with the joint test of
Lemma 4.3 (joint_test.py). If some choice passes, the graph is toroidal, and
copies are added while it still passes; every smaller vector is then
toroidal, and the clause "some atom has more copies" is added. If every
choice fails, copies are removed while it still fails; every larger vector
is then nontoroidal, and the clause "some atom of it has fewer copies" is
added. Each recorded vector is taken with its images under the twelve
symmetries of Psi_1. When the clauses become unsatisfiable, every vector of
the domain lies below a recorded toroidal graph or above a recorded
nontoroidal one; the final formula is written to coverage.cnf with a DRAT
proof of its unsatisfiability, coverage.drat.

Classification of each recorded nontoroidal graph, in order: isomorphic to
some T_i; a one-edge deletion or contraction outside Psi_1 after which every
choice of the joint test fails; a model of some T_i carried over from an
isomorphic graph classified earlier; a K5 minor; a model of some T_i found
by catalogue_models.search. The proof finds T48 and T58.

The covering keeps coverage.jsonl and the classification
candidate_outcomes.jsonl, so that a second start continues.
Run as: python generate_E4.py FOLDER
"""
import itertools as it
import json
import sys
import time
from collections import Counter
from functools import lru_cache
from pathlib import Path

import networkx as nx
from pysat.solvers import Solver

import common
import endpoint_domain as dom
from catalogue_models import search
from joint_test import joint_test, graph_record, CORE_VERTICES, CORE_EDGES
from minor_search import find_k5_minor, graph_invariant, first_isomorphism, PositiveModelCache


def clause(counts, direction):
    if direction == 'down':
        return [dom.variable(i, c + 1) for i, c in enumerate(counts) if c < dom.MAX_COUNTS[i]]
    assert direction == 'up'
    return [-dom.variable(i, c) for i, c in enumerate(counts) if c]


@lru_cache(maxsize=200000)
def passes(counts):
    return joint_test(dom.graph_from_counts(counts), certify=False, stop_on_pass=True)['passes']


def strengthen(counts, direction):
    result = list(counts)
    # Both properties are monotone. An increment that fails remains failing
    # after later increments; a deletion that passes remains passing after
    # later deletions. Thus one coordinate sweep is enough.
    for i, maximum in enumerate(dom.MAX_COUNTS):
        step = 1 if direction == 'down' else -1
        while (result[i] < maximum if step == 1 else result[i] > 0):
            candidate = result.copy()
            candidate[i] += step
            if not dom.legal(candidate):
                break
            if passes(tuple(candidate)) != (direction == 'down'):
                break
            result = candidate
    return tuple(result)


def symmetry_actions():
    """The permutations of the atoms induced by the twelve symmetries of Psi_1."""
    size = len(dom.MAX_COUNTS)
    basis = []
    for i in range(size):
        c = [0] * size
        c[i] = 1
        g = dom.graph_from_counts(c)
        for v in g:
            g.nodes[v]['fixed'] = v if v in CORE_VERTICES else ''
        basis.append(g)
    invariants = [graph_invariant(g) for g in basis]
    actions = []
    for permutation in it.permutations(range(3)):
        for swap in (False, True):
            mapping = {}
            for side in (0, 1):
                for i in range(3):
                    mapping['ace'[i] if side == 0 else 'bdf'[i]] = (
                        'ace' if (side ^ swap) == 0 else 'bdf')[permutation[i]]
            mapping.update(dict(zip('pqrz', 'rzpq' if swap else 'pqrz')))
            action = []
            for i, g in enumerate(basis):
                h = nx.relabel_nodes(g, mapping, copy=True)
                for v in h:
                    h.nodes[v]['fixed'] = v if v in CORE_VERTICES else ''
                matches = [j for j, target in enumerate(basis)
                           if invariants[j] == invariants[i]
                           and nx.is_isomorphic(h, target,
                                                node_match=nx.algorithms.isomorphism.categorical_node_match('fixed', ''))]
                assert len(matches) == 1, (i, mapping, matches)
                action.append(matches[0])
            assert sorted(action) == list(range(size))
            if action not in actions:
                actions.append(action)
    assert len(actions) == 12
    return actions


def orbit(counts, actions):
    result = set()
    for action in actions:
        transformed = [0] * len(counts)
        for i, j in enumerate(action):
            transformed[j] = counts[i]
        assert dom.legal(transformed)
        result.add(tuple(transformed))
    return sorted(result)


def read_jsonl(path):
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def cover(outdir):
    base_clauses, nvars = dom.base_cnf()
    clauses = [list(c) for c in base_clauses]
    coverage_path = outdir / 'coverage.jsonl'
    records = read_jsonl(coverage_path)
    seen = {(r['direction'], tuple(r['counts'])) for r in records}
    for record in records:
        assert record['clause'] == clause(record['counts'], record['direction'])
        clauses.append(record['clause'])
    actions = symmetry_actions()
    counts = Counter(r['direction'] for r in records)
    started = last = time.monotonic()
    with Solver(name='g3', bootstrap_with=clauses) as solver, coverage_path.open('a') as output:
        while solver.solve():
            proposed = tuple(dom.counts_from_model(solver.get_model()))
            assert dom.legal(proposed)
            direction = 'down' if passes(proposed) else 'up'
            strengthened = strengthen(proposed, direction)
            for vector in orbit(strengthened, actions):
                if (direction, vector) in seen:
                    continue
                g = dom.graph_from_counts(vector)
                test = joint_test(g, certify=True, stop_on_pass=direction == 'down')
                assert test['passes'] == (direction == 'down')
                block = clause(vector, direction)
                record = {'id': 'E4-%06d' % (len(records) + 1), 'direction': direction,
                          'counts': list(vector), 'graph': graph_record(g),
                          'joint_test': test, 'clause': block}
                output.write(json.dumps(record, sort_keys=True, separators=(',', ':')) + '\n')
                output.flush()
                records.append(record)
                seen.add((direction, vector))
                counts[direction] += 1
                clauses.append(block)
                solver.add_clause(block)
            if time.monotonic() - last > 10:
                print(json.dumps({'stage': 'covering', 'records': dict(counts),
                                  'seconds': round(time.monotonic() - started, 1)}), flush=True)
                last = time.monotonic()
    with (outdir / 'coverage.cnf').open('w') as stream:
        stream.write(f'p cnf {nvars} {len(clauses)}\n')
        for row in clauses:
            stream.write(' '.join(map(str, row)) + ' 0\n')
    # Only the final formula is proof-logged; the clauses added during the
    # covering are inputs of the proof, not learned clauses.
    with Solver(name='g3', bootstrap_with=clauses, with_proof=True) as proof_solver:
        assert proof_solver.solve() is False
        proof = proof_solver.get_proof()
    (outdir / 'coverage.drat').write_text('\n'.join(proof) + '\n')
    return records, {'records': dict(counts), 'clauses': len(clauses), 'variables': nvars}


def joint_minor(g):
    variable_edges = [e for e in sorted(tuple(sorted(e)) for e in g.edges()) if e not in CORE_EDGES]
    for kind in ('delete_edge', 'contract_edge'):
        for u, v in variable_edges:
            operation = {'type': kind, 'edge': [u, v]}
            if kind == 'delete_edge':
                h = g.copy()
                h.remove_edge(u, v)
            else:
                if u in CORE_VERTICES and v in CORE_VERTICES:
                    continue
                if v in CORE_VERTICES:
                    u, v = v, u
                h = nx.contracted_nodes(g, u, v, self_loops=False)
                operation['keep'] = u
            if not joint_test(h, certify=False, stop_on_pass=True)['passes']:
                return {'type': 'proper_joint_minor', 'operation': operation,
                        'graph': graph_record(h), 'joint_test': joint_test(h, certify=True)}
    return None


def classify_graph(g, catalog, positive_cache):
    invariant = graph_invariant(g)
    for name, h in catalog:
        if graph_invariant(h) != invariant:
            continue
        mapping = first_isomorphism(g, h)
        if mapping is not None:
            return {'type': 'exact_survivor', 'target': name,
                    'mapping': {str(v): str(w) for v, w in mapping.items()}}
    witness = joint_minor(g)
    if witness is not None:
        return witness
    witness = positive_cache.find(g)
    if witness is not None:
        return {'type': 'proper_obstruction_minor', 'target': witness['metadata']['target'],
                'branch_sets': witness['branch_sets'],
                'search': 'checked isomorphism transport of a positive catalogue model'}
    stats = {}
    model = find_k5_minor(g, trials=128, seed=0, exact_threshold=8, use_sat=True,
                          sat_conflicts=10000, sat_rounds=100, stats=stats)
    if model is not None:
        return {'type': 'k5_minor', 'branch_sets': model, 'search': stats}
    stats = {}
    model = search(g, catalog, max_contractions=5, stats=stats)
    if model is not None:
        return model
    return {'type': 'unresolved', 'search': stats}


def classify(outdir, records):
    path = outdir / 'candidate_outcomes.jsonl'
    previous = read_jsonl(path)
    existing = {r['id'] for r in previous}
    catalog = [(r['name'], nx.from_graph6_bytes(r['graph6'].encode())) for r in common.catalogue_records()]
    targets = dict(catalog)
    positive_cache = PositiveModelCache()
    for r in previous:
        if r['outcome']['type'] == 'proper_obstruction_minor':
            rec = next(x for x in records if x['id'] == r['id'])
            g = nx.Graph()
            g.add_nodes_from(rec['graph']['vertices'])
            g.add_edges_from(rec['graph']['edges'])
            positive_cache.add(g, targets[r['outcome']['target']], r['outcome']['branch_sets'],
                               metadata={'target': r['outcome']['target']})
    counts = Counter(r['outcome']['type'] for r in previous)
    started = last = time.monotonic()
    with path.open('a') as output:
        for record in records:
            if record['direction'] != 'up' or record['id'] in existing:
                continue
            g = nx.Graph()
            g.add_nodes_from(record['graph']['vertices'])
            g.add_edges_from(record['graph']['edges'])
            outcome = classify_graph(g, catalog, positive_cache)
            if outcome['type'] == 'proper_obstruction_minor':
                positive_cache.add(g, targets[outcome['target']], outcome['branch_sets'],
                                   metadata={'target': outcome['target']})
            result = {'id': record['id'], 'outcome': outcome}
            output.write(json.dumps(result, sort_keys=True, separators=(',', ':')) + '\n')
            output.flush()
            previous.append(result)
            existing.add(record['id'])
            counts[outcome['type']] += 1
            if time.monotonic() - last > 10:
                print(json.dumps({'stage': 'classification', 'classified': len(existing),
                                  'outcomes': dict(counts), 'seconds': round(time.monotonic() - started, 1)}),
                      flush=True)
                last = time.monotonic()
    return dict(counts), [r['id'] for r in previous if r['outcome']['type'] == 'unresolved']


if __name__ == '__main__':
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    records, cst = cover(out)
    print('covering:', cst, round(time.time() - t0), 's', flush=True)
    outcomes, unresolved = classify(out, records)
    summary = {'item': 'E4', 'where': 'Theorem 4.23', 'covering': cst, 'outcomes': outcomes,
               'unresolved': unresolved, 'seconds': round(time.time() - t0, 1)}
    (out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))
