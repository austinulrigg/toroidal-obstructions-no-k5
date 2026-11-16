"""Generate R12, the rotation table for H2 in Lemma 4.10.

Usage: python make_R12.py OUTPUT_FOLDER
Writes OUTPUT_FOLDER/rotations/R12.json. Uses only the standard library.
"""
import itertools
import json
import math
import sys
from collections import Counter
from pathlib import Path


HOW_TO_READ = (
    "For each graph: 'cyclic orders' lists, for every vertex, all cyclic orders of its neighbours, each starting "
    "with its alphabetically first neighbour (the other neighbours follow in every possible order, in lexicographic "
    "order). Rotation system number i (counting from 0) takes at each vertex the cyclic order whose position is the "
    "corresponding digit of i, written in the mixed radix given by the numbers of cyclic orders in the order of "
    "'vertices', the last vertex varying fastest (the order of Python's itertools.product). Facial walks are traced "
    "by the rule used throughout the supplement: after the dart u->v comes v->w, where w follows u in the cyclic "
    "order at v. The genus is (2 - V + E - F)/2. The attachment sequence of a facial walk lists the attachments it "
    "passes through, in walk order, one entry per passage; it is written from its lexicographically least cyclic "
    "shift, keeping the direction of the walk, and the faces of a row are sorted. Since many rotation systems give "
    "the same row, the distinct rows are listed once under 'distinct rows' and 'rows' gives, for every rotation "
    "system in order, the position of its row in that list. 'Up to cyclic shift and reversal' in the claims means "
    "the least sequence among all cyclic shifts of the sequence and of its reverse.")

PAPER_TEXT = r"""We consider all possible embeddings of $H_2$ and a closed regular neighborhood $\mathcal N$ of it. Enumerating the rotations of $H_2$ and tracing their boundary walks gives the following~\supp{R12}. If $\mathcal N$ has genus one, no boundary walk contains all four attachments, and every complementary region is then a disk by Lemma~\ref{lem:torus-complement}, a contradiction.

If $\mathcal N$ has genus zero, no boundary walk contains all four attachments, and exactly one pair contains all four between them. Their attachment sequences are $s_1,s_2$ and $t_1,t_2$, with each attachment appearing once. Lemma~\ref{lem:torus-complement} therefore puts $H_1'-A$ in the annulus bounded by these walks, with all four attachments on its boundary. By Lemma~\ref{lem:annulus-embedding}, $H_1'+E$ is planar, a contradiction."""


def canonical(sequence):
    sequence = tuple(sequence)
    return min(sequence[i:] + sequence[:i] for i in range(len(sequence))) if sequence else ()


def word(sequence):
    return ' '.join(sequence) if sequence else '-'


def generate():
    vertices = ['s1', 't1', 's2', 't2', 's3', 't3', 'p', 'q', 'r', 'z']
    attachments = ['s1', 't1', 's2', 't2']
    roots = set(attachments)
    edges = sorted({tuple(sorted((twin, root)))
                    for twins, triple in [(['p', 'q'], ['s1', 's2', 's3']),
                                          (['r', 'z'], ['t1', 't2', 't3'])]
                    for twin in twins for root in triple} | {('s3', 't3')})
    adjacency = {v: sorted(w if v == u else u for u, w in edges if v in (u, w))
                 for v in vertices}
    cyclic = {v: [[nb[0], *permutation] for permutation in itertools.permutations(nb[1:])]
              for v, nb in adjacency.items()}
    count = math.prod(math.factorial(len(adjacency[v]) - 1) for v in vertices)
    rows, distinct, row_index = [], [], {}
    genera, face_counts = Counter(), {}
    g1_full, g0_full = 0, 0
    pair_counts, pair_sequences, once_counts, signatures = Counter(), Counter(), Counter(), Counter()
    for choices in itertools.product(*(cyclic[v] for v in vertices)):
        rotation = dict(zip(vertices, choices))
        successor = {(u, v): (v, rotation[v][(rotation[v].index(u) + 1) % len(rotation[v])])
                     for v in vertices for u in adjacency[v]}
        unseen = set(successor)
        walks = []
        while unseen:
            first = current = min(unseen)
            walk = []
            while True:
                unseen.remove(current)
                walk.append(current[0])
                current = successor[current]
                if current == first:
                    break
            walks.append(walk)
        numerator = 2 - len(vertices) + len(edges) - len(walks)
        assert numerator >= 0 and numerator % 2 == 0
        genus = numerator // 2
        sequences = [[v for v in walk if v in roots] for walk in walks]
        row = (genus, tuple(sorted(word(canonical(s)) for s in sequences)))
        if row not in row_index:
            row_index[row] = len(distinct)
            distinct.append([genus, list(row[1])])
        rows.append(row_index[row])
        genera[genus] += 1
        face_counts.setdefault(genus, Counter())[len(walks)] += 1
        if genus == 1:
            g1_full += sum(roots <= set(s) for s in sequences)
        if genus == 0:
            g0_full += sum(roots <= set(s) for s in sequences)
            pairs = [(s, t) for s, t in itertools.combinations(sequences, 2) if roots <= set(s) | set(t)]
            pair_counts[len(pairs)] += 1
            for s, t in pairs:
                pair_sequences[' | '.join(sorted([word(canonical(s)), word(canonical(t))]))] += 1
                once_counts[all(s.count(v) + t.count(v) == 1 for v in attachments)] += 1
            signatures[' / '.join(row[1])] += 1
    assert len(rows) == count == 64
    assert genera == Counter({0: 16, 1: 32, 2: 16})
    assert g1_full == g0_full == 0
    assert pair_counts == Counter({1: 16})
    assert pair_sequences == Counter({'s1 s2 | t1 t2': 16})
    assert once_counts == Counter({True: 16})
    assert len(signatures) == 4 and set(signatures.values()) == {4}
    by_genus = {str(k): genera[k] for k in sorted(genera)}
    graph = dict(name='H2 = (J1 union J2) + s3t3', vertices=vertices, edges=edges,
                 attachments=attachments)
    graph.update({'cyclic orders': cyclic, 'number of rotation systems': count,
                  'rotation systems by genus': by_genus,
                  'row format': '[genus, attachment sequences of the facial walks]',
                  'distinct rows': distinct, 'rows': rows})
    claims = [
        dict(statement='H2 has 10 vertices, 13 edges and 64 rotation systems. A rotation with F facial walks has genus (5-F)/2.',
             computed={'vertices': 10, 'edges': 13, 'rotation systems': 64, 'by genus': by_genus,
                       'facial walk counts by genus': {str(g): {str(k): n for k, n in sorted(c.items())}
                                                      for g, c in sorted(face_counts.items())}}, holds=True),
        dict(statement='In genus one, no boundary walk contains all four attachments s1,t1,s2,t2.',
             computed={'genus-one rotations': genera[1], 'walks containing all four attachments': g1_full}, holds=True),
        dict(statement='In genus zero, no single boundary walk contains all four attachments. Exactly one pair contains all four between them, with sequences s1,s2 and t1,t2, each attachment appearing once.',
             computed={'genus-zero rotations': genera[0], 'walks containing all four attachments': g0_full,
                       'covering pairs per rotation': {str(k): n for k, n in sorted(pair_counts.items())},
                       'covering pair attachment sequences': dict(sorted(pair_sequences.items())),
                       'covering pairs with each attachment appearing once': once_counts[True]}, holds=True),
        dict(statement='In genus zero, the five attachment sequences are {si,tj}, {s1,s2}, {s(3-i)}, {t1,t2}, {t(3-j)} for i,j in {1,2}. Each of the four signatures occurs in four rotations.',
             computed={'boundary signatures': dict(sorted(signatures.items()))}, holds=True),
    ]
    return {'item': 'R12', 'where': 'Lemma 4.10 (dis:sole-two-rails), proof',
            'paper text': [{'label': 'dis:sole-two-rails', 'text': PAPER_TEXT}],
            'what': 'Every rotation system of H2, consisting of twins p,q on {s1,s2,s3}, twins r,z on {t1,t2,t3}, and the edge s3t3. Attachments: s1,t1,s2,t2; twins and s3,t3 omitted from attachment sequences. All rotations, including genus two, are recorded; only genera zero and one can occur for a regular neighbourhood inside the torus.',
            'how to read': HOW_TO_READ, 'graphs': [graph], 'claims': claims}


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python make_R12.py OUTPUT_FOLDER')
    target = Path(sys.argv[1]) / 'rotations' / 'R12.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(generate(), indent=1) + '\n')
    print(target)
