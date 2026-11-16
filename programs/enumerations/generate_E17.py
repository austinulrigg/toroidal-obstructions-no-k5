"""E17: every component on a triple of {a, b, c, d} is a single vertex (Lemma 5.28).

The program names the vertices with neighbourhood T by T itself: abc0,
abc1, ... on abc, abd0, ... on abd, acd0, ... on acd and bcd0, ... on bcd.
As in the proof, G is determined by the numbers of vertices with each of
the neighbourhoods {a,b,c}, {a,b,d}, {a,c,d}, {b,c,d}: at least two on each
of the first two, at most nine in total, and cd absent or present. These
are the 252 graphs of the lemma. Every count vector is taken (labeled), and
each candidate records its class up to interchanging a, b and c, d. Every
graph is classified (enumeration.classify). The proof finds T9, T10, T16.

Run as: python generate_E17.py FOLDER
"""
import itertools
import json
import sys
import time
from pathlib import Path

import networkx as nx

import enumeration as el
import psi2_runs as ext

TRIPLES = ("abc", "abd", "acd", "bcd")


def graph(counts, cd):
    G = nx.Graph()
    G.add_nodes_from("abcd")
    for t, n in zip(TRIPLES, counts):
        for i in range(n):
            G.add_edges_from(("%s%d" % (t, i), x) for x in t)
    if cd:
        G.add_edge("c", "d")
    return G


def canonical(counts, cd):
    """Representative of the orbit under a<->b and c<->d."""
    n1, n2, n3, n4 = counts
    orbit = [(n1, n2, n3, n4), (n1, n2, n4, n3), (n2, n1, n3, n4), (n2, n1, n4, n3)]
    return (min(orbit), cd)


if __name__ == "__main__":
    out = Path(sys.argv[1])
    t0 = time.time()
    parts = []
    for n1, n2, n3, n4 in itertools.product(range(2, 10), range(2, 10), range(8), range(8)):
        if n1 + n2 + n3 + n4 > 9:
            continue
        for cd in (False, True):
            P = graph((n1, n2, n3, n4), cd)
            parts.append((P, {"counts on abc, abd, acd, bcd": [n1, n2, n3, n4], "cd": cd,
                              "class up to a<->b, c<->d": [list(canonical((n1, n2, n3, n4), cd)[0]), cd]}))
    st = {"count vectors with cd absent or present": len(parts)}
    pay = {"vertices": "a,b,c,d; vertex 'Ti' is the i-th singleton on triple T",
           "domain": "n_abc>=2, n_abd>=2, n_acd>=0, n_bcd>=0, total<=9; cd absent or present"}
    s = ext.run_whole_item("E17", [("singletons", parts, st, pay, lambda P: P)], ["T9", "T10", "T16"], out,
                           structural=[], structural_payload={},
                           note="singleton components on the triples of {a,b,c,d}")
    d = json.load(open(out / "candidates.json"))
    classes = {}
    for r in d["candidates"]:
        o = r["outcome"]
        if o["outcome"] == "catalogue graph":
            key = tuple(map(tuple, [r["recipe"]["class up to a<->b, c<->d"][0]])) + (r["recipe"]["cd"],)
            classes.setdefault(o["catalogue"], set()).add(
                (tuple(r["recipe"]["counts on abc, abd, acd, bcd"]), r["recipe"]["cd"]))
    s["catalogue graphs with their count vectors"] = {k: sorted([list(c), cd] for c, cd in v)
                                                     for k, v in classes.items()}
    json.dump(s, open(out / "summary.json", "w"), indent=1)
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "survivors",
                                         "catalogue graphs with their count vectors",
                                         "agrees with the text", "seconds")}, indent=1))
