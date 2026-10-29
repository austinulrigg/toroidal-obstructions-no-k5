"""E10: a bridge with four attachments, (n_c, n_d) = (3,2) and H = B + ab + ac + bc planar (Lemma 5.9).

As in the proof, both H + ad and H + bd are nonplanar, and Lemma 2.7 makes
(H + ad + bd)/ab nonplanar; write x for the vertex formed from a and b (the
program calls it r). Construction 3.17 with attachments x, c, d and added
edges xc, xd lists the reduced bridges before x is separated. For each,
every edge xv is replaced by av or bv in every possible way, and a or b
without an interior neighbour gets one, again in every possible way. The
bridges are compared up to isomorphism fixing a, b, c, d. Each is
completed with three twins on abc and two on abd (the program calls them
p0, p1, p2 and r0, r1), and every resulting graph is classified
(enumeration.classify). The proof finds no obstruction here; T11, which
Proposition 5.10 gives for (n_c, n_d) = (3,2), comes from case B of E9.

Run as: python generate_E10.py FOLDER
"""
import itertools
import json
import sys

import networkx as nx
import enumeration as el

R = ("a", "b", "c", "d")


def lifts(Q):
    """Every lift of the part Q on roots r, c, d to roots a, b, c, d."""
    inner = sorted(v for v in Q if v not in ("r", "c", "d"))
    at_r = sorted(v for v in Q["r"] if v in inner)
    base = Q.copy()
    base.remove_node("r")
    base.add_nodes_from(R)
    for choice in itertools.product("ab", repeat=len(at_r)):
        P = base.copy()
        P.add_edges_from((x, v) for x, v in zip(choice, at_r))
        lonely = [x for x in "ab" if not any(u in inner for u in P[x])]
        for nbrs in itertools.product(inner, repeat=len(lonely)):
            H = P.copy()
            H.add_edges_from(zip(lonely, nbrs))
            yield H


def parts():
    spec = el.Spec(("r", "c", "d"), c_edges=[("r", "c"), ("r", "d")])
    quotients, st = el.generate_parts(spec, min_interior=1)
    buckets, out = {}, []
    lifted = 0
    for Q, recipe in quotients:
        for P in lifts(Q):
            lifted += 1
            key = el.rooted_key(P, R)
            if any(el.rooted_isomorphic(P, out[j][0], R) for j in buckets.get(key, [])):
                continue
            buckets.setdefault(key, []).append(len(out))
            out.append((el.relabel_interior(P, R), {"quotient part": el.graph_record(Q), "quotient recipe": recipe}))
    st.update({"quotient parts": len(quotients), "lifts": lifted, "parts up to rooted isomorphism": len(out)})
    return out, st


def assemble(P):
    G = nx.Graph()
    G.add_edges_from(("p%d" % i, x) for i in range(3) for x in "abc")
    G.add_edges_from(("r%d" % i, x) for i in range(2) for x in "abd")
    G.add_edges_from(P.edges())
    return G


if __name__ == "__main__":
    out = sys.argv[1]
    ps, st = parts()
    s = el.run_item("E10", [("lifts", None, assemble, {"parts": ps, "stats": st,
                                                "payload": {"quotient roots": ["r", "c", "d"],
                                                            "added edges": [["r", "c"], ["r", "d"]],
                                                            "lifted roots": list(R)}})],
                    [], out, note="G = lifted P + three twins on abc + two on abd")
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found", "seconds")}, indent=1))
    print(s["cases"])
