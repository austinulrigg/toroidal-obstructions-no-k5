"""E15: a cut {a, c, w} with b and d in different components, D missing a (Proposition 5.21).

The program's names: p, q are the twins on {a, c, w} (p', q' in the paper),
t0, t1 the twins on abc and t2, t3 those on abd; U is the paper's H_b, on
the attachments b, c, w, and r is the paper's x, the vertex formed by
contracting cw. As in the proof, D = {d} with the edge dw, U + cw is
planar, and Lemma 2.7 makes (U + cw + bc)/cw nonplanar.

Construction 3.17 with attachments b, r and the added edge br, which the
retained K3,3 subdivision must use. Each edge at r is assigned back to c
or to w in every possible way, and c or w without an interior neighbour
gets one, in every possible way. The bridges are compared up to
isomorphism fixing b, c, w. G = Psi_2 + p, q on acw + U + dw, and every G
is classified (enumeration.classify). Each candidate also records whether
U + cw is planar and U + cw + bw, U + cw + bc are nonplanar, the case of
the proof. The proof finds no obstruction here.

Run as: python generate_E15.py FOLDER
"""
import itertools
import json
import sys

import networkx as nx
import enumeration as el

R = ("b", "c", "w")


def lifts(Q):
    inner = sorted(v for v in Q if v not in ("b", "r"))
    at_r = sorted(v for v in Q["r"] if v in inner)
    base = Q.copy()
    base.remove_node("r")
    base.add_nodes_from(R)
    for choice in itertools.product("cw", repeat=len(at_r)):
        P = base.copy()
        P.add_edges_from((x, v) for x, v in zip(choice, at_r))
        lonely = [x for x in "cw" if not any(u in inner for u in P[x])]
        for nbrs in itertools.product(inner, repeat=len(lonely)):
            H = P.copy()
            H.add_edges_from(zip(lonely, nbrs))
            yield H


def parts():
    spec = el.Spec(("b", "r"), c_edges=[("b", "r")], required=[("b", "r")])
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


def assemble(U):
    G = nx.Graph()
    G.add_edges_from((t, x) for t in ("t0", "t1") for x in "abc")
    G.add_edges_from((t, x) for t in ("t2", "t3") for x in "abd")
    G.add_edges_from((p, x) for p in "pq" for x in "acw")
    G.add_edges_from(U.edges())
    G.add_edge("d", "w")
    return G


def case_flags(U):
    out = {}
    for name, E in (("U+cw", [("c", "w")]), ("U+cw+bw", [("c", "w"), ("b", "w")]),
                    ("U+cw+bc", [("c", "w"), ("b", "c")])):
        H = U.copy()
        H.add_edges_from(E)
        out[name + " planar"] = nx.check_planarity(H)[0]
    out["in the item's case"] = out["U+cw planar"] and not out["U+cw+bw planar"] and not out["U+cw+bc planar"]
    return out


if __name__ == "__main__":
    out = sys.argv[1]
    ps, st = parts()
    s = el.run_item("E15", [("lifts", None, assemble, {"parts": ps, "stats": st,
                                                "payload": {"quotient roots": ["b", "r"], "added edges": [["b", "r"]],
                                                            "required edges": [["b", "r"]], "lifted roots": list(R)}})],
                    [], out, note="G = Psi2 + p, q on acw + U + dw")
    data = json.load(open(out + "/candidates.json"))
    for r in data["candidates"]:
        U = nx.Graph()
        U.add_nodes_from(r["part"]["vertices"])
        U.add_edges_from(r["part"]["edges"])
        r["case flags"] = case_flags(U)
    json.dump(data, open(out + "/candidates.json", "w"))
    inc = [r for r in data["candidates"] if r["case flags"]["in the item's case"]]
    s["candidates in the item's case"] = len(inc)
    s["outcomes in the item's case"] = dict(__import__("collections").Counter(r["outcome"]["outcome"] for r in inc))
    json.dump(s, open(out + "/summary.json", "w"), indent=1)
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "candidates in the item's case", "outcomes in the item's case",
                                         "seconds")}, indent=1))
