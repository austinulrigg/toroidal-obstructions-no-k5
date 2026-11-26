"""E14: the case D = {d} with edge dw in an earlier proof of Proposition 5.27.

Here p, q are the twins on {a, c, w}, t0, t1 those on {a, b, c} and t2, t3
those on {a, b, d}. The bridge U has attachments b, c, w. Contracting cw
in U + cw + bc gives a nonplanar graph with b adjacent to the merged vertex r.

Construction 3.18 on b, r gives the candidates with br in the retained K3,3
subdivision. Assign each edge at r back to c or w. If either vertex has no
interior neighbor, try each possible neighbor there. Compare the bridges up
to isomorphism fixing b, c, w. For each one, classify the graph obtained by
adding Psi_2, the twins on {a, c, w} and dw. Record whether U + cw is planar
and both U + cw + bw and U + cw + bc are nonplanar.

    python generate_E14.py FOLDER"""
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
    s = el.run_item("E14", [("lifts", None, assemble, {"parts": ps, "stats": st,
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
