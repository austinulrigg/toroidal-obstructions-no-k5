"""E5: bridges of type (c) with at least two interior vertices (Proposition 4.25).

The names are those of the proof: a = s1, b = t1, c = s2, d = t3,
e = s3, f = t2, with twins p, q on {a, c, e} and r, z on {b, d, f}. The
bridge B contains P1 and has attachments a, b, c, d; P2 = cf and P3 = ed are
edges, so H_simp is the cycle c f d e c. The rest of G is H_ext, one of the
six sets of the proof (up to interchanging a, b and c, d): empty, {ad},
{h_a}, {ad, bc}, {ad, h_b}, {h_a, h_b}, where N(h_a) = acd and
N(h_b) = bcd (the program writes hA, hB).

B is a part of Construction 3.18 with attachments a, b, c, d and added
edges ac, bd and H_ext. The retained K3,3 subdivision must use h_a, h_b
with their edges to c and d, and ad, bc when present; it must use bd when
neither h_b nor bc occurs, and ac when neither h_a nor ad occurs. B has at
least two interior vertices. Every resulting graph G is classified
(enumeration.classify). The proof finds T21 and T63 only.

Run as: python generate_E5.py FOLDER
"""
import sys
import networkx as nx
import enumeration as el

R = ("a", "b", "c", "d")
HA, HB = ("a", "c", "d"), ("b", "c", "d")
ZS = [("Z0", [], {}, [("b", "d"), ("a", "c")]),
      ("Z_ad", [("a", "d")], {}, [("a", "d"), ("b", "d")]),
      ("Z_hA", [], {"hA": HA}, [("hA", "c"), ("hA", "d"), ("b", "d")]),
      ("Z_ad_bc", [("a", "d"), ("b", "c")], {}, [("a", "d"), ("b", "c")]),
      ("Z_ad_hB", [("a", "d")], {"hB": HB}, [("a", "d"), ("hB", "c"), ("hB", "d")]),
      ("Z_hA_hB", [], {"hA": HA, "hB": HB}, [("hA", "c"), ("hA", "d"), ("hB", "c"), ("hB", "d")])]


def exterior():
    G = nx.Graph()
    G.add_edges_from((v, s) for v in "pq" for s in "ace")   # S = {a, c, e}
    G.add_edges_from((v, t) for v in "rz" for t in "bfd")   # T = {b, f, d}
    G.add_edges_from([("c", "f"), ("e", "d")])              # R2 = s2t2, R3 = s3t3
    return G


def case(name, zedges, outside, req):
    spec = el.Spec(R, outside=outside, c_edges=[("a", "c"), ("b", "d")] + zedges, required=req)

    def assemble(P):
        G = exterior()
        G.add_edges_from(P.edges())
        G.add_edges_from(zedges)
        for h, nb in outside.items():
            G.add_edges_from((h, x) for x in nb)
        return G
    return (name, spec, assemble, {"min_interior": 2})


if __name__ == "__main__":
    out = sys.argv[1]
    s = el.run_item("E5", [case(*z) for z in ZS], ["T21", "T63"], out,
                    note="a=s1, b=t1, c=s2, d=t3, e=s3, f=t2; p,q on ace; r,z on bfd")
    import json
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["parts up to rooted isomorphism"], v["outcomes"])
