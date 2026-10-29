"""E19: a bridge on a triple whose three two-edge completions are nonplanar (Lemma 5.24).

The program's names: x1, x2 are the twins on abc and y1, y2 those on abd.
As in the proof, B is a bridge on T = {a, c, d} or {b, c, d} whose three
two-edge completions are nonplanar (the program names the bridge U), and
another bridge on one of the three other triples is contracted to a vertex
v. So G = Psi_2 + U + v with N(v) one of the three other triples: six
cases. The bridges U are listed by Construction 3.17 with attachments T:

  nonplanar U: no added edges;
  planar U:    one added edge xy, which the retained K3,3 subdivision must
               use, taken for each pair xy of T in turn; a part is kept
               when it is planar and its three two-edge completions are
               nonplanar, and the parts kept for the three pairs are
               compared up to isomorphism fixing T.

Every resulting graph is classified (enumeration.classify). The random
search for a K5 minor uses, for the k-th bridge (from 0) of the j-th
nonplanar case (from 0, in the order the program lists them), the seed
1753 + 876 j + k, with the vertices of the graph in the order of its
construction, and for the k-th form of the j-th planar case the seed
1 + 120 j + k, with the vertices in the order of its edges. These are the
settings that give the data folder; any settings give valid records. The proof finds a cut of
at most two vertices, a K5 minor or a proper nontoroidal minor for each of
the 5,934 graphs.

Run as: python generate_E19.py FOLDER
"""
import itertools
import json
import sys

import networkx as nx

import enumeration as el

TRIPLES = ("abc", "abd", "acd", "bcd")
ON = ("acd", "bcd")


def psi2():
    G = nx.Graph()
    for t, nb in {"x1": "abc", "x2": "abc", "y1": "abd", "y2": "abd"}.items():
        G.add_edges_from((t, v) for v in nb)
    return G


def two_edge_completions(P, T):
    x, y, t = T
    out = []
    for E in (((x, y), (x, t)), ((x, y), (y, t)), ((x, t), (y, t))):
        H = P.copy()
        H.add_edges_from(E)
        out.append(H)
    return out


def planar_forms(T):
    """The planar bridges on T with three nonplanar two-edge completions,
    up to isomorphism fixing T, with the counts for each added edge."""
    forms, buckets, stats = [], {}, {}
    for x, y in itertools.combinations(T, 2):
        spec = el.Spec(tuple(T), c_edges=[(x, y)], required=[(x, y)])
        parts, st = el.generate_parts(spec)
        kept = 0
        for P, recipe in parts:
            if not nx.check_planarity(P)[0]:
                continue
            if any(nx.check_planarity(H)[0] for H in two_edge_completions(P, T)):
                continue
            kept += 1
            key = el.rooted_key(P, tuple(T))
            if any(el.rooted_isomorphic(P, forms[i][0], tuple(T)) for i in buckets.get(key, [])):
                continue
            buckets.setdefault(key, []).append(len(forms))
            forms.append((P, dict(recipe, **{"added edge": [x, y]})))
        stats[x + y] = {"parts": len(parts), "planar with all three two-edge completions nonplanar": kept}
    return forms, stats


def make_assemble(S):
    def assemble(P):
        G = psi2()
        G.add_edges_from(P.edges())
        G.add_edges_from(("v", s) for s in S)
        return G
    return assemble


if __name__ == "__main__":
    out = sys.argv[1]
    nonplanar, planar, forms_info = {}, {}, {}
    for T in ON:
        spec = el.Spec(tuple(T))
        nonplanar[T] = el.generate_parts(spec) + (spec.payload(),)
        forms, stats = planar_forms(T)
        planar[T] = (forms, {"forms": len(forms)}, {"roots": list(T), "U": "planar, one added edge xy per pair"})
        forms_info[T] = {"count": len(forms), "by added edge": stats,
                         "forms": [{"index": i, "graph": el.graph_record(P), "recipe": r}
                                   for i, (P, r) in enumerate(forms)]}
    cases = []
    for T in ON:
        parts, st, payload = nonplanar[T]
        for S in TRIPLES:
            if S != T:
                cases.append(("U_%s_nonplanar_v_on_%s" % (T, S), None, make_assemble(S),
                              {"parts": parts, "stats": st, "payload": dict(payload, exterior="v on " + S)}))
    for T in ON:
        forms, st, payload = planar[T]
        for S in TRIPLES:
            if S != T:
                cases.append(("U_%s_planar_v_on_%s" % (T, S), None, make_assemble(S),
                              {"parts": forms, "stats": st, "payload": dict(payload, exterior="v on " + S)}))
    order = [c[0] for c in cases]

    def seed_of(name, k, n):
        j = order.index(name)
        if "_nonplanar_" in name:
            return 1753 + 876 * j + k
        return 1 + 120 * (j - 6) + k
    s = el.run_item("E19", cases, [], out, seed_of=seed_of, nodes_first=lambda name: "_nonplanar_" in name,
                    note="x1,x2 on abc; y1,y2 on abd; U on acd or bcd; v on another triple")
    with open(out + "/forms.json", "w") as f:
        json.dump(forms_info, f, indent=1)
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v.get("parts up to rooted isomorphism", v.get("forms")), v["outcomes"])
