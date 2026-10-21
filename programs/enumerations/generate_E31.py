"""E31: a bridge on {b, c, d} with all two-edge completions nonplanar (Lemma 6.17).

The program's names: p1, p2 are the twins on abc and q1, q2 those on ade.
As in the proof, G = Psi_3 + ce + B with B on {b, c, d} and B + bc
nonplanar, so B is a part of Construction 3.16 with attachments b, c, d
and the added edge bc. The symmetric case, G = Psi_3 + bd + B with B on
{c, d, e} and the added edge de, is run as well. Every graph is
classified, first by looking for a proper nontoroidal minor (multi_genus)
and then by enumeration.classify. Each catalogue graph found is checked
against the conditions G satisfies here, with certificates: a vertex of
degree three with two adjacent neighbours (Lemma 3.3(a)), two twin triples
meeting other than once, a choice of twin pairs leaving a component with
four or five neighbours (an earlier row of Theorem 6.1), and the planarity
of the three two-edge completions of B. The proof finds no obstruction
here: the one copy of T41 has B + bd + cd planar.

Run as: python generate_E31.py FOLDER
"""
import json
import sys

import enumeration as el

import psi3_runs as ext

CASES = [("U_on_bcd", ("b", "c", "d"), [("b", "c")], [("c", "e")],
          [("U+bc+bd", [("b", "c"), ("b", "d")]), ("U+bc+cd", [("b", "c"), ("c", "d")]),
           ("U+bd+cd", [("b", "d"), ("c", "d")])]),
         ("U_on_cde", ("c", "d", "e"), [("d", "e")], [("b", "d")],
          [("U+de+ce", [("d", "e"), ("c", "e")]), ("U+de+cd", [("d", "e"), ("c", "d")]),
           ("U+ce+cd", [("c", "e"), ("c", "d")])])]


def case(name, R, D, cross, tests):
    spec = el.Spec(R, c_edges=D, name=name)

    def assemble(P):
        G = ext.psi3()
        G.add_edges_from(cross)
        G.add_edges_from(P.edges())
        return G
    return (name, spec, assemble, {"min_interior": 1})


if __name__ == "__main__":
    out = sys.argv[1]
    s = ext.run_item_minor_first("E31", [case(*c) for c in CASES], [], out,
                                 note="Psi3 + ce + part U of Construction 3.16 with R = bcd, D = bc; "
                                      "and the mirror Psi3 + bd + U with R = cde, D = de")
    found = ext.annotate(out, {c[0]: (c[1], c[4]) for c in CASES}, want_row=3)
    s["catalogue graphs and their structural exclusions"] = found
    s["agrees with the text"] = s["complete"] and all(
        (x["structural exclusion"] or x["planar hypothesis tests"]) for xs in found.values() for x in xs)
    json.dump(s, open(out + "/summary.json", "w"), indent=1)
    print(json.dumps({k: s[k] for k in ("candidates", "outcomes", "catalogue graphs found",
                                         "catalogue graphs and their structural exclusions",
                                         "agrees with the text", "seconds")}, indent=1))
    for k, v in s["cases"].items():
        print(k, v["recipes"], v["parts up to rooted isomorphism"], v["outcomes"])
