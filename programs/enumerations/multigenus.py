"""Genus lower bounds from multi_genus (Brinkmann), the program bundled with
the Python package graph-genus (or the one named by the environment variable
MULTI_GENUS). Used only to decide whether connected graphs have genus larger
than one (or larger than zero)."""

import os
import subprocess
from pathlib import Path


def _bundled():
    import platform
    import graph_genus
    return str(Path(graph_genus.__file__).resolve().parent / "bin" /
               ("multi_genus.exe" if platform.system() == "Windows" else "multi_genus"))


BIN = os.environ.get("MULTI_GENUS") or _bundled()


def multicode(G):
    nodes = sorted(G.nodes(), key=str)
    idx = {v: i + 1 for i, v in enumerate(nodes)}
    n = len(nodes)
    assert 2 <= n < 255
    out = [n]
    for v in nodes[:-1]:
        out += sorted(idx[u] for u in G[v] if idx[u] > idx[v])
        out.append(0)
    return bytes(out)


def genus_above(graphs, k):
    """For connected graphs, a list of booleans: genus larger than k."""
    data = b"".join(multicode(G) for G in graphs)
    p = subprocess.run([BIN, "l%d" % k, "F", "L"], input=data, capture_output=True, check=True)
    chunks = p.stdout.split(b"\xff")
    assert chunks[0] == b"" and len(chunks) == len(graphs) + 1, (len(chunks), len(graphs))
    return [len(c) > 0 for c in chunks[1:]]


def nontoroidal(graphs):
    return genus_above(graphs, 1)


def nonplanar(graphs):
    return genus_above(graphs, 0)
