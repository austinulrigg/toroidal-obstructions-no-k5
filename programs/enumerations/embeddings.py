"""Torus embeddings (rotation systems of genus at most one) found with
multi_genus, used by E4 for its toroidal graphs.

multi_genus (the program bundled with the package graph-genus, or the one
named by the environment variable MULTI_GENUS) is called with the arguments
l1 w, so that it writes an embedding of genus at most one of each block. The
rotations are accepted only after checking that they list exactly the
neighbours of each vertex and that Euler's formula gives genus at most one;
the genus that multi_genus reports is not used. The genus is added over the
blocks, and the rotations of the blocks are joined at cut vertices one after
the other.
"""

from collections import Counter
import math
import subprocess
import time

import networkx as nx
from graph_genus import _core

import graphs as gr
import multigenus


def _checked_genus(G, rotation):
    if set(rotation) != set(G):
        raise RuntimeError("the rotation system has other vertices than the graph")
    for v in G:
        if Counter(rotation[v]) != Counter(G.neighbors(v)):
            raise RuntimeError("the rotation at a vertex does not list its neighbours")
    if G.number_of_edges() == 0:
        return 0
    faces = gr.faces_of_rotation(rotation)
    twice = 2 - len(G) + G.number_of_edges() - len(faces)
    if twice < 0 or twice % 2:
        raise RuntimeError("the rotation system gives an impossible Euler characteristic")
    return twice // 2


def _block_rotation(B, remaining):
    verts = sorted(B)
    n = len(verts)
    if n > 128 or B.number_of_edges() > 512:
        raise ValueError("a block too large for multi_genus")
    index = {v: i for i, v in enumerate(verts)}
    adjacency = [sorted(index[w] for w in B[v]) for v in verts]
    # multicode: the order, then for each vertex but the last its larger
    # neighbours (numbered from 1) followed by 0.
    data = [n]
    for i in range(n - 1):
        data.extend(j + 1 for j in adjacency[i] if j > i)
        data.append(0)
    proc = subprocess.run([multigenus.BIN, "l1", "w"], input=bytes(data), stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, timeout=remaining, check=False)
    if proc.returncode:
        raise RuntimeError("multi_genus exited with status %d: %s" % (
            proc.returncode, proc.stderr.decode(errors="replace")[-2000:]))
    if not proc.stdout:
        raise RuntimeError("multi_genus wrote no embedding of genus at most one")
    rot = _core._decode_multicode_rotation(proc.stdout)
    if len(rot) != n:
        raise RuntimeError("the embedding has another order than the block")
    encoded = [n]
    for neighbors in rot:
        encoded.extend(j + 1 for j in neighbors)
        encoded.append(0)
    if bytes(encoded) != proc.stdout:
        raise RuntimeError("multi_genus wrote more than one embedding")
    rotation = {verts[i]: [verts[j] for j in rot[i]] for i in range(n)}
    genus = _checked_genus(B, rotation)
    if genus > 1:
        raise RuntimeError("the embedding of a block has genus larger than one")
    return genus, rotation


def toroidal_embedding(G, seconds=120):
    """(genus, parts) for a graph of genus at most one: for each component
    (in the order of their smallest vertices) its vertices, genus and
    rotation system, as graphs.page_embedding returns them. Raises an
    exception if no embedding is found within the time limit."""
    if seconds is None or not math.isfinite(seconds) or seconds <= 0:
        raise ValueError("seconds must be a positive number")
    if G.is_directed() or G.is_multigraph() or nx.number_of_selfloops(G):
        raise ValueError("a simple graph is required")
    if len({str(v) for v in G}) != len(G):
        raise ValueError("two vertices have the same name as strings")
    deadline = time.monotonic() + seconds
    total, parts = 0, []
    for comp in sorted(nx.connected_components(G), key=min):
        H = G.subgraph(comp)
        verts = sorted(H)
        rotation = {v: [] for v in verts}
        component_genus = 0
        for block in sorted(nx.biconnected_components(H), key=min):
            B = H.subgraph(block)
            if len(B) == 2:
                u, v = sorted(B)
                rotation[u].append(v)
                rotation[v].append(u)
                continue
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise subprocess.TimeoutExpired([multigenus.BIN, "l1", "w"], seconds)
            genus, block_rot = _block_rotation(B, remaining)
            component_genus += genus
            for v in block:
                rotation[v].extend(block_rot[v])
        genus = _checked_genus(H, rotation)
        if genus != component_genus:
            raise RuntimeError("the joined rotations of the blocks have another genus")
        total += genus
        if total > 1:
            raise RuntimeError("the embeddings of the components have total genus larger than one")
        parts.append({"vertices": verts, "genus": genus, "rotation": {str(v): rotation[v] for v in verts}})
    return total, parts
