"""Shared helpers for the genus computations.

Graphs are dictionaries mapping each vertex to the set of its neighbours.
Graphs read from graph6 have the vertices 0, 1, ..., n-1 in graph6 order.

PAGE and multi_genus are the programs bundled with the Python package
graph-genus 0.1.3 (https://github.com/SanderGi/Genus). Their sources in that
release are PAGE/page.c and MultiGenus/multi_genus.c of the repository.
"""
import gzip
import json
import os
import platform
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------- the programs

def tool(name):
    """Path of a program bundled with graph-genus (or given by an environment
    variable MULTI_GENUS_BIN or PAGE_BIN)."""
    env = {'multi_genus': 'MULTI_GENUS_BIN', 'page': 'PAGE_BIN'}.get(name)
    if env and os.environ.get(env):
        return os.environ[env]
    import graph_genus
    exe = name + ('.exe' if platform.system() == 'Windows' else '')
    p = Path(graph_genus.__file__).resolve().parent / 'bin' / exe
    if not p.is_file():
        raise SystemExit('Could not find %s inside the graph-genus package. Install it with '
                         '"pip install graph-genus==0.1.3".' % exe)
    return str(p)


# ---------------------------------------------------------------- graphs

def from_graph6(s):
    s = s.strip()
    data = [ord(c) - 63 for c in s]
    if data[0] >= 63:
        raise ValueError('graphs with more than 62 vertices are not used here')
    n = data[0]
    bits = []
    for x in data[1:]:
        bits.extend((x >> k) & 1 for k in range(5, -1, -1))
    adj = {v: set() for v in range(n)}
    k = 0
    for j in range(1, n):
        for i in range(j):
            if bits[k]:
                adj[i].add(j)
                adj[j].add(i)
            k += 1
    return adj


def to_graph6(adj, order=None):
    order = sorted(adj) if order is None else list(order)
    n = len(order)
    assert n < 63
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if order[i] in adj[order[j]] else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [n + 63]
    for k in range(0, len(bits), 6):
        x = 0
        for b in bits[k:k + 6]:
            x = (x << 1) | b
        out.append(x + 63)
    return ''.join(chr(c) for c in out)


def edges(adj):
    return sorted((u, v) for u in adj for v in adj[u] if u < v)


def n_edges(adj):
    return sum(len(s) for s in adj.values()) // 2


def components(adj):
    seen, comps = set(), []
    for s in sorted(adj):
        if s in seen:
            continue
        comp, stack = [s], [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    comp.append(y)
                    stack.append(y)
        comps.append(sorted(comp))
    return comps


def induced(adj, vs):
    vs = set(vs)
    return {v: adj[v] & vs for v in vs}


def multicode(adj):
    """multicode of a graph with at least two and fewer than 255 vertices."""
    vs = sorted(adj)
    idx = {v: i + 1 for i, v in enumerate(vs)}
    assert 2 <= len(vs) < 255
    out = [len(vs)]
    for v in vs[:-1]:
        out += sorted(idx[u] for u in adj[v] if idx[u] > idx[v])
        out.append(0)
    return bytes(out)


def zero_indexed(adj):
    """The adjacency list format of graph_genus.embed."""
    vs = sorted(adj)
    idx = {v: i for i, v in enumerate(vs)}
    return [[idx[u] for u in sorted(adj[v])] for v in vs]


# ---------------------------------------------------------------- embeddings

def genus_of_rotation(rotation):
    """Orientable genus of a rotation system given as a list of cyclic
    neighbour lists (vertex i has the list rotation[i]); faces are traced
    and Euler's formula gives the genus. Raises ValueError if the lists do
    not form a rotation system of a connected graph."""
    n = len(rotation)
    succ = {}
    darts = 0
    for v, nb in enumerate(rotation):
        for i, u in enumerate(nb):
            succ[(u, v)] = (v, nb[(i + 1) % len(nb)])
            darts += 1
    for (u, v) in list(succ):
        if (v, u) not in succ:
            raise ValueError('edge listed at one end only')
    seen, faces = set(), 0
    for d in succ:
        if d in seen:
            continue
        faces += 1
        x = d
        while x not in seen:
            seen.add(x)
            x = succ[x]
    e = darts // 2
    two_g = 2 - (n - e + faces)
    if two_g < 0 or two_g % 2:
        raise ValueError('not an orientable embedding of a connected graph')
    return two_g // 2


def decode_rotations(data):
    """Split the output of multi_genus with the options w and L (a byte 255
    before each graph, then an embedding in multicode) into rotation lists."""
    chunks = data.split(b'\xff')
    assert chunks[0] == b''
    out = []
    for c in chunks[1:]:
        if not c:
            out.append(None)
            continue
        n = c[0]
        pos, rot = 1, []
        for _ in range(n):
            nb = []
            while c[pos] != 0:
                nb.append(c[pos] - 1)
                pos += 1
            pos += 1
            rot.append(nb)
        assert pos == len(c)
        out.append(rot)
    return out


# ---------------------------------------------------------------- multi_genus

def multi_genus_batch(graphs):
    """Genus of each connected graph by multi_genus (option w writes a
    minimum genus embedding of each graph, L separates them). Returns the
    list of genera (each the genus of the embedding written, checked by
    tracing its faces) and the running time of multi_genus in seconds."""
    import time
    data = b''.join(multicode(G) for G in graphs)
    t = time.perf_counter()
    p = subprocess.run([tool('multi_genus'), 'w', 'L'], input=data, capture_output=True, check=True)
    dt = time.perf_counter() - t
    rots = decode_rotations(p.stdout)
    assert len(rots) == len(graphs), 'multi_genus wrote %d embeddings for %d graphs' % (len(rots), len(graphs))
    out = []
    for G, rot in zip(graphs, rots):
        assert rot is not None
        vs = sorted(G)
        assert [sorted(vs[u] for u in nb) for nb in rot] == [sorted(G[v]) for v in vs], 'embedding of another graph'
        out.append(genus_of_rotation(rot))
    return out, dt


# ---------------------------------------------------------------- files

def read_list(path):
    """Rows of a whitespace-separated table, without the comment lines."""
    with open(path) as f:
        return [line.split() for line in f if line.strip() and not line.startswith('#')]


def load_json(path):
    path = str(path)
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as f:
        return json.load(f)


def save_json(path, obj, indent=None):
    path = str(path)
    text = json.dumps(obj, indent=indent)
    if path.endswith('.gz'):
        with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as f:
            f.write(text.encode())
    else:
        with open(path, 'w') as f:
            f.write(text + '\n')
