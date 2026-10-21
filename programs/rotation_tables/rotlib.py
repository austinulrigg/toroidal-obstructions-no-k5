"""Rotation tables: every rotation system of a small graph, its genus, and the
attachment sequence of each facial walk.

Conventions (stated in every table file):
- The neighbours of a vertex are taken in alphabetical order. Its cyclic orders
  are all orders that start with the alphabetically first neighbour, the others
  following in every possible order (lexicographic order of these
  permutations). A vertex of degree at most two has one cyclic order.
- Row i of a table is the rotation system obtained by writing i in the mixed
  radix given by the numbers of cyclic orders of the vertices, in the listed
  vertex order, the last vertex varying fastest (the order of
  itertools.product).
- A facial walk is traced by the rule used throughout the supplement: after
  the dart u->v comes v->w, where w follows u in the cyclic order at v.
- The attachment sequence of a facial walk lists the attachments it passes
  through, in walk order, one entry per passage. A row gives each sequence
  starting at its lexicographically least cyclic shift (keeping the direction
  of the walk) and lists the faces in sorted order; '-' is a facial walk that
  meets no attachment.
"""
import itertools


def cyclic_orders(nbrs):
    nbrs = sorted(nbrs)
    if len(nbrs) <= 2:
        return [nbrs]
    return [[nbrs[0]] + list(p) for p in itertools.permutations(nbrs[1:])]


def least_shift(seq):
    seq = tuple(seq)
    if not seq:
        return seq
    return min(seq[i:] + seq[:i] for i in range(len(seq)))


def canon(seq):
    """Up to cyclic shift and reversal."""
    seq = tuple(seq)
    return min(least_shift(seq), least_shift(seq[::-1]))


def fmt(seq, sep=''):
    return sep.join(seq) if seq else '-'


class RotationTable:
    def __init__(self, name, vertices, edges, attachments, sep=''):
        self.name = name
        self.vertices = list(vertices)
        assert len(set(self.vertices)) == len(self.vertices)
        self.adj = {v: set() for v in self.vertices}
        for u, w in edges:
            assert u != w and w not in self.adj[u], ('loop or parallel edge', u, w)
            self.adj[u].add(w)
            self.adj[w].add(u)
        self.edges = sorted(sorted(e) for e in edges)
        self.attachments = list(attachments)
        self.att = set(attachments)
        assert self.att <= set(self.vertices)
        self.sep = sep
        self.orders = {v: cyclic_orders(self.adj[v]) for v in self.vertices}
        self.radix = [len(self.orders[v]) for v in self.vertices]
        self.count = 1
        for r in self.radix:
            self.count *= r
        # darts
        self.darts = [(u, w) for u in self.vertices for w in sorted(self.adj[u])]
        self.dix = {d: i for i, d in enumerate(self.darts)}
        self.tail = [d[0] for d in self.darts]
        self.V = len(self.vertices)
        self.E = len(self.darts) // 2
        # for each vertex v and each cyclic order k: list of (index of dart u->v, index of next dart v->w)
        self.pieces = {}
        for v in self.vertices:
            pl = []
            for order in self.orders[v]:
                m = []
                for i, u in enumerate(order):
                    w = order[(i + 1) % len(order)]
                    m.append((self.dix[(u, v)], self.dix[(v, w)]))
                pl.append(m)
            self.pieces[v] = pl
        self.connected = self._connected()

    def _connected(self):
        seen = {self.vertices[0]}
        stack = [self.vertices[0]]
        while stack:
            x = stack.pop()
            for y in self.adj[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        return len(seen) == self.V

    def digits(self, i):
        out = []
        for r in reversed(self.radix):
            out.append(i % r)
            i //= r
        return out[::-1]

    def rotation(self, i):
        return {v: self.orders[v][k] for v, k in zip(self.vertices, self.digits(i))}

    def faces_of_digits(self, ds):
        """Faces as lists of dart indices, for the rotation with the given digits."""
        nxt = [0] * len(self.darts)
        for v, k in zip(self.vertices, ds):
            for a, b in self.pieces[v][k]:
                nxt[a] = b
        seen = [False] * len(self.darts)
        faces = []
        for s in range(len(self.darts)):
            if seen[s]:
                continue
            f = []
            d = s
            while not seen[d]:
                seen[d] = True
                f.append(d)
                d = nxt[d]
            assert d == s
            faces.append(f)
        return faces

    def genus(self, faces):
        two_g = 2 - self.V + self.E - len(faces)
        assert self.connected and two_g >= 0 and two_g % 2 == 0
        return two_g // 2

    def seq(self, dart_list):
        return tuple(self.tail[d] for d in dart_list if self.tail[d] in self.att)

    def row(self, faces):
        g = self.genus(faces)
        seqs = sorted(fmt(least_shift(self.seq(f)), self.sep) for f in faces)
        return [g, seqs]

    def all_digits(self):
        return itertools.product(*[range(r) for r in self.radix])

    def region_along(self, faces, s, t):
        """The region formed by the faces on the two sides of the edge st,
        glued along st. Returns ('disk', seq) or ('annulus', seq1, seq2)."""
        ds, dt = self.dix[(s, t)], self.dix[(t, s)]
        where = {}
        for fi, f in enumerate(faces):
            for j, d in enumerate(f):
                where[d] = (fi, j)
        (f1, j1), (f2, j2) = where[ds], where[dt]
        if f1 != f2:
            F1, F2 = faces[f1], faces[f2]
            w1 = F1[j1 + 1:] + F1[:j1]      # darts after s->t, up to the one before s->t
            w2 = F2[j2 + 1:] + F2[:j2]
            return ('disk', self.seq(w1 + w2))
        F = faces[f1]
        n = len(F)
        x = [F[(j1 + 1 + i) % n] for i in range((j2 - j1 - 1) % n)]   # after s->t, before t->s
        y = [F[(j2 + 1 + i) % n] for i in range((j1 - j2 - 1) % n)]   # after t->s, before s->t
        return ('annulus', self.seq(x), self.seq(y))

    def header(self):
        return {
            'name': self.name,
            'vertices': self.vertices,
            'edges': self.edges,
            'attachments': self.attachments,
            'cyclic orders': self.orders,
            'number of rotation systems': self.count,
        }


def region_key(reg, sep=''):
    if reg[0] == 'disk':
        return fmt(canon(reg[1]), sep)
    a, b = sorted([fmt(canon(reg[1]), sep), fmt(canon(reg[2]), sep)])
    return a + ' | ' + b + ' (annulus)'


def alternates(seq):
    """Two distinct labels u, v appearing cyclically as u, v, u, v."""
    labels = sorted(set(seq))
    for u, v in itertools.combinations(labels, 2):
        r = [x for x in seq if x in (u, v)]
        # collapse cyclically consecutive repeats
        c = []
        for x in r:
            if not c or c[-1] != x:
                c.append(x)
        if len(c) > 1 and c[0] == c[-1]:
            c.pop()
        if len(c) >= 4:
            return (u, v)
    return None


def combine_repeats(seq):
    """Merge cyclically consecutive equal appearances."""
    c = []
    for x in seq:
        if not c or c[-1] != x:
            c.append(x)
    while len(c) > 1 and c[0] == c[-1]:
        c.pop()
    return tuple(c)


def consecutive_edges(seq):
    """The graph H obtained by joining consecutive attachments, last to first."""
    n = len(seq)
    es = set()
    for i in range(n):
        u, v = seq[i], seq[(i + 1) % n]
        if u != v:
            es.add(tuple(sorted((u, v))))
    return es


def shape(seq):
    """Classify H (consecutive attachments joined) and whether seq is a facial walk
    of a plane embedding of H: 'cycle', 'triangle with a pendant edge',
    'three-edge star', 'path', or None."""
    if any(seq[i] == seq[(i + 1) % len(seq)] for i in range(len(seq))):
        return None
    es = consecutive_edges(seq)
    vs = sorted(set(seq))
    deg = {v: sum(v in e for e in es) for v in vs}
    n, m = len(vs), len(es)
    # connectedness is automatic (a closed walk)
    degs = sorted(deg.values())
    if m == n and degs == [2] * n:
        kind = 'cycle'
        walk_len = n
    elif n == 4 and m == 4 and degs == [1, 2, 2, 3]:
        kind = 'triangle with a pendant edge'
        walk_len = 5
    elif n == 4 and m == 3 and degs == [1, 1, 1, 3]:
        kind = 'three-edge star'
        walk_len = 6
    elif m == n - 1 and max(degs) <= 2:
        kind = 'path'
        walk_len = 2 * (n - 1)
    else:
        return None
    # a closed walk of H using each edge the right number of times: for these
    # graphs, the facial walk containing every vertex has length walk_len and
    # every closed walk of that length along H through all vertices with no
    # immediate repeat traverses the boundary of the outer face.
    if len(seq) != walk_len:
        return None
    return kind
