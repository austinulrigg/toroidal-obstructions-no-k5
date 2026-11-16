"""Covering a finite family of bridge assemblies by toroidal and nontoroidal
graphs (used for E4 and E12).

An assembly is a count vector n over a list of atoms (bridges, or edges, with
named attachments). Its graph is a fixed base graph plus n[i] copies of each
atom i. The family is given by bounds 0 <= n[i] <= maxima[i] and clauses over
the variables x(i,k) = [n[i] >= k].

The loop asks a SAT solver for a vector of the family not yet covered.
  - If its graph is toroidal, atoms are added one at a time while it stays
    toroidal, giving a maximal toroidal vector P. Every vector below P is
    toroidal (its graph is a subgraph of that of P), so the clause "some
    n[i] > P[i]" is added.
  - If its graph is nontoroidal, atoms are removed one at a time while it
    stays nontoroidal, giving a vector F whose one-atom deletions are all
    toroidal. Every vector above F is nontoroidal, so the clause "some
    n[i] < F[i]" is added.
When the solver reports that the clauses are unsatisfiable, every vector of
the family lies below some P or above some F. An obstruction in the family
is nontoroidal with every proper subgraph toroidal, so it equals one of the
vectors F. The final formula is written in DIMACS form, with a DRAT proof of
its unsatisfiability.

cover() is the loop above. parallel_cover() is the same loop with the
graphs tested by several worker processes: the coordinator alone owns the
SAT solver, each worker turns a vector into a cone (a vector P or F covering
it), and the state is saved every few seconds so that a second start
continues. Which vectors are proposed depends on the order in which the
workers finish, so two runs of parallel_cover() find different, equally
valid, families of vectors P and F.

Genus decisions come from multi_genus (multigenus.py and
enumeration.genus_above_one).
"""

from collections import OrderedDict
import gzip
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import queue
import shutil
import signal
import subprocess
import tempfile
import time
import traceback
from types import SimpleNamespace

from pysat.solvers import Solver

import enumeration as el

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parent.parent / "tools"


class Domain:
    def __init__(self, atoms, base, maxima, clauses_fn, name):
        """atoms: list of dicts with 'id', 'edges' (pairs; interior vertices
        start with '@'), 'attachments'. base: networkx graph. maxima: list of
        ints. clauses_fn(var) returns (clauses, top variable) over
        var(i, k) = the variable [n_i >= k]."""
        self.atoms, self.base, self.maxima, self.name = atoms, base, maxima, name
        self.var = {}
        v = 0
        for i, m in enumerate(maxima):
            for k in range(1, m + 1):
                v += 1
                self.var[(i, k)] = v
        self.nvars = v
        self.base_clauses = []
        for i, m in enumerate(maxima):
            for k in range(2, m + 1):
                self.base_clauses.append([-self.var[(i, k)], self.var[(i, k - 1)]])
        extra, top = clauses_fn(lambda i, k: self.var[(i, k)])
        self.base_clauses += extra
        self.nvars = max(self.nvars, top)

    def graph(self, n):
        """The base graph plus n[i] copies of atom i; interior vertex @x of
        copy j of atom i is named <atom id>#j@x."""
        G = self.base.copy()
        for i, c in enumerate(n):
            for j in range(c):
                tag = "%s#%d" % (self.atoms[i]["id"], j)
                G.add_edges_from(tuple(x if not x.startswith("@") else tag + x for x in e)
                                 for e in self.atoms[i]["edges"])
        return G

    def decode(self, model):
        pos = set(x for x in model if x > 0)
        return tuple(sum(1 for k in range(1, m + 1) if self.var[(i, k)] in pos)
                     for i, m in enumerate(self.maxima))

    def above_clause(self, P):
        """Excludes every vector <= P."""
        return [self.var[(i, P[i] + 1)] for i in range(len(P)) if P[i] < self.maxima[i]]

    def below_clause(self, F):
        """Excludes every vector >= F."""
        return [-self.var[(i, F[i])] for i in range(len(F)) if F[i] > 0]

    def nontoroidal(self, vectors):
        return el.genus_above_one([self.graph(n) for n in vectors])

    def small_enough(self, n, vertices=120, edges=500):
        G = self.graph(n)
        return G.number_of_nodes() <= vertices and G.number_of_edges() <= edges

    def failing_below(self, n):
        """For a vector too large to test: add its atoms one at a time until
        the graph is nontoroidal, and return that vector (or None if every
        stage is toroidal and the last stage is still too large)."""
        cur = [0] * len(n)
        for i in range(len(n)):
            for _ in range(n[i]):
                cur[i] += 1
                if not self.small_enough(cur):
                    return None
                if self.nontoroidal([tuple(cur)])[0]:
                    return tuple(cur)
        return None

    def extend(self, P, allowed=None):
        """A maximal toroidal vector above P (adding one atom at a time)."""
        P = list(P)
        while True:
            options = []
            for i in range(len(P)):
                if P[i] < self.maxima[i]:
                    Q = list(P)
                    Q[i] += 1
                    if allowed is None or allowed(Q):
                        options.append(tuple(Q))
            if not options:
                return tuple(P)
            flags = self.nontoroidal(options)
            good = [o for o, f in zip(options, flags) if not f]
            if not good:
                return tuple(P)
            P = list(good[0])

    def shrink(self, F):
        """A nontoroidal vector below F whose one-atom deletions are all
        toroidal (removing one atom at a time)."""
        F = list(F)
        while True:
            options = []
            for i in range(len(F)):
                if F[i] > 0:
                    Q = list(F)
                    Q[i] -= 1
                    options.append(tuple(Q))
            if not options:
                return tuple(F)
            flags = self.nontoroidal(options)
            bad = [o for o, f in zip(options, flags) if f]
            if not bad:
                return tuple(F)
            F = list(bad[0])


def image(n, perm):
    out = [0] * len(n)
    for i, c in enumerate(n):
        out[perm[i]] = c
    return tuple(out)


def cover(domain, log=None, solver_name="cadical153", allowed=None, max_iter=10 ** 7, symmetries=None,
          origins=None):
    """The covering loop with one process. Returns (passing, failing, cnf).

    With symmetries (atom permutations induced by relabelings of the base
    graph, so that every image graph is isomorphic to the original), each
    vector found is stored with all its distinct images, and origins, if
    given, receives for each stored vector (kind, index) the pair (index of
    the vector found, index of the symmetry)."""
    # The solver records a proof. The proof is not used (fresh_proof() writes
    # the final one), but recording it is part of the setting that gives the
    # data folder.
    s = Solver(name=solver_name, bootstrap_with=domain.base_clauses, with_proof=True)
    s.set_phases([-v for v in range(1, domain.nvars + 1)])
    passing, failing = [], []
    perms = symmetries or [list(range(len(domain.maxima)))]
    it = 0

    def store(kind, vec):
        lst = passing if kind == "P" else failing
        rep = len(lst)
        seen = set()
        for k, perm in enumerate(perms):
            im = image(vec, perm)
            if im in seen:
                continue
            seen.add(im)
            if origins is not None:
                origins[(kind, len(lst))] = (rep, k)
            lst.append(im)
            cl = domain.above_clause(im) if kind == "P" else domain.below_clause(im)
            assert cl, "empty covering clause"
            s.add_clause(cl)

    while s.solve():
        it += 1
        n = domain.decode(s.get_model())
        if not domain.small_enough(n):
            F0 = domain.failing_below(n)
            if F0 is None:
                raise RuntimeError("a vector too large to test")
            store("F", domain.shrink(F0))
        elif domain.nontoroidal([n])[0]:
            store("F", domain.shrink(n))
        else:
            store("P", domain.extend(n, allowed))
        if log and it % 200 == 0:
            log("iteration %d: %d toroidal, %d nontoroidal" % (it, len(passing), len(failing)))
        if it >= max_iter:
            raise RuntimeError("iteration limit")
    s.delete()
    cnf = list(domain.base_clauses) + [domain.above_clause(P) for P in passing] + \
        [domain.below_clause(F) for F in failing]
    return passing, failing, cnf


def read_dimacs(path):
    """(number of variables, clauses) of a DIMACS file, possibly
    gzip-compressed, with the clauses in the order of the file."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        lines = f.read().split("\n")
    head = lines[0].split()
    assert head[:2] == ["p", "cnf"]
    clauses, cur = [], []
    for line in lines[1:]:
        if not line.strip() or line.startswith("c"):
            continue
        for x in map(int, line.split()):
            if x == 0:
                clauses.append(cur)
                cur = []
            else:
                cur.append(x)
    assert not cur and len(clauses) == int(head[3])
    return int(head[2]), clauses


def write_dimacs(path, nvars, clauses):
    with open(path, "w") as f:
        f.write("p cnf %d %d\n" % (nvars, len(clauses)))
        for c in clauses:
            f.write(" ".join(map(str, c)) + " 0\n")


_membership = {}


def in_domain(domain, n):
    """Whether the vector satisfies the base clauses for some values of the
    auxiliary variables (a SAT call with the variables x(i,k) fixed)."""
    key = id(domain)
    if key not in _membership:
        _membership[key] = Solver(name="cadical153", bootstrap_with=domain.base_clauses)
    assum = [(domain.var[(i, k)] if k <= c else -domain.var[(i, k)])
             for i, c in enumerate(n) for k in range(1, domain.maxima[i] + 1)]
    return _membership[key].solve(assumptions=assum)


def find_drat_trim():
    """drat-trim given by the environment variable DRAT_TRIM, or on the PATH,
    or tools/drat-trim of the repository, compiled from tools/drat-trim.c the
    first time it is needed if a C compiler is available; None otherwise."""
    exe = os.environ.get("DRAT_TRIM") or shutil.which("drat-trim")
    if exe:
        return exe
    local = TOOLS / "drat-trim"
    if local.is_file() and os.access(local, os.X_OK):
        return str(local)
    src = TOOLS / "drat-trim.c"
    cc = shutil.which("cc") or shutil.which("clang") or shutil.which("gcc")
    if src.is_file() and cc:
        r = subprocess.run([cc, "-O2", "-o", str(local), str(src)], capture_output=True, text=True)
        if r.returncode == 0 and local.is_file():
            return str(local)
    return None


def fresh_proof(cnf_path, drat_path):
    """Solves the final formula from scratch with Glucose 4 and writes its
    DRAT proof. Returns True when drat-trim verifies the proof, False when
    it does not, and None when drat-trim is not available."""
    from pysat.formula import CNF
    s = Solver(name="glucose4", bootstrap_with=CNF(from_file=str(cnf_path)).clauses, with_proof=True)
    if s.solve():
        raise RuntimeError("the final formula is satisfiable: the covering is not complete")
    with open(drat_path, "w") as f:
        for line in s.get_proof():
            f.write(line.rstrip("\n") + "\n")
    s.delete()
    exe = find_drat_trim()
    if exe is None:
        return None
    chk = subprocess.run([exe, str(cnf_path), str(drat_path)], capture_output=True, text=True)
    return "s VERIFIED" in chk.stdout.splitlines()


# ---------------------------------------------------------------- parallel version

def _sha_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _domain_payload(domain, symmetries):
    return {
        "name": domain.name,
        "atoms": domain.atoms,
        "maxima": list(domain.maxima),
        "base_vertices": sorted(domain.base.nodes(), key=str),
        "base_edges": sorted((sorted(e, key=str) for e in domain.base.edges()), key=lambda e: tuple(map(str, e))),
        "base_clauses": domain.base_clauses,
        "threshold_variables": [[i, k, v] for (i, k), v in sorted(domain.var.items())],
        "nvars": domain.nvars,
        "symmetries": symmetries,
    }


def _fingerprint(payload):
    """The family, the programs that decide genus, and the multi_genus binary:
    a saved state is resumed only when all of these are unchanged."""
    import multigenus
    canonical = dict(payload, base_clauses=[list(c) for c in sorted(
        tuple(sorted(c, key=lambda lit: (abs(lit), lit))) for c in payload["base_clauses"])])
    details = {"family": canonical,
               "programs": {name: _sha_file(HERE / name)
                            for name in ("covering.py", "enumeration.py", "graphs.py", "multigenus.py")},
               "multi_genus": _sha_file(Path(multigenus.BIN).resolve())}
    return hashlib.sha256(json.dumps(details, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _atomic_pickle(path, state):
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            pickle.dump(state, stream, protocol=pickle.HIGHEST_PROTOCOL)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


class _GenusOracle:
    """Genus decisions of one worker, with a bounded cache keyed by the
    smallest image of the vector under the symmetries."""

    def __init__(self, payload, cache_size):
        import networkx as nx
        import multigenus
        self.nx, self.mg = nx, multigenus
        base = nx.Graph()
        base.add_nodes_from(payload["base_vertices"])
        base.add_edges_from(payload["base_edges"])
        self.domain = SimpleNamespace(atoms=payload["atoms"], maxima=payload["maxima"], base=base)
        self.maxima = tuple(payload["maxima"])
        self.cache = OrderedDict()
        self.cache_size = max(0, int(cache_size))
        self.inverse_permutations = []
        for permutation in payload["symmetries"]:
            inverse = [0] * len(permutation)
            for i, j in enumerate(permutation):
                inverse[j] = i
            self.inverse_permutations.append(tuple(inverse))

    def graph(self, vector):
        return Domain.graph(self.domain, vector)

    def key(self, vector):
        return min(tuple(vector[i] for i in inverse) for inverse in self.inverse_permutations)

    def __call__(self, vectors):
        keys = [self.key(vector) for vector in vectors]
        missing, seen, answers = [], set(), {}
        for key in keys:
            if key in self.cache:
                answers[key] = self.cache[key]
                self.cache.move_to_end(key)
            elif key not in seen:
                seen.add(key)
                missing.append(key)
        connected, disconnected = [], []
        for key in missing:
            G = self.graph(key)
            if len(G) >= 255:
                raise RuntimeError("a graph with 255 or more vertices: multi_genus cannot decide it")
            if G.number_of_edges() == 0:
                answers[key] = False
            elif self.nx.is_connected(G):
                connected.append((key, G))
            else:
                disconnected.append((key, G))
        if connected:
            flags = self.mg.nontoroidal([G for _, G in connected])
            answers.update((key, bool(flag)) for (key, _), flag in zip(connected, flags))
        if disconnected:
            # Two components of genus one give genus two, so the genus is
            # added over the components (enumeration.genus_above_one).
            flags = el.genus_above_one([G for _, G in disconnected])
            answers.update((key, bool(flag)) for (key, _), flag in zip(disconnected, flags))
        for key in missing:
            if self.cache_size:
                self.cache[key] = answers[key]
                self.cache.move_to_end(key)
                while len(self.cache) > self.cache_size:
                    self.cache.popitem(last=False)
        return [answers[key] for key in keys]

    def small_enough(self, vector):
        G = self.graph(vector)
        return len(G) <= 120 and G.number_of_edges() <= 500

    def process(self, seed):
        """The cone of a proposed vector: (kind, vector), kind 'P' for a
        maximal toroidal vector above it, 'F' for a minimal nontoroidal one
        below it."""
        seed = tuple(seed)
        initial = seed
        if not self.small_enough(seed):
            cur = [0] * len(seed)
            found = False
            for i, count in enumerate(seed):
                for _ in range(count):
                    cur[i] += 1
                    if not self.small_enough(cur):
                        raise RuntimeError("no small nontoroidal vector below a vector too large to test")
                    if self([cur])[0]:
                        initial, found = tuple(cur), True
                        break
                if found:
                    break
            if not found:
                raise RuntimeError("no small nontoroidal vector below a vector too large to test")
            kind = "F"
        else:
            kind = "F" if self([seed])[0] else "P"
        cur = list(initial)
        if kind == "P":
            # An increment that makes the graph nontoroidal keeps it so after
            # later increments, so one pass over the atoms gives a maximal vector.
            for i, maximum in enumerate(self.maxima):
                while cur[i] < maximum:
                    cur[i] += 1
                    if self([cur])[0]:
                        cur[i] -= 1
                        break
        else:
            # A deletion that leaves a toroidal graph stays toroidal after
            # further deletions, so one pass gives a minimal vector.
            for i in range(len(cur)):
                while cur[i] > 0:
                    cur[i] -= 1
                    if not self([cur])[0]:
                        cur[i] += 1
                        break
        return {"seed": seed, "kind": kind, "vector": tuple(cur)}


def _worker_main(index, payload, cache_size, tasks, results):
    try:
        if os.name == "posix":
            # multi_genus processes started by this worker share its process
            # group, so stopping the worker also stops them.
            os.setsid()
        oracle = _GenusOracle(payload, cache_size)
        results.put({"event": "ready", "worker": index})
        while True:
            task = tasks.get()
            if task is None:
                return
            job, seed = task
            try:
                record = oracle.process(seed)
                results.put({"event": "result", "worker": index, "job": job, **record})
            except BaseException as exc:
                results.put({"event": "error", "worker": index, "job": job, "seed": seed,
                             "error": repr(exc), "traceback": traceback.format_exc()})
                return
    except BaseException as exc:
        results.put({"event": "fatal", "worker": index, "error": repr(exc), "traceback": traceback.format_exc()})


def _stop_workers(processes, task_queues):
    for tasks in task_queues:
        try:
            tasks.put_nowait(None)
        except queue.Full:
            pass
    for process in processes:
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                if process.is_alive():
                    process.terminate()
        elif process.is_alive():
            process.terminate()
    for process in processes:
        process.join(timeout=1)
        if process.is_alive():
            process.kill()
            process.join(timeout=1)
    for tasks in task_queues:
        tasks.cancel_join_thread()
        tasks.close()


def _exact_vector_block(domain, vector):
    clause = []
    for i, (count, maximum) in enumerate(zip(vector, domain.maxima)):
        if count:
            clause.append(-domain.var[i, count])
        if count < maximum:
            clause.append(domain.var[i, count + 1])
    return clause


def parallel_cover(domain, output_dir, workers=2, symmetries=None, log=None, max_seconds=None,
                   checkpoint_seconds=10, cache_size=20000):
    """The covering loop with worker processes. Returns a dict with
    'complete', 'P', 'F', 'origins' and 'cnf'.

    The coordinator proposes a vector of the family not yet covered and
    blocks it (with its images) by a temporary clause until a worker returns
    its cone; the temporary clauses are never part of the final formula.
    Every returned cone is checked to contain its proposed vector, so the
    temporary clauses are implied by the final ones, and the loop ends only
    when the solver finds no vector and no worker is busy. The state is
    saved in cover_checkpoint.pkl every checkpoint_seconds seconds; a second
    start with the same family and programs continues from it. max_seconds
    stops this start early, with the state saved."""
    symmetries = [list(p) for p in (symmetries or [list(range(len(domain.maxima)))])]
    identity = list(range(len(domain.maxima)))
    if symmetries[0] != identity:
        raise ValueError("the identity must be the first symmetry")
    for permutation in symmetries:
        if sorted(permutation) != identity or any(domain.maxima[i] != domain.maxima[j]
                                                  for i, j in enumerate(permutation)):
            raise ValueError("a symmetry is not a permutation of the atoms preserving the maxima")
    payload = _domain_payload(domain, symmetries)
    fingerprint = _fingerprint(payload)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "cover_checkpoint.pkl"
    started = time.monotonic()
    solver, results = None, None
    processes, task_queues = [], []
    pending = {}
    if path.exists():
        with path.open("rb") as stream:
            state = pickle.load(stream)
        if state.get("fingerprint") != fingerprint:
            raise RuntimeError("%s was made with another family or other programs; move it away to start "
                               "afresh" % path)
        if log:
            log("continuing: %d toroidal, %d nontoroidal vectors, %d cones" %
                (len(state["P"]), len(state["F"]), state["iterations"]))
    else:
        state = {"fingerprint": fingerprint, "P": [], "F": [], "origins": {}, "iterations": 0,
                 "complete": False, "next_job": 0}
    P, F, origins = state["P"], state["F"], state["origins"]
    indices = {"P": {tuple(v): i for i, v in enumerate(P)}, "F": {tuple(v): i for i, v in enumerate(F)}}

    def final_cnf():
        return list(domain.base_clauses) + [domain.above_clause(v) for v in P] + [domain.below_clause(v) for v in F]

    def save():
        _atomic_pickle(path, state)

    def returned():
        return dict(state, cnf=final_cnf())

    if state["complete"]:
        return returned()
    try:
        solver = Solver(name="cadical153", bootstrap_with=final_cnf())
        solver.set_phases([-v for v in range(1, domain.nvars + 1)])
        context = mp.get_context("spawn")
        results = context.Queue()
        for index in range(workers):
            tasks = context.Queue(maxsize=1)
            process = context.Process(target=_worker_main, args=(index, payload, cache_size, tasks, results))
            process.start()
            task_queues.append(tasks)
            processes.append(process)
        idle, ready, blocked = set(), set(), set()
        exhausted = False
        last_checkpoint = time.monotonic()

        def store(kind, vector):
            vectors = P if kind == "P" else F
            vector = tuple(vector)
            rep = indices[kind][vector] if vector in indices[kind] else len(vectors)
            for k, permutation in enumerate(symmetries):
                im = image(vector, permutation)
                if im in indices[kind]:
                    continue
                indices[kind][im] = len(vectors)
                origins[kind, len(vectors)] = (rep, k)
                vectors.append(im)
                solver.add_clause(domain.above_clause(im) if kind == "P" else domain.below_clause(im))

        while True:
            now = time.monotonic()
            if max_seconds is not None and now - started >= max_seconds:
                break
            if now - last_checkpoint >= checkpoint_seconds:
                save()
                last_checkpoint = time.monotonic()
                if log:
                    log("%d toroidal, %d nontoroidal vectors, %d cones, %d proposed vectors in progress" %
                        (len(P), len(F), state["iterations"], len(pending)))
            while idle and not exhausted:
                if not solver.solve():
                    exhausted = True
                    break
                seed = tuple(domain.decode(solver.get_model()))
                if seed in blocked:
                    raise RuntimeError("the solver proposed a blocked vector")
                worker = min(idle)
                idle.remove(worker)
                job = state["next_job"]
                state["next_job"] += 1
                pending[job] = {"seed": seed, "worker": worker}
                for permutation in symmetries:
                    im = image(seed, permutation)
                    if im not in blocked:
                        solver.add_clause(_exact_vector_block(domain, im))
                        blocked.add(im)
                task_queues[worker].put((job, seed))
            if not pending and exhausted:
                state["complete"] = True
                break
            try:
                event = results.get(timeout=.25)
            except queue.Empty:
                dead = [i for i, process in enumerate(processes) if process.exitcode is not None]
                if dead:
                    try:
                        event = results.get(timeout=.25)
                    except queue.Empty:
                        raise RuntimeError("worker processes %s stopped" % dead)
                else:
                    continue
            worker = event["worker"]
            if event["event"] == "ready":
                ready.add(worker)
                idle.add(worker)
                continue
            if event["event"] in ("error", "fatal"):
                raise RuntimeError("worker %d: %s\n%s" % (worker, event["error"], event["traceback"]))
            request = pending.get(event["job"])
            if request is None or request["worker"] != worker or tuple(event["seed"]) != request["seed"]:
                raise RuntimeError("a worker returned a result for another vector")
            vector, seed, kind = tuple(event["vector"]), request["seed"], event["kind"]
            covered = (kind == "P" and all(a <= b for a, b in zip(seed, vector))) or (
                kind == "F" and all(a >= b for a, b in zip(seed, vector)))
            if not covered:
                raise RuntimeError("a returned cone does not contain its proposed vector")
            store(kind, vector)
            pending.pop(event["job"])
            idle.add(worker)
            state["iterations"] += 1
        save()
        return returned()
    finally:
        _stop_workers(processes, task_queues)
        if results is not None:
            results.cancel_join_thread()
            results.close()
        if solver is not None:
            solver.delete()
