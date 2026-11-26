"""Compute the genera of T1–T74 with PAGE and multi_genus and of the
enumeration graphs with multi_genus. Write the results and timings to results/.

    python run_genus.py
"""
import argparse
import csv
import gzip
import multiprocessing as mp
import os
import platform
import subprocess
import sys
import time
from collections import Counter, defaultdict

from common import (HERE, components, from_graph6, genus_of_rotation, induced, load_json, multi_genus_batch,
                    n_edges, read_list, save_json, tool, zero_indexed)

RESULTS = HERE / 'results'


def log(who, msg):
    print('[%s %s] %s' % (time.strftime('%H:%M:%S'), who, msg), flush=True)


def embed_timed(adj, algorithm):
    """Genus by one program, component by component; returns (genus, seconds)."""
    import graph_genus as gg
    total, secs = 0, 0.0
    for comp in components(adj):
        H = induced(adj, comp)
        if n_edges(H) == 0:
            continue
        A = zero_indexed(H)
        t = time.perf_counter()
        g, rot = gg.embed(A, algorithm=algorithm)
        secs += time.perf_counter() - t
        if genus_of_rotation(rot) != g:
            raise AssertionError('%s reported genus %d, but its embedding has another genus' % (algorithm, g))
        if [sorted(nb) for nb in rot] != [sorted(nb) for nb in A]:
            raise AssertionError('%s returned an embedding of another graph' % algorithm)
        total += g
    return total, secs


def obstructions():
    return [(r[0], from_graph6(r[1])) for r in read_list(HERE / 'obstructions.txt')]


def enumeration_graphs():
    return [(int(r[0]), r[1]) for r in read_list(HERE / 'enumeration_graphs.txt')]


def items_of(uses):
    """enumeration -> sorted list of the indices of the graphs it uses."""
    by_item = defaultdict(set)
    for g in uses['graphs']:
        for u in g['uses']:
            by_item[u['item']].add(g['index'])
    return {k: sorted(v) for k, v in sorted(by_item.items(), key=lambda kv: int(kv[0][1:]))}


def multi_genus_worker(out_path):
    res = {'obstructions': {}, 'enumerations': {}, 'enumeration graph genus': {}}
    for name, adj in obstructions():
        g, s = embed_timed(adj, 'multi_genus')
        res['obstructions'][name] = {'genus': g, 'seconds': s}
    log('multi_genus', 'T1-T74 done')
    graphs = {i: from_graph6(s) for i, s in enumeration_graphs()}
    for item, idx in items_of(load_json(HERE / 'enumeration_graph_uses.json.gz')).items():
        gen, secs = multi_genus_batch([graphs[i] for i in idx])
        for i, g in zip(idx, gen):
            old = res['enumeration graph genus'].setdefault(i, g)
            assert old == g, 'two runs of multi_genus disagree on graph %d' % i
        res['enumerations'][item] = {'graphs': len(idx), 'seconds': secs}
        log('multi_genus', 'graphs of %s: %d graphs in %.2f s, genus %s'
            % (item, len(idx), secs, dict(sorted(Counter(gen).items()))))
    save_json(out_path, res)


def page_worker(out_path, threads):
    os.environ['THREADS'] = str(threads)
    res = {'obstructions': {}, 'threads': threads}
    for name, adj in obstructions():
        g, s = embed_timed(adj, 'page')
        res['obstructions'][name] = {'genus': g, 'seconds': s}
    log('PAGE', 'T1-T74 done')
    save_json(out_path, res)


def machine():
    info = {'system': platform.platform(), 'machine': platform.machine(), 'python': sys.version.split()[0],
            'cores': os.cpu_count()}
    if platform.system() == 'Darwin':
        try:
            info['processor'] = subprocess.run(['sysctl', '-n', 'machdep.cpu.brand_string'], capture_output=True,
                                               text=True).stdout.strip()
        except OSError:
            pass
    else:
        info['processor'] = platform.processor()
    import graph_genus
    info['graph-genus'] = graph_genus.__version__
    info['programs'] = 'page and multi_genus as bundled with graph-genus %s' % graph_genus.__version__
    tool('multi_genus'), tool('page')          # stop early if either is missing
    return info


def write_tables(t1, t2, mg, pg, info, wall):
    RESULTS.mkdir(exist_ok=True)
    by_item = items_of(load_json(HERE / 'enumeration_graph_uses.json.gz'))

    rows1 = []
    for name, adj in t1:
        p, m = pg['obstructions'][name], mg['obstructions'][name]
        rows1.append([name, len(adj), n_edges(adj), p['genus'], '%.3f' % p['seconds'], m['genus'],
                      '%.3f' % m['seconds']])
    with open(RESULTS / 'obstructions.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['graph', 'vertices', 'edges', 'genus (PAGE)', 'PAGE time (s)', 'genus (multi_genus)',
                    'MG time (s)'])
        w.writerows(rows1)

    gen = {int(k): v for k, v in mg['enumeration graph genus'].items()}
    rows2 = []
    for item, idx in by_item.items():
        gs = [gen[i] for i in idx]
        rows2.append([item, len(idx), min(gs), max(gs), '%.3f' % mg['enumerations'][item]['seconds']])
    with open(RESULTS / 'enumeration_graphs.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['enumeration', 'graphs', 'smallest genus', 'largest genus', 'MG time (s)'])
        w.writerows(rows2)
    with gzip.open(RESULTS / 'enumeration_graphs_genus.txt.gz', 'wt') as f:
        f.write('# index genus (multi_genus)\n')
        for i, _ in t2:
            f.write('%d %d\n' % (i, gen[i]))

    allg = Counter(gen[i] for i, _ in t2)
    L = ['# Genus computations', '',
         'PAGE and `multi_genus` ran concurrently on %s (%d cores) with Python %s and graph-genus %s. '
         'The run took %.0f seconds with PAGE set to `THREADS=%d`. '
         '[Machine and software versions](run_info.json).'
         % (info.get('processor') or info['machine'], info['cores'], info['python'], info['graph-genus'],
            wall, pg['threads']), '',
         '## The 74 obstructions', '',
         'Both programs give genus two for each graph. Times include program startup. '
         'For $T_{25}$ the genus and time are summed over its two components.', '',
         '| Graph | Vertices | Edges | Genus | PAGE (s) | `multi_genus` (s) |',
         '|---|---|---|---|---|---|']
    for r in rows1:
        g = str(r[3]) if r[3] == r[5] else '%s (PAGE) / %s (multi_genus)' % (r[3], r[5])
        L.append('| %s | %d | %d | %s | %s | %s |' % (r[0], r[1], r[2], g, r[4], r[6]))
    L += ['', '## Nontoroidal graphs used in the enumerations', '',
          'Each row counts graphs up to isomorphism within that enumeration. A graph used in several '
          'enumerations appears in each count. The time is for one run of `multi_genus` on all graphs in the row.', '',
          '| Enumeration | Graphs | Smallest genus | Largest genus | `multi_genus` (s) |',
          '|---|---|---|---|---|']
    for r in rows2:
        L.append('| %s | %d | %d | %d | %s |' % tuple(r))
    distribution = ' and '.join(f'{n:,} have genus {g}' for g, n in sorted(allg.items()))
    L += ['', f'Of the {len(t2):,} graphs up to isomorphism, {distribution}. '
          '[Genus by graph index](enumeration_graphs_genus.txt.gz).']
    (RESULTS / 'TABLES.md').write_text('\n'.join(L) + '\n')
    save_json(RESULTS / 'run_info.json', dict(info, seconds=round(wall, 1), PAGE_threads=pg['threads']), indent=1)


def main():
    ap = argparse.ArgumentParser(description='The genus computations of the paper, with PAGE and multi_genus.')
    ap.add_argument('--cores', type=int, default=os.cpu_count(),
                    help='cores to use: one for multi_genus, the rest for PAGE (default: all)')
    a = ap.parse_args()
    threads = max(1, a.cores - 1)
    info = machine()
    log('run', 'multi_genus on 1 core and PAGE with %d thread%s on %s'
        % (threads, '' if threads == 1 else 's', info.get('processor') or 'this machine'))
    t0 = time.time()
    t1, t2 = obstructions(), enumeration_graphs()
    RESULTS.mkdir(exist_ok=True)
    ctx = mp.get_context('spawn')
    pm = ctx.Process(target=multi_genus_worker, args=(str(RESULTS / 'multi_genus_raw.json'),))
    pp = ctx.Process(target=page_worker, args=(str(RESULTS / 'page_raw.json'), threads))
    pm.start()
    pp.start()
    pm.join()
    pp.join()
    if pm.exitcode or pp.exitcode:
        raise SystemExit('a genus program failed; see the messages above')
    mg = load_json(RESULTS / 'multi_genus_raw.json')
    pg = load_json(RESULTS / 'page_raw.json')
    wall = time.time() - t0
    for name, _ in t1:
        assert mg['obstructions'][name]['genus'] == 2 and pg['obstructions'][name]['genus'] == 2, name
    assert len(mg['enumeration graph genus']) == len(t2) and all(g >= 2 for g in mg['enumeration graph genus'].values())
    write_tables(t1, t2, mg, pg, info, wall)
    os.remove(RESULTS / 'multi_genus_raw.json')
    os.remove(RESULTS / 'page_raw.json')
    log('run', 'done in %.0f s. T1-T74 have genus 2 by both programs. All %d nontoroidal graphs of the '
        'enumerations have genus at least 2 by multi_genus. The tables are in %s' % (wall, len(t2), RESULTS / 'TABLES.md'))


if __name__ == '__main__':
    main()
