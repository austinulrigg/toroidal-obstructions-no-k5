"""Generate the enumerations E1-E30, or the named items.

    python make_enumerations.py OUTPUT_FOLDER [E1 E2 ...] [--E4-covering-of-data]

Each generator writes its item to OUTPUT_FOLDER/enumerations. E1-E29 keep
working records under OUTPUT_FOLDER/work and skip completed computations.
With --E4-covering-of-data, E4 uses the published covering and recomputes
its DRAT proof, embeddings and classifications. Without the option, E4
searches for a covering."""
import os
import subprocess
import sys
import time

import common
import package

HERE = os.path.dirname(os.path.abspath(__file__))
ITEMS = ['E%d' % k for k in range(1, 31)]


def finished(work):
    """Whether the program of the item has finished in its work folder."""
    return os.path.exists(os.path.join(work, 'summary.json'))


def run(label, work, extra=()):
    prog = os.path.join(HERE, 'generate_%s.py' % label)
    steps = [[work, 'with_bd'], [work, 'without_bd'], [work, 'summarize']] if label == 'E20' else [[work] + list(extra)]
    for args in steps:
        subprocess.run([sys.executable, prog] + args, cwd=HERE, check=True)


PACKAGE = {
    'E6': lambda w, o: package.package_bridges('E6', w, o),
    'E7': lambda w, o: package.package_bridges('E7', w, o),
    'E18': package.package_E18,
    'E11': package.package_E11,
    'E1': package.package_E1,
    'E2': package.package_E2,
    'E10': package.package_E10,
    'E13': package.package_E13,
    'E3': package.package_E3,
    'E4': package.package_E4,
    'E12': package.package_E12,
}


def package_item(label, work, out):
    if label in PACKAGE:
        return PACKAGE[label](work, out)
    if label in package.WHOLE:
        s = package.package_whole(label, work, out)
        if label == 'E29':
            package.package_E29_bases(work, out)
        return s
    raise SystemExit('no packaging for %s' % label)


def main():
    args = [a for a in sys.argv[1:] if a != '--E4-covering-of-data']
    reuse = len(args) < len(sys.argv) - 1
    if not args:
        raise SystemExit(__doc__)
    out = os.path.abspath(args[0])
    items = args[1:] or ITEMS
    for label in items:
        if label not in ITEMS:
            raise SystemExit('unknown item %s' % label)
    print('data folder (for O1):', common.DATA)
    print('processes:', common.workers())
    for label in items:
        if label == 'E30':
            destination = os.path.join(out, 'enumerations', label)
            run(label, out)
            print('E30: written to %s' % destination, flush=True)
            continue
        work = os.path.join(out, 'work', label)
        t0 = time.time()
        if finished(work):
            print('%s: the program has already finished in %s' % (label, work), flush=True)
        else:
            print('%s: running generate_%s.py' % (label, label), flush=True)
            extra = []
            if label == 'E4' and reuse:
                extra = ['--use-covering-of', os.path.join(common.DATA, 'enumerations', 'E4')]
            run(label, work, extra)
        s = package_item(label, work, out)
        print('%s: written to %s (%.0f s)' % (label, os.path.join(out, 'enumerations', label), time.time() - t0),
              flush=True)
        if isinstance(s, dict) and 'verdicts' in s:
            print('   ', s['verdicts'], flush=True)


if __name__ == '__main__':
    main()
