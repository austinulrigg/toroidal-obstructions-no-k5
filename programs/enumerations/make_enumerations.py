"""Run the enumerations and write them in the layout of the data folder.

    python make_enumerations.py OUTPUT_FOLDER [E1 E2 ...] [--E4-covering-of-data]

For each item named (all of E1, ..., E29 if none is named) this runs its
program generate_E<k>.py, which writes its records to OUTPUT_FOLDER/work/E<k>,
and then writes the files of the item to OUTPUT_FOLDER/enumerations/E<k>
(package.py). A program that has already finished in OUTPUT_FOLDER/work/E<k>
is not run again; delete that folder to run it again. An interrupted program
that keeps a file progress.jsonl continues where it stopped (E3, E4 and E12
keep their own files for this).

E4 runs its covering with several processes, and which graphs it records
depends on the order in which they finish; with --E4-covering-of-data,
generate_E4.py takes the covering recorded in the data folder instead and
redoes everything after it, which reproduces the data folder.

The output can be compared with the data folder by
    python ../compare_files.py ../../supplemental_data OUTPUT_FOLDER enumerations
"""
import os
import subprocess
import sys
import time

import common
import package

HERE = os.path.dirname(os.path.abspath(__file__))
ITEMS = ['E%d' % k for k in range(1, 30)]


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
