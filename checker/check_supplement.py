"""Check every witness of the supplemental data, using only the files of the
data folder.

    python check_supplement.py ../supplemental_data
    python check_supplement.py ../supplemental_data --items O2,E11,R8

The program first compares every file of the folder with its SHA-256 checksum
in SHA256SUMS.txt. It then reads each item, rebuilds every graph from its
record, and checks every witness by its kind (witnesses.py), together with the
structure of each item that its module describes. It also checks that each
item's summary agrees with its records, and the DRAT proofs of E4, E5 and E13
with drat-trim (found on the PATH, or given by the environment variable
DRAT_TRIM; without it those proofs are reported as not checked).

It prints one line per item and ends with "all checks passed" or with the
list of failures. Exit status 0 means that every check passed.
"""
import argparse
import hashlib
import importlib
import os
import sys
import time
import traceback


# item -> module that checks it (each module has a function check(root, item) returning a short report)
ITEMS = {}
for k in ('O1', 'O2', 'O3', 'O4', 'TO1', 'TO2'):
    ITEMS[k] = 'items_obstructions'
for i in range(1, 22):
    ITEMS['W%d' % i] = 'items_witnesses'
for i in range(1, 12):
    ITEMS['R%d' % i] = 'items_rotations'
for i in (1, 2, 3, 4, 7, 8, 12):
    ITEMS['E%d' % i] = 'items_bridges'
for i in (5, 13):
    ITEMS['E%d' % i] = 'items_coverings'
for i in (6, 9, 10, 11, 14, 15, 16, 17, 18, 19):
    ITEMS['E%d' % i] = 'items_enumerations_a'
for i in range(20, 35):
    ITEMS['E%d' % i] = 'items_enumerations_b'


def order(k):
    pre = k.rstrip('0123456789')
    return (['O', 'W', 'R', 'E', 'TO'].index(pre), int(k[len(pre):]))


def check_checksums(root):
    """Every file listed in SHA256SUMS.txt is present with that checksum, and
    no other file is present (files whose names start with a dot, which some
    systems add to folders, are ignored)."""
    listing = os.path.join(root, 'SHA256SUMS.txt')
    if not os.path.exists(listing):
        return None
    want = {}
    with open(listing) as f:
        for line in f:
            h, p = line.rstrip('\n').split('  ', 1)
            want[p] = h
    bad = []
    for p, h in sorted(want.items()):
        full = os.path.join(root, *p.split('/'))
        if not os.path.isfile(full):
            bad.append('missing: ' + p)
            continue
        d = hashlib.sha256()
        with open(full, 'rb') as f:
            for b in iter(lambda: f.read(1 << 20), b''):
                d.update(b)
        if d.hexdigest() != h:
            bad.append('changed: ' + p)
    for dp, dn, fn in os.walk(root):
        dn[:] = [x for x in dn if not x.startswith('.')]
        for x in fn:
            if x.startswith('.'):
                continue
            p = os.path.relpath(os.path.join(dp, x), root).replace(os.sep, '/')
            if p != 'SHA256SUMS.txt' and p not in want:
                bad.append('not listed: ' + p)
    return len(want), bad


def main():
    if not __debug__:
        raise SystemExit('Run this checker without the option -O: its checks use assert statements.')
    ap = argparse.ArgumentParser(description='Check every witness of the supplemental data.')
    ap.add_argument('root', help='the supplemental data folder')
    ap.add_argument('--items', help='comma-separated items to check (default: all)')
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    if not os.path.isdir(os.path.join(root, 'enumerations')):
        raise SystemExit('%s does not look like the supplemental data folder' % root)
    items = sorted(ITEMS, key=order) if not a.items else [x.strip() for x in a.items.split(',')]
    for item in items:
        if item not in ITEMS:
            raise SystemExit('unknown item %s' % item)
    failures = []
    t_all = time.time()
    t = time.time()
    res = check_checksums(root)
    if res is None:
        print('files  no SHA256SUMS.txt, checksums not compared', flush=True)
    else:
        n, bad = res
        if bad:
            failures.append(('files', '; '.join(bad[:10]) + (' ...' if len(bad) > 10 else '')))
            print('files  FAIL  %6.1f s  %d problems: %s' % (time.time() - t, len(bad), '; '.join(bad[:10])), flush=True)
        else:
            print('files  ok    %6.1f s  all %d files match SHA256SUMS.txt' % (time.time() - t, n), flush=True)
    for item in items:
        mod = importlib.import_module(ITEMS[item])
        t = time.time()
        try:
            rep = mod.check(root, item)
            print('%-5s  ok    %6.1f s  %s' % (item, time.time() - t, rep), flush=True)
        except AssertionError as e:          # CheckFailed is an AssertionError too
            failures.append((item, str(e) or repr(e)))
            print('%-5s  FAIL  %6.1f s  %s' % (item, time.time() - t, str(e) or repr(e)), flush=True)
            traceback.print_exc()
        except Exception as e:
            failures.append((item, 'error: %r' % e))
            print('%-5s  ERROR %6.1f s' % (item, time.time() - t), flush=True)
            traceback.print_exc()
    print('%d items in %.0f s' % (len(items), time.time() - t_all))
    if failures:
        print('FAILURES:')
        for item, msg in failures:
            print('  %s: %s' % (item, msg))
        sys.exit(1)
    print('all checks passed')


if __name__ == '__main__':
    main()
