"""Compare regenerated files with the corresponding supplemental data.

    python compare_files.py DATA_FOLDER OUTPUT_FOLDER [SUBFOLDER ...]

Compare files byte for byte, then compare the contents of JSON files that
differ (including compressed JSON). Key order, formatting and values under
the key "seconds" do not affect the result. Report differences confined to
rotation systems separately. Missing, extra or otherwise different files
give a nonzero exit status."""
import gzip
import json
import os
import sys


def read_json(p):
    with (gzip.open(p, 'rt', encoding='utf-8') if p.endswith('.gz') else open(p, encoding='utf-8')) as f:
        return json.load(f)


def differences(a, b, path, found):
    """Appends to found every place (a path into the JSON content) where a
    and b differ. A rotation system (the value of a key "rotation") that
    differs counts as one place."""
    if type(a) != type(b):
        found.append(path or '(top)')
    elif isinstance(a, dict):
        if set(a) != set(b):
            found.append((path or '(top)') + ' (keys ' + str(sorted(set(a) ^ set(b))) + ')')
            return
        for k in a:
            if k == 'seconds':
                continue
            if k == 'rotation':
                if a[k] != b[k]:
                    found.append(path + '/rotation')
                continue
            differences(a[k], b[k], path + '/' + str(k), found)
    elif isinstance(a, list):
        if len(a) != len(b):
            found.append((path or '(top)') + ' (lengths %d and %d)' % (len(a), len(b)))
            return
        for i, (x, y) in enumerate(zip(a, b)):
            differences(x, y, path + '[%d]' % i, found)
    elif a != b:
        found.append(path or '(top)')


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    data, out = sys.argv[1], sys.argv[2]
    subs = sys.argv[3:] or ['']
    same = fmt = 0
    bad, rot = [], []
    for sub in subs:
        top = os.path.join(out, sub)
        for dp, dn, fn in os.walk(top):
            dn.sort()
            for f in sorted(fn):
                p = os.path.join(dp, f)
                rel = os.path.relpath(p, out)
                q = os.path.join(data, rel)
                if not os.path.exists(q):
                    bad.append('%s: not in the data folder' % rel)
                    continue
                if open(p, 'rb').read() == open(q, 'rb').read():
                    same += 1
                    continue
                if rel.endswith('.json') or rel.endswith('.json.gz'):
                    found = []
                    differences(read_json(p), read_json(q), '', found)
                    if not found:
                        fmt += 1
                    elif all(x.endswith('/rotation') for x in found):
                        rot.append('%s: differs in %d rotation system%s only (e.g. %s)'
                                   % (rel, len(found), '' if len(found) == 1 else 's', found[0]))
                    else:
                        others = [x for x in found if not x.endswith('/rotation')]
                        bad.append('%s: differs at %d place%s (e.g. %s)'
                                   % (rel, len(others), '' if len(others) == 1 else 's', ', '.join(others[:3])))
                else:
                    bad.append('%s: differs' % rel)
        for dp, dn, fn in os.walk(os.path.join(data, sub)):
            for f in fn:
                if f.startswith('.'):
                    continue
                rel = os.path.relpath(os.path.join(dp, f), data)
                if not os.path.exists(os.path.join(out, rel)):
                    bad.append('%s: in the data folder, not in the output' % rel)
    print('%d files identical, %d with the same content in another format' % (same, fmt))
    for r in rot:
        print('  ' + r)
    for b in bad:
        print('  ' + b)
    if bad:
        print('THE OUTPUT DIFFERS FROM THE DATA FOLDER')
    elif rot:
        print('the output agrees with the data folder except for rotation systems of toroidal graphs')
    else:
        print('the output agrees with the data folder')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
