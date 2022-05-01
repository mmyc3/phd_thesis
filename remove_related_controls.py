#!/resources/tools/apps/software/lang/Python/3.7.4-GCCcore-8.3.0/bin/python
# Code written by Catalin Voinescu, UCL

import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Remove related controls, keep cases')
    parser.add_argument('cases', type=argparse.FileType('r'),
        help='File of case IDs, tab-separated, two columns: FID and ID')
    parser.add_argument('controls', type=argparse.FileType('r'),
        help='File of control IDs, tab-separated, two columns: FID and ID')
    parser.add_argument('king', type=argparse.FileType('r'), default=sys.stdin,
        help='File created by KING')
    parser.add_argument('keep', type=argparse.FileType('w'), default=sys.stdout,
        help='Output. File of IDs to be kept, tab-separated, two columns: FID and ID')
    parser.add_argument('--verbose', '-v', action='store_true',
        help='List removed controls')

    args = parser.parse_args()

    cases = []
    for line in args.cases:
        line = line.strip()
        if not line:
            continue
        fid, id = line.split()
        cases.append((fid, id))

    controls = []
    for line in args.controls:
        line = line.strip()
        if not line:
            continue
        fid, id = line.split()
        controls.append((fid, id))

    case_set = set(cases)
    control_set = set(controls)
    remove_set = set()

    header = True
    for line in args.king:
        line = line.strip()
        if not line:
            continue
        if header and line.startswith('FID1\tID1\t'):
            header = False
            continue

        fid1, id1, fid2, id2, rest = line.split('\t', 4)
        key1 = (fid1, id1)
        key2 = (fid2, id2)

        if key1 in remove_set or key2 in remove_set:
            continue  # one of the two has already been removed, no conflict

        if key1 in case_set:
            remove_set.add(key2)
        else:
            remove_set.add(key1)

    if args.verbose:
        for (fid, id) in controls:
            if (fid, id) in remove_set:
                print('Removed control', fid, id, file=sys.stderr)

    for (fid, id) in cases:
        if (fid, id) in remove_set:
            print('*** REMOVED CASE ***', fid, id, file=sys.stderr)

    print('Got:', len(cases), 'cases and', len(controls), 'controls', file=sys.stderr)
    print('Removed:', len(case_set & remove_set), 'cases,',
        len(control_set & remove_set), 'controls.', file=sys.stderr)

    # Remove the samples. We do it here, rather than subtract from sets,
    # to keep the same order in the original case and control files.
    cases = [key for key in cases if key not in remove_set]
    controls = [key for key in controls if key not in remove_set]

    for (fid, id) in cases:
        print(fid + '\t' + id, file=args.keep)
    for (fid, id) in controls:
        print(fid + '\t' + id, file=args.keep)


if __name__ == "__main__":
    main()
