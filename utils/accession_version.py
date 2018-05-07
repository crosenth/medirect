#!/usr/bin/env python3
"""
parse accession.version from ncbi feature table csv id column
"""

import argparse
import csv
import sys


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        'csv',
        type=argparse.FileType('r'),
        nargs='?',
        default=sys.stdin,
        help='csv file from mefetch with column id')
    p.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='list of all version downloaded')

    args = p.parse_args()
    reader = csv.DictReader(args.csv)
    writer = csv.DictWriter(args.out, reader.fieldnames)
    writer.writeheader()
    for r in reader:
        r['id'] = r['id'].split('|')[1]
        writer.writerow(r)


if __name__ == '__main__':
    main()
