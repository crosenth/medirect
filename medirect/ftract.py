# This file is part of medirect
#
#    medirect is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    medirect is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with medirect.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import csv
import medirect
import re
import sys


class Ftract(medirect.MEDirect):

    def add_arguments(self, parser):
        parser.add_argument(
            'table',
            nargs='?',
            default=sys.stdin,
            type=argparse.FileType('r'),
            help="An ncbi feature table in text format")
        parser.add_argument(
            '-feature', '--feature',
            action='append',
            dest='features',
            help='parse only specific record features in feature '
                 'table columns via string pattern '
                 'feature_key:qualifier_key:qualifier_value.'
                 'ex: rRNA:product:16S')
        return parser

    def filter_features(self, records, features):
        """
        http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html
        Parsing a five column, tab delimited file with a fasta file
        header line.
        example -
        > seqid
        line 1
        column 1: seq_start
        column 2: seq_stop
        column 3: feature_key
        line 2
        column 4: qualifier_key
        column 5: qualifier_value
        """

        # parse features for columns 3-5
        if features:
            features = [f.split(':') for f in features]
            for f in features:
                # make sure features are valid
                if len(f) != 3:
                    msg = str(f) + ' is not a valid feature argument'
                    raise argparse.ArgumentTypeError(msg)
            features = zip(*features)
            features = ['' if '' in f else f for f in features]
            features = [[x if x else '.' for x in f] for f in features]
            features = ['|'.join(f) for f in features]
        else:
            features = ['', '', '']

        # create column patterns
        column1 = '\D?(?P<seq_start>\d+)'
        column2 = '\D?(?P<seq_stop>\d+)'
        column3, column4, column5 = features
        column3 = '.*?(?P<feature_key>{})+.*?'.format(column3)
        column4 = '.*?(?P<qualifier_key>{})+.*?'.format(column4)
        column5 = '.*?(?P<qualifier_value>{})+.*?'.format(column5)

        # three line types, lines 1 and 2 are tab delimited
        seqid_line = re.compile(
            '^>Feature (?P<seqid>.*)', re.IGNORECASE)
        line1 = re.compile(
            '^{}\t{}\t{}'.format(column1, column2, column3), re.IGNORECASE)
        line2 = re.compile(
            '^\t\t\t{}\t{}'.format(column4, column5), re.IGNORECASE)

        # iterate and match lines
        seqid, seq_start, seq_stop = None, None, None
        yielded = False
        for line in records:
            match = re.search(seqid_line, line)
            if match:
                seqid = match.group('seqid')
                continue

            match = re.search(line2, line)
            if match and not yielded:
                if seq_stop < seq_start:
                    # swap coordinates and switch strands
                    yield seqid, seq_stop, seq_start, '2'
                else:
                    yield seqid, seq_start, seq_stop, '1'
                yielded = True
                continue

            match = re.search(line1, line)
            if match:
                seq_start = int(match.group('seq_start'))
                seq_stop = int(match.group('seq_stop'))
                yielded = False
                continue

    def main(self, args, *other_args):
        out = csv.writer(args.out)
        out.writerow(['id', 'seq_start', 'seq_stop', 'strand'])
        for f in self.filter_features(args.table, args.features):
            out.writerow(f)


def run():
    Ftract()
