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
import logging
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
            help='Parse only specific record features in feature '
                 'table columns via string pattern '
                 'feature_key:qualifier_key:qualifier_value.'
                 'ex: rRNA:product:16S')
        parser.add_argument(
            '-full-format', '--full-format',
            action='store_true',
            help=('If specified, also output the '
                  'feature, qualifier, and qualifier value'))
        parser.add_argument(
            '-min-length', '--min-length',
            type=int,
            help='Minimum sequence length to parse feature')
        parser.add_argument(
            '-on-error', '--on-error',
            choices=['halt', 'continue'],
            default='continue',
            help='If an exception is encountered in a feature table '
                 'line encountered halt or continue? [%(default)s]')
        return parser

    def coordinates(self, start, stop, strand='1'):
        if stop < start:
            start, stop, strand = stop, start, '2'
        return start, stop, strand

    def filter_features(
            self, records, features, on_error, full_format, min_length):
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
            # any empty feature overrides all other features
            features = ['' if '' in f else f for f in features]
            features = [set(x.lower() for x in f) for f in features]
            features = ['|'.join(f) for f in features]
        else:
            features = ['', '', '']

        feat_key_pattern, qual_key_pattern, qual_val_pattern = features

        # general patterns for all 5 feature table columns
        column1 = '[<>]?(?P<seq_start>\d+)'
        column2 = '[<>]?(?P<seq_stop>\d+)'
        column3 = '(?P<feature_key>{})'.format(feat_key_pattern or '.*')
        column4 = '(?P<qualifier_key>{})'.format(qual_key_pattern or '.*')
        column5 = '(?P<qualifier_value>{})'.format(qual_val_pattern or '.*')

        # general line1 coordinates pattern
        line1 = re.compile('^{}\t{}'.format(column1, column2))

        # these three patterns look for features passed by user
        seqid_line = re.compile('^>Feature (?P<seqid>\S*)', re.IGNORECASE)
        coordinates_feature_key = re.compile(
            '^{}\t{}\t?{}'.format(column1, column2, column3), re.IGNORECASE)
        qualifiers = re.compile(
            '^\t\t\t{}\t{}'.format(column4, column5), re.IGNORECASE)

        # iterate and match lines
        seqid, seq_start, seq_stop = None, None, None
        feature, qual, qual_value = None, None, None
        for line in records:
            '''
            Match for line2 first because that is the most
            common line in a feature table.
            '''
            if line.startswith('\t'):
                if seq_start and seq_stop and re.search(qualifiers, line):
                    if full_format:
                        qual_match = qualifiers.search(line)
                        qual = qual_match.group('qualifier_key')
                        qual_value = qual_match.group('qualifier_value')
                        yield(
                            seqid,
                            *self.coordinates(seq_start, seq_stop),
                            feature,
                            qual,
                            qual_value)
                    else:
                        yield (seqid, *self.coordinates(seq_start, seq_stop))
                    qual, qual_value = None, None
            elif re.search(line1, line):
                match = re.search(coordinates_feature_key, line)
                feature, qual, qual_value = None, None, None
                if match:
                    seq_start = int(match.group('seq_start'))
                    seq_stop = int(match.group('seq_stop'))
                    if (min_length is None or
                            abs(seq_stop - seq_start) >= min_length):
                        feature = match.group('feature_key')
                        if not (qual_key_pattern or qual_val_pattern):
                            # immediately yield if no qualifier patterns
                            yield (
                                seqid,
                                *self.coordinates(seq_start, seq_stop))
                            seq_start, seq_stop = None, None
                    else:
                        seq_start, seq_stop = None, None
                else:
                    seq_start, seq_stop = None, None
            elif line.startswith('>Feature'):
                match = re.search(seqid_line, line)
                seqid = match.group('seqid')
                seq_start, seq_stop = None, None
                feature, qual, qual_value = None, None, None
            else:
                msg = str(seqid) + ' contains invalid line: ' + line
                logging.error(msg)
                if on_error == 'halt':
                    raise ValueError(msg)

    def main(self, args, *other_args):
        self.setup_logging()
        out = csv.writer(args.out)
        if args.full_format:
            out.writerow(['id', 'seq_start', 'seq_stop', 'strand',
                          'feature', 'qual', 'qual_value'])
        else:
            out.writerow(['id', 'seq_start', 'seq_stop', 'strand'])
        # remove any blank lines
        table = (line for line in args.table if line.strip())
        ff = self.filter_features(
            table,
            args.features,
            args.on_error,
            args.full_format,
            args.min_length)
        for f in ff:
            out.writerow(f)


def run():
    Ftract()
