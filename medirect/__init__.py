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
#    along with medirect.  If not, see <http://www.gnu.org/licenses/>
'''
main medirect object, all commands must extend this class
'''
import argparse
import logging
import sys


class MEDirect:
    TOOL = 'medirect'

    def __init__(self, testing=None):
        if testing is None:
            parser = self.add_arguments(self.arg_parser())
            args, other_args = parser.parse_known_args()
            self.log = args.log
            self.verbosity = args.verbosity
            self.main(args, *other_args)

    def add_arguments(self, parser):
        '''
        define argparse arguments here and return the arg_parser
        '''
        raise NotImplementedError(
            'add_arguments must be implemented when extending MEDirect')

    def arg_parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            '-log', '--log',
            metavar='FILE',
            help='Send logging to a file')
        parser.add_argument(
            '-v', '--verbose',
            action='count',
            dest='verbosity',
            default=1,
            help='Increase verbosity of screen output '
                 '(eg, -v is verbose, -vv more so)')
        parser.add_argument(
            '-out', '--out',
            metavar='',
            default=sys.stdout,
            type=argparse.FileType('w'),
            help='Output location [sys.stdout]')
        return parser

    def main(self, args, *other_args):
        raise NotImplementedError(
            'main must be implemented when extending MEDirect')

    def setup_logging(self):
        """
        setup global logging
        """
        if self.log:
            logfile = open(self.log, 'a')
        else:
            logfile = sys.stderr

        loglevel = {
            0: logging.ERROR,
            1: logging.WARNING,
            2: logging.INFO,
            3: logging.DEBUG,
        }.get(self.verbosity, logging.DEBUG)

        log_format = ('%(asctime)s %(levelname)s %(lineno)s %(message)s')

        logging.basicConfig(
            stream=logfile,
            format=log_format,
            level=loglevel,
            datefmt='%Y-%m-%d %H:%M:%S')
