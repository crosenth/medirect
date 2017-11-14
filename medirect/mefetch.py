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

import csv
import functools
import itertools
import logging
import medirect
import multiprocessing
import os
import retrying
import sys

from Bio import Entrez
from html.parser import HTMLParser


class MEFetch(medirect.MEDirect):
    RETMAX = 10000

    def add_arguments(self, parser):
        parser.add_argument(
            '-email', '--email',
            required=True,
            help="The name of the run")
        parser.add_argument(
            '-debug', '--debug',
            action='store_true',
            help='no multiprocessing')
        parser.add_argument(
            '-in-order', '--in-order',
            action='store_true',
            help=('Order of results determined by ncbi. Note, '
                  'processor utilization can become bottlenecked '
                  'waiting for slower ncbi request return. [%(default)s]'))

        input_group = parser.add_argument_group(title='input')
        input_group.add_argument(
            '-id', '--id',
            metavar='',
            default=sys.stdin,
            help='File or comma separated list of ids to fetch [stdin]')
        input_group.add_argument(
            '-csv', '--csv',
            action='store_true',
            help='If the input is in csv format')

        input_group.add_argument(
            '-format', '--format',
            metavar='',
            help='Format of record or report')
        input_group.add_argument(
            '-mode', '--mode',
            metavar='',
            help='text, xml, asn.1, json')

        proc_group = parser.add_argument_group(title='processing')
        retry_group = proc_group.add_mutually_exclusive_group()
        retry_group.add_argument(
            '-max-retry', '--max-retry',
            metavar='INT',
            type=lambda x: None if x == '-1' else int(x),
            default=10,
            help=('Max number of retries after consecutive '
                  'http exceptions [%(default)s].  Use -1 for '
                  'continuous retrying.'))
        proc_group.add_argument(
            '-proc', '--proc',
            default=3,
            metavar='INT',
            type=int,
            help='Number of available processors [%(default)s]')
        proc_group.add_argument(
            '-retmax', '--retmax',
            metavar='INT',
            default=self.RETMAX,
            type=int,
            help=('number of records returned per request '
                  'or chunksize [%(default)s]'))
        proc_group.add_argument(
            '-retry', '--retry',
            metavar='MILLISECONDS',
            type=int,
            default=60000,
            help=('Number of milliseconds to wait '
                  'between -max-retry(ies) [%(default)s]'))
        return parser

    def main(self, args, *unknown_args):
        Entrez.email = args.email

        # zip unknown args to be passed as general efetch args
        unknown_args = [a.strip('-') for a in unknown_args]
        base_args = dict(zip(unknown_args, unknown_args[1:]))

        # piped in csv
        if args.id is sys.stdin and args.csv:
            retmax = 1
            ids = (i for i in args.id if not i.startswith('#'))
            chunks = csv.DictReader(ids)
        # piped in edirect xml file
        elif args.id is sys.stdin:
            retmax = min(self.RETMAX, args.retmax)  # 10k is ncbi max
            edirect = parse_edirect(args.id.read())
            count = int(edirect['count'])
            base_args = dict(base_args, **edirect)
            chunks = ({'retstart': r} for r in range(0, count, retmax))
        # -id csv file
        elif os.path.isfile(args.id) and args.csv:
            retmax = 1
            ids = (i for i in open(args.id) if not i.startswith('#'))
            chunks = csv.DictReader(ids)
        # -id txt file of ids
        elif os.path.isfile(args.id):
            retmax = min(self.RETMAX, args.retmax)  # 10k is ncbi max
            ids = (i for i in open(args.id) if not i.startswith('#'))
            ids = (i.strip() for i in ids)
            chunks = ({'id': i} for i in chunker(ids, retmax))
        # comma separated list of ids
        else:
            retmax = min(self.RETMAX, args.retmax)  # 10k is ncbi max
            ids = (i.strip() for i in args.id.split(','))
            chunks = ({'id': i} for i in chunker(ids, retmax))

        # add efetch retmax argument
        base_args.update(retmax=retmax)

        # convert these efetch arguments to Entrez.efetch arguments
        if args.format:
            base_args.update(rettype=args.format)
        if args.mode:
            base_args.update(retmode=args.mode)

        efetches = functools.partial(
            efetch, args.retry, args.max_retry, **base_args)

        pool = multiprocessing.Pool(processes=args.proc)

        if args.in_order:
            results = pool.imap(efetches, chunks)
        elif args.debug:
            results = map(efetches, chunks)
        else:
            results = pool.imap_unordered(efetches, chunks)

        for r in results:
            # remove blank lines and append a single newline
            r = '\n'.join(l for l in r.split('\n') if l.strip()) + '\n'
            args.out.write(r)


def liststr(l):
    """
    Concise string representation of iteratable
    """

    s = None
    if hasattr(l, '__iter__') and not isinstance(l, str):
        if len(l) == 0:
            s = 'None'
        elif len(l) == 1:
            s = str(l[0])
        else:
            s = '{}...{}'.format(l[0], l[-1])
    else:
        s = l
    return s


def edirect_pprint(**kwds):
    """
    Output Entrez commands in command line format:
    http://www.ncbi.nlm.nih.gov/books/NBK179288/
    """

    entrez = dict(id='-id', db='-db', rettype='-format', retmode='-mode',
                  strand='-strand', complexity='-complexity',
                  seq_start='-seq_start', seq_stop='-seq_stop')
    keys = set(entrez.keys()).intersection(list(kwds.keys()))
    info = ['efetch']
    info += ['{} {}'.format(entrez[k], kwds[k]) for k in keys]
    return ' '.join(info)


def efetch(retry, max_retry, chunks, **args):
    """
    Wrap Entrez.efetch with some http exception retrying
    """

    def print_retry_message(exception):
        """
        http exceptions: http://www.w3.org/Protocols/rfc2616/rfc2616-sec10.html
        retrying api: https://pypi.python.org/pypi/retrying
        """
        seconds = float(retry) / 1000
        msg = '{}, retrying in {} seconds... {} max retry(ies)'.format(
            exception, seconds, max_retry or 'no')
        logging.error(msg)
        return True

    @retrying.retry(
        retry_on_exception=print_retry_message,
        wait_fixed=retry,
        stop_max_attempt_number=max_retry)
    def rfetch(chunk, **args):
        pprint_chunk = dict((k, liststr(v)) for k, v in chunk.items())
        logging.info(edirect_pprint(**dict(pprint_chunk, **args)))
        args.update(**chunk)
        db = args.pop('db')
        return Entrez.efetch(db, **args).read()

    return rfetch(chunks, **args)


def chunker(iterable, n):
    """
    Continuously chunk an iterator n items at a time
    """

    while True:
        chunk = list(itertools.islice(iterable, n))
        if chunk:
            yield chunk
        else:
            return


def parse_edirect(text):
    class EsearchParser(HTMLParser):
        def __init__(self):
            HTMLParser.__init__(self)
            self.data = {}

        def handle_starttag(self, tag, attrs):
            if tag == 'querykey':
                self.starttag = 'query_key'
            else:
                self.starttag = tag

        def handle_data(self, data):
            data = data.strip()
            if data:
                self.data[self.starttag] = data

        def get_data(self):
            return self.data
    parser = EsearchParser()
    parser.feed(text)
    return parser.get_data()


def run():
    MEFetch()
