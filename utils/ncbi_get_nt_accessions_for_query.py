#!/usr/bin/env python3
import argparse
import subprocess
import sys

"""
    Use etools and metools to search NCBI NT for accessions matching a query. Some useful queries to start with:
    16s rRNA for bacteria with some classified name: "16s[All Fields] AND rRNA[Feature Key] AND Bacteria[Organism] AND 500 : 99999999999[Sequence Length] NOT(environmental samples[Organism] OR unclassified Bacteria[Organism])"
    16s rRNA for type-organism bacteria: "16s[All Fields] AND rRNA[Feature Key] AND Bacteria[Organism] AND 500 : 99999999999[Sequence Length] AND sequence_from_type[Filter]"
    16s rRNA for TM7: "16s[All Fields] AND rRNA[Feature Key] AND Bacteria[Organism] AND 500 : 99999999999[Sequence Length]  AND Candidatus Saccharibacteria[Organism]"
"""


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # ins
    parser.add_argument(
        '--query',
        type=str,
        required=True,
        help='Query to search for in the NCBI NT database')

    parser.add_argument(
        '--email',
        type=str,
        required=True,
        help='email to pass along to NCBI as the requestor')

    parser.add_argument(
        '--retry_max',
        type=int,
        default=10,
        help='Maximum number of retries when attempting to retrieve an accession.')

    parser.add_argument(
        '--retry_delay',
        type=int,
        default=60000,
        help='Miliseconds to delay between retries.')

    parser.add_argument(
        '--ncbi_concurrent_connections',
        type=int,
        default=3,
        help='Number of concurrent connections to NCBI when retrieving.')


    # Out
    parser.add_argument(
        '--log',
        type=str,
        help='(Optional) Logfile location')

    parser.add_argument(
        '--out',
        type=argparse.FileType('a'),
        default=sys.stdout,
        help='Accessions IDs that met criteria. Default stdout')

    return parser

if __name__ == "__main__":

    args = build_parser().parse_args()
    # First make our search using esearch, saving the stdout
    search_process = subprocess.Popen(
        [
            'esearch',
            '-db', 'nucleotide',
            '-query', args.query,

        ],
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    # Communicate returns the stdout / stderr as strings. Fine for the short record pointing to the search result
    search_stdout, search_stderr = search_process.communicate()

    fetch_command=[
            'mefetch',
            '-vv',
            '-email', args.email,
            '-mode', 'text',
            '-format', 'acc',
            '-mat-retry', str(args.retry_max),
            '-retry', str(args.retry_delay),
            '-proc', str(args.ncbi_concurrent_connections),
        ]
    if args.log:
        fetch_command+=['-log', args.log]

    fetch_process = subprocess.Popen(
        fetch_command,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    # Write our search result to the stdin of the fetch process
    fetch_process.stdin.write(search_stdout)
    # Close it to signal to the fetch process we are done
    fetch_process.stdin.close()
    # Stream directly to the output to avoid memory issues given the large size of the possible records
    args.out.write(fetch_process.stdout.read())

