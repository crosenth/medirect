#!/usr/bin/env python3
"""
Extract Genbank records into csv and fasta formats
"""
import argparse
import csv
import re
import sys

from Bio import SeqIO

# columns of output files
ANNOTATION_COLS = ['seqname', 'version', 'accession', 'name',
                   'description', 'tax_id', 'modified_date', 'download_date',
                   'version_num', 'source', 'keywords', 'organism', 'length',
                   'ambig_count', 'strain', 'mol_type', 'isolate',
                   'isolation_source', 'seq_start', 'seq_stop']

PUBMED_COLS = ['pubmed_id', 'version', 'accession']

REFERENCE_COLS = ['pubmed_id', 'title', 'authors',
                  'journal', 'consrtm', 'comment']

REFSEQ_INFO_COLS = ['seqname', 'accession', 'gi', 'seq_start', 'seq_stop']

# https://www.ncbi.nlm.nih.gov/Sequin/acc.html
REFSEQ = '[A-Z]{2}_\w+'
ACCESSION = '[A-Z]+\d+'
COORDINATES = ':(?P<seq_start>\d+)-(?P<seq_stop>\d+)'
REFSEQ_SOURCE = re.compile(
    'REFSEQ.*?(?P<accession>{REFSEQ}|{ACCESSION})({COORDINATES})?'.format(
        REFSEQ=REFSEQ, ACCESSION=ACCESSION, COORDINATES=COORDINATES),
    re.DOTALL)

GI_SOURCE = re.compile('gi:(?P<gi>\d+)')

ACGT = frozenset('ACGT')


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # ins
    parser.add_argument(
        'genbank',
        type=argparse.FileType('r'),
        default=sys.stdin,
        nargs='?',
        help='genbank record(s)')

    parser.add_argument(
        'download_date',
        help='date records were extracted were downloaded %%d-%%b-%%Y')

    # outs
    parser.add_argument(
        'fasta',
        type=argparse.FileType('w'),
        help='sequences')
    parser.add_argument(
        'annotations',
        type=dictwriter(ANNOTATION_COLS),
        help=str(ANNOTATION_COLS))
    parser.add_argument(
        'pubmed_info',
        type=dictwriter(PUBMED_COLS),
        help=str(PUBMED_COLS))
    parser.add_argument(
        'references',
        type=dictwriter(REFERENCE_COLS),
        help=str(REFERENCE_COLS))
    parser.add_argument(
        'refseq_info',
        type=dictwriter(REFSEQ_INFO_COLS),
        help=str(REFSEQ_INFO_COLS))

    return parser


def dictwriter(columns):
    def writer(filename):
        out = csv.DictWriter(
            open(filename, 'w'),
            fieldnames=columns,
            extrasaction='ignore')
        out.writeheader()
        return out
    return writer


def main():
    args = build_parser().parse_args()

    references = []

    for i, g in enumerate(SeqIO.parse(args.genbank, 'genbank')):
        sys.stderr.write('processing record ' + str(i) + '\r')

        try:
            record = parse_record(g)
            record.update({'download_date': args.download_date})

            refs = parse_references(g)
            for r in refs:
                r.update({
                    'version': record['version'],
                    'accession': record['accession']})
            references.extend(refs)

            args.annotations.writerow(record)

            # refseqs
            if '_' in record['accession']:
                row = parse_refseq_source(g)
                row.update({'seqname': record['seqname']})
                args.refseq_info.writerow(row)

            args.fasta.write('>{}\n{}\n'.format(record['seqname'], g.seq))
        except Exception as e:
            print(g)
            raise(e)

    # deduplicate and write references and pubmed_info
    ref_set = set(tuple(r[c] for c in REFERENCE_COLS) for r in references)
    args.references.writerows([dict(zip(REFERENCE_COLS, r)) for r in ref_set])

    pub_set = set(tuple(r[c] for c in PUBMED_COLS) for r in references)
    args.pubmed_info.writerows([dict(zip(PUBMED_COLS, r)) for r in pub_set])


def parse_coordinates(record):
    accessions = record.annotations['accessions']
    if 'REGION:' in accessions:
        coordinates = accessions[accessions.index('REGION:') + 1]
        matches = re.findall('\d+', coordinates)
        if len(matches) == 1:
            # some records are strange...
            seq_start, seq_stop = matches[0], len(record.seq)
        else:
            seq_start, seq_stop = matches
    else:
        seq_start, seq_stop = 1, len(record.seq)
    return seq_start, seq_stop


def parse_record(record):
    version, accession, version_num = parse_version(record)

    source = next(i for i in record.features if i.type == 'source')
    quals = source.qualifiers

    # get tax_id
    if 'db_xref' in quals:
        for i in quals.get('db_xref', []):
            if i.startswith('taxon:'):
                tax_id = i[6:]
                break
    # occasionally, tax_ids are missing
    else:
        tax_id = ''

    seq_start, seq_stop = parse_coordinates(record)
    if all([seq_start, seq_stop]):
        seq_id = '{}_{}_{}'.format(accession, seq_start, seq_stop)
    else:
        seq_id = record.id

    info = dict(accession=accession,
                ambig_count=sum(1 for b in record.seq if b not in ACGT),
                modified_date=record.annotations['date'],
                description=record.description,
                keywords=';'.join(record.annotations.get('keywords', [])),
                length=len(record),
                name=record.name,
                organism=record.annotations['organism'],
                seq_start=seq_start,
                seq_stop=seq_stop,
                seqname=seq_id,
                source=record.annotations['source'],
                mol_type=';'.join(quals.get('mol_type', '')),
                strain=';'.join(quals.get('strain', '')),
                isolate=';'.join(quals.get('isolate', '')),
                isolation_source=';'.join(quals.get('isolation_source', '')),
                tax_id=tax_id,
                version=version,
                version_num=version_num)

    return info


def parse_references(record):
    """
    Parse reference annotations that have a pubmed_id
    """
    references = []
    if 'references' in record.annotations:
        refs = [r for r in record.annotations['references'] if r.pubmed_id]
        for r in refs:
            references.append(
                dict(title=r.title,
                     authors=r.authors,
                     comment=r.comment,
                     consrtm=r.consrtm,
                     journal=r.journal,
                     pubmed_id=r.pubmed_id))
    return references


def parse_refseq_source(record):
    acc = re.search(REFSEQ_SOURCE, record.annotations['comment'])
    gi = re.search(GI_SOURCE, record.annotations['comment'])

    if not acc:
        raise ValueError('Cannot parse record')

    result = {}
    if acc:
        result.update(acc.groupdict())
    if gi:
        result.update(gi.groupdict())

    return result


def parse_version(record):
    """
    Return the accession and version of a Bio.SeqRecord.SeqRecord
    """
    annotations = record.annotations
    accession = annotations.get('accessions', [''])[0]
    if accession:
        if 'sequence_version' in annotations:
            version_num = str(annotations.get('sequence_version'))
        elif record.id.startswith(accession + '.'):
            version_num = record.id.split('.', 1)[1]
        else:
            version_num = '1'
        version = accession + '.' + version_num
    else:
        version = ''
    return version, accession, version_num


if __name__ == '__main__':
    main()
