=====================================
medirect: multiprocessed ncbi edirect
=====================================

medirect - a multiprocessed utility for retrieving records and parsing feature tables from ncbi

.. contents:: Table of Contents

authors
=======

* `Chris Rosenthal <crosenth@gmail.com>`_

about
=====

As a bioinformatician I build a lot of nucleic bacteria reference databases.  I 
created this package to help do that quickly.  The utilities in this package 
(mefetch and ffetch) can be slotted in along with other ncbi utilities and 
follow the same edirect
`documenation <https://www.ncbi.nlm.nih.gov/books/NBK25501/>`_,
`guidelines and requirements <https://www.ncbi.nlm.nih.gov/books/NBK25497/#_chapter2_Usage_Guidelines_and_Requiremen_>`_
and
`usage policies <https://www.ncbi.nlm.nih.gov/home/about/policies.shtml>`_.

For large data requests I highlight two points from the usage policy:

* Run retrieval scripts on weekends or between 9 pm and 5 am Eastern Time weekdays for any series of more than 100 requests.
* Make no more than 3 requests every 1 second.

Some additional documentation for using ffetch:

* `feature tables <http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html>`_

edirect `ftp downloads <https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/>`_

dependencies
============

* Python 3.x
* `biopython <https://pypi.python.org/pypi/biopython>`_ >= 1.68
* `retrying <https://pypi.python.org/pypi/retrying>`_ >= 1.3.3

installation
============

medirect can be installed in two ways:

For regular users::

  % pip3 install medirect

For developers::

  % pip3 install git://github.com/crosenth/medirect.git
  # or
  % git clone git://github.com/crosenth/medirect.git 
  % cd medirect
  % python3 setup.py install

examples
========

The mefetch executable works exactly like edirect
`efetch <https://www.ncbi.nlm.nih.gov/books/NBK179288/efetch>`_ with a an 
additional multiprocessing argument ``-proc`` and a few more features.

By allowing additional processes to download records The ``-proc`` argument 
allows a linear download speed increase downloading large datasets.

For example::

  % esearch -db nucleotide -query 'Rhizobium' | time mefetch -email user@ema.il -mode text -format acc -proc 1 > accessions.txt
  0.53s user 0.11s system 0% cpu 11:43.11 total

Which is equivalent to ncbi efetch::

  % esearch -db nucleotide -query 'Rhizobium' | time efetch -mode text -format acc > accessions.txt
  0.53s user 0.11s system 0% cpu 12:47.54 total

Adding another processor ``-proc 2``::

  % esearch -db nucleotide -query 'Rhizobium' | time mefetch -email user@ema.il -proc 2 -mode text -format acc > accessions.txt
  0.46s user 0.08s system 0% cpu 5:17.51 total

And another ``-proc 3`` (default)::

  % esearch -db nucleotide -query 'Rhizobium' | time mefetch -email user@ema.il -proc 3 -mode text -format acc > accessions.txt
  0.35s user 0.10s system 0% cpu 2:57.01 total

And ``-proc 4`` (see usage policy)::

  % esearch -db nucleotide -query 'Rhizobium' | time mefetch -email user@ema.il -proc 4 -mode text -format acc > accessions.txt
  0.35s user 0.08s system 0% cpu 1:40.54 total

Results can be returned in the same order as efetch using the ``-in-order``
argument.  Otherwise, the order will be determined by how fast ncbi returns
results per process.

The ``-retmax`` argument (or chunksize) determines the number of results 
returned per ``-proc``.  By default, it is set to the 10,000 max records per
`documentation <https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_>`_.
Setting the ``-retmax`` to higher than 10,000 will automatically be set back
down to 10,000.

By default the ``-id`` reads stdin xml output from ``esearch``.  The ``-id`` 
argument can also take input as a comma delimited list of ids or text
file of ids.  When coupled with the ``-csv`` argument the input can be a csv
file with additional argument columns.  This is useful for bulk downloads with 
different positional arguments.

ftract allows csv output of different features from ncbi
`feature tables <http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html>`_.
The required ``-feature`` argument is comma separated
feature_key:qualifier_key:qualifier_value
::

  % mefetch -id KN150849 -db nucleotide -email user@ema.il -format ft | ftract --feature rrna::16s
  id,seq_start,seq_stop,strand
  KN150849.1,594136,595654,2
  KN150849.1,807985,809503,2
  KN150849.1,2227751,2229271,1

And pipe this back into mefetch to download these three regions in genbank format::

  % mefetch -id KN150849 -db nucleotide -email user@ema.il -format ft | ftract --feature rrna:product:16s | mefetch -db nucleotide -email crosenth@uw.edu -csv -format gb


And finally, return all the Burkholderia gladioli 16s rrna products in fasta format like this::

  % esearch -query 'Burkholderia gladioli' -db 'nucleotide' | mefetch -email user@ema.il -format ft | ftract --feature rrna:product:16s | mefetch -db nucleotide -email user@ema.il -csv -format fasta

issues
======

Please use the Issue Tracker(s) available on Github or Bitbucket to report any bugs
or feature requests.  For all other inquiries email `Chris Rosenthal <crosenth@gmail.com>`_.

license
=======

Copyright (c) 2016 Chris Rosenthal

Released under the `GPLv3 <http://www.gnu.org/copyleft/gpl.html>`_ License
