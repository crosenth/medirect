=====================================
medirect: ncbi edirect multiprocessed
=====================================

medirect - a multiprocessed utility for retrieving records and parsing feature tables from ncbi

.. contents:: Table of Contents

authors
=======

* Chris Rosenthal

about
=====

I build nucleotide bacteria reference databases for a living.  I built this 
package to do that quickly.  The utilities in this package can be slotted
in along with other ncbi utilities and should follow the same 
`guidelines and requirements <https://www.ncbi.nlm.nih.gov/books/NBK25497/#_chapter2_Usage_Guidelines_and_Requiremen_>`_
and 
`usage policies <https://www.ncbi.nlm.nih.gov/home/about/policies.shtml>`_.

For large data requests I am highlighting two points from the usage policy:

* Run retrieval scripts on weekends or between 9 pm and 5 am Eastern Time weekdays for any series of more than 100 requests.
* Make no more than 3 requests every 1 second.

Some additional documentation:

* `feature tables <http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html>`_
* `ncbi E-utilities <https://www.ncbi.nlm.nih.gov/books/NBK25501/>`_

dependencies
============

* Python 3.x
* `biopython <https://pypi.python.org/pypi/biopython>`_ >= 1.68
* `retrying <https://pypi.python.org/pypi/retrying>`_ >= 1.3.3

installation
============

medirect can be installed in two ways:

For regular users:
::

  % pip install medirect

For developers::

  % git clone git://github.com/crosenth/medirect.git 
  % cd medirect
  % python setup.py install
  # or
  % pip install -U .

license
=======

Copyright (c) 2016 Chris Rosenthal

Released under the `GPLv3 <http://www.gnu.org/copyleft/gpl.html>`_ License
