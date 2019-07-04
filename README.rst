========
pyprecag
========

A suite of tools for Precision Agriculture data analysis

.. image:: https://travis-ci.com/CSIRO-Precision-Agriculture/pyprecag.svg?branch=master
   :target: https://travis-ci.com/CSIRO-Precision-Agriculture/pyprecag

* Homepage: https://github.com/CSIRO-Precision-Agriculture/pyprecag
* Documentation: https://CSIRO-Precision-Agriculture.github.io/pyprecag_docs
* `Installation <https://csiro-precision-agriculture.github.io/pyprecag_docs/installation.html#installation>`__
* Version: 0.3.0
* License: `CSIRO Open Source Software License Agreement <https://github.com/CSIRO-Precision-Agriculture/pyprecag/blob/master/LICENSE>`__

pyprecag is currently supported only on Python version 2.7

.. include:: ./CITATION.rst

Dependencies
------------

* Python version 2.7
* `GDAL <https://www.gdal.org/>`_
* `Fiona <https://github.com/Toblerity/Fiona>`_
* `Rasterio <https://github.com/mapbox/rasterio>`_
* `VESPER <https://sydney.edu.au/agriculture/pal/software/vesper.shtml>`_ (for Kriging)

.. include:: ./docs/installation.rst

Development
===========

There is a makefile in the project root with targets for the most common
development operations such as lint checks, running unit tests, building the
documentation, and building installer packages. `tox` does not have a target,
as `make tox` is more typing than `tox`.

Run make with no target to see the list of targets::

    $ make

`Bumpversion <https://pypi.python.org/pypi/bumpversion>`_ is used to manage the
package version numbers. This ensures that the version number is correctly
incremented in all required files. Please see the bumpversion documentation for
usage instructions, and do not edit the version strings directly.

Version numbers follow the `Semantic versioning guidelines <semver.org>`_.

.. include:: ./contributing.rst


