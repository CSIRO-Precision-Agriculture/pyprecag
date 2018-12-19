===============================
pyprecag
===============================

A suite of tools for Precision Agriculture

* Free software: CSIRO Open Source Software License Agreement
* Homepage: https://github.com/CSIRO-Precision-Agriculture/pyprecag
* Documentation: https://pyprecag.readthedocs.org
* Version: 0.0.4

pyprecag supports Python version 2.7


.. include:: ./installation.rst

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

Dependancies
============

* Python version 2.7
* `GDAL <https://www.gdal.org/>`_
* `Fiona <https://github.com/Toblerity/Fiona>`_
* `Rasterio <https://github.com/mapbox/rasterio>`_
* `Vesper <https://sydney.edu.au/agriculture/pal/software/vesper.shtml>`_ (for Kriging)
