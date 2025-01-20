========
pyprecag
========
A suite of tools for Precision Agriculture data analysis

* Homepage: https://github.com/CSIRO-Precision-Agriculture/pyprecag
* Documentation: https://CSIRO-Precision-Agriculture.github.io/pyprecag_docs
* `Installation <https://csiro-precision-agriculture.github.io/pyprecag_docs/installation.html#installation>`__
* Version: 0.4.3.rc3
* License: `CSIRO Open Source Software License Agreement <https://github.com/CSIRO-Precision-Agriculture/pyprecag/blob/master/LICENSE>`__

pyprecag is currently supported for Python version 2.7 and 3.7

.. include:: ./CITATION.rst

Dependencies
------------

* Python version 2.7, 3.7
* `GDAL <https://www.gdal.org/>`_
* `Fiona <https://fiona.readthedocs.io/en/stable/`_
* `Rasterio <https://rasterio.readthedocs.io/en/stable/>`_
* `Geopandas <https://geopandas.org/en/stable/`_
* `VESPER <https://precision-agriculture.sydney.edu.au/resources/software/>`_ (for Kriging)

.. include:: ./docs/installation.rst

Development
===========

There is a makefile in the project root with targets for the most common
development operations such as lint checks, running unit tests, building the
documentation, and building installer packages. `tox` does not have a target,
as `make tox` is more typing than `tox`.

Run make with no target to see the list of targets:

.. code-block:: bash
    $ make

`Bump-my-version <https://callowayproject.github.io/bump-my-version/>`_ is used to manage the
package version numbers. This ensures that the version number is correctly
incremented in all required files. Please see the bumpversion documentation for
usage instructions, and do not edit the version strings directly.

Version numbers follow the `Semantic versioning guidelines <semver.org>`_.

.. include:: ./contributing.rst


