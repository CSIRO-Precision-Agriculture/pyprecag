Installation
============

pyprecag is available through the Python Packaging Index and can be installed with pip, however some dependencies require additional steps to install properly.

It is recommended to install pyprecag in a virtual environment so that the dependencies do not cause issues with your system-level Python installation.

VESPER Kriging is only supported on Windows platforms with the `VESPER <https://sydney.edu.au/agriculture/pal/software/vesper.shtml>`_ software installed.

Install via pip::

    pip install pyprecag

Linux
-----

The only dependency that causes issues is `GDAL <https://www.gdal.org/>`_ . The Python package is available from `PyPI <https://pypi.org/project/GDAL/>`_ .
However, the `libgdal-dev` dependencies are required, and the location of the header files needs to be specified when installing. These libraries are available via  `UbuntuGIS <https://wiki.ubuntu.com/UbuntuGIS>`_ and other avenues.

On Debian systems, this process should work.
Add the unstable release of UbuntuGIS, get and install packages with::

    sudo apt-get install software-properties-common
    sudo apt-add-repository ppa:ubuntugis/ubuntugis-unstable
    sudo apt-get update
    sudo apt-get install libgdal-dev gdal-bin

If the above doesn't work try removing all gdal related files from the `/etc/apt/sources.list.d` direcory and trying again (it's possible old files are interfering with install attempts).

Once the `libgdal-dev` libraries and `gdal-bin` are installed, pip install GDAL while specifying the version and the location of the header files with::

    pip install GDAL==$(gdal-config --version) --global-option=build_ext --global-option="-I/usr/include/gdal"

Then pyprecag should install without issue.

Docker
------

The Dockerfile included in the examples directory builds a container with Pyprecag installed in an Ubuntu 18.04 image. Build the container with::

    docker build -f docs/Dockerfile -t <tag> .

The tests can be run inside this container like this::

    docker run <container id or tag> make test

Windows
-------

The dependencies `GDAL <https://www.gdal.org/>`_, `Fiona <https://github.com/Toblerity/Fiona>`_ , and `Rasterio <https://github.com/mapbox/rasterio>`_ , all cause issues when installing on Windows.

The easiest way to install all of these is with the Windows Binaries provided by Christoph Gohlke at https://www.lfd.uci.edu/~gohlke/pythonlibs/. Download the Python 2.7 version for your system architecture then install them with::

    pip install /path/to/<downloaded_file.whl>

It is recommended to install GDAL first, as it is a requirements for both Fiona and Rasterio. GDAL might also return an error and require you to install Visual C++.

They are also available from the conda-forge channel in conda. pyprecag is not currently available to install with conda, however can be pip installed in a conda environment and use the conda-installed versions of GDAL, Fiona and Rasterio.

Once those three dependencies are installed, pyprecag should install without issue.
