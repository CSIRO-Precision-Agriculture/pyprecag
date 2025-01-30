# pyprecag

A suite of tools for Precision Agriculture data analysis

-   Homepage: <https://github.com/CSIRO-Precision-Agriculture/pyprecag>
-   Documentation:
    <https://CSIRO-Precision-Agriculture.github.io/pyprecag_docs>
-   [Installation](https://csiro-precision-agriculture.github.io/pyprecag_docs/installation.html#installation)
-   Version: 0.4.3
-   License: [GNU General Public Licence version 3
    (GPLv3)](https://github.com/CSIRO-Precision-Agriculture/pyprecag/blob/master/LICENSE)

Precision Agriculture Tools (PAT) Plugin for QGIS contains tools for
processing and analysing precision agriculture data which uses the
pyPrecAg python module. For more infomation visit
<https://github.com/CSIRO-Precision-Agriculture/PAT_QGIS_Plugin>

## Citation

If you use pyprecag for any published work, please cite using the
reference below:

*Ratcliff, Christina; Gobbett, David; & Bramley, Rob (2025): pyprecag -
a Python package for the analysis of precision agriculture data. CSIRO.
v5. Software. http://hdl.handle.net/102.100.100/78353*

## Dependencies

-   Python version 3.7+
-   [GDAL](https://www.gdal.org/)
-   [Fiona](https://fiona.readthedocs.io/en/stable/)
-   [Rasterio](https://rasterio.readthedocs.io/en/stable/)
-   [Geopandas](https://geopandas.org/en/stable/)
-   [VESPER](https://precision-agriculture.sydney.edu.au/resources/software/)
    (for Kriging)

### Installation

pyprecag is available through the Python Packaging Index and can be
installed with pip, however some dependencies require additional steps
to install properly.

It is recommended to install pyprecag in a virtual environment so that
the dependencies do not cause issues with your system-level Python
installation.

VESPER Kriging is only supported on Windows platforms with the
[VESPER](https://sydney.edu.au/agriculture/pal/software/vesper.shtml)
software installed.

Install via pip:

``` console
pip install pyprecag
```

## Linux

The only dependency that causes issues is [GDAL](https://www.gdal.org/)
. The Python package is available from
[PyPI](https://pypi.org/project/GDAL/) . However, the
[libgdal-dev]{.title-ref} dependencies are required, and the location of
the header files needs to be specified when installing. These libraries
are available via [UbuntuGIS](https://wiki.ubuntu.com/UbuntuGIS) and
other avenues.

On Debian systems, this process should work. Add the unstable release of
UbuntuGIS, get and install packages with:

``` console
sudo apt-get install software-properties-common
sudo apt-add-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get update
sudo apt-get install libgdal-dev gdal-bin
```

If the above doesn\'t work try removing all gdal related files from the
[/etc/apt/sources.list.d]{.title-ref} direcory and trying again (it\'s
possible old files are interfering with install attempts).

Once the [libgdal-dev]{.title-ref} libraries and [gdal-bin]{.title-ref}
are installed, pip install GDAL while specifying the version and the
location of the header files with:

``` console
pip install GDAL==$(gdal-config --version) --global-option=build_ext --global-option="-I/usr/include/gdal"
```

Then pyprecag should install without issue.

## Docker

The Dockerfile included in the examples directory builds a container
with Pyprecag installed in an Ubuntu 18.04 image. Build the container
with:

``` console
docker build -f docs/Dockerfile -t <tag> .
```

The tests can be run inside this container like this:

``` console
docker run <container id or tag> make test
```

## Windows

The dependencies [GDAL](https://www.gdal.org/),
[Fiona](https://github.com/Toblerity/Fiona) , and
[Rasterio](https://github.com/mapbox/rasterio) , all cause issues when
installing on Windows.

The easiest way to install all of these is with the Windows Binaries
provided by Christoph Gohlke at
<https://www.lfd.uci.edu/~gohlke/pythonlibs/>. Download the Python 2.7
version for your system architecture then install them with:

``` console
pip install /path/to/<downloaded_file.whl>
```

It is recommended to install GDAL first, as it is a requirements for
both Fiona and Rasterio. GDAL might also return an error and require you
to install Visual C++.

They are also available from the conda-forge channel in conda. pyprecag
is not currently available to install with conda, however can be pip
installed in a conda environment and use the conda-installed versions of
GDAL, Fiona and Rasterio.

Once those three dependencies are installed, pyprecag should install
without issue.

### Development

There is a makefile in the project root with targets for the most common
development operations such as lint checks, running unit tests, building
the documentation, and building installer packages. [tox]{.title-ref}
does not have a target, as [make tox]{.title-ref} is more typing than
[tox]{.title-ref}.

Run make with no target to see the list of targets:

``` bash
$ make
```

[Bump-my-version](https://callowayproject.github.io/bump-my-version/) is
used to manage the package version numbers. This ensures that the
version number is correctly incremented in all required files. Please
see the bumpversion documentation for usage instructions, and do not
edit the version strings directly.

Version numbers follow the [Semantic versioning guidelines](semver.org).

### Contributing

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

## Bug reports

When [reporting a
bug](https://github.com/CSIRO-Precision-Agriculture/pyprecag/issues)
please include:

> -   Your operating system name and version.
> -   Any details about your local setup that might be helpful in
>     troubleshooting.
> -   The pyprecag package version.
> -   Detailed steps to reproduce the bug.

## Documentation improvements

pyprecag could always use more documentation, whether as part of the
official pyprecag docs, in docstrings, or even on the web in blog posts,
articles, and such.

:::: note
::: title
Note
:::

This project uses Google-style docstrings. Contributed code should
follow the same conventions. For examples, please see the [Napoleon
examples](http://sphinxcontrib-napoleon.readthedocs.org/en/latest/example_google.html),
or the [Google Python Style
Guide](https://github.com/google/styleguide/blob/gh-pages/pyguide.md).
::::

## Feature requests and feedback

The best way to send feedback is to [file an
issue](https://github.com/CSIRO-Precision-Agriculture/pyprecag/issues)

If you are proposing a feature:

-   Explain in detail how it would work.
-   Keep the scope as narrow as possible, to make it easier to
    implement.
-   Remember that this is a volunteer-driven project, and that
    contributions are welcome :)

Or, implement the feature yourself and submit a pull request.

## Development

To set up pyprecag for local development:

1.  Fork the [pyprecag]{.title-ref} repo on GitHub.
2.  Clone your fork locally:

``` bash
$ git clone git@github.com:your_name_here/pyprecag.git
```

3.  Create a branch for local development:

``` bash
git checkout -b name-of-your-bugfix-or-feature
Now you can make your changes locally.
```

4.  When you\'re done making changes, run all the tests, doc builder and
    pylint checks using the project makefile:

``` bash
make clean lint test docs
```

5.  Commit your changes and push your branch to GitHub:

``` console
git add .
git commit -m "Your detailed description of your changes."
git push origin name-of-your-bugfix-or-feature
```

6.  Submit a pull request through the GitHub website.

## Pull Request Guidelines

If you need some code review or feedback while you\'re developing the
code just make the pull request.

For merging, you should:

1.  Include passing tests.
2.  Update documentation when there\'s new API, functionality etc.
3.  Add a note to `CHANGELOG.rst` about the changes.
4.  Add yourself to `AUTHORS.rst`.
