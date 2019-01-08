FROM ubuntu:18.04

RUN     apt-get update && \
        apt-get -y install software-properties-common python-pip && \
        apt-add-repository ppa:ubuntugis/ubuntugis-unstable && \
        apt-get update && \
        apt-get -y install libgdal-dev gdal-bin

RUN     pip install --user pip && \
        hash -r pip && \
        pip install --upgrade pip tox
RUN     pip install GDAL==$(gdal-config --version) --global-option=build_ext --global-option="-I/usr/include/gdal"

COPY    README.rst setup.cfg setup.py MANIFEST.in Makefile /
COPY    docs /docs
COPY    pyprecag /pyprecag

RUN     pip install -e.[dev,test]
