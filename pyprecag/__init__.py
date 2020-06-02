# -*- coding: utf-8 -*-
""" pyprecag

    A suite of tools for Precision Agriculture
"""

import logging
import os
import sys
import subprocess
import tempfile
from . import config

__author__ = 'Christina Ratcliff',
__email__ = 'Christina.Ratcliff@csiro.au',
__version__ = '0.3.1'

TEMPDIR = os.path.join(tempfile.gettempdir(), 'PrecisionAg')

if not os.path.exists(TEMPDIR):
    os.mkdir(TEMPDIR)

LOGGER = logging.getLogger('pyprecag')
LOGGER.addHandler(logging.NullHandler())
# DEBUG = config.get_debug_mode() # LOGGER.isEnabledFor(logging.DEBUG)

# check to see if GDAL_DATA is defined
if not os.environ.get('GDAL_DATA', None):
    try:
        # If gdal-config is available find the datadir (UNIX-systems)
        process = subprocess.Popen(['gdal-config','--datadir'],
            stdout=subprocess.PIPE)
        output, error = process.communicate()
        gdal_data = output.split('\n')[0]
        assert os.path.isfile(os.path.join(gdal_data, 'gcs.csv'))
    except:
        # Otherwise, search for gdal folder relative to the python executable
        env_root = os.path.split(sys.executable)[0]
        gdal_data = os.path.join(env_root, 'Library', 'share', 'gdal')
        while not os.path.isfile(os.path.join(gdal_data, 'gcs.csv')):
            gdal_data = os.path.split(gdal_data)[0]
        assert os.path.isfile(
            os.path.join(gdal_data, 'gcs.csv')
        ), 'Could not find GDAL_DATA directory'
    LOGGER.warn(
        'Environment variable GDAL_DATA does not exist. Setting to {}'.format(
            gdal_data)
    )
    os.environ['GDAL_DATA'] = gdal_data

if not os.environ.get('PROJ_LIB', None):
    try:
        # Find relative to pyproj installed path
        import pyproj.datadir
        proj_lib = pyproj.datadir.pyproj_datadir
        assert os.path.isfile(os.path.join(proj_lib, 'epsg'))
    except:
        # Otherwise, find relative to gdal_data
        gdal_data = os.path.normpath(os.environ['GDAL_DATA'])
        proj_lib = os.path.join(gdal_data.split(os.sep + 'gdal')[0], 'proj')
        if not os.path.exists(proj_lib):
            proj_lib = os.path.split(proj_lib)[0]
        assert os.path.isfile(
            os.path.join(proj_lib, 'epsg')
        ), 'Could not find PROJ_LIB directory'
    LOGGER.warn(
        'Environment variable PROJ_LIB does not exist. Setting to {}'.format(
            proj_lib)
    )
    os.environ['PROJ_LIB'] = proj_lib
