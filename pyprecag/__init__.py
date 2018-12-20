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
__version__ = '0.0.4'

# Consider using tempfile.mkdtemp('PrecisionAg_') which will create a unique
# name and create the folder
TEMPDIR = os.path.join(tempfile.gettempdir(), 'PrecisionAg')

if not os.path.exists(TEMPDIR):
    os.mkdir(TEMPDIR)

LOGGER = logging.getLogger('pyprecag')
LOGGER.addHandler(logging.NullHandler())
DEBUG = config.get_config_key('debug_mode') # LOGGER.isEnabledFor(logging.DEBUG)

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

vesper_exe = config.get_config_key('vesperEXE')
if vesper_exe is None or vesper_exe == '' or not os.path.exists(vesper_exe):
    # Update if Vesper is installed.
    if os.path.exists(r'C:/Program Files (x86)/Vesper/Vesper1.6.exe'):
        config.set_config_key('vesperEXE', r'C:/Program Files (x86)/Vesper/Vesper1.6.exe')
    else:  # Otherwise report it.
        if not os.path.exists(vesper_exe):
            LOGGER.warning(
                'Vesper*.exe at "{}" does not exist. Please install and configure to allow for kriging to occur'.format(
                    vesper_exe))
        else:
            LOGGER.warning(
                'Vesper*.exe not found. Please install and configure to allow for kriging to occur')
#else:
    #LOGGER.info('Found Vesper*.exe at {}'.format(vesper_exe))

# Other applications ie QGIS can be used to set this value so force retrival from file by delete the variable.
del vesper_exe

# __gdal_version__ = get_gdal_release_name().decode('utf-8')    #from Fiona
# __gdal_version__ = gdal_version()    #From RasterIO

#
# # Set default logging handler to avoid "No handler found" warnings.
# try:  # Python 2.7+
#     from logging import NullHandler
# except ImportError:
#     class NullHandler(logging.Handler):
#         def emit(self, record):
#             pass
#
# def excepthook(*args):
#     """Catch any uncaught exception gets logged when assert/raise is used, not just AssertionError:
#     source: http: // code.activestate.com / recipes / 577074 - logging - asserts /
#       Args:
#           *args (): The stack trace, message to log
#       """
#     logging.error('', exc_info=args)
#
# sys.excepthook = excepthook
