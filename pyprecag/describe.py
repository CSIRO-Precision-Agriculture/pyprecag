import csv
import datetime
import difflib
import logging
import os
import re
import time
import warnings
from collections import OrderedDict, defaultdict
from operator import itemgetter

import chardet
import fiona
import geopandas
import numpy as np
import pandas as pd
from geopandas import GeoDataFrame
from osgeo import gdal

from . import crs as pyprecag_crs
from . import TEMPDIR, config

try:
    from pandas.errors import ParserWarning  # 0.20+
except:
    from pandas.io.common import ParserWarning

from unidecode import unidecode

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())  # Handle logging, no logging has been configured
# DEBUG = config.get_debug_mode()  # LOGGER.isEnabledFor(logging.DEBUG)
# LOGGER.setLevel("DEBUG")


class VectorDescribe:
    def __init__(self, input_data):
        """Get a description for a vector file

        Args:
            input_data (str):The input file
        """
        self.source = None
        self.file_encoding = None
        self.column_properties = None
        self.geometry_type = None
        self.is_mz_aware = None
        self.feature_count = None
        self.extent = None
        self.crs = pyprecag_crs.crs()

        self.source = input_data

        # Read the file and populate
        self.describeFile()

        return

    def open_geo_dataframe(self):
        """Create geopandas from file"""
        return GeoDataFrame.from_file(self.source, encoding=self.file_encoding)

    def describeFile(self):
        """
        Describe a vector File and set class properties
        """

        # Use this open so as to not hold open and use up memory
        gdf = GeoDataFrame.from_file(self.source, encoding=self.file_encoding)

        with fiona.open(self.source) as fio_coll:
            self.crs.getFromWKT(fio_coll.crs_wkt)

        # No idea why this works but it does so use it.
        exec ('rawstring = "{}"'.format(repr(','.join(gdf.columns.values))))
        result = chardet.detect(rawstring)
        self.file_encoding = result['encoding']

        self.feature_count = len(gdf)
        self.extent = list(gdf.total_bounds)

        # find the first element containing multi string, otherwise just use the first element from list.
        self.geometry_type = next((eaString for eaString in set(gdf.geom_type) if 'MULTI' in eaString.upper()),
                                  gdf.geom_type[0])

        self.is_mz_aware = gdf.geometry[0].has_z
        self.column_properties = get_column_properties(gdf)

        del gdf

    def get_column_names(self):
        return self.column_properties.keys()

    def get_alias_column_names(self):
        return [val['alias'] for key, val in self.column_properties.items()]

    def get_column_types(self):
        return [val['type'] for key, val in self.column_properties.items()]

    def get_shapefile_names(self):
        return [val['shapefile'] for key, val in self.column_properties.items()]


class CsvDescribe:
    def __init__(self, csv_filename):
        """A description of key elements relating to a comma or tab delimited text file.

        Args:
            csv_filename (str): a comma or tab delimited text

        """
        if not os.path.exists(csv_filename):
            raise IOError("Invalid path: {}".format(csv_filename))

        self.source = csv_filename
        self.dataframe = None
        self.file_encoding = 'ascii'
        self.dialect = None
        self.has_column_header = True
        self.row_count = None
        self.column_properties = None

        self.describe_file()
        self.column_count = len(self.column_properties)
        return

    def set_pandas_dataframe(self, pandas_dataframe):
        self.dataframe = pandas_dataframe

    def get_pandas_dataframe_fromfile(self):
        self.set_pandas_dataframe(self.open_pandas_dataframe())

    def open_pandas_dataframe(self, **kwargs):

        # suppress warnings for conflicting dialects using pandas.read_csv
        try:
            warnings.simplefilter(action='ignore', category=ParserWarning)
        except:
            pass

        if self.has_column_header:
            pdf = pd.read_csv(self.source, dialect=self.dialect, encoding=self.file_encoding, **kwargs)
        else:
            pdf = pd.read_csv(self.source, dialect=self.dialect, prefix='Column', header=None,
                              encoding=self.file_encoding, **kwargs)

        return pdf

    def get_column_names(self):
        return self.column_properties.keys()

    def get_alias_column_names(self):
        return [val['alias'] for key, val in self.column_properties.items()]

    def get_shapefile_column_names(self):
        return [val['shapefile'] for key, val in self.column_properties.items()]

    def get_column_types(self):
        return [val['type'] for key, val in self.column_properties.items()]

    def describe_file(self):
        """Describe a CSV File and set class properties
        """
        with open(self.source, 'r') as f:
            # sniff into 10KB of the file to check its dialect
            # this will sort out the delimiter and quote character.
            self.dialect = csv.Sniffer().sniff(f.read(10 * 1024))
            f.seek(0)  # reset read to start of file

            # read header based on the 10k of file.
            header = csv.Sniffer().has_header(f.read(10 * 1024))
            f.seek(0)  # reset read to start of file
            if not header:
                warnings.warn("The CSV file doesn't appear to contain column headers")
                self.has_column_header = False

            f.seek(0)  # reset read to start of file

        detector = chardet.UniversalDetector()
        with open(self.source, 'rb') as eaop:
            for line in eaop.readlines(100):
                detector.feed(line)
                if detector.done:
                    break
            detector.close()
        self.file_encoding = detector.result['encoding']

        pandas_df = self.open_pandas_dataframe()
        self.row_count = len(pandas_df)
        # store a dictionary of original and alias names along with column types. In most cases, objects types will
        # be strings this will enable lookups if necessary
        self.column_properties = get_column_properties(pandas_df)
        return


def get_esri_shapefile_schema(inputGeoDataFrame):
    """Construct an esri compatible schema for use with fiona.
        - remaps to fiona dtypes
        - Adheres to ESRI column naming standards  - 10 alpha numeric characters including '_' underscore.

    Args:
        inputGeoDataFrame (geopandas.GeoDataframe):

    Returns (dict): A Fiona compatible Dictionary

    """
    # construct the schema using geopandas
    schema = geopandas.io.file.infer_schema(inputGeoDataFrame)

    # Edit it to ESRI Shapefile Standards
    properties = OrderedDict([
        (re.sub('[^A-Za-z0-9_]+', '', name)[:10], fld_type) for name, fld_type in schema['properties'].iteritems()
    ])
    schema['properties'] = properties
    return schema


def save_geopandas_tofile(inputGeoDataFrame, output_filename, overwrite=True, file_encoding='ascii'):
    """Save a geodataframe to file.
        - adds functionality to asses and rename columns to ESRI compatible 10 alpha-numeric characters.
        - Maps lists and boolean column types to string.

    Args:
        inputGeoDataFrame (geopandas.geodataframe.GeoDataFrame): The Geodataframe to save
        output_filename (str): The output filename
        overwrite (bool):  Overwrite Existing file
        file_encoding (str): encoding type for output file.

    """
    if not isinstance(inputGeoDataFrame, GeoDataFrame):
        raise TypeError('Invalid Type : inputGeodataFrame')

    # if out_shapefilename doesn't include a path then add tempdir as well as overwriting it
    if output_filename is not None and not os.path.isabs(output_filename):
        output_filename = os.path.join(TEMPDIR, output_filename)
        overwrite = True

    if os.path.exists(output_filename) and not overwrite:
        raise IOError('Output file ({}) already exists, and overwrite is false'.format(output_filename))

    if os.path.splitext(output_filename)[-1] != '.shp':
        raise NotImplementedError('Currently only support shapefiles.... ')

    step_time = time.time()
    driver = 'ESRI Shapefile'
    if driver == 'ESRI Shapefile':

        inputGeoDataFrame = inputGeoDataFrame.copy()
        fldProp = get_column_properties(inputGeoDataFrame)

        # get a list of either bool or list columns and convert to string.
        fix_cols = [(key, val['type']) for key, val in fldProp.items() if val['type'] in ['bool', 'list']]
        fix_cols += [(key, val['dtype']) for key, val in fldProp.items() if 'datetime' in val['dtype'].lower()]

        # Convert them to Strings
        for col, col_type in fix_cols:
            LOGGER.info('Converting column {} datatype from {} to str'.format(col, col_type))
            if col_type == 'list':
                inputGeoDataFrame[col] = inputGeoDataFrame[col].apply(lambda x: ",".join(map(str, x)))
            else:
                inputGeoDataFrame[col] = inputGeoDataFrame[col].astype(str)

        # rename columns to alias names. columns must be listed in the same order
        inputGeoDataFrame.columns = [val['shapefile'] for key, val in fldProp.items()]

        '''Saving to file sometimes throws an error similar to
        CPLE_AppDefined in Value xxxx of field Timestamp of feature xxxx not successfully written. Possibly due to too
        larger number with respect to field width. This is a known GDAL Error. The following two lines will hide this
        from the user but may hide other message.
        https://gis.stackexchange.com/a/68042 is also an option that works
        '''

        gdal.UseExceptions()
        gdal.PushErrorHandler('CPLQuietErrorHandler')
        if file_encoding == 'ascii':
            inputGeoDataFrame.to_file(output_filename, driver=driver)
        else:
            inputGeoDataFrame.to_file(output_filename, driver=driver, encoding=file_encoding)

    if config.get_debug_mode():
        LOGGER.info('{:<30} {:<15} {dur}'.format('Saved to file',output_filename,
                                              dur=datetime.timedelta(seconds=time.time() - step_time)))


def get_dataframe_encoding(dataframe):
    exec ('rawstring = "{}"'.format(repr(','.join(dataframe.columns))))
    result = chardet.detect(rawstring)
    return result['encoding']


def get_column_properties(dataframe):
    """ Get a dictionary representing Column Properties for a pandas dataframe or a geopandas geodataframe.
       Includes:
            alias - removes spaces and replaces unicode chars with a sensible string using unidecode. ie oC to degC
            shapefile - An ESRI compatible 10 char alpha-numeric (excludes '-' & '_') column name.
            type - The fiona compatible column type
            dtype - The Pandas/Geopandas compatible column type.

        At Present it does not store the column width precision etc.
    Args:
        dataframe ([pandas.core.frame.DataFrame or geopandas.geodataframe.GeoDataFrame]):

    Returns:
        collections.OrderedDict: Representing properties of a column.

    TODO: Consider converting to a dictionary class with the get_shapefile_column_names etc
           see: https://stackoverflow.com/questions/1305532/convert-python-dict-to-object

    """
    column_desc = OrderedDict()

    for col, _type in zip(dataframe.columns, dataframe.dtypes):
        if col.lower() == 'geometry' or _type.name == 'geometry':
            fldtype = 'geometry'
        elif _type.name == 'object':
            fldtype = type(dataframe.iloc[0][col]).__name__
            if fldtype == 'unicode':
                fldtype = 'str'
        else:
            fldtype = type(np.zeros(1, _type).item()).__name__
            if fldtype == 'long':
                fldtype = 'int'

        if isinstance(col, unicode):
            aliasFld = unidecode(unicode(col))
        else:
            aliasFld = col

        # create a shapefile valid name 10 alpha numeric and underscore characters.
        # to keep underscore('_'), addit after the 9
        shpFld = re.sub('[^A-Za-z0-9_-]+', '', col)[:10]

        column_desc[col] = {'alias': aliasFld.replace(' ', ''),
                            'shapefile': shpFld,
                            'type': fldtype,
                            'dtype': str(dataframe[col].dtype)}
    return column_desc


def predictCoordinateColumnNames(column_names):
    """ Get the Longitude/easting and latitude/northing columns from a list of column_names

    Args:
        column_names (List[str]): A list of column names

    Returns:
        List[str]: [xColumn,yColumn]   Best matched column names

    """
    x_column = None
    y_column = None
    for eaVal in ['y', 'x']:

        valList = []
        for eaFld in config.get_config_key('geoCSV')['{}Coordinate_ColumnName'.format(eaVal)]:
            seqMatchDict = defaultdict(dict)

            # get a list of close matches by comparing known values to column_names
            close_matches = difflib.get_close_matches(eaFld.upper(), map(lambda x: x.upper(), column_names))
            if len(close_matches) > 0:
                # For each close match, calculate the match ratio and select the largest value
                for guess in column_names:
                    '''save the results to a dictionary, key is columnname and value is the match ratio. the ratio is
                    calculated on the occurrence of letters in the string in any order.
                     ie matching HEADING to EASTING has a ration of 0.714 because there are multiple similar matches'''
                    seqMatchDict[guess] = difflib.SequenceMatcher(None, eaFld.upper(), guess.upper(), True).ratio()

                # create short list of matches and ratios
                valList.append(max(seqMatchDict.iteritems(), key=lambda x: x[1]))

        # select the largest ratio as the best match
        if len(valList) > 0:
            best_match = max(valList, key=itemgetter(1))[0]
            exec ('{}_column = "{}"'.format(eaVal, best_match))

    LOGGER.debug('GeoCSV Columns:     x = {}, y = {}'.format(x_column, y_column))
    return [x_column, y_column]
