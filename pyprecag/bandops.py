import os
import collections
import numpy as np
from numpy import ma
import rasterio


class BandMapping(collections.MutableMapping, dict):
    """A dictionary used to manage band types and band numbers.

    If it has not been set it will have the default value of 0

    The list of keys is confined to those in __defaults, values must be integers.

    Attributes:
        __defaults = values to use as defaults
    """
    __defaults = {'red': 0, 'green': 0, 'blue': 0, 'infrared': 0, 'rededge': 0, 'mask': 0}

    def __init__(self, *args, **kwargs):
        self.update(**self.__defaults)
        self.update(*args, **kwargs)

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        """ check if key is from a predefined list and the value is an integer"""
        allowed_keys = self.__defaults.keys()

        if key not in allowed_keys:
            raise AttributeError(
                'BandMapping has no attribute {}. Allowed keys are {}'.format(key, ', '.join(self.__defaults.keys())))

        if not isinstance(value, int):
            raise ValueError('{v} is not an integer'.format(v=value))

        dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        dict.__setitem__(self, key, self.__defaults[key])

    def __iter__(self):
        return dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __contains__(self, x):
        return dict.__contains__(self, x)

    def allocated_bands(self):
        """Get a list of bands numbers already allocated to a band type ie not zero."""
        mapped_bands = dict(set(self.items()) - set(self.__defaults.items()))
        return sorted(mapped_bands.values())


class CalculateIndices(object):
    """Functions used for validating and calculating image indices.

    Args:
        **kwargs (): Components of the BandMapping object
    """

    def __init__(self, **kwargs):
        self.band_map = BandMapping()
        self.band_map.update(**kwargs)

        # store the equation functions in a dictionary, then can retrieve them by a string.
        self.equation = {'NDVI': self.ndvi,
                         'PCD': self.pcd,
                         'NDRE': self.ndre,
                         'GNDVI': self.gndvi,
                         'CHLRE': self.chlre}

        self.__image = None

    def valid_indices(self):
        """Generate a list of indices supported by the band mapping"""

        ind_list = []

        # the following will work when a single band integer or a numpy array is loaded against each band.
        if np.nansum(self.band_map['infrared']) > 0 and np.nansum(self.band_map['red']) > 0:
            ind_list.append('NDVI')
            ind_list.append('PCD')

        if np.nansum(self.band_map['infrared']) > 0 and np.nansum(self.band_map['green']) > 0:
            ind_list.append('GNDVI')

        if np.nansum(self.band_map['infrared']) > 0 and np.nansum(self.band_map['rededge']) > 0:
            ind_list.append('NDRE')
            ind_list.append('CHLRE')

        return ind_list

    def calculate(self, index_name, raster, src_nodata=None, dest_nodata=-9999):
        """ Using a given raster, calculate an image index. Valid indices include.
            NDVI - Normalised difference vegetation index
            PCD - Plant cell density index
            GNDVI - Green normalised difference vegetation index
            CHLRE - Chlorophyll red-edge index
            NDRE - Normalised difference red-edge index

        Args:
            index_name (str): The name of the index to calculate options include NDVI, PCD, GNDVI, NDRE and CHLRE
            raster (str): The input raster. This can be a filename, rasterio.io.Memoryfile or opened rasterio dataset
            src_nodata (int): The nodata value of the image. This will only be used if the input raster has a value of None
            dest_nodata (int): The value to use as the output no data value. If np.nan is required use None
        Returns:
            numpy.ndarray: The result of the index calculation
        """
        if raster is None:
            raise ValueError('Image file required')
        elif isinstance(raster, rasterio.io.MemoryFile):
            pass
        elif isinstance(raster, rasterio.DatasetReader):
            pass
        elif not os.path.exists(raster):
            raise ValueError('Image file does not exist')

        if src_nodata is not None and not isinstance(src_nodata, (int, float, long)):
            raise ValueError('src_nodata should be numeric (int, float or long)')

        if not dest_nodata and not isinstance(dest_nodata, (int, float, long)):
            raise ValueError('dst_nodata should be numeric (int, float or long)')

        self.__image = raster
        self.__src_nodata = src_nodata

        with np.errstate(divide='ignore', invalid='ignore'):
            result = self.equation[index_name.upper()]()

        if dest_nodata is not None:
            result[np.isnan(result)] = dest_nodata

        del self.__image, self.__src_nodata
        return result

    def _get_band_from_image(self, band_num):
        """extract a band from an image an apply masking"""
        if isinstance(self.__image, rasterio.DatasetReader):
            raster = self.__image
        elif isinstance(self.__image, rasterio.io.MemoryFile):
            raster = self.__image.open()
        else:
            raster = rasterio.open(os.path.normpath(self.__image))

        mask_band = 1
        if self.band_map['mask'] > 0:
            mask_band = self.band_map['mask']

        if raster.nodata is not None:
            mask = raster.read(mask_band, masked=True).mask
        elif self.__src_nodata is not None:
            mask = ma.masked_values(raster.read(mask_band, masked=True), self.__src_nodata).mask
        else:
            mask = np.zeros(shape=(raster.height, raster.width), dtype=bool)

        band = np.where(~mask, raster.read(band_num), np.nan).astype('float32')

        if not isinstance(self.__image, rasterio.DatasetReader):
            raster.close()

        return band

    def ndvi(self):
        """Calculate a normalised difference vegetation index. Requires Red and InfraRed Bands"""
        if self.band_map['infrared'] == 0:
            raise ValueError('InfraRed band not specified')

        if self.band_map['red'] == 0:
            raise ValueError('Red band not specified')

        red = self._get_band_from_image(self.band_map['red'])
        ir = self._get_band_from_image(self.band_map['infrared'])
        result = np.true_divide(ir - red, ir + red)
        return result.astype('float32')

    def gndvi(self):
        """Calculate a normalised difference vegetation index. Requires Green and InfraRed Bands"""
        if np.nansum(self.band_map['green']) == 0:
            raise ValueError('Green band not specified')

        if np.nansum(self.band_map['infrared']) == 0:
            raise ValueError('InfraRed band not specified')

        green = self._get_band_from_image(self.band_map['green'])
        ir = self._get_band_from_image(self.band_map['infrared'])
        result = np.true_divide(ir - green, ir + green)
        return result.astype('float32')

    def ndre(self):
        """Calculate a normalised difference red-edge index. Requires Red-edge and InfraRed Bands"""
        if np.nansum(self.band_map['infrared']) == 0:
            raise ValueError('InfraRed band not specified')

        if np.nansum(self.band_map['rededge']) == 0:
            raise ValueError('Red Edge band not specified')

        re = self._get_band_from_image(self.band_map['rededge'])
        ir = self._get_band_from_image(self.band_map['infrared'])
        result = np.true_divide(ir - re, ir + re)
        return result.astype('float32')

    def pcd(self):
        """Calculate a plant cell density index. Requires Red and InfraRed Bands"""

        if np.nansum(self.band_map['infrared']) == 0:
            raise ValueError('InfraRed band not specified')

        if np.nansum(self.band_map['red']) == 0:
            raise ValueError('Red band not specified')

        red = self._get_band_from_image(self.band_map['red'])
        ir = self._get_band_from_image(self.band_map['infrared'])
        result = np.true_divide(ir, red)
        return result.astype('float32')

    def chlre(self):
        """Calculate a Chlorophyll Red-edge index. Requires Red-edge and InfraRed Bands"""
        if np.nansum(self.band_map['infrared']) == 0:
            raise ValueError('InfraRed band not specified')

        if np.nansum(self.band_map['rededge']) == 0:
            raise ValueError('Red Edge band not specified')

        re = self._get_band_from_image(self.band_map['rededge'])
        ir = self._get_band_from_image(self.band_map['infrared'])
        result = np.true_divide(ir, re) - 1
        return result.astype('float32')
