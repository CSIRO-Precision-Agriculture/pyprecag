import inspect
import logging

import os
import time
import warnings
from datetime import timedelta
from tempfile import NamedTemporaryFile

from osgeo import gdal
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio import shutil as rio_shutil
from rasterio.warp import calculate_default_transform, reproject
from rasterio.windows import get_data_window, from_bounds, intersection

from scipy.ndimage import generic_filter
from . import crs as pyprecag_crs
from .bandops import CalculateIndices
from . import config, TEMPDIR

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())


def create_raster_transform(bounds, pixel_size, snap_extent_to_pixel=True, buffer_by_pixels=0):
    """Create parameters required for creating a new raster file based on a known extent and pixel
     size.

    snap_extent_to_pixel can be used to ensure the bounding coordinates are a divisible of the
     pixel size.

    Args:
        bounds (float, float, float, float): the bounding box coordinates (xmin, ymin, xmax, ymax)
        pixel_size (int, float): the required pixel size
        snap_extent_to_pixel (bool): round the extent coordinates to be divisible by the pixel size
        buffer_by_pixels (int): The number of pixels to buffer the input bounding box by.

    Returns:
        affine.Affine: the Rasterio Transformation object
        int: the width or number of columns
        int: the height or number of rows.
        [float, float, float, float]: new bounding box coordinates updated if snapping was used.
    """

    if not isinstance(bounds, (tuple, list)):
        raise TypeError('Bounds should be a tuple or list for xmin, ymin, xmax, ymax')

    if not isinstance(pixel_size, (int, long, float)):
        raise TypeError('pixel_size must be numeric number.')
    if not isinstance(snap_extent_to_pixel, (int, bool)):
        raise TypeError('snap_extent_to_pixel must be boolean.')

    # adjust values to an area slightly larger than the bounds
    if buffer_by_pixels > 0:
        bounds = (bounds[0] - (pixel_size * buffer_by_pixels),
                  bounds[1] - (pixel_size * buffer_by_pixels),
                  bounds[2] + (pixel_size * buffer_by_pixels),
                  bounds[3] + (pixel_size * buffer_by_pixels))

    # We may want to snap the output grids to a multiple of the grid size, allowing
    # adjacent blocks to align nicely.
    if snap_extent_to_pixel:
        x_min, y_min, x_max, y_max = raster_snap_extent(*bounds, pixel_size=pixel_size)
    else:
        x_min, y_min, x_max, y_max = bounds

    width = int((x_max - x_min) / pixel_size)  # columns
    height = int((y_max - y_min) / pixel_size)  # rows

    # create an affine transformation matrix to associate the array to the coordinates.
    from rasterio.transform import from_origin
    transform = from_origin(x_min, y_max, pixel_size, pixel_size)

    return transform, width, height, (x_min, y_min, x_max, y_max)


def raster_snap_extent(x_min, y_min, x_max, y_max, pixel_size):
    """Calculate a new raster extent where bounding coordinates are a divisible of the pixel size.

   The LL will be rounded down, and the UR will be rounded up.

    Using this function will ensure that rasters have the same pixel origin, when they are created
    which in the long run, will save resampling when multiple raster with the same pixel size
    are compared.

    Args:
        x_min (float): x_min coordinate will be rounded down to the nearest pixel edge
        y_min (float): y_min coordinate will be rounded down to the nearest pixel edge
        x_max (float): x_max coordinate will be rounded up to the nearest pixel edge
        y_max (float): y_max coordinate will be rounded up to the nearest pixel edge
        pixel_size (float): The pixel size representing the divisor

    Returns:
        List[float]]: Representing the Bounding box (xmin, ymin, xmax, ymax)

    """

    for argCheck in [('x_min', x_min), ('y_min', y_min), ('x_max', x_max), ('y_max', y_max),
                     ('pixel_size', pixel_size)]:
        if not isinstance(argCheck[1], (int, long, float)):
            raise TypeError('{} must be a floating number.'.format(argCheck[0]))

    b_box = [x_min - (x_min % pixel_size),  # calc new xMin
             y_min - (y_min % pixel_size),  # calc new yMin
             (x_max + pixel_size) - ((x_max + pixel_size) % pixel_size),  # calc new xMax
             (y_max + pixel_size) - ((y_max + pixel_size) % pixel_size)]  # calc new yMax

    return b_box


def rescale(raster, min_value, max_value, band_num=1, ignore_nodata=True):
    """ Rescale a single band between a set number of values.

    If ignore_nodata is used, then the selected band will be opened as a numpy masked array
    and the specified nodata values will be excluded

    It returns the calculated single band and can be written to file using
    rasterio.open(os.path.normpath(),'w')

    Args:
        raster (rasterio.io.DatasetReader): Raster file opened by rasterio.open(os.path.normpath())
        min_value (int): The lower/min value to use during rescaling
        max_value (int): The Upper/max value to use during rescaling
        band_num (int): The band number to apply rescaling too.
        ignore_nodata (bool): Ignore nodata values during rescaling.
                    If False, the nodata pixels and values will be used during the calculation
                    If True, band is read as a masked array and nodata pixels and values excluded

    Returns:
        numpy.ndarray: A single band as a numpy array.
        or
        numpy.ma.core.MaskedArray:  Single band numpy array with nodata being stored in the mask
    """

    if not isinstance(raster, rasterio.DatasetReader):
        raise TypeError("Input should be a rasterio.DatasetReader created using "
                        "rasterio.open(os.path.normpath())")

    for argCheck in [('min_value', min_value), ('max_value', max_value)]:
        if not isinstance(argCheck[1], (int, long, float)):
            raise TypeError('{} must be numeric.'.format(argCheck[0]))

    band = raster.read(band_num, masked=ignore_nodata)

    # formula: https://gis.stackexchange.com/a/28563
    # use np.nanXXX to create consistent results.
    band = band.astype(np.float64)

    rescaled = ((band - np.nanmin(band)) * (max_value - min_value)
                / (np.nanmax(band) - np.nanmin(band)) + min_value)

    # pick and assign the most appropriate dtype for the result
    rescaled = rescaled.astype(np.dtype(rasterio.dtypes.get_minimum_dtype(rescaled)))

    return rescaled


def normalise(raster, band_num=1, ignore_nodata=True):
    """Normalise a single band by adjusting to a mean of zero and standard deviation of 1

    If ignore_nodata is used, then the selected band will be opened as a numpy masked array and
    the specified nodata values will be excluded from calculations

    Returns calculated single band that can be written to file using
     rasterio.open(os.path.normpath(),'w')

    Args:
        raster (rasterio.io.DatasetReader): Raster file opened by rasterio.open(os.path.normpath())
        band_num (int):       The band number to apply rescaling too.
        ignore_nodata (bool): Ignore nodata values during rescaling.
                              If False, nodata pixels and values will be used during the calculation
                              If True, band will be read as a masked array and nodata pixels
                                  and values excluded.
    Returns:
        numpy.ndarray: A single band as a numpy array.
        or
        numpy.ma.core.MaskedArray:    A single band as a numpy array with nodata being stored in
                                   the mask
    """

    if not isinstance(raster, rasterio.DatasetReader):
        raise TypeError("Input should be a rasterio.DatasetReader created using "
                        "rasterio.open(os.path.normpath())")

    band = raster.read(band_num, masked=ignore_nodata)

    # convert to float64 for accuracy
    band = band.astype(np.float64)

    normalised = (band - np.nanmean(band)) / np.nanstd(band)

    normalised = normalised.astype(np.dtype(rasterio.dtypes.get_minimum_dtype(normalised)))

    return normalised


def nancv(x):
    """ A function used with scipy.ndimage.generic_filter to calculate the coefficient
    variant of pixels/values excluding nan (nodata) values. It can be used in conjunction with
     focal_statistics.

    example using a 3x3: generic_filter(band,nancv, mode='constant', cval=np.nan, size=3)

    A variation of https://stackoverflow.com/a/14060024
    """
    return np.true_divide(np.nanstd(x), np.nanmean(x))


def pixelcount(x):
    """ A function used with scipy.ndimage.generic_filter to count the number of real values/pixels
        (ie not nan) when applying a NxN filter.

        A Count of 0 will be replace by np.nan. It can be used in conjunction with focal statistics.
        A variation of https://stackoverflow.com/a/14060024"""
    val = sum(~np.isnan(x))
    if val == 0:
        return np.nan
    else:
        return val


def focal_statistics(raster, band_num=1, ignore_nodata=True, size=3, function=np.nanmean,
                     clip_to_mask=False, out_colname=None):
    """Derives for each pixel a statistic of the values with the specified neighbourhood.

    Currently the neighbourhood is restricted to square neighbourhoods ie 3x3 size.

    Any numpy statistical functions are supported along with custom functions.

    Nodata values are converted to np.nan and will be excluded from the statistical calculation.
    Nodata pixels may be assigned a value if at least one pixel in the neighbourhood has a valid
    value. To remove/mask these values from the final output, set the clip_to_mask setting to True.

    Using a size of 1 returns the selected band with the converted nodata values. No statistical
    functions are applied.

    An string out_colname is returned and can be used as a filename or column name during future
    analysis. If None, it is derived from the input raster, size and statistical function used.
          For single band inputs   <stat><size>x<size>_<raster name>
             eg.   mean3x3_area1_yield    apply a mean 3x3 filter for raster area1_yield

          For multi band inputs   <function><size>x<size>bd<band_num>_<raster name>
             eg.   mean3x3b3_area2       apply a mean 3x3 filter for band 3 of the raster area2

    Source: https://stackoverflow.com/a/30853116/9567306
    https://stackoverflow.com/questions/46953448/local-mean-filter-of-a-numpy-array-with-missing-data/47052791#47052791

    Args:
        raster (rasterio.io.DatasetReader): Raster file opened using
                                 rasterio.open(os.path.normpath())
        band_num (int):       The band number to apply focal statistics to.
        ignore_nodata (bool): If true, the nodata value of the raster will be converted
                              to np.nan and excluded from statistical calculations.
        size (int):           Size of the neighbourhood filter used for statistics calculations.
                              Currently restricted to a square neighbourhood ie 3x3, 5x5 etc.
        function (function):  a functions to apply to the raster. These can include numpy functions
                              like np.nanmean or custom ones.
        clip_to_mask (bool):  If true, remove values assigned to nodata pixels
        out_colname (str):    An output string used to describe the filter result and can be used
                              as a column or filename If NONE, then it will be derived.
    Returns:
        numpy.ndarray:        A 1D numpy array of double (float32) values
        str:                  a string representation of the inputs

    """

    if not isinstance(raster, rasterio.DatasetReader):
        raise TypeError("Input should be a rasterio.DatasetReader created using rasterio.open()")

    if not isinstance(size, int) or size % 2 == 0:
        raise TypeError("Size should be an odd number integer greater than one. Only "
                        "Square Filters are supported.")

    if not isinstance(ignore_nodata, bool):
        raise TypeError('{} should be a boolean.'.format(ignore_nodata))

    start_time = time.time()
    col_name = []
    mask = raster.read_masks(band_num)
    if ignore_nodata:
        # change nodata values to np.nan
        band = np.where(mask, raster.read(band_num), np.nan)
    else:
        band = raster.read(band_num)

    # convert to float64 for accuracy
    band = band.astype(np.float64)
    if size > 1:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            filtered = generic_filter(band, function, mode='constant', cval=np.nan, size=size)

        col_name += [function.func_name.replace('nan', ''), '{0}x{0}'.format(size)]
    else:
        filtered = band
        col_name += ['pixel']

    if raster.count > 1:  # the number of bands in the raster
        col_name += ['bd{}'.format(band_num)]

    if clip_to_mask:
        # reapply the mask to remove values assigned to nodata pixels
        filtered = np.where(mask, filtered, np.nan)

    title = os.path.splitext(os.path.basename(raster.name))[0]
    if out_colname is None or out_colname == '':
        out_colname = '{}_{}'.format(''.join(col_name), title)

    if config.get_debug_mode():
        LOGGER.info('{:50}  {dur:17} min: {:>.4f} max: {:>.4f}'.format(
            out_colname, np.nanmin(filtered), np.nanmax(filtered),
            dur=timedelta(seconds=time.time() - start_time)))

    return filtered.astype(np.float32), out_colname


def calculate_image_indices(image_file, band_map, out_image_file, indices=[], out_nodata=-9999):
    """Creates a multi-band image where each band represents a calculated index.

    Rasterio's band tags are used to document which band represents which index.

    The band mapping matches a band number to a band type ie band 3 is the Red band to enable the
    index calculation to occur. It also identifies a band where the nodata value removes the
    non-vine or bare earth signal. This nodata is assigned to the output image.

     Indices currently supported are:
        NDVI - Normalised difference vegetation index
        PCD - Plant cell density index
        GNDVI - Green normalised difference vegetation index
        CHLRE - Chlorophyll red-edge index
        NDRE - Normalised difference red-edge index

    Args:
        image_file (str): image_file (str): the input image file
        band_map (pyprecag.bandops.BandMapping): a Band mapping matching a band number to band type.
        out_image_file (str): the name and location of the output image file.
        indices (List[str]): indices (List[str]): The list of indices to calculate.
        out_nodata (int): the value to use in the output image for the nodata

    Returns:
        None:
    """

    if not isinstance(indices, list):
        raise TypeError('indices must be a list of indices to calculate')

    start_time = time.time()
    meta = rasterio.open(image_file).meta.copy()
    meta.update({'driver': 'GTiff',
                 'dtype': rasterio.float32,
                 'count': len(indices),
                 'nodata': out_nodata})

    ci = CalculateIndices(**band_map)

    with rasterio.open(out_image_file, 'w', compress="NONE", **meta) as dest:
        # Calculate each index and write them as a separate band.
        for i, eaIndex in enumerate(indices, 1):
            index_arr = ci.calculate(eaIndex, image_file)
            dest.write(index_arr.astype(rasterio.float32), i)
            dest.update_tags(i, name=eaIndex)
            del index_arr

    LOGGER.info(
        '{:<30} {:>10}   {:<15} {dur}'.format('Indices Calculate for Image', '', ', '.join(indices),
                                              dur=timedelta(seconds=time.time() - start_time)))


def reproject_image(image_file, out_imagefile, out_epsg, band_nums=[],
                    image_epsg=0, image_nodata=None, resampling=Resampling.nearest):
    """Reproject selected image bands from one coordinate system to another.

    "image_epsg" and "image_nodata" can be used to set the coordinate system and image nodata
    values when they are not already specified within the image file.

    Args:
        image_file (str): An input image path and name
        out_imagefile (str): An output image path and name
        out_epsg (int): The epsg number representing output coordinate system
        band_nums (List): The bands to reproject. If list is empty, all bands will be used
        image_epsg (int): the epsg number for the image. Only used if not provided by the image.
        image_nodata (int): the nodata value for the image. Only used if not provided by the image.
        resampling (rasterio.enums.Resampling): The resampling technique to use during reprojection.

    Returns:
        None:
    """
    start_time = time.time()

    if out_epsg > 0:
        dst_crs = rasterio.crs.CRS.from_epsg(out_epsg)

    with rasterio.open(os.path.normpath(image_file)) as src:
        meta = src.meta.copy()
        meta.update({'driver': 'GTiff'})

        if len(band_nums) == 0:
            band_nums = range(1, src.count + 1)
        meta.update({'count': len(band_nums)})

        if src.crs is None:
            if image_epsg is None or image_epsg == 0:
                raise ValueError('Input coordinate system required - image_file does not '
                                 'contain a coordinate system, and in_epsg is 0')
            else:
                src_crs = rasterio.crs.CRS.from_epsg(image_epsg)
                meta.update({'crs': src_crs})
        else:
            image_epsg = pyprecag_crs.getCRSfromRasterFile(image_file).epsg_number
            src_crs = src.crs

        if image_nodata is not None or image_nodata != meta['nodata']:
            meta.update({'nodata': image_nodata})

        # create a raster in memory by copying input data or reprojecting ------------------------
        if out_epsg == image_epsg:
            with rasterio.open(out_imagefile, 'w', **meta) as dest:
                for i, ea_band in enumerate(band_nums, 1):
                    dest.write(src.read(ea_band), i)
                    dest.update_tags(i, **src.tags(ea_band))

                dest.update_tags(**src.tags())

            LOGGER.info('{:<30} {:>10}   {:<15} {dur}'.format(
                'Processed Image', '', 'CRS: {} To  {}, nodata: {} To {}'.format(
                    image_epsg, out_epsg, src.nodata, image_nodata),
                dur=timedelta(seconds=time.time() - start_time)))

        else:
            # calculate the new affine transform for to use for reprojection.
            transform, width, height = calculate_default_transform(src_crs, dst_crs, src.width,
                                                                   src.height, *src.bounds)

            # update the parameters for the output file
            meta.update({'crs': dst_crs, 'transform': transform, 'width': width, 'height': height})

            with rasterio.open(out_imagefile, 'w', **meta) as dest:
                for i in band_nums:
                    reproject(source=rasterio.band(src, i),
                              destination=rasterio.band(dest, i),
                              src_transform=src.transform,
                              src_crs=src.crs,
                              dst_transform=transform,
                              dst_crs=dst_crs,
                              resampling=resampling)

                    # image statistics have changed so copy only the other tags
                    cleaned_tags = dict([(key, val) for key, val in src.tags(i).iteritems()
                                         if not key.upper().startswith('STATISTIC')])
                    if len(cleaned_tags) > 0:
                        dest.update_tags(i, **cleaned_tags)

            LOGGER.info('{:<30} {:>10}   {:<15} {dur}'.format(
                'Reproject Image', '', 'From {} To  {}'.format(image_epsg, out_epsg),
                dur=timedelta(seconds=time.time() - start_time)))

            del transform


def stack_and_clip_rasters(raster_files, use_common=True, output_tif=None):
    """Combine multiple single band files into a single multi band raster where one band represents
    one file. The filename will be saved in the bands tag. Optionally a minimum common area mask
    can be applied.

    All input rasters MUST be of the same coordinate system and pixel size.

    Args:
        raster_files (List[str]): list of raster files to process
        use_common (bool):  if true input rasters will be masked to the minimum common data area
        output_tif (str):  output tif name. if omitted will be created in TEMPDIR

    Returns:
        str(): The output tif file
        list(): A list of the band tags. (ie order of allocation to bands)
    """
    # if output_tif doesn't include a path then add tempdir as well as overwriting it
    if output_tif is not None and not os.path.isabs(output_tif):
        output_tif = os.path.join(TEMPDIR, output_tif)

    if output_tif is None or config.get_debug_mode():
        # get a unique name, it will create, open the file and delete after saving the variable
        with NamedTemporaryFile(prefix='{}_'.format(
                inspect.getframeinfo(inspect.currentframe())[2]),
                suffix='.tif', dir=TEMPDIR) as new_file:
            output_tif = new_file.name
        del new_file

    if not isinstance(raster_files, list):
        raise TypeError('Invalid Type: raster_files should be a list')

    if len(raster_files) == 0:
        raise TypeError('Invalid Type: Empty list of raster files')

    not_exists = [my_file for my_file in raster_files if not os.path.exists(my_file)]
    if len(not_exists) > 0:
        raise IOError('raster_files: {} raster file(s) do '
                      'not exist\n\t({})'.format(len(not_exists), '\n\t'.join(not_exists)))

    check_pixelsize = []
    check_crs = []

    for ea_raster in raster_files:
        with rasterio.open(ea_raster) as src:
            if src.crs is None:
                check_crs.append(ea_raster)
            if src.res not in check_pixelsize:
                check_pixelsize.append(src.res)

    if len(check_pixelsize) == 1:
        resolution = check_pixelsize[0]
    else:
        raise TypeError("raster_files are of different pixel sizes - {}".format(list(set(check_pixelsize))))

    if len(check_crs) > 0:
        raise TypeError("{} raster(s) don't have coordinates "
                        "systems assigned \n\t{}".format(len(check_crs),'\n\t'.join(check_crs)))

    start_time = time.time()
    step_time = time.time()

    # Make sure ALL masks are saved inside the TIFF file not as a sidecar.
    gdal.SetConfigOption('GDAL_TIFF_INTERNAL_MASK', 'YES')

    try:
        # Find minimum overlapping extent of all rasters ------------------------------------
        for i, ea_raster in enumerate(raster_files, start=1):
            with rasterio.open(ea_raster) as src:

                band1 = src.read(1, masked=True)

                # get minimum data extent
                data_window = get_data_window(band1)

                # find the intersection of all the windows
                if i == 1:
                    min_window = data_window
                else:
                    # create a new window using the last coordinates based on this image in case
                    # extents/pixel origins are different
                    min_img_window = from_bounds(*min_bbox, transform=src.transform)

                    # find the intersection of the windows.
                    min_window = intersection(min_img_window, data_window).round_lengths('ceil')

                # convert the window co coordinates
                min_bbox = src.window_bounds(min_window)

                del min_window
        if config.get_debug_mode():
            LOGGER.info('{:<30} {:<15} {dur} {}'.format(
                'Found Common Extent', '', min_bbox,
                dur=timedelta(seconds=time.time() - step_time)))

    except rasterio.errors.WindowError as e:
        # reword 'windows do not intersect' error message
        if not e.args:
            e.args = ('',)
        e.args = ("Rasters (Images) do not overlap",)

        raise  # re-raise current exception

    # Create the metadata for the overlapping extent image. ---------------------------------------
    transform, width, height, bbox = create_raster_transform(min_bbox, resolution[0],
                                                             buffer_by_pixels=5)

    # update the new image metadata ---------------------------------------------------------------
    kwargs = rasterio.open(raster_files[0]).meta.copy()
    kwargs.update({'count': len(raster_files), 'nodata': -9999,
                   'height': height, 'width': width,
                   'transform': transform})

    ''' combine individual rasters into one multi-band tiff using reproject will 
           fit each raster to the same pixel origin
           will crop to the min bbox
           standardise the nodata values
           apply nodata values to entire image including internal blocks'''

    with rasterio.open(output_tif, 'w', **kwargs) as dst:
        band_list = []
        for i, ea_raster in enumerate(raster_files, start=1):
            image = os.path.basename(os.path.splitext(ea_raster)[0])
            band_list.append(image)

            with rasterio.open(ea_raster) as src:
                reproject(source=rasterio.band(src, 1),
                          destination=rasterio.band(dst, i),
                          src_transform=src.transform,
                          dst_transform=transform,
                          resampling=Resampling.nearest)

            dst.update_tags(i, **{'name': image})

        dst.descriptions = band_list
    if config.get_debug_mode():
        LOGGER.info('{:<30} {:<15} {dur}'.format('Rasters Combined', '',
                                                 dur=timedelta(seconds=time.time() - step_time)))

    # find the common data area -------------------------------------------------------
    if use_common:
        with rasterio.open(output_tif, 'r+') as src:
            # find common area across all bands as there maybe internal nodata values in some bands.
            mask = []

            # loop through all all the masks
            for ea_mask in src.read_masks():
                if len(mask) == 0:
                    mask = ea_mask
                else:
                    mask = mask & ea_mask

            # change mask values to  0 = nodata, 1 = valid data
            mask[mask == 255] = 1

            # apply mask to all bands of file
            src.write_mask(mask)

        if config.get_debug_mode():
            # write mask to new file
            tmp_file = os.path.join(TEMPDIR, output_tif.replace('.tif', '_commonarea.tif'))

            kwargs_tmp = kwargs.copy()
            kwargs_tmp.update({'count': 1, 'nodata': 0,
                               'dtype': rasterio.dtypes.get_minimum_dtype(mask)})

            with rasterio.open(tmp_file, 'w', **kwargs_tmp) as tmp_dst:
                tmp_dst.write(mask, 1)

        LOGGER.info('{:<30} {:<15} {dur}'.format('Rasters Combined and clipped', '',
                                                 dur=timedelta(seconds=time.time() - start_time)))

    gdal.SetConfigOption('GDAL_TIFF_INTERNAL_MASK', None)

    return output_tif, band_list


def update_band_statistics_GDAL(raster_file):
    """ update band statistics using GDAL. Not yet built into rasterio """
    ds = gdal.Open(raster_file, gdal.GA_Update)
    for i in range(ds.RasterCount):
        ds.GetRasterBand(i + 1).ComputeStatistics(0)
    ds = band = None  # save, close
    return


def save_in_memory_raster_to_file(memory_raster, out_image):
    """Save a rasterio memory file as a TIF to disk

    if out_image does not contain a path, it will be save to the Temp Dir.

    Args:
        memory_raster (rasterio.io.MemoryFile):
        out_image (str): the tif image name.

    Returns:
        str: the image path and name.
    """

    if not isinstance(memory_raster, rasterio.io.MemoryFile):
        raise TypeError('Input raster is not a raster.io.MemoryFile')

    if os.path.splitext(out_image)[-1].lower() != '.tif':
        raise ValueError('File Extension of {} is not supported. Please change to '
                         '.tif (GeoTiff)'.format(os.path.splitext(out_image)[-1]))

    if out_image is not None and not os.path.isabs(out_image):
        out_image = os.path.join(TEMPDIR, out_image)

    start_time = time.time()

    with memory_raster.open() as src:
        with rasterio.open(out_image, 'w', tfw='YES', **src.meta.copy()) as dest:
            for i in src.indexes:
                dest.write(src.read(i), i)

                cleaned_tags = dict([(key, val) for key, val in src.tags(i).iteritems()
                                     if not key.upper().startswith('STATISTIC')])

                if len(cleaned_tags) > 0:
                    dest.update_tags(i, **cleaned_tags)

    if config.get_debug_mode():
        LOGGER.info('{:<30} {:<15} {dur}'.format('Saved to file', out_image,
                                                 dur=timedelta(seconds=time.time() - start_time)))

    return out_image
