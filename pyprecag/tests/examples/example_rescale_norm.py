import logging
import os
import tempfile
import time

import numpy as np
import rasterio
from osgeo import gdal

from pyprecag.tests.make_dummy_data import make_dummy_tif_files
from pyprecag.raster_ops import normalise, rescale

os.environ['GDAL_DATA']=r'C:/Miniconda3/envs/pvtools27/Library/share/gdal'
os.environ['PROJ_LIB']=r'C:/Miniconda3/envs/pvtools27/Library/share'

DEBUG = True

pyFile = os.path.basename(__file__)

out_dir = tempfile.gettempdir()
out_dir = r'C:\data\temp'
out_dir = os.path.join(out_dir, os.path.splitext(pyFile)[0])

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

log_file = '{0}.log'.format(pyFile.replace('.py', time.strftime('_%Y-%m-%d', time.localtime())))
# logging.captureWarnings(False)
# logging.basicConfig(level=logging.DEBUG, format="%(levelname)s - %(name)s : %(message)s", filename=os.path.join(out_dir, log_file), filemode='w')
# https://stackoverflow.com/questions/9321741/printing-to-screen-and-writing-to-a-file-at-the-same-time
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
# set a format which is simpler for console use
formatter = logging.Formatter("%(message)s")
# tell the handler to use this format
console.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger('').addHandler(console)

fileList = make_dummy_tif_files(out_dir)

start_time = time.time()
for ea in fileList:
    file_time = time.time()
    raster = ea

    outNorm = os.path.join(out_dir, os.path.basename(raster).replace('.tif', '_pyprecag-normalise.tif'))
    outMulti = os.path.join(out_dir, os.path.basename(raster).replace('.tif', '_pyprecag-Multiband.tif'))

    outRescale1 = os.path.join(out_dir, os.path.basename(raster).replace('.tif', '_pyprecag-rescale0-1.tif'))
    outRescale255 = os.path.join(out_dir, os.path.basename(raster).replace('.tif', '_pyprecag-rescale0-255.tif'))

    print('{lne}\nRescale-Normalise \t {}'.format(ea,lne='=' * 100))

    gdalRaster = gdal.Open(raster)
    wktproj = gdalRaster.GetProjectionRef()
    del gdalRaster

    with rasterio.open(os.path.normpath(raster)) as src:
        meta = src.meta.copy()
        meta['crs'] = wktproj
        meta['count'] = 1
        band = src.read(1, masked=True)
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('dType:', src.dtypes, 'BandCount:', src.count, 'NoData Val:',
                                                   src.nodata, ''))

        band = band.astype(np.float64)
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Pixels:-', 'Total:', band.size, 'Valid:', band.count(), 'Nodata:',
                                                     band.size - band.count()))
        print('{:<30}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}'.format('type', 'min', 'max', 'mean', 'sum',
                                                                               'std-ddof0', 'std-ddof1'))
        print('{:<30}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:<20}\t{:8} {}'.format('Orig',
                                                                                        np.nanmin(band),
                                                                                        np.nanmax(band),
                                                                                        np.nanmean(band),
                                                                                        np.nansum(band),
                                                                                        np.nanstd(band, ddof=0),
                                                                                        np.nanstd(band, ddof=1),
                                                                                        'dtype', band.dtype))
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++
        rescaled1 = rescale(src, 0, 1, ignore_nodata=True)
        meta['dtype'] = rescaled1.dtype
        with rasterio.open(os.path.normpath(outRescale1), 'w', **meta) as out:
            out.write_band(1, rescaled1)

        # update band statistics using GDAL. Not yet built into rasterio
        ds = gdal.Open(outRescale1, gdal.GA_Update)
        for i in range(ds.RasterCount):
            ds.GetRasterBand(i + 1).ComputeStatistics(0)
        ds = band = None  # save, close

        with rasterio.open(os.path.normpath(outRescale1)) as out:
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('dType:', out.dtypes, 'BandCount:', out.count, 'NoData Val:',
                                                      out.nodata, outRescale1))

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++
        rescaled255 = rescale(src, 0, 255)
        meta['dtype'] = rescaled255.dtype
        with rasterio.open(os.path.normpath(outRescale255), 'w', **meta) as out:
            out.write_band(1, rescaled255)

        # update band statistics using GDAL. Not yet built into rasterio
        ds = gdal.Open(outRescale255, gdal.GA_Update)
        for i in range(ds.RasterCount):
            ds.GetRasterBand(i + 1).ComputeStatistics(0)
        ds = band = None  # save, close

        with rasterio.open(os.path.normpath(outRescale255)) as out:
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('dType:', out.dtypes, 'BandCount:', out.count, 'NoData Val:',
                                                      out.nodata, outRescale255))

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++
        normalised = normalise(src, ignore_nodata=True)
        meta['dtype'] = normalised.dtype
        with rasterio.open(os.path.normpath(outNorm), 'w', **meta) as out:
            out.write_band(1, normalised)

        # update band statistics using GDAL. Not yet built into rasterio
        ds = gdal.Open(outNorm, gdal.GA_Update)
        for i in range(ds.RasterCount):
            ds.GetRasterBand(i + 1).ComputeStatistics(0)
        ds = band = None  # save, close

        with rasterio.open(os.path.normpath(outNorm)) as out:
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('dType:', out.dtypes, 'BandCount:', out.count, 'NoData Val:',
                                                      out.nodata, outNorm))

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++
        meta['count'] = 3
        meta['dtype'] = normalised.dtype

        with rasterio.open(os.path.normpath(outMulti), 'w', **meta) as out:
            out.write_band(1, normalised)
            out.write_band(2, rescaled1)
            out.write_band(3, rescaled255)

        # update band statistics using GDAL. Not yet built into rasterio
        ds = gdal.Open(outMulti, gdal.GA_Update)
        for i in range(ds.RasterCount):
            ds.GetRasterBand(i + 1).ComputeStatistics(0)
        ds = band = None  # save, close

        with rasterio.open(os.path.normpath(outMulti)) as out:
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('dType:', out.dtypes, 'BandCount:', out.count, 'NoData Val:',
                                                      out.nodata, outMulti))
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++
print "\n\n***********\n{} is Complete\n".format(pyFile)
