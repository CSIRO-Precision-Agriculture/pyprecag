import logging
import os

from osgeo import gdal
import rasterio
import tempfile
import time
import numpy as np

from pyprecag.tests import make_dummy_data
import pyprecag.crs
from pyprecag.describe import VectorDescribe
from pyprecag.processing import extractPixelStatisticsForPoints, randomPixelSelection

DEBUG = True
pyFile = os.path.basename(__file__)
# check to see if GDAL_DATA is defined
if not os.environ.get('GDAL_DATA', None):
    UserWarning('Environment variable GDAL_DATA does not exist. Setting to {}'.format(
        r'C:/Miniconda3/envs/pvtools27/Library/share/gdal'))
    os.environ['GDAL_DATA'] = r'C:/Miniconda3/envs/pvtools27/Library/share/gdal'

if not os.environ.get('PROJ_LIB', None):
    UserWarning('Environment variable PROJ_LIB does not exist. Setting to {}'.format(
        r'C:/Miniconda3/envs/pvtools27/Library/share'))
    os.environ['PROJ_LIB'] = r'C:/Miniconda3/envs/pvtools27/Library/share'

out_dir = tempfile.gettempdir()
out_dir = r'C:\data\temp'
out_dir = os.path.join(out_dir, os.path.splitext(pyFile)[0])

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

log_file = '{0}.log'.format(pyFile.replace('.py', time.strftime('_%Y-%m-%d', time.localtime())))
logging.captureWarnings(False)
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(name)s : %(message)s", filename=os.path.join(out_dir, log_file), filemode='w')
# https://stackoverflow.com/questions/9321741/printing-to-screen-and-writing-to-a-file-at-the-same-time
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# set a format which is simpler for console use
formatter = logging.Formatter("%(message)s")
# tell the handler to use this format
console.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger('').addHandler(console)

fileList = make_dummy_data.make_dummy_tif_files(out_dir)

# Need to get the wktproj of the raster from gdal NOT rasterio.
# RasterIO works from the proj4 string NOT the wkt string so aussie zones details gets lost.
rast_crs = pyprecag.crs.getCRSfromRasterFile(fileList[0])

ptcount=50
with rasterio.open(os.path.normpath(fileList[0])) as src:
    rand_pts_gdf,rand_pts_crs = randomPixelSelection(src, rast_crs, ptcount )

start_time = time.time()
_ = extractPixelStatisticsForPoints(rand_pts_gdf, rand_pts_crs, fileList, function_list=[np.nanmean, np.nanstd],
                                    size_list=[1, 3, 7],
                                    output_csvfile=os.path.join(out_dir,'random_points_{}.csv'.format(ptcount)))

print(os.path.join(out_dir,'random_points_{}.csv'.format(ptcount)))
print("\n\n***********\n{} is Complete\n".format(pyFile))
