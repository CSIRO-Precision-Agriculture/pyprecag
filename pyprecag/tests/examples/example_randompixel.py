import datetime
import logging
import os
import tempfile
import time

import rasterio
from osgeo import gdal

from pyprecag.tests.make_dummy_data import make_dummy_tif_files

import pyprecag.crs
from pyprecag.processing import randomPixelSelection

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
#logging.captureWarnings(False)
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(name)s : %(message)s", filename=os.path.join(out_dir, log_file), filemode='w')
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

startTime = time.time()
for ea in fileList:
    loopTime = time.time()
    logging.info ('{lne}\nrandomPixelSelection \t {}'.format(ea,lne='=' * 100))

    out_shp = os.path.join(out_dir, os.path.basename(ea).replace('.tif', '_randpts.shp'))
    rast_crs = pyprecag.crs.getCRSfromRasterFile(ea)

    with rasterio.open(os.path.normpath(ea)) as raster:
        rand_gdf, _ = randomPixelSelection(raster, rast_crs, 50, out_shapefile=out_shp)

    logging.info('{:<30} {dur:<15} {}'.format('Complete', out_shp,
                                       dur=datetime.timedelta(seconds=time.time() - loopTime)))

logging.info('{}'.format('_' * 100))
logging.info('{:<30} {dur:<15} {}'.format('Completed All {} Files'.format(len(fileList)), '',
                                          dur=datetime.timedelta(seconds=time.time() - startTime)))

print "\n\n***********\n{} is Complete\n".format(pyFile)
