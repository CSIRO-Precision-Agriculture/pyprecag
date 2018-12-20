import datetime
import logging
import os,sys
import tempfile
import time
import rasterio
import geopandas
from pyprecag.bandops import BandMapping, CalculateIndices
from pyprecag.processing import calc_indices_for_block, resample_bands_to_block

os.environ['GDAL_DATA']=r'C:/Miniconda3/envs/pvtools27/Library/share/gdal'
os.environ['PROJ_LIB']=r'C:/Miniconda3/envs/pvtools27/Library/share'
from osgeo import gdal

global DEBUG
DEBUG = True

pyFile = os.path.basename(__file__)

out_dir = tempfile.gettempdir()
out_dir = r'C:\data\temp'
out_dir = os.path.join(out_dir, os.path.splitext(pyFile)[0])

if not os.path.exists(out_dir): os.mkdir(out_dir)

log_file = '{0}.log'.format(pyFile.replace('.py', time.strftime('_%Y-%m-%d_%H%M', time.localtime())))
logging.captureWarnings(False)
logging.basicConfig(level=logging.INFO, format="%(message)s", filename=os.path.join(out_dir, log_file))   #, filemode='a')
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

logging.info('Python: {} {}'.format(sys.executable, sys.version))
logging.info('GDAL Version: {}'.format(gdal.__version__))
logging.info('RasterIO Version: {}'.format(rasterio.__version__ ))
logging.info('GeoPandas Version: {}'.format(geopandas.__version__))

"""Full options: ie reproject, calculate indices, process for individual polygons. """
fileList = [(r'..\data\area1_rgbi_jan_50cm_84sutm54.tif', 32754, 2.5, BandMapping(infrared=4,rededge=1,mask=5), 0,
                r'..\data\area1_blocklayout_94mga54.shp','')
            , (r'..\data\area1_rgbi_jan_50cm_84sutm54.tif', 28354, 2, BandMapping(infrared=4, rededge=1), 0,'', '')
            , (r'..\data\area1_rgbi_jan_50cm_84sutm54.bil', 32754, 2, BandMapping( green=2, infrared=4,mask=5), 0,
                r'..\data\area1_blocklayout_94mga54.shp','NAME')
            , (r'..\data\area1_rgbi_jan_50cm_84sutm54.bil', 28354, 2, BandMapping(green=2, infrared=4, rededge=1, mask=5), 0,
                r'..\data\area1_blocklayout_94mga54.shp', 'NAME')
            , (r'..\data\area1_rgbi_jan_50cm_84sutm54.tif', 32754,2, BandMapping(red=3, infrared=4), 0,
                r'..\data\PolyMZ_wgs84_MixedPartFieldsTypes.shp','part_type')
            ]

startTime = time.time()
for i,ea in enumerate(fileList, start=1):
    loopTime = time.time()

    image,image_epsg,pix_size,bm,nodata_val,poly,grpfield = ea
    logging.info('\n{lne}\n{}, {} \n{}\n pixel: {}, nodata: {} \n grpby: {}, poly: {} .'.format(image,image_epsg,bm,pix_size,nodata_val,grpfield,poly, lne='=' * 100))

    mydir = os.path.join(out_dir, 'it{}'.format(i))
    if not os.path.exists(mydir): os.mkdir(mydir)
    logging.info('\n***CalculateImageIndices {} of {}\n'.format(i, len(fileList)))
    indices = CalculateIndices(**bm).valid_indices()
    files = calc_indices_for_block(image, pix_size, bm, mydir, indices, image_epsg, nodata_val, poly, grpfield, out_epsg=28354)
    logging.info ('{:<30} {:>10}   {dur:<15}, {}'.format('Indices Complete!!',indices,files, dur=datetime.timedelta(seconds=time.time() - loopTime)))

    logging.info('\n***Resample Band {} of {}'.format(i, len(fileList)))
    files = resample_bands_to_block( image, pix_size, mydir, [6], image_epsg, nodata_val, poly, grpfield, out_epsg=28354)
    logging.info('{:<30} {:>10}   {dur:<15}, {}'.format('Resampling Complete !!', i, files, dur=datetime.timedelta(seconds=time.time() - loopTime)))

print "\n\n***********\n{} is Complete\n".format(pyFile)
