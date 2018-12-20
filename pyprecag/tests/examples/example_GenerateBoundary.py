import datetime
import glob
import logging
import operator
import os
import tempfile
import time

from pyprecag import convert
from pyprecag import processing

DEBUG = True
pyFile = os.path.basename(__file__)

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
console.setLevel(logging.DEBUG)
# set a format which is simpler for console use
formatter = logging.Formatter("%(message)s")
# tell the handler to use this format
console.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger('').addHandler(console)

fileList = [r'../data/area1_yield_file_ascii_wgs84.csv',
            r'../data/area2_yield_file_ISO-8859-1.csv'
           ]

startTime = time.time()
for ea in fileList:
    loopTime = time.time()
    logging.info('{lne}\ncreatePolygonFromPointTrail\t {}'.format(ea,lne='=' * 100))

    outPtsFile = os.path.join(out_dir, os.path.splitext(os.path.basename(ea))[0] + '_points.shp')
    outPolyFile = os.path.join(out_dir, os.path.splitext(os.path.basename(ea))[0] + '_poly.shp')

    gdfPts, gdfCrs = convert.convertCsvToPoints(ea, outPtsFile, coord_columns_EPSG=4326, out_EPSG=-1)

    processing.createPolygonFromPointTrail(gdfPts,gdfCrs, outPolyFile,
                                           thin_dist_m=2.5,
                                           aggregate_dist_m=25,
                                           buffer_dist_m=10,
                                           shrink_dist_m=3)

    logging.info('{:<30} {dur:<15} {}'.format('File Complete!!', outPolyFile,
                                       dur=datetime.timedelta(seconds=time.time() - loopTime)))

logging.info('{}'.format('_' * 100))
logging.info('{:<30} {dur:<15} {}'.format('Completed All {} Files'.format(len(fileList)), '',
                                          dur=datetime.timedelta(seconds=time.time() - startTime)))
print "\n\n***********\n{} is Complete\n".format(pyFile)
