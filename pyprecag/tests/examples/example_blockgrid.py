import datetime
import logging
import os
import tempfile
import time

from pyprecag import processing

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

fileList = [r'../data/area1_onebox_94mga54.shp',
            r'../data/area2_onebox_94mga54.shp'
           ]

startTime = time.time()
for ea in fileList:
    loopTime = time.time()
    logging.info('{lne}\nBlockGrid \t {}'.format(ea,lne='=' * 100))

    outputRaster = os.path.join(out_dir, os.path.basename(ea.replace('.shp', '.tif')))
    outputVesperFile = os.path.join(out_dir, os.path.basename(ea.replace('.shp', '_v.txt')))
    processing.BlockGrid(in_shapefilename=ea,
                         pixel_size=2.5,
                         out_rasterfilename=outputRaster,
                         out_vesperfilename=outputVesperFile,
                         snap=True,
                         overwrite=True)
    logging.info('{:<30} {dur:<15} {}'.format('Complete', outputRaster,
                                       dur=datetime.timedelta(seconds=time.time() - loopTime)))

logging.info('{}'.format('_' * 100))
logging.info('{:<30} {dur:<15} {}'.format('Completed All {} Files'.format(len(fileList)), '',
                                          dur=datetime.timedelta(seconds=time.time() - startTime)))

print "\n\n***********\n{} is Complete\n".format(pyFile)
