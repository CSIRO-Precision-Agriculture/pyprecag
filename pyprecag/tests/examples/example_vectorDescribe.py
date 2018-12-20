import datetime
import glob
import logging
import os
import tempfile
import time

from pyprecag.describe import VectorDescribe

DEBUG = True

pyFile = os.path.basename(__file__)

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

fileList = glob.glob(os.path.join(r'../data/*.shp'))

startTime = time.time()
for ea in fileList:
    loopTime = time.time()
    print('{lne}\nvectorDescribe \t {}'.format(ea, lne='=' * 100))

    vectDescObj = VectorDescribe(ea)
    print('\tFile Encoding: {}'.format(vectDescObj.file_encoding))
    print('\tFeature Count: {}'.format(vectDescObj.feature_count))
    print('\tGeometry Type: {}'.format(vectDescObj.geometry_type))
    print('\tZM Aware: {}'.format(vectDescObj.is_mz_aware))
    print('\tColumn Count: {}'.format(len(vectDescObj.column_properties)))
    # print('\tOGRDriver: {}'.format(vectDescObj.ogrDriver.name))
    print('\tExtent: {}'.format(vectDescObj.extent))

    print('\tColumn Names: {}'.format(vectDescObj.get_column_names()))
    print('\tAlias Names: {}'.format(vectDescObj.get_alias_column_names()))
    print('\tShapefile Names: {}'.format(vectDescObj.get_column_names()))
    print('\tColumn Types: {}'.format(vectDescObj.get_column_types()))
    print('\n\tColumn Properties:')
    print(vectDescObj.column_properties)

    if vectDescObj.crs.epsg is not None:
        print('\n\tCoordinate System Properties:')
        print('\t\tEPSG: {}   Predicted: {}'.format(vectDescObj.crs.epsg, vectDescObj.crs.epsg_predicted))
        print('\t\tIs Projected?: {}'.format(bool(vectDescObj.crs.srs.IsProjected())))
        print('\t\tSRS: {}'.format(str(vectDescObj.crs.srs).replace('\n', '\n\t\t')))


print('{}'.format('_' * 100))
print('{:<30} {dur:<15} {}'.format('Completed All {} Files'.format(len(fileList)), '',
                                          dur=datetime.timedelta(seconds=time.time() - startTime)))
print "\n\n***********\n{} is Complete\n".format(pyFile)
