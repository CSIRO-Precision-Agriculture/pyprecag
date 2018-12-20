import glob
import os
import tempfile
import time
import datetime

import logging
from pyprecag.describe import CsvDescribe, predictCoordinateColumnNames

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

fileList = glob.glob(os.path.join(r'../data/*.csv'))

startTime = time.time()
for ea in fileList:
    loopTime = time.time()
    print('{lne}\ncsvDescribe \t {}'.format(ea,lne='=' * 100))

    csvDescObj = CsvDescribe(ea)

    print('\tDialect: \tdelimiter={}\tquotechar={}\tlineterminator={}'.format(
        repr(csvDescObj.dialect.delimiter), repr(csvDescObj.dialect.quotechar),
        repr(csvDescObj.dialect.lineterminator)))

    print('\tRow Count: {}'.format(csvDescObj.row_count))
    print('\tColumn Count: {}'.format(csvDescObj.column_count))
    print('\tFile Encoding: {}'.format(csvDescObj.file_encoding))
    if csvDescObj.has_column_header:
        print('\tCoordinate Columns (XY): {} '.format(predictCoordinateColumnNames(csvDescObj.get_column_names())))

    print('\tColumn Names: {}'.format(csvDescObj.get_column_names()))
    print('\tAlias Names: {}'.format(csvDescObj.get_alias_column_names()))
    print('\tColumn Types: {}'.format(csvDescObj.get_column_types()))
    print('\tColumn Properties:')
    # general.displayTableDictionary(csvDescObj.column_properties, '\t\t')


print('{}'.format('_' * 100))
print('{:<30} {dur:<15} {}'.format('Completed All {} Files'.format(len(fileList)), '',
                                          dur=datetime.timedelta(seconds=time.time() - startTime)))

print "\n\n***********\n{} is Complete\n".format(pyFile)
