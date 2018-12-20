import datetime
import logging
import os
import tempfile
import time

from pyprecag import convert, vector_ops

DEBUG = True
pyFile = os.path.basename(__file__)

out_dir = tempfile.gettempdir()
out_dir = r'C:\data\temp'
out_dir = os.path.join(out_dir, os.path.splitext(pyFile)[0])

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

log_file = '{0}.log'.format(pyFile.replace('.py', time.strftime('_%Y-%m-%d', time.localtime())))
logging.captureWarnings(False)
logging.basicConfig(level=logging.INFO, format="%(message)s", filename=os.path.join(out_dir, log_file), filemode='w')
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

fileList = [(4326, 28354, r'../data/area1_yield_file_ascii_wgs84.csv'),
            (4326, 28354, r"../data/area2_yield_file_ISO-8859-1.csv")
            # , (4326, -1, r'C:\data\_projects\PVTools_data\TestData\National_Paddocks\Graham_Spackman\David_Storey_Stranraer\2015\2015_Sorghum_David Storey Stranraer.csv',)
            # , (4326, -1, r'C:\data\_projects\PVTools_data\TestData\National_Paddocks\Hilary_Wittwer\Steve_Lyneham_Glenlivit_7\2013\2013_Barley_Steve_Lyneham_Glenlivit_7.csv')
            # , (4326, -1, r'C:\data\_projects\PVTools_data\TestData\National_Paddocks\Chris_Minehan\Ben_Beck_P3\2014\2014_Canola_Ben_Beck_P3.csv')
            # , (4326, -1, r'C:\data\_projects\PVTools_data\TestData\National_Paddocks\Luke_Marquis\Paul_and_Ainsley_Foulds_Yippee_Downs\2015\2015_Wheat_Paul_and_Ainsley_Foulds_Yippee_Downs.csv')
           ]

# NOTE: out_EPSG= of -1 will calculate it from the data

startTime = time.time()
for ea in fileList:
    logging.info('{lne}\nthinPoints \t {}'.format(ea, lne='=' * 100))

    in_epsg, out_epsg, csv_file = ea

    output_shp = os.path.join(out_dir, os.path.splitext(os.path.basename(csv_file))[0] + '_thinned.shp')

    loopTime = time.time()
    tabHeader = '\n{:<30} {:>10}   {:<15} {}'.format('Task', 'FeatCount', 'Duration H:M:S', 'Comments')
    logging.info(tabHeader)
    logging.info('-' * len(tabHeader))
    gdfPoints, gdfPtsCrs = convert.convertCsvToPoints(csv_file, coord_columns_EPSG=in_epsg, out_EPSG=out_epsg)

    logging.info('{:<30} {:>10,}   {:<15}'.format('Convert To Points and Reprj', len(gdfPoints),
                                                 datetime.timedelta(seconds=time.time() - loopTime)))
    loopTime = time.time()

    vector_ops.thin_point_by_distance(gdfPoints, gdfPtsCrs, 2.5, output_shp)
    logging.info ('{:<30} {:>10,}   {dur:<15}'.format('File Complete!!', len(gdfPoints),
                                            dur=datetime.timedelta(seconds=time.time() - loopTime)))
logging.info('{}'.format('_' * 100))
logging.info('{:<30} {dur:<15} {}'.format('Completed All {} Files'.format(len(fileList)), '',
                                          dur=datetime.timedelta(seconds=time.time() - startTime)))

print "\n\n***********\n{} is Complete\n".format(pyFile)
