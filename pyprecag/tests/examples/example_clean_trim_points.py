import datetime
import logging
import os
import tempfile
import time

from pyprecag import convert
from pyprecag.processing import cleanTrimPoints

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

fileList = [(4326, 28354, "Yield", r'../data/area1_yield_file_ascii_wgs84.csv',
             r'../data/area1_onebox_94mga54.shp'),
            (4326, 28354, 'Yld Mass(Dry)(tonne/ha)',
             r"../data/area2_yield_file_ISO-8859-1.csv",
             r'../data/area2_onebox_94mga54.shp')
           ]

for ea in fileList:
    startTime = time.time()
    logging.info ('{lne}\ncleanTrimPoints \t {}'.format(ea, lne='=' * 100))

    in_epsg, out_epsg, process_column, csv_file, poly = ea

    output_csv = os.path.join(out_dir, os.path.splitext(os.path.basename(csv_file))[0] + '_trimmed.csv')
    output_shp = os.path.join(out_dir, os.path.splitext(os.path.basename(csv_file))[0] + '_trimmed.shp')
    out_remove_shp = os.path.join(out_dir, os.path.splitext(os.path.basename(csv_file))[0] + '_removed.shp')

    loopTime = time.time()
    tabHeader = '\n{:<30} {:>10}   {:<15} {}'.format('Task', 'FeatCount', 'Duration H:M:S', 'Comments')
    logging.info(tabHeader)
    logging.info('-' * len(tabHeader))
    gdfPoints, gdfPtsCrs = convert.convertCsvToPoints(csv_file, coord_columns_EPSG=in_epsg, out_EPSG=out_epsg)

    logging.info('{:<30} {:>10,}   {dur:<15}'.format('Convert To Points', len(gdfPoints),
                                                 dur=datetime.timedelta(seconds=time.time() - loopTime)))

    output, _ = cleanTrimPoints(gdfPoints, gdfPtsCrs, process_column, output_csv, out_keep_shapefile=output_shp,
                             out_removed_shapefile=out_remove_shp, boundary_polyfile=poly, thin_dist_m=1)

    logging.info( '{:<30}\t{:>10,}\tDuration (H:M:S): {dur:<15}\t'.format('File Complete !!', len(output), '',
                                                              dur=datetime.timedelta(seconds=time.time() - startTime)))

    print "\n\n***********\n{} is Complete\n".format(pyFile)
