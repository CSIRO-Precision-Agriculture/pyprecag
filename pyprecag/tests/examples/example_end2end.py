import datetime
import os
import tempfile

import numpy as np
from osgeo import gdal
import rasterio

import time
import logging
from pyprecag import convert, processing, kriging_ops, crs, describe ,raster_ops
from pyprecag.bandops import BandMapping, CalculateIndices
from pyprecag.processing import calc_indices_for_block, resample_bands_to_block

global DEBUG
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

fileList = [([], 4326, 30, 2.5, 'Yld Mass(Dry)(tonne/ha)',r'../data/area2_yield_file_ISO-8859-1.csv',
             r'../data/area2_onebox_94mga54.shp','','')
            , ([], 4326, 10, 2.5, 'Yield', r'../data/area1_yield_file_ascii_wgs84.csv',
             r'../data/PolyMZ_wgs84_MixedPartFieldsTypes.shp', 'part_type', r'..\data\area1_rgbi_jan_50cm_84sutm54.tif')
           ]

startTime = time.time()
for ea in fileList:
    coord_columns, csvepsg, block_size, pixel_size, data_column, csvFile, polyFile,grpby,tifFile = ea
    logging.info('{}\nProcessing {}'.format('_' * 100, os.path.basename(csvFile)))
    fileTime = time.time()
    taskTime = time.time()
    tabHeader = '\n{:<30} {:<15} {}'.format('Task', 'Duration H:M:S', 'Comments')
    logging.info(tabHeader)
    logging.info('-' * len(tabHeader))

    fileSubName = os.path.join(out_dir, os.path.splitext(os.path.basename(csvFile))[0])

    gdfPts, ptsCrs = convert.convertCsvToPoints(csvFile, coord_columns_EPSG=csvepsg, coord_columns=coord_columns,
                                                out_EPSG=-1)

    # this will allow you to create the polygon for the blockgrid step, but not use it for clipping in the clean/trim step
    polyShp = polyFile

    taskTime = time.time()
    if polyFile == '' or polyFile is None:
        polyShp = fileSubName + '_bnd.shp'
        processing.createPolygonFromPointTrail(gdfPts, ptsCrs, polyShp,
                                               thin_dist_m=1,
                                               aggregate_dist_m=25,
                                               buffer_dist_m=10,
                                               shrink_dist_m=3)

        taskTime = time.time()

    descPoly = describe.VectorDescribe(polyShp)
    if ptsCrs.epsg != descPoly.crs.epsg:
        gdfPoly = descPoly.open_geo_dataframe()
        gdfPoly.to_crs(epsg=ptsCrs.epsg_number, inplace=True)
        logging.info('Projecting Polygon to {}'.format(ptsCrs.epsg))
        gdfPoly.to_file(fileSubName + '_projpoly.shp', driver='ESRI Shapefile')
        polyShp = fileSubName + '_projpoly.shp'

    processing.BlockGrid(in_shapefilename=polyShp,
                         pixel_size=pixel_size,
                         out_rasterfilename=fileSubName + '_block.tif',
                         out_vesperfilename=fileSubName + '_block_v.txt',
                         snap=True,
                         overwrite=True)

    out_csvFile = fileSubName + '_trimmed.csv'
    outGDF, outCRS = processing.cleanTrimPoints(gdfPts, ptsCrs, data_column, out_csvFile,
                                                boundary_polyfile=polyFile, thin_dist_m=2.5)

    bat_file, ctrl_file = kriging_ops.prepareForVesperKrig(outGDF, data_column, fileSubName + '_block_v.txt', out_dir,
                                                           block_size=block_size, coord_columns=[])

    kriging_ops.run_vesper(ctrl_file)
    taskTime = time.time()

    pred, _, _ = kriging_ops.vesperTextToRaster(ctrl_file, outCRS.epsg_number)
    fileList = [pred]
    rast_crs = crs.getCRSfromRasterFile(pred)

    taskTime = time.time()
    with rasterio.open(os.path.normpath(pred)) as src:
        out_meta = src.meta.copy()
        out_meta['crs'] = rast_crs.crs_wkt
        rescaled = raster_ops.rescale(src, 0, 1)
        rescaled2 = raster_ops.rescale(src, 0, 255)
        normalised = raster_ops.normalise(src)

        rand_pts_gdf, rand_pts_crs = processing.randomPixelSelection(src, rast_crs, 50, fileSubName + '_randompts.shp')

    out_meta['count'] = 1
    for ea in [(rescaled,'_rescaled0-1.tif'),(rescaled2,'_rescaled0-255.tif'),(normalised,'_normalised.tif')]:
        fileList += [fileSubName + ea[-1]]
        with rasterio.open(os.path.normpath(fileList[-1]), 'w', **out_meta) as dst:
            dst.write_band(1, ea[0])

    out_meta['count'] = 3
    with rasterio.open(os.path.normpath(fileSubName + '_multi_normstd.tif'), 'w', **out_meta) as dst:
        dst.write_band(1, rescaled)
        dst.write_band(2, rescaled2)
        dst.write_band(3, normalised)

    del rescaled, rescaled2, normalised

    logging.info('{:<30}\t{dur:<15}'.format('Norm/Std', dur=datetime.timedelta(seconds=time.time() - taskTime)))
    taskTime = time.time()

    if tifFile != '':
        bm = BandMapping(red=3, infrared=4, rededge=1, mask=5)
        indices = CalculateIndices(**bm).valid_indices()
        files = calc_indices_for_block(tifFile, pixel_size, bm, out_dir, indices, image_nodata=0, image_epsg=32754,
                                       polygon_shapefile=polyFile, out_epsg=28354)
        fileList +=files

        logging.info('{:<30}\t{dur:<15}'.format('calc_indices_for_block', dur=datetime.timedelta(seconds=time.time() - taskTime)))
        taskTime = time.time()

        files = resample_bands_to_block(tifFile, pixel_size, out_dir,band_nums=[6],image_nodata=0, image_epsg=32754,
                                      polygon_shapefile=polyFile, out_epsg=28354)

        fileList += files


        logging.info('{:<30}\t{dur:<15}'.format('resample_bands_to_block', dur=datetime.timedelta(seconds=time.time() - taskTime)))
        taskTime = time.time()

    _ = processing.extractPixelStatisticsForPoints(rand_pts_gdf, rand_pts_crs, fileList, function_list=[np.nanmean, np.nanstd],
                                        size_list=[1, 3, 7],
                                        output_csvfile=fileSubName + '_gridextract.csv')

    logging.info('{:<30}\t{dur:<15}\t{}'.format('*File Complete', os.path.basename(csvFile),
                                                dur=datetime.timedelta(seconds=time.time() - fileTime)))
logging.info('{}'.format('_' * 100))
logging.info('{:<30} {dur:<15} {}'.format('Completed All {} Files'.format(len(fileList)), 'Includes Vesper Kriging Processing',
                                          dur=datetime.timedelta(seconds=time.time() - startTime)))


print "\n\n***********\n{} is Complete\n".format(pyFile)
