import logging
import os
import platform
import shutil
import tempfile
import unittest

import pandas as pd
import numpy as np
import rasterio
import time

from pyprecag import convert, crs
from pyprecag.bandops import CalculateIndices, BandMapping
from pyprecag.describe import VectorDescribe, CsvDescribe, predictCoordinateColumnNames
from pyprecag.processing import clean_trim_points, create_polygon_from_point_trail, block_grid, random_pixel_selection, \
    calc_indices_for_block, resample_bands_to_block, extract_pixel_statistics_for_points, kmeans_clustering
from pyprecag.raster_ops import rescale, normalise
try:
    from pyprecag.kriging_ops import prepare_for_vesper_krige, vesper_text_to_raster, run_vesper
except ImportError as e:
    if "Vesper" not in e.message:
        raise e

pyFile = os.path.basename(__file__)

TmpDir = tempfile.gettempdir()
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])

this_dir = os.path.abspath(os.path.dirname(__file__))

logging.captureWarnings(True)
logging.basicConfig(level=logging.INFO, format="%(message)s")

fileCSV = os.path.realpath(this_dir + "/data/area1_yield_ascii_wgs84.csv")
fileBox = os.path.realpath(this_dir + "/data/area1_onebox_94mga54.shp")
fileBoxes = os.path.realpath(this_dir + "/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp")
fileImage = os.path.realpath(this_dir + "/data/area1_rgbi_jan_50cm_84sutm54.tif")

epsg = 28354


class test_End2End(unittest.TestCase):
    gridextract_files=[]
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_End2End, cls).setUpClass()
        if os.path.exists(TmpDir):
            print 'Folder Exists.. Deleting {}'.format(TmpDir)
            shutil.rmtree(TmpDir)

        os.mkdir(TmpDir)

        global testFailed
        testFailed = False

    @classmethod
    def tearDownClass(cls):
        if not testFailed:
            print ('Tests Passed .. Deleting {}'.format(TmpDir))
            shutil.rmtree(TmpDir)

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def run(self, result=None):
        global testFailed
        unittest.TestCase.run(self, result) # call superclass run method
        if len(result.failures) > 0 or len(result.errors) > 0:
            testFailed = True

    def test01_csvDescribe_ASCII(self):
        csvDesc = CsvDescribe(fileCSV)

        self.assertEqual(csvDesc.file_encoding, 'ascii')
        self.assertEqual(csvDesc.row_count, 13756)
        self.assertEqual(csvDesc.column_count, 24)
        self.assertEqual(predictCoordinateColumnNames(csvDesc.get_column_names()), ['Lon', 'Lat'])
        self.assertTrue(csvDesc.has_column_header)

    def test02_createPolygonFromPointTrail(self):
        global filePoints,filePoly
        filePoints = os.path.join(TmpDir, os.path.splitext(os.path.basename(fileCSV))[0] + '_points.shp')
        filePoly = os.path.join(TmpDir, os.path.splitext(os.path.basename(fileCSV))[0] + '_poly.shp')

        if not os.path.exists(filePoly):
            gdfPts, gdfCrs = convert.convert_csv_to_points(fileCSV, filePoints, coord_columns_epsg=4326, out_epsg=epsg)

            create_polygon_from_point_trail(gdfPts, gdfCrs, filePoly,
                                            thin_dist_m=2.5,
                                            aggregate_dist_m=25,
                                            buffer_dist_m=7,
                                            shrink_dist_m=3)

        self.assertTrue(os.path.exists(filePoly), True)

    def test03_vectorDescribe(self):
        vDesc = VectorDescribe(filePoly)
        self.assertEqual(vDesc.crs.epsg_number, epsg)
        self.assertFalse(vDesc.is_mz_aware)
        self.assertEqual(vDesc.geometry_type, 'Polygon')

        vDesc = VectorDescribe(filePoints)
        self.assertEqual(vDesc.crs.epsg_number, epsg)
        self.assertFalse(vDesc.is_mz_aware)
        self.assertEqual(vDesc.geometry_type, 'Point')

        vDesc = VectorDescribe(fileBox)
        self.assertEqual(vDesc.crs.epsg_number, epsg)
        self.assertFalse(vDesc.is_mz_aware)
        self.assertEqual(vDesc.geometry_type, 'Polygon')

    def test04_blockGrid(self):

        global fileBlockTif,fileBlockTxt
        fileBlockTif = os.path.join(TmpDir, os.path.splitext(os.path.basename(fileCSV))[0] + '_block.tif')
        fileBlockTxt = os.path.join(TmpDir, os.path.splitext(os.path.basename(fileCSV))[0] + '_block_v.txt')

        if not os.path.exists(fileBlockTif):
            block_grid(in_shapefilename=fileBox,
                       pixel_size=2.5,
                       out_rasterfilename=fileBlockTif,
                       out_vesperfilename=fileBlockTxt,
                       snap=True,
                       overwrite=True)

        vDesc = VectorDescribe(fileBox)

        self.assertTrue(os.path.exists(fileBlockTif))
        self.assertTrue(os.path.exists(fileBlockTxt))

        with rasterio.open(os.path.normpath(fileBlockTif)) as dataset:
            self.assertEqual(dataset.count,1)
            self.assertEqual(dataset.width,55)
            self.assertEqual(dataset.height,55)
            self.assertEqual(dataset.nodatavals,(-9999.0,))
            self.assertEqual(dataset.dtypes, ('int16',))

            print(vDesc.crs.crs_wkt)
            print('SRS:\t{}'.format(dataset.crs))
        print('Temp Files in {}'.format(TmpDir))

    def test05_cleanTrimPoints(self):
        global fileTrimmed,data_col
        fileTrimmed = os.path.join(TmpDir, os.path.splitext(os.path.basename(fileCSV))[0] + '_normtrimmed.csv')
        file_shp = os.path.join(TmpDir, os.path.splitext(os.path.basename(fileCSV))[0] + '_normtrimmed.shp')
        file_removed = os.path.join(TmpDir, os.path.splitext(os.path.basename(fileCSV))[0] + '_remove.shp')

        data_col = r'Yield'

        gdfPoints, gdfPtsCrs = convert.convert_csv_to_points(fileCSV, coord_columns_epsg=4326, out_epsg=epsg)
        gdfOut, crsOut = clean_trim_points(gdfPoints, gdfPtsCrs, data_col, fileTrimmed, out_keep_shapefile=file_shp,
                                           out_removed_shapefile=file_removed, boundary_polyfile=fileBox,
                                           thin_dist_m=2.5)

        self.assertTrue(os.path.exists(fileTrimmed))
        self.assertTrue(os.path.exists(file_shp))
        self.assertTrue(os.path.exists(file_removed))
        self.assertEqual(gdfOut.crs, {'init': 'epsg:{}'.format(epsg), 'no_defs': True})
        self.assertEqual(len(gdfOut), 648)
        self.assertIn('nrm_' + data_col, gdfOut.columns)
        self.assertIn('Easting',gdfOut.columns)
        self.assertIn('Northing', gdfOut.columns)
        self.assertIn('EN_EPSG', gdfOut.columns)


    @unittest.skipIf(
        platform.system() != 'Windows',
        'Vesper only present on Windows'
    )
    def test06_prepareForVesperKrig(self):
        descCSV = CsvDescribe(fileTrimmed)
        dfCSV = descCSV.open_pandas_dataframe()

        global fileControl
        subFile = os.path.splitext(os.path.basename(fileCSV))[0]
        fileControl = subFile + '_control_' + data_col + '.txt'

        if not os.path.exists(fileControl):
            bat_file, fileControl = prepare_for_vesper_krige(dfCSV, data_col, fileBlockTxt, TmpDir, block_size=30,
                                                             control_textfile=fileControl, coord_columns=[], epsg=epsg)

            self.assertTrue(os.path.exists(bat_file))
            self.assertTrue(os.path.exists(fileControl))

        self.assertTrue(os.path.exists( os.path.join(TmpDir,r'Vesper',subFile + '_vesperdata_' + data_col + '.csv')))
        dfCSV = pd.read_csv(os.path.join(TmpDir,r'Vesper',subFile + '_vesperdata_' + data_col + '.csv'))

        x_column, y_column = predictCoordinateColumnNames(dfCSV.columns)
        self.assertEqual(x_column.upper(), 'EASTING')
        self.assertEqual(y_column.upper(), 'NORTHING')

        print('Running Vesper, Please wait....')
        run_vesper(fileControl)


    @unittest.skipIf(
        platform.system() != 'Windows',
        'Vesper only present on Windows'
    )
    def test07_vesperTextToRaster(self):
        global out_PredTif
        out_PredTif, out_SETif, out_CITxt = vesper_text_to_raster(fileControl, epsg)
        for eaFile in [out_PredTif, out_SETif, out_CITxt]:
            self.assertTrue(os.path.exists(eaFile))

        with rasterio.open(os.path.normpath(out_PredTif)) as dataset:
            self.assertEqual(dataset.count, 1)
            self.assertEqual(dataset.width, 54)
            self.assertEqual(dataset.height, 53)
            self.assertEqual(dataset.nodatavals, (-9999.0,))
            self.assertEqual(dataset.dtypes, ('float32',))
            self.assertEqual(dataset.crs,  rasterio.crs.CRS.from_epsg(28354))

        with rasterio.open(os.path.normpath(out_SETif)) as dataset:
            self.assertEqual(dataset.count, 1)
            self.assertEqual(dataset.width, 54)
            self.assertEqual(dataset.height, 53)
            self.assertEqual(dataset.nodatavals, (-9999.0,))
            self.assertEqual(dataset.dtypes, ('float32',))
            self.assertEqual(dataset.crs, rasterio.crs.CRS.from_epsg(28354))

    def test08_randomPixelSelection(self):
        global rand_gdf, rand_crs
        out_randompts = os.path.join(TmpDir, os.path.basename(fileBlockTif).replace('.tif', '_randpts.shp'))
        rast_crs = crs.getCRSfromRasterFile(fileBlockTif)

        with rasterio.open(os.path.normpath(fileBlockTif)) as raster:
            rand_gdf, rand_crs = random_pixel_selection(raster, rast_crs, 50, out_shapefile=out_randompts)
        self.assertEqual(len(rand_gdf), 50)
        self.assertTrue(os.path.exists(out_randompts))
        self.assertEqual(rast_crs.epsg, rand_gdf.crs)

    def test09_calcImageIndices_allopts(self):
        out_fold = os.path.join(TmpDir,'calcindex_allopts')
        if not os.path.exists(out_fold): os.mkdir(out_fold)
        bm = BandMapping(green=2, infrared=4, rededge=1, mask=5)
        indices = CalculateIndices(**bm).valid_indices()
        files = calc_indices_for_block(fileImage,2.5,bm,out_fold,indices,image_nodata=0, image_epsg=32754,
                                       polygon_shapefile=fileBox,out_epsg=28354)

        self.gridextract_files += files

        self.assertEqual(len(files), 3)

        with rasterio.open(files[0]) as src:
            self.assertEqual(src.nodata, -9999)
            self.assertEqual(src.crs, rasterio.crs.CRS.from_epsg(28354))
            self.assertEqual(src.meta['dtype'], 'float32')
            self.assertEqual(src.res, (2.5,2.5))
            self.assertEqual(src.count,1)

    def test10_resampleBands2Block_allopts(self):
        out_fold = os.path.join(TmpDir,'resamp2block_allopts')
        if not os.path.exists(out_fold): os.mkdir(out_fold)

        files = resample_bands_to_block(fileImage,2.5,out_fold,band_nums=[6],image_nodata=0, image_epsg=32754,
                                        polygon_shapefile=fileBox,out_epsg=28354)

        self.gridextract_files +=files
        self.assertEqual(len(files),1)

        with rasterio.open(files[0]) as src:
            self.assertEqual(src.nodata,0.0)
            self.assertEqual(src.crs, rasterio.crs.CRS.from_epsg(28354))
            self.assertEqual(src.meta['dtype'], 'float32')
            self.assertEqual(src.res, (2.5,2.5))
            self.assertEqual(src.count,1)

    def test11_rescaleNormaliseRaster(self):
        in_file = self.gridextract_files[-1]
        rast_crs = crs.getCRSfromRasterFile(in_file)

        with rasterio.open(os.path.normpath(in_file)) as src:
            rescaled = rescale(src, 0, 255)
            rescaled2 = rescale(src, 0, 5)
            norm = normalise(src)
            out_meta = src.meta.copy()

        out_meta['crs'] = rast_crs.crs_wkt
        out_meta['count'] = 1  # contains only one band
        out_meta['dtype'] = np.float32

        out_rescale = os.path.join(TmpDir, os.path.basename(in_file).replace('.tif', '_rescale0-255.tif'))
        with rasterio.open(os.path.normpath(out_rescale), 'w', **out_meta) as out:
            out.write_band(1, rescaled)

        out_normalised = os.path.join(TmpDir, os.path.basename(in_file).replace('.tif', '_normalised.tif'))
        with rasterio.open(os.path.normpath(out_normalised), 'w', **out_meta) as out:
            out.write_band(1, rescaled2)

        self.assertAlmostEqual(float(np.nanmax(norm)), 2.0000722408294678, 4)
        self.assertAlmostEqual(float(np.nanmin(norm)), -2.266947031021118, 4)

        self.assertEqual(np.nanmin(rescaled), 0)
        self.assertEqual(np.nanmax(rescaled), 255)

        self.assertEqual(np.nanmin(rescaled2), 0)
        self.assertEqual(np.nanmax(rescaled2), 5)

    def test12_kmeansCluster(self):
        out_img = os.path.join(TmpDir, 'kmeans-cluster_3cluster_3rasters.tif')

        with self.assertRaises(TypeError) as msg:
            outDF = kmeans_clustering(self.gridextract_files + [fileImage] ,out_img)
            self.assertEqual("Pixel Sizes Don't Match - [(0.5, 0.5), (2.5, 2.5)]", str(msg.exception))

        with self.assertRaises(TypeError) as msg:
            outDF = kmeans_clustering(self.gridextract_files + [os.path.realpath(this_dir + '/data/area1_onebox_NDRE_250cm.tif')] ,out_img)
            self.assertEqual("1 raster(s) don't have coordinates systems assigned", str(msg.exception))

        outDF = kmeans_clustering(self.gridextract_files, out_img)

        self.assertTrue(os.path.exists(out_img))
        self.assertTrue(os.path.exists(out_img.replace('.tif', '_statistics.csv')))
        self.assertEqual(3, len(outDF['zone'].unique()))

        with rasterio.open(out_img) as src:
            self.assertEqual(1, src.count)
            self.assertEqual(src.crs.to_string(), '+init=epsg:28354')
            self.assertEqual(0, src.nodata)
            band1 = src.read(1, masked=True)
            self.assertItemsEqual(np.array([0, 1, 2, 3]), np.unique(band1.data))

    def test99_gridExtract(self):
        out_fold = os.path.join(TmpDir, 'gridextract')
        if not os.path.exists(out_fold): os.mkdir(out_fold)
        global rand_gdf, rand_crs
        stats_gdf, stats_crs = extract_pixel_statistics_for_points(rand_gdf, rand_crs, self.gridextract_files, function_list=[np.nanmean, np.nanstd],
                                                                   size_list=[1, 3], output_csvfile=os.path.join(out_fold, 'grid_extract.csv'))

        self.assertEqual(len(stats_gdf.columns), 17)
        self.assertEqual(stats_gdf['std3x3_Band6_250cm'].isna().sum(), 0)
        self.assertEqual(stats_gdf['EPSG'].unique()[0], 28354)
