import platform
import shutil
import tempfile
import unittest
from pathlib import Path

from pyprecag.tests import setup_folder, KEEP_TEST_OUTPUTS
import geopandas as gpd
from pyprecag import convert, crs
from pyprecag.bandops import CalculateIndices, BandMapping
from pyprecag.describe import CsvDescribe, predictCoordinateColumnNames
from pyprecag.processing import *
from pyprecag.raster_ops import rescale, normalise
from pyprecag.kriging_ops import prepare_for_vesper_krige, vesper_text_to_raster, run_vesper, VesperControl

logging.captureWarnings(True)
logging.basicConfig(level=logging.INFO, format="%(message)s")

PY_FILE = os.path.basename(__file__)
THIS_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
TEMP_FOLD = os.path.join(tempfile.gettempdir(), os.path.splitext(PY_FILE)[0])

FILE_CSV = os.path.realpath(os.path.join(THIS_DIR, "area1_yield_ascii_wgs84.csv"))
FILE_BOX = os.path.realpath(os.path.join(THIS_DIR, "area1_onebox_94mga54.shp"))
FILE_BOXES = os.path.realpath(os.path.join(THIS_DIR, "PolyMZ_wgs84_MixedPartFieldsTypes.shp"))

RASTERS_DIR = os.path.realpath(os.path.join(THIS_DIR, "rasters"))
FILE_IMAGE = os.path.realpath(os.path.join(RASTERS_DIR, "area1_rgbi_jan_50cm_84sutm54.tif"))

EPSG = 28354

class TestEnd2End(unittest.TestCase):
    failedTests = []
    gridextract_files = []

    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestEnd2End, cls).setUpClass()

        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            #shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def run(self, result=None):
        unittest.TestCase.run(self, result)  # call superclass run method
        if len(result.failures) > 0 or len(result.errors) > 0:
            self.failedTests.append(self._testMethodName)

    def test01_csvDescribe_ASCII(self):
        csv_desc = CsvDescribe(FILE_CSV)

        self.assertEqual('ascii', csv_desc.file_encoding)
        self.assertEqual(14756, csv_desc.row_count)
        self.assertEqual(24, csv_desc.column_count )
        self.assertEqual(['Lon', 'Lat'], predictCoordinateColumnNames(csv_desc.get_column_names()))
        self.assertTrue(csv_desc.has_column_header)

    def test02_createPolygonFromPointTrail(self):
        out_fold = setup_folder(self.TmpDir, new_folder=self._testMethodName)

        global filePoints, filePoly
        filePoints = os.path.join(out_fold, os.path.splitext(os.path.basename(FILE_CSV))[0] + '_points.shp')
        filePoly = os.path.join(out_fold, os.path.splitext(os.path.basename(FILE_CSV))[0] + '_poly.shp')

        if not os.path.exists(filePoly):
            gdf_points, gdf_pts_crs = convert.convert_csv_to_points(FILE_CSV, filePoints,
                                                                    coord_columns_epsg=4326,
                                                                    out_epsg=EPSG)

            result = create_polygon_from_point_trail(gdf_points, None, filePoly,
                                            thin_dist_m=2.5,
                                            aggregate_dist_m=25,
                                            buffer_dist_m=7,
                                            shrink_dist_m=3)

        self.assertTrue(True, os.path.exists(filePoly))

    def test03_vectorDescribe(self):
        vect_desc = VectorDescribe(filePoly)

        self.assertFalse(vect_desc.is_mz_aware)
        self.assertEqual('Polygon', vect_desc.geometry_type )

        vect_desc = VectorDescribe(filePoints)

        self.assertFalse(vect_desc.is_mz_aware)
        self.assertEqual('Point', vect_desc.geometry_type )

        vect_desc = VectorDescribe(FILE_BOX)

        self.assertFalse(vect_desc.is_mz_aware)
        self.assertEqual('Polygon', vect_desc.geometry_type )

    def test04_blockGrid(self):

        global fileBlockTif, file_block_txt
        fileBlockTif = os.path.join(self.TmpDir, os.path.splitext(os.path.basename(FILE_CSV))[0] + '_block.tif')
        file_block_txt = os.path.join(self.TmpDir, os.path.splitext(os.path.basename(FILE_CSV))[0] + '_block_v.txt')

        vect_desc = VectorDescribe(FILE_BOX)

        if not os.path.exists(fileBlockTif):
            block_grid(in_shapefilename=FILE_BOX,
                       pixel_size=2.5,
                       out_rasterfilename=fileBlockTif,
                       out_vesperfilename=file_block_txt,
                       out_epsg=vect_desc.crs.epsg_number,
                       snap=True,
                       overwrite=True)

        vect_desc = VectorDescribe(FILE_BOX)

        self.assertTrue(os.path.exists(fileBlockTif))
        self.assertTrue(os.path.exists(file_block_txt))

        with rasterio.open(os.path.normpath(fileBlockTif)) as src:
            self.assertEqual(1, src.count, 'Incorrect band count')
            self.assertEqual(95, src.width, 'Incorrect image width')
            self.assertEqual(95, src.height, 'Incorrect image height')
            self.assertEqual((-9999.0,), src.nodatavals, 'Incorrect image nodata value')
            self.assertEqual(('int16',), src.dtypes, 'Incorrect data type')

            print(vect_desc.crs.crs_wkt)
            print('SRS:\t{}'.format(src.crs))
        print('Temp Files in {}'.format(self.TmpDir))

    def test05_cleanTrimPoints(self):
        global fileTrimmed, data_col
        file_ptstoShp = os.path.join(self.TmpDir, os.path.splitext(os.path.basename(FILE_CSV))[0] + 'Pts2shp.shp')

        fileTrimmed = os.path.join(self.TmpDir, os.path.splitext(os.path.basename(FILE_CSV))[0] + '_normtrimmed.csv')
        file_shp = os.path.join(self.TmpDir, os.path.splitext(os.path.basename(FILE_CSV))[0] + '_normtrimmed.shp')
        file_removed = os.path.join(self.TmpDir, os.path.splitext(os.path.basename(FILE_CSV))[0] + '_remove.shp')

        data_col = r'Yield'

        gdf_points, _ = convert.convert_csv_to_points(FILE_CSV,
                                                               out_shapefilename=file_ptstoShp,
                                                               coord_columns_epsg=4326,
                                                               out_epsg=EPSG)

        gdf_out, _ = clean_trim_points(gdf_points, None, data_col, fileTrimmed,
                                             out_keep_shapefile=file_shp,
                                             out_removed_shapefile=file_removed,
                                             boundary_polyfile=FILE_BOX,
                                             thin_dist_m=2.5)

        self.assertTrue(os.path.exists(fileTrimmed))
        self.assertTrue(os.path.exists(file_shp))
        self.assertTrue(os.path.exists(file_removed))

        tmp = gpd.read_file(file_removed)
        self.assertEqual(tmp.crs, gdf_points.crs)

        out_stats = tmp.groupby(by='filter').agg(count=pd.NamedAgg(column='filter', aggfunc='count'))
        tmp = pd.DataFrame.from_records(data=[{'filter': '01 null/missing data', 'count': 1000},
                                              {'filter': '02 Duplicate XY', 'count': 931},
                                              {'filter': '03 clip', 'count': 12211},
                                              {'filter': '04 <= zero', 'count': 62},
                                              {'filter': '05 3 std iter 1', 'count': 3},
                                              {'filter': '06 3 std iter 2', 'count': 4},
                                              {'filter': '07 3 std iter 3', 'count': 1},
                                              {'filter': '09 pointXY (2.5m)', 'count': 2}], index='filter')

        pd.testing.assert_frame_equal(tmp, out_stats)

        self.assertEqual(EPSG, gdf_out.crs.to_epsg())
        self.assertEqual(542, len(gdf_out))
        self.assertIn('nrm_' + data_col, gdf_out.columns)
        self.assertIn('Easting', gdf_out.columns)
        self.assertIn('Northing', gdf_out.columns)
        self.assertIn('EN_EPSG', gdf_out.columns)

    @unittest.skipIf(platform.system() != 'Windows', 'Vesper only present on Windows')
    def test06_prepareForVesperKrig(self):
        csv_desc = CsvDescribe(fileTrimmed)
        df_csv = csv_desc.open_pandas_dataframe()

        global file_control
        sub_file = os.path.splitext(os.path.basename(FILE_CSV))[0]
        file_control = sub_file + '_control_' + data_col + '.txt'

        vc = VesperControl()
        vc.update(xside=30, yside=30)

        if not os.path.exists(file_control):
            bat_file, file_control = prepare_for_vesper_krige(df_csv, data_col, file_block_txt,
                                                              self.TmpDir,
                                                              control_textfile=file_control,
                                                              coord_columns=[], epsg=EPSG,
                                                              control_options=vc)

            self.assertTrue(os.path.exists(bat_file))
            self.assertTrue(os.path.exists(file_control))

        self.assertTrue(os.path.exists(os.path.join(self.TmpDir, r'Vesper',
                                                    sub_file + '_vesperdata_' + data_col + '.csv')))
        df_csv = pd.read_csv(os.path.join(self.TmpDir, r'Vesper',
                                          sub_file + '_vesperdata_' + data_col + '.csv'))

        x_column, y_column = predictCoordinateColumnNames(df_csv.columns)
        self.assertEqual('EASTING', x_column.upper())
        self.assertEqual('NORTHING', y_column.upper())

        print('Running Vesper, Please wait....')
        run_vesper(file_control)

    @unittest.skipIf(platform.system() != 'Windows','VESPER only present on Windows')
    def test07_vesperTextToRaster(self):
        global out_predtif
        out_predtif, out_setif, out_citxt = vesper_text_to_raster(file_control, EPSG)
        for eaFile in [out_predtif, out_setif, out_citxt]:
            self.assertTrue(os.path.exists(eaFile))

        with open(out_citxt) as test_file:
            test_lines = test_file.readlines()[:2]
            SE = test_lines[0].split(':')[-1].strip()
            CI95 = test_lines[1].split(':')[-1].strip()

        with rasterio.open(os.path.normpath(out_predtif)) as src:
            self.assertEqual(1, src.count, 'Incorrect band count')
            self.assertEqual(94, src.width, 'Incorrect image width')
            self.assertEqual(93, src.height, 'Incorrect image height')
            self.assertEqual((-9999.0,), src.nodatavals, 'Incorrect image nodata value')
            self.assertEqual(('float32',), src.dtypes, 'Incorrect data type')

            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), src.crs, 'Incorrect coordinate system')
            self.assertTrue(src.tags()['PAT_MedianPredSE'])
            self.assertTrue(src.tags()['PAT_95ConfLevel'])
            self.assertEqual(src.tags()['PAT_MedianPredSE'], SE)
            self.assertEqual(src.tags()['PAT_95ConfLevel'], CI95)

        with rasterio.open(os.path.normpath(out_setif)) as src:
            self.assertEqual(1, src.count, 'Incorrect band count')
            self.assertEqual(94, src.width, 'Incorrect image width')
            self.assertEqual(93, src.height, 'Incorrect image height')
            self.assertEqual((-9999.0,), src.nodatavals, 'Incorrect image nodata value')
            self.assertEqual(('float32',), src.dtypes, 'Incorrect data type')
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), src.crs, 'Incorrect coordinate system')

    def test08_randomPixelSelection(self):
        global rand_gdf, rand_crs
        out_randompts = os.path.join(self.TmpDir,
                                     os.path.basename(fileBlockTif).replace('.tif', '_randpts.shp'))
        rast_crs = crs.getCRSfromRasterFile(fileBlockTif)

        with rasterio.open(os.path.normpath(fileBlockTif)) as raster:
            rand_gdf, rand_crs = random_pixel_selection(raster, rast_crs,
                                                        50, out_shapefile=out_randompts)
        self.assertEqual(len(rand_gdf), 50)
        self.assertTrue(os.path.exists(out_randompts))
        self.assertEqual(rast_crs.epsg, rand_gdf.crs)

    def test09_calcImageIndices_allopts(self):

        out_fold = setup_folder(self.TmpDir, new_folder=self._testMethodName)
        bm = BandMapping(green=2, infrared=4, rededge=1, mask=5)
        indices = CalculateIndices(**bm).valid_indices()
        files = calc_indices_for_block(FILE_IMAGE, 2.5, bm, out_fold, indices, image_nodata=0,
                                       image_epsg=32754, polygon_shapefile=FILE_BOX, out_epsg=28354)

        self.gridextract_files += files

        self.assertEqual(len(files), 3)

        with rasterio.open(files[0]) as src:
            self.assertEqual(src.nodata, -9999)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), src.crs, 'Incorrect coordinate system ')
            self.assertEqual(src.meta['dtype'], 'float32')
            self.assertEqual(src.res, (2.5, 2.5))
            self.assertEqual(1, src.count, 'Incorrect band count')

    def test10_resampleBands2Block_allopts(self):
        out_fold = setup_folder(self.TmpDir, new_folder=self._testMethodName)

        files = resample_bands_to_block(FILE_IMAGE, 2.5, out_fold, band_nums=[6],
                                        image_nodata=0, image_epsg=32754,
                                        polygon_shapefile=FILE_BOX, out_epsg=28354)

        self.gridextract_files += files
        self.assertEqual(len(files), 1)

        with rasterio.open(files[0]) as src:
            self.assertEqual(src.nodata, 0.0)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), src.crs, 'Incorrect coordinate system ')
            self.assertEqual(src.meta['dtype'], 'float32')
            self.assertEqual(src.res, (2.5, 2.5))
            self.assertEqual(1, src.count, 'Incorrect band count')

    def test11_rescaleNormaliseRaster(self):
        out_fold = setup_folder(self.TmpDir, new_folder=self._testMethodName)
        in_file = self.gridextract_files[-1]

        rast_crs = crs.getCRSfromRasterFile(in_file)

        with rasterio.open(os.path.normpath(in_file)) as src:
            rasterio.shutil.copy(src,os.path.join(out_fold,"INPUT_Band6_250cm.tif"))

            rescaled0_255 = rescale(src, 0, 255)
            rescaled0_5 = rescale(src, 0, 5)
            norm = normalise(src)
            out_meta = src.meta.copy()

        out_meta['crs'] = rast_crs.crs_wkt
        out_meta['count'] = 1  # contains only one band
        out_meta['dtype'] = np.float32

        out_rescale255 = os.path.join(out_fold, os.path.basename(in_file).replace('.tif', '_rescale0-255.tif'))
        with rasterio.open(os.path.normpath(out_rescale255), 'w', **out_meta) as out:
            out.write_band(1, rescaled0_255)

        self.assertEqual(0, np.nanmin(rescaled0_255),
                         'Raster Min for ' + str(Path(out_rescale255).relative_to(TEMP_FOLD)))
        self.assertEqual(255, np.nanmax(rescaled0_255),
                         'Raster Max for ' + str(Path(out_rescale255).relative_to(TEMP_FOLD)))

        out_rescale5 = os.path.join(out_fold, os.path.basename(in_file).replace('.tif', '_rescale0-5.tif'))
        with rasterio.open(os.path.normpath(out_rescale5), 'w', **out_meta) as out:
            out.write_band(1, rescaled0_5)

        self.assertEqual(0, np.nanmin(rescaled0_5),
                         'Raster Min for ' + str(Path(out_rescale5).relative_to(TEMP_FOLD)))
        self.assertEqual(5, np.nanmax(rescaled0_5),
                         'Raster Max for ' + str(Path(out_rescale5).relative_to(TEMP_FOLD)))

        out_normalised = os.path.join(out_fold, os.path.basename(in_file).replace('.tif', '_normalised.tif'))
        with rasterio.open(os.path.normpath(out_normalised), 'w', **out_meta) as out:
            out.write_band(1, norm)

        self.assertAlmostEqual(-2.245925, float(np.nanmin(norm)), 2,
                               'Raster Min for ' + str(Path(out_normalised).relative_to(TEMP_FOLD)))

        self.assertAlmostEqual(2.000072, float(np.nanmax(norm)), 2,
                               'Raster Max for ' + str(Path(out_normalised).relative_to(TEMP_FOLD)))


    def test12_kmeansCluster(self):
        out_img = os.path.join(setup_folder(self.TmpDir, new_folder=self._testMethodName),
                               'kmeans-cluster_3cluster_3rasters.tif')

        with self.assertRaises(TypeError) as msg:
            _ = kmeans_clustering(self.gridextract_files + [FILE_IMAGE], out_img)
        self.assertIn("raster_files are of different pixel sizes", str(msg.exception))

        with self.assertRaises(TypeError) as msg:
            _ = kmeans_clustering(self.gridextract_files +
                                  [os.path.realpath(os.path.join(RASTERS_DIR, 'area1_onebox_NDRE_250cm.tif'))],
                                  out_img)

        self.assertIn("1 raster(s) don't have coordinates systems assigned", str(msg.exception))

        out_df = kmeans_clustering(self.gridextract_files, out_img)

        self.assertTrue(os.path.exists(out_img))
        self.assertTrue(os.path.exists(out_img.replace('.tif', '_statistics.csv')))
        self.assertEqual(5, len(out_df['zone'].unique()))

        with rasterio.open(out_img) as src:
            self.assertEqual(1, src.count)
            self.assertEqual(src.crs.to_epsg(), 28354)
            self.assertEqual(0, src.nodata)
            band1 = src.read(1, masked=True)

            six.assertCountEqual(self, np.array([0, 1, 2, 3]), np.unique(band1.data))

    def test99_gridExtract(self):
        out_fold = setup_folder(self.TmpDir, new_folder=self._testMethodName)

        global rand_gdf, rand_crs
        stats_gdf, stats_crs = extract_pixel_statistics_for_points(
            rand_gdf, rand_crs, self.gridextract_files, function_list=[np.nanmean, np.nanstd],
            size_list=[1, 3], output_csvfile=os.path.join(out_fold, 'grid_extract.csv'))

        self.assertEqual(len(stats_gdf.columns), 17)
        self.assertEqual(stats_gdf['std3x3_Band6_250cm'].isna().sum(), 0)
        self.assertEqual(stats_gdf['EPSG'].unique()[0], 28354)
