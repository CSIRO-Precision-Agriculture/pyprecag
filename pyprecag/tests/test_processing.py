import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
import geopandas as gpd
from pandas import testing as pdts

from pyprecag.tests import make_dummy_tif_files, setup_folder, KEEP_TEST_OUTPUTS, warn_with_traceback

from pyprecag.bandops import BandMapping, CalculateIndices
from pyprecag.crs import crs
from pyprecag import convert, raster_ops
from pyprecag.convert import add_point_geometry_to_dataframe

from pyprecag.processing import *

PY_FILE = os.path.basename(__file__)
TEMP_FOLD = os.path.join(tempfile.gettempdir(), os.path.splitext(PY_FILE)[0])
THIS_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')

logging.captureWarnings(True)
logging.basicConfig(level=logging.INFO, format="%(message)s")

config.set_debug_mode(False)


class Test_BlockGrid(unittest.TestCase):
    failedTests = []

    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(Test_BlockGrid, cls).setUpClass()

        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD, new_folder=__class__.__name__)
        cls.singletif, cls.multitif = make_dummy_tif_files(cls.TmpDir)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def run(self, result=None):
        unittest.TestCase.run(self, result)  # call superclass run method
        if self.id() in result.failed_tests or len(result.errors) > 0:
            self.failedTests.append(self._testMethodName)
        else:
            if os.path.exists(self.test_outdir) and not KEEP_TEST_OUTPUTS:
                shutil.rmtree(self.test_outdir)

    def test_BlockGrid(self):
        poly = os.path.realpath(os.path.join(THIS_DIR, "area2_onebox_94mga54.shp"))

        file_sub_name = os.path.join(self.test_outdir, os.path.splitext(os.path.basename(poly))[0])
        vect_desc = VectorDescribe(poly)
        output_files = block_grid(in_shapefilename=poly,
                                  pixel_size=5,
                                  out_rasterfilename=file_sub_name + '_block.tif',
                                  out_vesperfilename=file_sub_name + '_block_v.txt',
                                  out_epsg=vect_desc.crs.epsg_number,
                                  snap=True,
                                  overwrite=True)

        self.assertTrue(os.path.exists(file_sub_name + '_block.tif'))
        self.assertTrue(os.path.exists(file_sub_name + '_block_v.txt', ))
        self.assertTrue(len(output_files), 1)

        with rasterio.open(os.path.normpath(file_sub_name + '_block.tif')) as src:
            self.assertEqual(1, src.count, 'Incorrect band count')
            self.assertEqual(87, src.width, 'Incorrect image width')
            self.assertEqual(66, src.height, 'Incorrect image height')
            self.assertEqual((-9999.0,), src.nodatavals, 'Incorrect image nodata value')
            self.assertEqual(('int16',), src.dtypes, 'Incorrect data type')
            self.assertEqual(28354, src.crs.to_epsg(), 'Incorrect EPSG')

    def test_BlockGrid_GrpBy(self):
        poly = os.path.realpath(os.path.join(THIS_DIR, "PolyMZ_wgs84_MixedPartFieldsTypes.shp"))

        file_sub_name = os.path.join(self.test_outdir, os.path.splitext(os.path.basename(poly))[0])

        output_files = block_grid(in_shapefilename=poly,
                                  pixel_size=5,
                                  out_rasterfilename=file_sub_name + '_block.tif',
                                  out_vesperfilename=file_sub_name + '_block_v.txt',
                                  out_epsg=28354,
                                  groupby='part_type',
                                  snap=True,
                                  overwrite=True)

        self.assertTrue(2, len(output_files))
        self.assertTrue(os.path.exists(file_sub_name + '_block_MultiPart.tif'))
        self.assertTrue(os.path.exists(file_sub_name + '_block_SinglePart_v.txt', ))

        with rasterio.open(os.path.normpath(file_sub_name + '_block_MultiPart.tif')) as src:
            self.assertEqual(1, src.count, 'Incorrect band count')
            self.assertEqual(99, src.width, 'Incorrect image width')
            self.assertEqual(68, src.height, 'Incorrect image height')
            self.assertEqual((-9999.0,), src.nodatavals, 'Incorrect image nodata value')
            self.assertEqual(('int16',), src.dtypes, 'Incorrect data type')
            self.assertEqual(28354, src.crs.to_epsg(), 'Incorrect EPSG')


class Test_CleanTrim(unittest.TestCase):
    failedTests = []

    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(Test_CleanTrim, cls).setUpClass()

        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD, new_folder=__class__.__name__)
        cls.singletif, cls.multitif = make_dummy_tif_files(cls.TmpDir)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def run(self, result=None):
        unittest.TestCase.run(self, result)  # call superclass run method
        if self.id() in result.failed_tests or len(result.errors) > 0:
            self.failedTests.append(self._testMethodName)
        else:
            if os.path.exists(self.test_outdir) and not KEEP_TEST_OUTPUTS:
                shutil.rmtree(self.test_outdir)

    def test_cleanTrimPoints_area1(self):
        in_csv = os.path.join(THIS_DIR, "area1_yield_ascii_wgs84.csv")
        in_poly = os.path.join(THIS_DIR, "area1_onebox_94mga54.shp")
        out_csv = os.path.join(self.test_outdir, os.path.basename(in_csv))
        out_shp = os.path.join(self.test_outdir, os.path.basename(in_csv).replace('.csv', '.shp'))
        out_rm_shp = os.path.join(self.test_outdir, os.path.basename(in_csv).replace('.csv', '_remove.shp'))

        gdf_points, gdf_pts_crs = convert.convert_csv_to_points(in_csv, coord_columns_epsg=4326,
                                                                out_epsg=28354)
        out_gdf, out_crs = clean_trim_points(gdf_points, None, 'Yield',
                                             out_csv, out_keep_shapefile=out_shp,
                                             out_removed_shapefile=out_rm_shp,
                                             boundary_polyfile=in_poly, thin_dist_m=2.5)

        self.assertIsInstance(out_gdf, GeoDataFrame)
        self.assertTrue(os.path.exists(out_csv))

        self.assertEqual(28354, out_gdf.crs.to_epsg())  # {'init': 'EPSG:28354', 'no_defs': True})
        self.assertIn('EN_EPSG', out_gdf.columns)
        self.assertEqual(542, len(out_gdf))

        tmp = pd.DataFrame.from_records(data=[{'filter': '01 null/missing data', 'count': 1000},
                                              {'filter': '02 Duplicate XY', 'count': 931},
                                              {'filter': '03 clip', 'count': 12211},
                                              {'filter': '04 <= zero', 'count': 62},
                                              {'filter': '05 3 std iter 1', 'count': 3},
                                              {'filter': '06 3 std iter 2', 'count': 4},
                                              {'filter': '07 3 std iter 3', 'count': 1},
                                              {'filter': '09 pointXY (2.5m)', 'count': 2}], index='filter')

        res_df = gpd.read_file(out_rm_shp)

        self.assertEqual(14214, len(res_df))

        res_stats = res_df.value_counts('filter', sort=False).to_frame('count')
        pdts.assert_frame_equal(tmp, res_stats)

    def test_cleanTrimPoints_area2(self):
        in_csv = os.path.join(THIS_DIR, "area2_yield_ISO-8859-1.csv")
        in_poly = os.path.join(THIS_DIR, "area2_onebox_94mga54.shp")
        out_csv = os.path.join(self.test_outdir, os.path.basename(in_csv))
        out_shp = os.path.join(self.test_outdir, os.path.basename(in_csv).replace('.csv', '.shp'))
        out_rm_shp = os.path.join(self.test_outdir, os.path.basename(in_csv).replace('.csv', '_remove.shp'))

        gdf_points, gdf_pts_crs = convert.convert_csv_to_points(in_csv, coord_columns_epsg=4326, out_epsg=28354)
        out_gdf, out_crs = clean_trim_points(gdf_points, gdf_pts_crs, 'Yld Mass(Dry)(tonne/ha)',
                                             out_csv, out_keep_shapefile=out_shp,
                                             out_removed_shapefile=out_rm_shp,
                                             boundary_polyfile=in_poly, thin_dist_m=2.5)

        self.assertIsInstance(out_gdf, GeoDataFrame)
        self.assertTrue(os.path.exists(out_csv))

        self.assertEqual(28354, out_gdf.crs.to_epsg())
        self.assertEqual(554, len(out_gdf))
        self.assertIn('EN_EPSG', out_gdf.columns)

        res_df = gpd.read_file(out_rm_shp)
        self.assertEqual(989, len(res_df))

        tmp = pd.DataFrame.from_records(data=[{'filter': '01 clip', 'count': 951},
                                              {'filter': '03 pointXY (2.5m)', 'count': 38}],
                                        index='filter')

        res_stats = res_df.value_counts('filter', sort=False).to_frame('count')
        pdts.assert_frame_equal(tmp, res_stats)


class Test_Processing(unittest.TestCase):
    failedTests = []

    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(Test_Processing, cls).setUpClass()

        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD, new_folder=__class__.__name__)
        cls.singletif, cls.multitif = make_dummy_tif_files(cls.TmpDir)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def run(self, result=None):
        unittest.TestCase.run(self, result)  # call superclass run method
        if self.id() in result.failed_tests or len(result.errors) > 0:
            self.failedTests.append(self._testMethodName)
        else:
            if os.path.exists(self.test_outdir) and not KEEP_TEST_OUTPUTS:
                shutil.rmtree(self.test_outdir)

    def test_createPolygonFromPointTrail(self):
        in_csv = os.path.join(THIS_DIR, "area2_yield_ISO-8859-1.csv")

        out_polyfile = os.path.join(self.test_outdir, os.path.splitext(os.path.basename(in_csv))[0] + '_poly.shp')

        gdf_points, gdf_pts_crs = convert.convert_csv_to_points(in_csv, None,
                                                                coord_columns_epsg=4326, out_epsg=28354)

        create_polygon_from_point_trail(gdf_points, None, out_polyfile,
                                        thin_dist_m=2.5,
                                        aggregate_dist_m=25,
                                        buffer_dist_m=10,
                                        shrink_dist_m=3)

        self.assertTrue(os.path.exists(out_polyfile))

        vect_desc = VectorDescribe(out_polyfile)
        self.assertEqual(28354, vect_desc.crs.epsg_number)
        self.assertFalse(vect_desc.is_mz_aware)
        self.assertEqual('Polygon', vect_desc.geometry_type)
        self.assertEqual(28354, result.crs.to_epsg())

    def test_randomPixelSelection(self):
        raster_file = self.singletif
        out_shp = os.path.join(self.test_outdir, os.path.basename(raster_file).replace('.tif', '_randpts.shp'))
        rast_crs = pyprecag_crs.getCRSfromRasterFile(raster_file)

        with rasterio.open(os.path.normpath(raster_file)) as raster:
            rand_gdf, rand_crs = random_pixel_selection(raster, rast_crs, 50, out_shapefile=out_shp)

        self.assertEqual(len(rand_gdf), 50)
        self.assertTrue(os.path.exists(out_shp))
        self.assertEqual(rand_crs, rast_crs)

    def test_PersistorAllYears(self):
        raster_files = glob.glob(os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'Year*.tif')))
        raster_files += [os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))]

        out_img = os.path.join(self.test_outdir, 'persistor_allyears.tif')

        # check for raised pixel size error
        with self.assertRaises(TypeError) as msg:
            persistor_all_years(raster_files, out_img, True, 10)
        self.assertIn("raster_files are of different pixel sizes", str(msg.exception))
        raster_files.remove(os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif')))

        persistor_all_years(raster_files, out_img, True, 10)
        self.assertTrue(os.path.exists(out_img))

        src_img = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'persistor_allyears.tif'))
        with rasterio.open(src_img) as src, rasterio.open(os.path.normpath(out_img)) as test:
            self.assertEqual(src.profile, test.profile)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), test.crs)

            np.testing.assert_array_equal(src.read(), test.read())

            band1 = test.read(1, masked=True)
            six.assertCountEqual(self, [-9999, 0, 1, 2, 3], np.unique(band1.data).tolist())

    def test_PersistorTargetProb(self):
        raster_files = glob.glob(os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'Year*.tif')))

        out_img = os.path.join(self.test_outdir, 'persistor_targetprob.tif')

        persistor_target_probability(raster_files, 10, 75,
                                     raster_files, -10, 75, out_img)

        self.assertTrue(os.path.exists(out_img))
        src_img = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'persistor_targetprob.tif'))
        with rasterio.open(src_img) as src, rasterio.open(os.path.normpath(out_img)) as test:
            self.assertEqual(src.profile, test.profile)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), test.crs)

            np.testing.assert_array_equal(src.read(), test.read())

            band1 = test.read(1, masked=True)
            six.assertCountEqual(self, [-9999, -1, 0, 1], np.unique(band1.data).tolist())


class TestKMeansCluster(unittest.TestCase):
    failedTests = []

    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestKMeansCluster, cls).setUpClass()

        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD, new_folder=__class__.__name__)

        cls.singletif, cls.multitif = make_dummy_tif_files(cls.TmpDir)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def run(self, result=None):

        unittest.TestCase.run(self, result)  # call superclass run method
        if self.id() in result.failed_tests or len(result.errors) > 0:
            # self.failedTests += [ea.split('.')[-1] for ea in result.failed_tests]
            self.failedTests.append(self._testMethodName)
        else:
            if os.path.exists(self.test_outdir) and not KEEP_TEST_OUTPUTS:
                shutil.rmtree(self.test_outdir)

    def test_kmeansCluster(self):
        raster_files = glob.glob(os.path.realpath(os.path.join(THIS_DIR, 'rasters', '*one*.tif')))
        raster_files += [os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))]

        out_img = os.path.join(self.test_outdir, 'kmeans-cluster_3cluster_3rasters.tif')

        with self.assertRaises(TypeError) as msg:
            _ = kmeans_clustering(raster_files, out_img)

        self.assertIn('raster_files are of different pixel sizes', str(msg.exception))
        raster_files.remove(os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif')))

        with self.assertRaises(TypeError) as msg:
            _ = kmeans_clustering(raster_files, out_img)

        self.assertIn("1 raster(s) don't have coordinates systems assigned", str(msg.exception))
        raster_files.remove(os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_onebox_NDRE_250cm.tif')))

        out_df = kmeans_clustering(raster_files, out_img)

        self.assertTrue(os.path.exists(out_img))
        self.assertTrue(os.path.exists(out_img.replace('.tif', '_statistics.csv')))
        self.assertEqual(5, len(out_df['zone'].unique()))

        with rasterio.open(out_img) as src:
            self.assertEqual(1, src.count)
            if hasattr(src.crs, 'to_proj4'):
                self.assertEqual(src.crs.to_proj4().lower(), '+init=epsg:28354')
            else:
                self.assertEqual(src.crs.to_string().lower(), '+init=epsg:28354')
            self.assertEqual(0, src.nodata)
            band1 = src.read(1, masked=True)
            six.assertCountEqual(self, np.array([0, 1, 2, 3]), np.unique(band1.data))

    def test_kmeansCluster_alias(self):
        raster_files = glob.glob(os.path.realpath(os.path.join(THIS_DIR, 'rasters', '*one*.tif')))
        raster_files += [os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))]

        raster_files = {ea: os.path.basename(ea).split('_')[2] for ea in raster_files}

        out_img = os.path.join(self.test_outdir, 'kmeans-cluster_3cluster_3rasters.tif')

        with self.assertRaises(TypeError) as msg:
            _ = kmeans_clustering(raster_files, out_img)

        self.assertIn('raster_files are of different pixel sizes', str(msg.exception))
        del raster_files[os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))]

        with self.assertRaises(TypeError) as msg:
            _ = kmeans_clustering(raster_files, out_img)

        self.assertIn("1 raster(s) don't have coordinates systems assigned", str(msg.exception))
        del raster_files[os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_onebox_NDRE_250cm.tif'))]

        out_df = kmeans_clustering(raster_files, out_img)

        out_alias = [ea.split(',')[0] for ea in out_df.columns if ',' in ea]

        # compare contents of lists disregarding order use assertCountEqual
        self.assertCountEqual(raster_files.values(), set(out_alias), 'Incorrect Alias')

        self.assertTrue(os.path.exists(out_img))
        self.assertTrue(os.path.exists(out_img.replace('.tif', '_statistics.csv')))
        self.assertEqual(5, len(out_df['zone'].unique()))

        with rasterio.open(out_img) as src:
            self.assertEqual(1, src.count)
            if hasattr(src.crs, 'to_proj4'):
                self.assertEqual(src.crs.to_proj4().lower(), '+init=epsg:28354')
            else:
                self.assertEqual(src.crs.to_string().lower(), '+init=epsg:28354')
            self.assertEqual(0, src.nodata)
            band1 = src.read(1, masked=True)
            six.assertCountEqual(self, np.array([0, 1, 2, 3]), np.unique(band1.data))


class TestStripTrials(unittest.TestCase):
    failedTests = []

    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestStripTrials, cls).setUpClass()

        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD, new_folder=__class__.__name__)

        cls.singletif, cls.multitif = make_dummy_tif_files(cls.TmpDir)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def run(self, result=None):

        unittest.TestCase.run(self, result)  # call superclass run method
        if self.id() in result.failed_tests or len(result.errors) > 0:
            # self.failedTests += [ea.split('.')[-1] for ea in result.failed_tests]
            self.failedTests.append(self._testMethodName)
        else:
            if os.path.exists(self.test_outdir) and not KEEP_TEST_OUTPUTS:
                shutil.rmtree(self.test_outdir)

    def test_CreateStripTreatmentPoints(self):
        import pandas
        import geopandas
        print('Pandas Version is ', pandas.__version__)
        print('GeoPandas Version is ', geopandas.__version__)

        out_points_file = os.path.join(self.test_outdir, 'testCreateStripTreatment_Points.shp')
        out_lines_file = os.path.join(self.test_outdir, 'testCreateStripTreatment_Lines.shp')
        data = [(10, LineString([(740800, 6169700, 5), (741269, 6169700, 5)])),
                (20, LineString([(741000, 6169800, 6), (741003, 6170012, 6), (741003.5, 6170012.5, 6)])),
                (30, LineString([(741308, 6169916, 7), (741250, 6169950, 8), (741155, 6169974, 8), (741100, 6170000, 8)])),
                (30, LineString([(741413, 6169853, 7), (741372, 6169874, 7), (741345, 6169899, 7), (741308, 6169916, 7)])),
                (50, LineString([(740800, 6169912, 8), (740900, 6170094, 8)]))]

        in_lines_gdf = gpd.GeoDataFrame(pd.DataFrame.from_records(data, columns=['TrialID', 'geometry']),
                                        geometry='geometry', crs=28354)

        if config.get_debug_mode():
            in_lines_gdf['length'] = in_lines_gdf.length
            save_geopandas_tofile(in_lines_gdf, out_lines_file.replace('_Lines', '_Input_Lines'), overwrite=True)

        gdf_lines_crs = crs()
        gdf_lines_crs.getFromEPSG(28354)

        out_points_gdf, out_crs, out_lines_gdf = \
            create_points_along_line(in_lines_gdf, gdf_lines_crs, 7, 25,
                                     out_points_shapefile=out_points_file,
                                     out_lines_shapefile=out_lines_file)

        # Test Points Output-----------------------------------------------------------
        self.assertIsInstance(out_points_gdf, GeoDataFrame)

        stats = out_points_gdf.groupby(by='TrialID', dropna=False).agg(count=pd.NamedAgg(column='TrialID', aggfunc='count'))
        pd.testing.assert_frame_equal(pd.DataFrame.from_records([{'TrialID': 0, 'count': 198},
                                                                 {'TrialID': 1, 'count': 87},
                                                                 {'TrialID': 2, 'count': 90},
                                                                 {'TrialID': 3, 'count': 147}], index='TrialID'), stats)

        self.assertEqual(28354, out_crs.epsg_number)
        self.assertTrue(os.path.exists(out_points_file))
        self.assertEqual({'Point'}, set(out_points_gdf.geom_type))
        self.assertEqual(['E Strip', 'N Strip', 'NE Strip', 'NW Strip', 'S Strip', 'SE Strip',
                          'SW Strip', 'Strip', 'W Strip'],
                         sorted(list(out_points_gdf['Strip_Name'].unique())))

        self.assertEqual(out_points_gdf.crs, out_lines_gdf.crs)

        # Test Line Output-----------------------------------------------------------
        self.assertIsInstance(out_lines_gdf, GeoDataFrame)
        self.assertEqual(12, len(out_lines_gdf))

        self.assertTrue(os.path.exists(out_lines_file))
        self.assertEqual({'LineString'}, set(out_lines_gdf.geom_type))
        self.assertEqual(['E Strip', 'N Strip', 'NE Strip', 'NW Strip', 'S Strip', 'SE Strip',
                          'SW Strip', 'Strip', 'W Strip'],
                         sorted(list(out_lines_gdf['Strip_Name'].unique())))

        del out_points_gdf, out_lines_gdf, out_crs

        lines_crs = pyprecag_crs.crs()
        lines_crs.getFromEPSG(28353)
        lines = GeoSeries([LineString(((740811.9268002876, 6169890.435495311),
                                       (740949.7570626185, 6170048.754039882),
                                       (741145.3270294399, 6170037.578613207),
                                       (741309.2332873486, 6169946.312628689),
                                       (741318.5461429111, 6169847.596359722)))])
        in_lines_gdf = GeoDataFrame({'geometry': lines, 'block': ['test']}, crs=28353)

        # without saving to file on a single multi vertex line
        out_points_gdf, out_crs, out_lines_gdf = create_points_along_line(in_lines_gdf,
                                                                          lines_crs, 31, 25)

        self.assertIsInstance(out_points_gdf, GeoDataFrame)
        self.assertEqual(66, len(out_points_gdf), "Row count doesn't match")

        self.assertEqual(lines_crs.epsg_number, out_crs.epsg_number)
        self.assertEqual(out_points_gdf.crs, out_lines_gdf.crs)

        self.assertEqual({'Point'}, set(out_points_gdf.geom_type))
        self.assertEqual(['N Strip', 'S Strip', 'Strip'],
                         sorted(list(out_points_gdf['Strip_Name'].unique())))

        self.assertEqual(3, len(out_lines_gdf))
        self.assertEqual(['N Strip', 'S Strip', 'Strip'],
                         sorted(list(out_lines_gdf['Strip_Name'].unique())))

    def test_TTestAnalysis(self):
        columns = ['FID', 'TrialID', 'Strip_Name', 'PointID', 'DistOnLine', 'geometry']
        data = [[30, 1, "Strip", 1, 2.5, Point(300169, 6181357)],
                [31, 1, "E Strip", 1, 2.5, Point(300189, 6181363)],
                [32, 1, "W Strip", 1, 2.5, Point(300150, 6181352)],
                [33, 1, "Strip", 2, 22.5, Point(300174, 6181338)],
                [34, 1, "E Strip", 2, 22.5, Point(300194, 6181343)],
                [35, 1, "W Strip", 2, 22.5, Point(300155, 6181333)],
                [36, 1, "Strip", 3, 42.5, Point(300180, 6181319)],
                [37, 1, "E Strip", 3, 42.5, Point(300199, 6181324)],
                [38, 1, "W Strip", 3, 42.5, Point(300160, 6181314)],
                [39, 1, "Strip", 4, 62.5, Point(300185, 6181299)],
                [40, 1, "E Strip", 4, 62.5, Point(300204, 6181305)],
                [41, 1, "W Strip", 4, 62.5, Point(300166, 6181294)],
                [42, 1, "Strip", 5, 82.5, Point(300190, 6181280)],
                [43, 1, "E Strip", 5, 82.5, Point(300209, 6181285)],
                [44, 1, "W Strip", 5, 82.5, Point(300171, 6181275)],
                [45, 1, "Strip", 6, 102.5, Point(300195, 6181261)],
                [46, 1, "E Strip", 6, 102.5, Point(300215, 6181266)],
                [47, 1, "W Strip", 6, 102.5, Point(300176, 6181256)],
                [48, 1, "Strip", 7, 122.5, Point(300200, 6181242)],
                [49, 1, "E Strip", 7, 122.5, Point(300220, 6181247)],
                [50, 1, "W Strip", 7, 122.5, Point(300181, 6181236)],
                [99, 3, "Strip", 4, 67.6, Point(300937, 6181897)],
                [100, 3, "SE Strip", 4, 67.6, Point(300949, 6181881)],
                [102, 3, "Strip", 5, 87.6, Point(300921, 6181884)],
                [103, 3, "SE Strip", 5, 87.6, Point(300933, 6181869)],
                [105, 3, "Strip", 6, 107.6, Point(300905, 6181872)],
                [106, 3, "SE Strip", 6, 107.6, Point(300917, 6181856)],
                [108, 3, "Strip", 7, 127.6, Point(300889, 6181860)],
                [109, 3, "SE Strip", 7, 127.6, Point(300901, 6181844)],
                [111, 3, "Strip", 8, 147.6, Point(300873, 6181848)],
                [112, 3, "SE Strip", 8, 147.6, Point(300886, 6181832)],
                [114, 3, "Strip", 9, 167.6, Point(300858, 6181836)],
                [115, 3, "SE Strip", 9, 167.6, Point(300870, 6181820)],
                [117, 3, "Strip", 10, 187.6, Point(300842, 6181823)],
                [118, 3, "SE Strip", 10, 187.6, Point(300854, 6181807)],
                [120, 3, "Strip", 11, 207.6, Point(300826, 6181811)],
                [121, 3, "SE Strip", 11, 207.6, Point(300838, 6181795)],
                [123, 3, "Strip", 12, 227.6, Point(300810, 6181799)],
                [124, 3, "SE Strip", 12, 227.6, Point(300822, 6181783)],
                [126, 3, "Strip", 13, 247.6, Point(300794, 6181787)],
                [127, 3, "SE Strip", 13, 247.6, Point(300806, 6181771)]]

        gdf_points = GeoDataFrame(data, columns=columns, geometry='geometry', crs=28354)
        crs_points = crs()
        crs_points.getFromEPSG(28354)

        values_rast = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'Yield_withStrip_PRED_2m.tif'))
        control_raster = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'Yield_withoutStrip_PRED_2m.tif'))
        zone_raster = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'k-means_3clusters_3rasters_2m.tif'))

        result = ttest_analysis(gdf_points, crs_points,
                                values_rast, self.test_outdir,
                                zone_raster=zone_raster,
                                control_raster=control_raster)

        self.assertEqual(13, len(result.columns))
        file_csv = glob.glob(os.path.join(self.test_outdir, 'Trial*.csv'))
        file_graph = glob.glob(os.path.join(self.test_outdir, 'Trial*graph.png'))
        file_map = glob.glob(os.path.join(self.test_outdir, 'Trial*map.png'))

        self.assertEqual(4, len(file_csv), "CSV File count doesn't match")
        self.assertEqual(4, len(file_graph), "Graph File count doesn't match")
        self.assertEqual(2, len(file_map), "Map File count doesn't match")


class TestExtractRasterStatisticsForPoints(unittest.TestCase):
    failedTests = []

    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestExtractRasterStatisticsForPoints, cls).setUpClass()

        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD, new_folder=__class__.__name__)

        cls.singletif, cls.multitif = make_dummy_tif_files(cls.TmpDir)

        df = pd.read_csv(os.path.realpath(os.path.join(THIS_DIR, 'area1_dummypoints_94mga54.csv')))
        cls.gdf_points, cls.crs_points = add_point_geometry_to_dataframe(df, ['X', 'Y'], 28354)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def run(self, result=None):
        unittest.TestCase.run(self, result)  # call superclass run method
        if self.id() in result.failed_tests or len(result.errors) > 0:
            # self.failedTests += [ea.split('.')[-1] for ea in result.failed_tests]
            self.failedTests.append(self._testMethodName)
        else:
            if os.path.exists(self.test_outdir) and not KEEP_TEST_OUTPUTS:
                shutil.rmtree(self.test_outdir)

    def test_SingleBand_ptsMGA(self):

        raster_file = self.singletif

        out_csv = os.path.join(self.test_outdir, os.path.basename(raster_file).replace('.tif', '_b1grdext.csv'))
        rast_crs = pyprecag_crs.getCRSfromRasterFile(raster_file)

        with rasterio.open(os.path.normpath(raster_file)) as raster:
            out_gdf, _ =  extract_pixel_statistics_for_points(self.gdf_points, None, [raster_file],
                                                              function_list=[np.nanmean, raster_ops.nancv],
                                                              size_list=[1, 3, 7], output_csvfile=out_csv)

        self.assertTrue(os.path.exists(out_csv), "Output csv file doesn't exist")
        self.assertEqual(len(out_gdf), len(self.gdf_points), "Row count doesn't match")

        self.assertEqual(len(self.gdf_points.columns) + 8, len(out_gdf.columns), "Column count doesn't match")
        self.assertEqual(rast_crs.epsg_number, out_gdf.crs.to_epsg(), "Coordinate system doesn't match raster")

    def test_SingleBand_ptsWGS84(self):

        # reproject points to wgs84
        ptcrs = crs()
        ptcrs.getFromEPSG(4326)
        self.gdf_points.to_crs(epsg=4326, inplace=True)

        raster_file = self.singletif

        out_csv = os.path.join(self.test_outdir, os.path.basename(raster_file).replace('.tif', '_b1grdextwgs84.csv'))
        rast_crs = pyprecag_crs.getCRSfromRasterFile(raster_file)

        out_gdf, _ = extract_pixel_statistics_for_points(self.gdf_points, ptcrs, [raster_file],
                                                               function_list=[np.nanmean, raster_ops.nancv],
                                                               size_list=[1, 3, 7], output_csvfile=out_csv)

        self.assertTrue(os.path.exists(out_csv), "Output csv file doesn't exist")
        self.assertEqual(len(self.gdf_points), len(out_gdf), "Row count doesn't match")
        self.assertEqual(len(self.gdf_points.columns) + 8, len(out_gdf.columns), "Column count doesn't match")
        self.assertEqual(ptcrs.epsg_number, out_gdf.crs.to_epsg(), "Coordinate system doesn't match GDF")
        self.assertEqual(2, out_gdf['mean7x7_dummy_singleband_94mga54'].isnull().sum(),
                         'There should be 2 nodata values')


class TestCalculateImageIndices(unittest.TestCase):
    failedTests = []

    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestCalculateImageIndices, cls).setUpClass()

        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD, new_folder=__class__.__name__)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def run(self, result=None):
        unittest.TestCase.run(self, result)  # call superclass run method
        if self.id() in result.failed_tests or len(result.errors) > 0:
            # self.failedTests += [ea.split('.')[-1] for ea in result.failed_tests]
            self.failedTests.append(self._testMethodName)
        else:
            if os.path.exists(self.test_outdir) and not KEEP_TEST_OUTPUTS:
                shutil.rmtree(self.test_outdir)

    def test_allOptions(self):

        # '''All Options includes:
        #     Use a non-vine mask.
        #     Original image nodata is None so set to 0
        #     Reproject Image
        #     Use Shapefile AND groupby field
        # '''

        image_file = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))
        poly_shapefile = os.path.realpath(os.path.join(THIS_DIR, 'PolyMZ_wgs84_MixedPartFieldsTypes.shp'))
        indices = CalculateIndices(red=3, infrared=4, mask=5).valid_indices()
        files = calc_indices_for_block(image_file, 2,
                                       BandMapping(red=3, green=2, infrared=4, rededge=1, mask=5),
                                       self.test_outdir, indices=indices, image_epsg=32754, image_nodata=0,
                                       polygon_shapefile=poly_shapefile,
                                       groupby='part_type', out_epsg=28354)

        dst_crs = rasterio.crs.CRS.from_epsg(28354)
        self.assertEqual(len(files), 4)
        with rasterio.open(files[0]) as src:
            self.assertEqual(1, src.count, 'Incorrect band count')
            self.assertEqual(src.crs.wkt, dst_crs.wkt, 'Incorrect coordinate system')
            self.assertEqual(src.tags(1)['name'], indices[0])
            self.assertEqual(-9999, src.nodata, 'Incorrect nodata value')
            self.assertEqual('float32', src.meta['dtype'], 'Incorrect data type')

            check_df = pd.DataFrame([[62, 77, 300725.0, 6181571.0, -9999],
                                     [67, 38, 300647.0, 6181561.0, 0.2031]],
                                    columns=['row', 'col', 'x', 'y', 'expected'])

            # get values using coordinates
            check_df["actual_xy"] = [x[0] for x in src.sample(check_df[['x', 'y']].values)]

            # derive coords from row/col
            # check_df[['xn', 'yn']] = check_df.apply(lambda r: pd.Series(src.xy(r['row'], r['col'])), axis=1)
            # check_df["actual_xyn"] = [x[0] for x in src.sample(check_df[['xn', 'yn']].values)]

            # get values using row/col
            data = src.read(1)
            check_df["actual_rc"] = check_df.apply(lambda r: data[int(r['row']), int(r['col'])], axis=1)

            check_df.to_csv(os.path.join(self.test_outdir, "{}_PixelVals.csv".format(self._testMethodName)))
            with np.printoptions(precision=6, suppress=True):
                np.testing.assert_almost_equal(check_df['expected'].tolist(), check_df['actual_rc'].tolist(), 4,
                                               'Incorrect Pixel Values by Row/Col - {}'.format(
                                                   str(Path(src.name).relative_to(self.TmpDir))))

                np.testing.assert_almost_equal(check_df['expected'].tolist(), check_df['actual_xy'].tolist(), 4,
                                               'Incorrect Pixel Values by XY - {}'.format(
                                                   str(Path(src.name).relative_to(self.TmpDir))))

    def test_dontApplyNonVineMask(self):

        image_file = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))
        poly_shapefile = os.path.realpath(os.path.join(THIS_DIR, 'PolyMZ_wgs84_MixedPartFieldsTypes.shp'))
        indices = CalculateIndices(red=3, infrared=4, mask=5).valid_indices()
        files = calc_indices_for_block(image_file, 2,
                                       BandMapping(red=3, green=2, infrared=4, rededge=1),
                                       self.test_outdir, indices=indices, image_epsg=32754, image_nodata=0,
                                       polygon_shapefile=poly_shapefile,
                                       groupby='part_type', out_epsg=28354)

        dst_crs = rasterio.crs.CRS.from_epsg(28354)
        self.assertEqual(4, len(files))
        with rasterio.open(files[0]) as src:
            self.assertEqual(1, src.count, 'Incorrect band count')
            self.assertEqual(src.crs.wkt, dst_crs.wkt, 'Incorrect coordinate system')
            self.assertEqual(src.tags(1)['name'], indices[0])
            self.assertEqual(-9999, src.nodata, 'Incorrect nodata value')

            self.assertEqual('float32', src.meta['dtype'], 'Incorrect data type')

            check_df = pd.DataFrame([[62, 77, 300725.0, 6181571.0, -9999.0],
                                     [67, 38, 300647.0, 6181561.0, 0.0225]],
                                    columns=['row', 'col', 'x', 'y', 'expected'])

            # get values using coordinates
            check_df["actual_xy"] = [x[0] for x in src.sample(check_df[['x', 'y']].values)]

            # derive coords from row/col
            # check_df[['xn', 'yn']] = check_df.apply(lambda r: pd.Series(src.xy(r['row'], r['col'])), axis=1)
            # check_df["actual_xyn"] = [x[0] for x in src.sample(check_df[['xn', 'yn']].values)]

            # get values using row/col
            data = src.read(1)
            check_df["actual_rc"] = check_df.apply(lambda r: data[int(r['row']), int(r['col'])], axis=1)

            check_df.to_csv(os.path.join(self.test_outdir, "{}_PixelVals.csv".format(self._testMethodName)))

            with np.printoptions(precision=6, suppress=True):
                np.testing.assert_almost_equal(check_df['expected'].tolist(),
                                               check_df['actual_rc'].tolist(), 4,
                                               'Incorrect Pixel Values by Row/Col - {}'.format(
                                                   str(Path(src.name).relative_to(self.TmpDir))))

                np.testing.assert_almost_equal(check_df['expected'].tolist(),
                                               check_df['actual_xy'].tolist(), 4,
                                               'Incorrect Pixel Values by XY - {}'.format(
                                                   str(Path(src.name).relative_to(self.TmpDir))))

    def test_noShapefile(self):
        # Use Full Image......
        #     No Shapfile,
        #     No Non-Vine mask

        image_file = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))
        indices = CalculateIndices(red=3, infrared=4, mask=5).valid_indices()
        files = calc_indices_for_block(image_file, 2,
                                       BandMapping(red=3, green=2, infrared=4, rededge=1),
                                       self.test_outdir, indices=indices, image_epsg=32754, image_nodata=0,
                                       out_epsg=28354)

        dst_crs = rasterio.crs.CRS.from_epsg(28354)
        self.assertEqual(len(files), len(indices))
        with rasterio.open(files[0]) as src:
            self.assertEqual(1, src.count, 'Incorrect band count')
            self.assertEqual(src.crs.wkt, dst_crs.wkt, 'Incorrect Coordinate System')
            self.assertEqual(src.tags(1)['name'], indices[0])
            self.assertEqual(-9999, src.nodata, 'Incorrect nodata value')
            self.assertEqual('float32', src.meta['dtype'], 'Incorrect data type')

            check_df = pd.DataFrame([[114, 278, 300725.0, 6181571.0, -0.0506],
                                     [119, 239, 300647.0, 6181561.0, 0.0225],
                                     [180, 356, 300881.342, 6181439.444, -9999.0]],
                                    columns=['row', 'col', 'x', 'y', 'expected'])

            # get values using coordinates
            check_df["actual_xy"] = [x[0] for x in src.sample(check_df[['x', 'y']].values)]

            # derive coords from row/col
            # check_df[['xn', 'yn']] = check_df.apply(lambda r: pd.Series(src.xy(r['row'], r['col'])), axis=1)
            # check_df["actual_xyn"] = [x[0] for x in src.sample(check_df[['xn', 'yn']].values)]

            # get values using row/col
            data = src.read(1)
            check_df["actual_rc"] = check_df.apply(lambda r: data[int(r['row']), int(r['col'])], axis=1)

            check_df.to_csv(os.path.join(self.test_outdir, "{}_PixelVals.csv".format(self._testMethodName)))

            with np.printoptions(precision=6, suppress=True):
                np.testing.assert_almost_equal(check_df['expected'].tolist(),
                                               check_df['actual_rc'].tolist(), 4,
                                               'Incorrect Pixel Values by Row/Col - {}'.format(
                                                   str(Path(src.name).relative_to(self.TmpDir))))

                np.testing.assert_almost_equal(check_df['expected'].tolist(),
                                               check_df['actual_xy'].tolist(), 4,
                                               'Incorrect Pixel Values by XY - {}'.format(
                                                   str(Path(src.name).relative_to(self.TmpDir))))

    def test_noGroupby(self):

        bm = BandMapping(green=2, infrared=4, rededge=1, mask=5)

        image_file = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))
        poly_shapefile = os.path.realpath(os.path.join(THIS_DIR, 'PolyMZ_wgs84_MixedPartFieldsTypes.shp'))
        indices = CalculateIndices(**bm).valid_indices()
        files = calc_indices_for_block(image_file, 2,
                                       BandMapping(red=3, green=2, infrared=4, rededge=1),
                                       self.test_outdir, indices=indices, image_epsg=32754, image_nodata=0,
                                       polygon_shapefile=poly_shapefile,
                                       out_epsg=28354)

        self.assertEqual(len(indices), 3)
        self.assertEqual(len(files), len(indices))
        with rasterio.open(files[0]) as src:
            self.assertEqual(-9999, src.nodata, 'Incorrect nodata value')
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), src.crs, 'Incorrect coordinate system ')
            self.assertEqual(('float32',), src.dtypes, 'Incorrect data type')
            self.assertEqual(src.res, (2.0, 2.0))
            self.assertEqual(1, src.count, 'Incorrect band count')


class TestResampleToBlock(unittest.TestCase):
    failedTests = []

    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestResampleToBlock, cls).setUpClass()
        cls.TmpDir = setup_folder(base_folder=TEMP_FOLD, new_folder=__class__.__name__)

    @classmethod
    def tearDownClass(cls):
        if len(cls.failedTests) == 0 and not KEEP_TEST_OUTPUTS:
            print('Tests Passed .. Deleting {}'.format(TEMP_FOLD))
            shutil.rmtree(TEMP_FOLD)

    def setUp(self):
        self.startTime = time.time()
        self.test_outdir = setup_folder(self.TmpDir, new_folder=self._testMethodName)

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def run(self, result=None):

        unittest.TestCase.run(self, result)  # call superclass run method

        if self.id() in result.failed_tests or len(result.errors) > 0:
            self.failedTests.append(self._testMethodName)
        else:
            if os.path.exists(self.test_outdir) and not KEEP_TEST_OUTPUTS:
                shutil.rmtree(self.test_outdir)

    def test_allOptions(self):

        # All Options includes:
        #     Use a non-vine mask.
        #     Original image nodata is None so set to 0
        #     Reproject Image
        #     Use Shapefile AND groupby field

        image_file = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))
        poly_shapefile = os.path.realpath(os.path.join(THIS_DIR, 'PolyMZ_wgs84_MixedPartFieldsTypes.shp'))
        files = resample_bands_to_block(image_file, 2,
                                        self.test_outdir, band_nums=[6], image_epsg=32754, image_nodata=0,
                                        polygon_shapefile=poly_shapefile,
                                        groupby='part_type', out_epsg=28354)

        dst_crs = rasterio.crs.CRS.from_epsg(28354)
        self.assertEqual(len(files), 2)

        with rasterio.open(files[0]) as src:
            self.assertEqual(src.count, 1, 'Incorrect band count')
            self.assertEqual(src.crs.wkt, dst_crs.wkt, 'Incorrect coordinate system ')
            self.assertEqual(src.nodata, 0, 'Incorrect image nodata value')
            self.assertEqual('float32', src.meta['dtype'], 'Incorrect data type')

            check_df = pd.DataFrame([[61, 46, 300663.0, 6181573.0, 0.0],
                                     [22, 52, 300675.0, 6181651.0, 899.3889]],
                                    columns=['row', 'col', 'x', 'y', 'expected'])

            # get values using coordinates
            check_df["actual_xy"] = [x[0] for x in src.sample(check_df[['x', 'y']].values)]

            # derive coords from row/col
            # check_df[['xn', 'yn']] = check_df.apply(lambda r: pd.Series(src.xy(r['row'], r['col'])), axis=1)
            # check_df["actual_xyn"] = [x[0] for x in src.sample(check_df[['xn', 'yn']].values)]

            # get values using row/col
            data = src.read(1)
            check_df["actual_rc"] = check_df.apply(lambda r: data[int(r['row']), int(r['col'])], axis=1)

            check_df.to_csv(os.path.join(self.test_outdir, "{}_PixelVals.csv".format(self._testMethodName)))

            with np.printoptions(precision=6, suppress=True):
                np.testing.assert_almost_equal(check_df['expected'].tolist(),
                                               check_df['actual_rc'].tolist(), 4,
                                               'Incorrect Pixel Values by Row/Col - {}'.format(
                                                   str(Path(src.name).relative_to(self.TmpDir))))

                np.testing.assert_almost_equal(check_df['expected'].tolist(),
                                               check_df['actual_xy'].tolist(), 4,
                                               'Incorrect Pixel Values by XY - {}'.format(
                                                   str(Path(src.name).relative_to(self.TmpDir))))

    def test_noShapefile(self):
        # Use Full Image......
        #     No Shapfile,
        #     No Non-Vine mask

        image_file = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))
        files = resample_bands_to_block(image_file, 2, self.test_outdir, band_nums=[6], image_epsg=32754,
                                        image_nodata=0, out_epsg=28354)

        dst_crs = rasterio.crs.CRS.from_epsg(28354)

        with rasterio.open(files[0]) as src:
            self.assertEqual(1, src.count, 'Incorrect band count')
            self.assertEqual(src.crs.wkt, dst_crs.wkt, 'Incorrect coordinate system ')
            self.assertEqual(src.nodata, 0, 'Incorrect image nodata value')
            self.assertEqual('float32', src.meta['dtype'], 'Incorrect data type')

            check_df = pd.DataFrame([[114, 278, 300725.0, 6181571.0, 0.0],
                                     [119, 239, 300647.0, 6181561.0, 916.85],
                                     [180, 356, 300881.342, 6181439.444, 0.0]],
                                    columns=['row', 'col', 'x', 'y', 'expected'])

            # get values using coordinates
            check_df["actual_xy"] = [x[0] for x in src.sample(check_df[['x', 'y']].values)]

            # derive coords from row/col
            # check_df[['xn', 'yn']] = check_df.apply(lambda r: pd.Series(src.xy(r['row'], r['col'])), axis=1)
            # check_df["actual_xyn"] = [x[0] for x in src.sample(check_df[['xn', 'yn']].values)]

            # get values using row/col
            data = src.read(1)
            check_df["actual_rc"] = check_df.apply(lambda r: data[int(r['row']), int(r['col'])], axis=1)

            check_df.to_csv(os.path.join(self.test_outdir, "{}_PixelVals.csv".format(self._testMethodName)))

            with np.printoptions(precision=6, suppress=True):
                np.testing.assert_almost_equal(check_df['expected'].tolist(),
                                               check_df['actual_rc'].tolist(), 4,
                                               'Incorrect Pixel Values by Row/Col - {}'.format(
                                                   str(Path(src.name).relative_to(self.TmpDir))))

                np.testing.assert_almost_equal(check_df['expected'].tolist(),
                                               check_df['actual_xy'].tolist(), 4,
                                               'Incorrect Pixel Values by XY - {}'.format(
                                                   str(Path(src.name).relative_to(self.TmpDir))))

    def test_noGroupby(self):

        image_file = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))
        poly_shapefile = os.path.realpath(os.path.join(THIS_DIR, 'PolyMZ_wgs84_MixedPartFieldsTypes.shp'))
        files = resample_bands_to_block(image_file, 2, self.test_outdir, band_nums=[6], image_epsg=32754,
                                        image_nodata=0, polygon_shapefile=poly_shapefile,
                                        out_epsg=28354)

        self.assertEqual(1, len(files))

        with rasterio.open(files[0]) as src:
            self.assertEqual(0.0, src.nodata, 'Incorrect Image Nodata value')
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), src.crs, 'Incorrect coordinate system ')
            self.assertEqual(('float32',), src.dtypes, 'Incorrect data type')
            self.assertEqual((2.0, 2.0), src.res, 'Incorrect pixel size')
            self.assertEqual(1, src.count, 'Incorrect band count')

    def test_nonStandardNoDataNotSet(self):

        # change input image nodata to 7777 in the image but leave the nodata as none

        new_image = os.path.join(self.test_outdir, 'geotif_32754_no-nodata-7777.tif')
        crs_sutm54 = rasterio.crs.CRS.from_epsg(32754)

        image_file = os.path.realpath(os.path.join(THIS_DIR, 'rasters', 'area1_rgbi_jan_50cm_84sutm54.tif'))
        with rasterio.open(image_file) as src:
            meta = src.meta.copy()
            del meta['crs']

            with rasterio.open(new_image, 'w', crs=crs_sutm54, **meta) as dest:
                for i in range(1, src.count + 1):
                    data = src.read(i)
                    np.putmask(data, data == 0, 7777)
                    dest.write(data, i)

        poly_shapefile = os.path.realpath(os.path.join(THIS_DIR, 'PolyMZ_wgs84_MixedPartFieldsTypes.shp'))
        files = resample_bands_to_block(new_image, 2.5, self.test_outdir, band_nums=[6], image_epsg=32754,
                                        image_nodata=7777, polygon_shapefile=poly_shapefile,
                                        out_epsg=28354)

        self.assertEqual(len(files), 1)

        with rasterio.open(files[0]) as src:
            self.assertEqual(7777, src.nodata, 'Incorrect image nodata value')
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), src.crs, 'Incorrect coordinate system ')
            self.assertEqual(('float32',), src.dtypes, 'Incorrect data type')
            self.assertEqual((2.5, 2.5), src.res, 'Incorrect pixel size')
            self.assertEqual(1, src.count, 'Incorrect band count')
