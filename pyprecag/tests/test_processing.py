import shutil
import tempfile
import unittest

from geopandas import GeoSeries, GeoDataFrame

from pyprecag.bandops import BandMapping, CalculateIndices
from pyprecag.crs import crs
from pyprecag.tests import make_dummy_data
from pyprecag import convert, raster_ops
from pyprecag.processing import *

pyFile = os.path.basename(__file__)
TmpDir = tempfile.gettempdir()
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])

this_dir = os.path.abspath(os.path.dirname(__file__))

logging.captureWarnings(False)
logging.basicConfig(level=logging.INFO, format="%(message)s")


class TestProcessing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestProcessing, cls).setUpClass()

        if os.path.exists(TmpDir):
            print 'Folder Exists.. Deleting {}'.format(TmpDir)
            shutil.rmtree(TmpDir)

        os.mkdir(TmpDir)
        cls.singletif, cls.multitif = make_dummy_data.make_dummy_tif_files(TmpDir)
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
        unittest.TestCase.run(self, result)  # call superclass run method
        if len(result.failures) > 0 or len(result.errors) > 0:
            testFailed = True

    def test_BlockGrid(self):
        poly = os.path.realpath(this_dir + "/data/area2_onebox_94mga54.shp")

        file_sub_name = os.path.join(TmpDir, os.path.splitext(os.path.basename(poly))[0])

        block_grid(in_shapefilename=poly,
                   pixel_size=5,
                   out_rasterfilename=file_sub_name + '_block.tif',
                   out_vesperfilename=file_sub_name + '_block_v.txt',
                   snap=True,
                   overwrite=True)

        self.assertTrue(os.path.exists(file_sub_name + '_block.tif'))
        self.assertTrue(os.path.exists(file_sub_name + '_block_v.txt', ))

        with rasterio.open(os.path.normpath(file_sub_name + '_block.tif')) as dataset:
            self.assertEqual(dataset.count, 1)
            self.assertEqual(dataset.width, 47)
            self.assertEqual(dataset.height, 26)
            self.assertEqual(dataset.nodatavals, (-9999.0,))
            self.assertEqual(dataset.dtypes, ('int16',))

    def test_cleanTrimPoints(self):
        in_csv = os.path.join(this_dir + "/data/area2_yield_ISO-8859-1.csv")
        in_poly = os.path.join(this_dir + "/data/area2_onebox_94mga54.shp")
        out_csv = os.path.join(TmpDir, os.path.basename(in_csv))
        out_shp = os.path.join(TmpDir, os.path.basename(in_csv).replace('.csv', '.shp'))
        out_rm_shp = os.path.join(TmpDir, os.path.basename(in_csv).replace('.csv', '_remove.shp'))

        gdf_points, gdf_pts_crs = convert.convert_csv_to_points(in_csv, coord_columns_epsg=4326,
                                                                out_epsg=28354)
        out_gdf, out_crs = clean_trim_points(gdf_points, gdf_pts_crs, 'Yld Mass(Dry)(tonne/ha)',
                                             out_csv, out_keep_shapefile=out_shp,
                                             out_removed_shapefile=out_rm_shp,
                                             boundary_polyfile=in_poly, thin_dist_m=2.5)

        self.assertIsInstance(out_gdf, GeoDataFrame)
        self.assertTrue(os.path.exists(out_csv))
        self.assertTrue(gdf_pts_crs, out_crs)
        self.assertEqual(out_gdf.crs, {'init': 'epsg:28354', 'no_defs': True})
        self.assertEqual(len(out_gdf), 554)
        self.assertIn('EN_EPSG', out_gdf.columns)

    def test_createPolygonFromPointTrail(self):
        in_csv = os.path.join(this_dir + "/data/area2_yield_ISO-8859-1.csv")

        out_polyfile = os.path.join(TmpDir, os.path.splitext(os.path.basename(in_csv))[0] + '_poly.shp')

        gdf_points, gdf_pts_crs = convert.convert_csv_to_points(
            in_csv, None, coord_columns_epsg=4326, out_epsg=28354)

        create_polygon_from_point_trail(gdf_points, gdf_pts_crs, out_polyfile,
                                        thin_dist_m=2.5,
                                        aggregate_dist_m=25,
                                        buffer_dist_m=10,
                                        shrink_dist_m=3)

        self.assertTrue(os.path.exists(out_polyfile), True)

        vect_desc = VectorDescribe(out_polyfile)
        self.assertEqual(vect_desc.crs.epsg_number, 28354)
        self.assertFalse(vect_desc.is_mz_aware)
        self.assertEqual(vect_desc.geometry_type, 'Polygon')

    def test_randomPixelSelection(self):
        raster_file = self.singletif
        out_shp = os.path.join(
            TmpDir, os.path.basename(raster_file).replace('.tif', '_randpts.shp'))
        rast_crs = pyprecag_crs.getCRSfromRasterFile(raster_file)

        with rasterio.open(os.path.normpath(raster_file)) as raster:
            rand_gdf, rand_crs = random_pixel_selection(raster, rast_crs, 50, out_shapefile=out_shp)

        self.assertEqual(len(rand_gdf), 50)
        self.assertTrue(os.path.exists(out_shp))
        self.assertEqual(rand_crs, rast_crs)

    def test_kmeansCluster(self):
        raster_files = glob.glob(os.path.realpath(this_dir + '/data/rasters/*one*.tif'))
        raster_files += [os.path.realpath(this_dir +
                                          '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif')]

        out_img = os.path.join(TmpDir, 'kmeans-cluster_3cluster_3rasters.tif')

        with self.assertRaises(TypeError) as msg:
            _ = kmeans_clustering(raster_files, out_img)

        self.assertIn('raster_files are of different pixel sizes', str(msg.exception))
        raster_files.remove(os.path.realpath(this_dir + '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif'))

        with self.assertRaises(TypeError) as msg:
            _ = kmeans_clustering(raster_files, out_img)

        self.assertIn("1 raster(s) don't have coordinates systems assigned", str(msg.exception))
        raster_files.remove(os.path.realpath(this_dir + '/data/rasters/area1_onebox_NDRE_250cm.tif'))

        out_df = kmeans_clustering(raster_files, out_img)

        self.assertTrue(os.path.exists(out_img))
        self.assertTrue(os.path.exists(out_img.replace('.tif', '_statistics.csv')))
        self.assertEqual(3, len(out_df['zone'].unique()))

        with rasterio.open(out_img) as src:
            self.assertEqual(1, src.count)
            self.assertIn(src.crs.to_string(), ['EPSG:28354','+init=epsg:28354'])

            self.assertEqual(0, src.nodata)
            band1 = src.read(1, masked=True)
            self.assertItemsEqual(np.array([0, 1, 2, 3]), np.unique(band1.data))

    def test_PersistorAllYears(self):
        raster_files = glob.glob(os.path.realpath(this_dir + '/data/rasters/Year*.tif'))
        raster_files += [os.path.realpath(this_dir +
                                          '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif') ]

        out_img = os.path.join(TmpDir, 'persistor_allyears.tif')

        # check for raised pixel size error
        with self.assertRaises(TypeError) as msg:
            persistor_all_years(raster_files,out_img, True, 10)
        self.assertIn("raster_files are of different pixel sizes", str(msg.exception))
        raster_files.remove(os.path.realpath(this_dir + '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif'))

        persistor_all_years(raster_files, out_img, True, 10)
        self.assertTrue(os.path.exists(out_img))

        src_img = os.path.realpath(this_dir + '/data/rasters/persistor_allyears.tif')
        with rasterio.open(src_img) as src, \
            rasterio.open(os.path.normpath(out_img)) as test:
            self.assertEqual(src.profile, test.profile)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), test.crs)

            np.testing.assert_array_equal(src.read(), test.read())

            band1 = test.read(1, masked=True)
            self.assertItemsEqual([-9999, 0, 1, 2, 3], np.unique(band1.data).tolist())

    def test_PersistorTargetProb(self):
        raster_files = glob.glob(os.path.realpath(this_dir + '/data/rasters/Year*.tif'))

        out_img = os.path.join(TmpDir, 'persistor_targetprob.tif')

        persistor_target_probability(raster_files, 10, 75,
                                     raster_files, -10, 75, out_img)

        self.assertTrue(os.path.exists(out_img))
        src_img = os.path.realpath(this_dir + '/data/rasters/persistor_targetprob.tif')
        with rasterio.open(src_img) as src, \
                rasterio.open(os.path.normpath(out_img)) as test:
            self.assertEqual(src.profile, test.profile)
            self.assertEqual(rasterio.crs.CRS.from_epsg(28354), test.crs)

            np.testing.assert_array_equal(src.read(), test.read())

            band1 = test.read(1, masked=True)
            self.assertItemsEqual([-9999, -1, 0, 1], np.unique(band1.data).tolist())


class TestStripTrials(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestStripTrials, cls).setUpClass()

        if os.path.exists(TmpDir):
            print 'Folder Exists.. Deleting {}'.format(TmpDir)
            shutil.rmtree(TmpDir)

        os.mkdir(TmpDir)
        cls.singletif, cls.multitif = make_dummy_data.make_dummy_tif_files(TmpDir)
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
        unittest.TestCase.run(self, result)  # call superclass run method
        if len(result.failures) > 0 or len(result.errors) > 0:
            testFailed = True

    def test_CreateStripTreatmentPoints(self):
        import pandas
        import geopandas
        print ('Pandas Version is ', pandas.__version__)
        print ('GeoPandas Version is ', geopandas.__version__)

        out_points_file = os.path.join(TmpDir, 'testCreateStripTreatment_Points.shp')
        out_lines_file = os.path.join(TmpDir, 'testCreateStripTreatment_Lines.shp')

        in_lines_gdf = GeoDataFrame(
            {'geometry': [LineString([(740800, 6169700, 5), (741269, 6169700, 5)]),
                          LineString([(741000, 6169800, 6), (741003, 6170012, 6)]),
                          LineString([(741400, 6169800, 7), (741355, 6169780, 7),
                                      (741300, 6169800, 7), (741250, 6169950, 7),
                                      (741150, 6169950, 7), (741100, 6170000, 7)]),
                          LineString([(740800, 6169912, 8), (740900, 6170094, 8)])]
                , 'TrialID': [1, 2, 3, 4]}, crs={'init': 'epsg:28354'})

        gdf_lines_crs = crs()
        gdf_lines_crs.getFromEPSG(28354)

        out_points_gdf, out_crs, out_lines_gdf = \
            create_points_along_line(in_lines_gdf, gdf_lines_crs, 7, 25,
                                     out_points_shapefile=out_points_file,
                                     out_lines_shapefile=out_lines_file)

        # Test Points Output-----------------------------------------------------------
        self.assertIsInstance(out_points_gdf, GeoDataFrame)
        self.assertEqual(573, len(out_points_gdf))
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
        in_lines_gdf = GeoDataFrame({'geometry': lines, 'block': ['test']})

        # without saving to file on a single multi vertex line
        out_points_gdf, out_crs, out_lines_gdf = create_points_along_line(in_lines_gdf,
                                                                          lines_crs, 31, 25)

        self.assertIsInstance(out_points_gdf, GeoDataFrame)
        self.assertEqual(69, len(out_points_gdf))

        self.assertEqual(lines_crs.epsg_number, out_crs.epsg_number)
        self.assertEqual(out_points_gdf.crs, out_lines_gdf.crs)

        self.assertEqual({'Point'}, set(out_points_gdf.geom_type))
        self.assertEqual(['N Strip', 'S Strip', 'Strip'],
                         sorted(list(out_points_gdf['Strip_Name'].unique())))

        self.assertEqual(3, len(out_lines_gdf))
        self.assertEqual(['N Strip', 'S Strip', 'Strip'],
                         (list(out_lines_gdf['Strip_Name'].unique())))

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

        gdf_points = GeoDataFrame(data, columns=columns, geometry='geometry',
                                  crs={'init': 'epsg:28354'})
        crs_points = crs()
        crs_points.getFromEPSG(28354)

        values_rast = os.path.realpath(this_dir + "/data/rasters/Yield_withStrip_PRED_2m.tif")
        control_raster = os.path.realpath(this_dir + "/data/rasters/Yield_withoutStrip_PRED_2m.tif")
        zone_raster = os.path.realpath(this_dir + "/data/rasters/k-means_3clusters_3rasters_2m.tif")

        result = ttest_analysis(gdf_points, crs_points,
                                values_rast, TmpDir,
                                zone_raster=zone_raster,
                                control_raster=control_raster)

        self.assertEqual(13, len(result.columns))
        file_csv = glob.glob(os.path.join(TmpDir, '*.csv'))
        file_graph = glob.glob(os.path.join(TmpDir, '*graph.png'))
        file_map = glob.glob(os.path.join(TmpDir, '*map.png'))

        self.assertEqual(4, len(file_csv))
        self.assertEqual(4, len(file_graph))
        self.assertEqual(2, len(file_map))


class TestExtractRasterStatisticsForPoints(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestExtractRasterStatisticsForPoints, cls).setUpClass()

        if os.path.exists(TmpDir):
            print 'Folder Exists.. Deleting {}'.format(TmpDir)
            shutil.rmtree(TmpDir)

        os.mkdir(TmpDir)
        cls.singletif, cls.multitif = make_dummy_data.make_dummy_tif_files(TmpDir)
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
        unittest.TestCase.run(self, result)  # call superclass run method
        if len(result.failures) > 0 or len(result.errors) > 0:
            testFailed = True

    def test_SingleBand_MGA(self):

        raster_file = self.singletif
        out_csv = os.path.join(TmpDir,
                               os.path.basename(raster_file).replace('.tif', '_b1grdext.csv'))
        rast_crs = pyprecag_crs.getCRSfromRasterFile(raster_file)

        with rasterio.open(os.path.normpath(raster_file)) as raster:
            pts_gdf, pts_crs = random_pixel_selection(raster, rast_crs, 50)

            out_gdf, out_crs = \
                extract_pixel_statistics_for_points(pts_gdf, pts_crs, [raster_file],
                                                    function_list=[np.nanmean, raster_ops.nancv],
                                                    size_list=[1, 3, 7], output_csvfile=out_csv)

        self.assertTrue(os.path.exists(out_csv))
        self.assertTrue(len(out_gdf), len(pts_gdf))
        self.assertTrue(len(pts_gdf.columns), len(out_gdf.columns) + 5)
        self.assertTrue(out_crs, rast_crs)

    def test_SingleBand_WGS84(self):
        raster_file = self.singletif
        _ = pyprecag_crs.getCRSfromRasterFile(raster_file)

        raster_file = self.singletif
        out_csv = os.path.join(TmpDir,
                               os.path.basename(raster_file).replace('.tif', '_b1grdextwgs84.csv'))
        rast_crs = pyprecag_crs.getCRSfromRasterFile(raster_file)

        with rasterio.open(os.path.normpath(raster_file)) as raster:
            pts_gdf, pts_crs = random_pixel_selection(raster, rast_crs, 50)

        pts_crs.getFromEPSG(4326)
        pts_gdf.to_crs(epsg=4326, inplace=True)

        out_gdf, out_crs = \
            extract_pixel_statistics_for_points(pts_gdf, pts_crs, [raster_file],
                                                function_list=[np.nanmean, raster_ops.nancv],
                                                size_list=[1, 3, 7], output_csvfile=out_csv)

        self.assertTrue(os.path.exists(out_csv))
        self.assertTrue(len(out_gdf), len(pts_gdf))
        self.assertTrue(len(pts_gdf.columns), len(out_gdf.columns) + 5)
        self.assertTrue(out_crs, rast_crs)
        self.assertTrue(pts_crs, out_crs)
        self.assertEqual(out_gdf['mean7x7_test_singleband_94mga54'].isnull().sum(), 0)


class TestCalculateImageIndices(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestCalculateImageIndices, cls).setUpClass()
        if os.path.exists(TmpDir):
            print 'Folder Exists.. Deleting {}'.format(TmpDir)
            shutil.rmtree(TmpDir)

        os.makedirs(TmpDir)

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
        unittest.TestCase.run(self, result)  # call superclass run method
        if len(result.failures) > 0 or len(result.errors) > 0:
            testFailed = True

    def test_allOptions(self):

        """ All Options includes:
            Use a non-vine mask.
            Original image nodata is None so set to 0
            Reproject Image
            Use Shapefile AND groupby field
        """

        out_dir = os.path.join(TmpDir, 'TestCalculateImageIndices', 'all-opts')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif')
        poly_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        indices = CalculateIndices(red=3, infrared=4, mask=5).valid_indices()
        files = calc_indices_for_block(image_file, 2,
                                       BandMapping(red=3, green=2, infrared=4, rededge=1, mask=5),
                                       out_dir, indices=indices, image_epsg=32754, image_nodata=0,
                                       polygon_shapefile=poly_shapefile,
                                       groupby='part_type', out_epsg=28354)

        dst_crs = rasterio.crs.CRS.from_epsg(28354)
        self.assertEqual(len(files), 4)
        with rasterio.open(files[0]) as src:
            self.assertEqual(src.count, 1)
            self.assertEqual(src.crs.wkt, dst_crs.wkt)
            self.assertEqual(src.tags(1)['name'], indices[0])
            self.assertEqual(src.nodata, -9999)
            self.assertEqual(src.meta['dtype'], 'float32')

            # coords 300725.0, 6181571.0
            self.assertEqual(src.read(1)[47, 62], -9999)
            # coords (300647.0, 6181561.0)
            self.assertAlmostEqual(src.read(1)[52, 23], 0.20820355, 4)

    def test_dontApplyNonVineMask(self):
        out_dir = os.path.join(TmpDir, 'TestCalculateImageIndices', 'no-nonvine')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif')
        poly_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        indices = CalculateIndices(red=3, infrared=4, mask=5).valid_indices()
        files = calc_indices_for_block(image_file, 2,
                                       BandMapping(red=3, green=2, infrared=4, rededge=1),
                                       out_dir, indices=indices, image_epsg=32754, image_nodata=0,
                                       polygon_shapefile=poly_shapefile,
                                       groupby='part_type', out_epsg=28354)

        dst_crs = rasterio.crs.CRS.from_epsg(28354)
        self.assertEqual(len(files), 4)
        with rasterio.open(files[0]) as src:
            self.assertEqual(src.count, 1)
            self.assertEqual(src.crs.wkt, dst_crs.wkt)
            self.assertEqual(src.tags(1)['name'], indices[0])
            self.assertEqual(src.nodata, -9999)
            self.assertEqual(src.meta['dtype'], 'float32')

            # coords 300725.0, 6181571.0
            self.assertEqual(src.read(1)[47, 62], -9999)
            # coords (300647.0, 6181561.0)
            self.assertAlmostEqual(src.read(1)[52, 23], 0.070868947, 4)

    def test_noShapefile(self):

        """ Use Full Image......
            No Shapfile,
            No Non-Vine mask
        """
        out_dir = os.path.join(TmpDir, 'TestCalculateImageIndices', 'no-shapefile')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif')
        indices = CalculateIndices(red=3, infrared=4, mask=5).valid_indices()
        files = calc_indices_for_block(image_file, 2,
                                       BandMapping(red=3, green=2, infrared=4, rededge=1),
                                       out_dir, indices=indices, image_epsg=32754, image_nodata=0,
                                       out_epsg=28354)

        dst_crs = rasterio.crs.CRS.from_epsg(28354)
        self.assertEqual(len(files), len(indices))
        with rasterio.open(files[0]) as src:
            self.assertEqual(src.count, 1)
            self.assertEqual(src.crs.wkt, dst_crs.wkt)
            self.assertEqual(src.tags(1)['name'], indices[0])
            self.assertEqual(src.nodata, -9999)
            self.assertEqual(src.meta['dtype'], 'float32')

            # test for values at coords
            row, col = src.index(300725, 6181571)
            self.assertAlmostEqual(src.read(1)[int(row), int(col)], -0.05087604, 4)

            row, col = src.index(300647.0, 6181561.0)
            self.assertAlmostEqual(src.read(1)[int(row), int(col)], 0.02232674, 4)

            row, col = src.index(300881.342, 6181439.444)
            self.assertEqual(src.read(1)[int(row), int(col)], -9999)

    def test_noGroupby(self):
        out_dir = os.path.join(TmpDir, 'TestCalculateImageIndices', 'no-groupby')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        bm = BandMapping(green=2, infrared=4, rededge=1, mask=5)

        image_file = os.path.realpath(this_dir + '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif')
        poly_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        indices = CalculateIndices(**bm).valid_indices()
        files = calc_indices_for_block(image_file, 2,
                                       BandMapping(red=3, green=2, infrared=4, rededge=1),
                                       out_dir, indices=indices, image_epsg=32754, image_nodata=0,
                                       polygon_shapefile=poly_shapefile,
                                       out_epsg=28354)

        self.assertEqual(len(indices), 3)
        self.assertEqual(len(files), len(indices))
        with rasterio.open(files[0]) as src:
            self.assertEqual(src.nodata, -9999)
            self.assertEqual(src.crs, rasterio.crs.CRS.from_epsg(28354))
            self.assertEqual(src.dtypes, ('float32',))
            self.assertEqual(src.res, (2.0, 2.0))
            self.assertEqual(src.count, 1)


class TestResampleToBlock(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(TestResampleToBlock, cls).setUpClass()

        if os.path.exists(TmpDir):
            print 'Folder Exists.. Deleting {}'.format(TmpDir)
            shutil.rmtree(TmpDir)

        os.makedirs(TmpDir)

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
        unittest.TestCase.run(self, result)  # call superclass run method
        if len(result.failures) > 0 or len(result.errors) > 0:
            testFailed = True

    def test_allOptions(self):

        """ All Options includes:
            Use a non-vine mask.
            Original image nodata is None so set to 0
            Reproject Image
            Use Shapefile AND groupby field
        """

        out_dir = os.path.join(TmpDir, 'TestResampleToBlock', 'all-opts')

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif')
        poly_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        files = resample_bands_to_block(image_file, 2,
                                        out_dir, band_nums=[6], image_epsg=32754, image_nodata=0,
                                        polygon_shapefile=poly_shapefile,
                                        groupby='part_type', out_epsg=28354)

        dst_crs = rasterio.crs.CRS.from_epsg(28354)
        self.assertEqual(len(files), 2)
        with rasterio.open(files[0]) as src:
            self.assertEqual(src.count, 1)
            self.assertEqual(src.crs.wkt, dst_crs.wkt)
            self.assertEqual(src.nodata, 0)
            self.assertEqual(src.meta['dtype'], 'float32')

            # coords 300725.0, 6181571.0
            self.assertEqual(src.read(1)[47, 62], 0)
            # coords (300647.0, 6181561.0)
            self.assertAlmostEqual(src.read(1)[52, 23], 855.23999, 4)

    def test_noShapefile(self):

        """ Use Full Image......
            No Shapfile,
            No Non-Vine mask
        """
        out_dir = os.path.join(TmpDir, 'TestResampleToBlock', 'no-shapefile')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif')
        files = resample_bands_to_block(image_file, 2, out_dir, band_nums=[6], image_epsg=32754,
                                        image_nodata=0, out_epsg=28354)

        dst_crs = rasterio.crs.CRS.from_epsg(28354)

        with rasterio.open(files[0]) as src:
            self.assertEqual(src.count, 1)
            self.assertEqual(src.crs.wkt, dst_crs.wkt)
            self.assertEqual(src.nodata, 0)
            self.assertEqual(src.meta['dtype'], 'float32')

            # test for values at coords
            row, col = src.index(300725, 6181571)
            self.assertAlmostEqual(src.read(1)[int(row), int(col)], 0, 4)

            row, col = src.index(300647.0, 6181561.0)
            self.assertAlmostEqual(src.read(1)[int(row), int(col)], 917.34998, 4)

            row, col = src.index(300881.342, 6181439.444)
            self.assertEqual(src.read(1)[int(row), int(col)], 0)

    def test_noGroupby(self):
        out_dir = os.path.join(TmpDir, 'TestResampleToBlock', 'no-groupby')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif')
        poly_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        files = resample_bands_to_block(image_file, 2, out_dir, band_nums=[6], image_epsg=32754,
                                        image_nodata=0, polygon_shapefile=poly_shapefile,
                                        out_epsg=28354)

        self.assertEqual(len(files), 1)

        with rasterio.open(files[0]) as src:
            self.assertEqual(src.nodata, 0.0)
            self.assertEqual(src.crs, rasterio.crs.CRS.from_epsg(28354))
            self.assertEqual(src.dtypes, ('float32',))
            self.assertEqual(src.res, (2.0, 2.0))
            self.assertEqual(src.count, 1)

    def test_nonStandardNoDataNotSet(self):
        """ change input image nodata to 7777 in the image but leave the nodata as none."""

        out_dir = os.path.join(TmpDir, 'TestResampleToBlock', 'nonstandard-nodata-notset')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # change nodata to 7777 but don't set it in the output
        new_image = os.path.join(out_dir, 'geotif_32754_no-nodata-7777.tif')
        crs_sutm54 = rasterio.crs.CRS.from_epsg(32754)

        image_file = os.path.realpath(this_dir + '/data/rasters/area1_rgbi_jan_50cm_84sutm54.tif')
        with rasterio.open(image_file) as src:
            meta = src.meta.copy()
            del meta['crs']

            with rasterio.open(new_image, 'w', crs=crs_sutm54, **meta) as dest:
                for i in range(1, src.count + 1):
                    data = src.read(i)
                    np.putmask(data, data == 0, 7777)
                    dest.write(data, i)

        poly_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        files = resample_bands_to_block(new_image, 2.5, out_dir, band_nums=[6], image_epsg=32754,
                                        image_nodata=7777, polygon_shapefile=poly_shapefile,
                                        out_epsg=28354)

        self.assertEqual(len(files), 1)

        with rasterio.open(files[0]) as src:
            self.assertEqual(src.nodata, 7777)
            self.assertEqual(src.crs, rasterio.crs.CRS.from_epsg(28354))
            self.assertEqual(src.dtypes, ('float32',))
            self.assertEqual(src.res, (2.5, 2.5))
            self.assertEqual(src.count, 1)
