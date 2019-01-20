
import shutil
import tempfile
import unittest

from pyprecag.bandops import BandMapping, CalculateIndices
from pyprecag.tests import make_dummy_data
from pyprecag import convert, raster_ops
from pyprecag.processing import *

pyFile = os.path.basename(__file__)
TmpDir = tempfile.gettempdir()
TmpDir = os.path.join(TmpDir, os.path.splitext(pyFile)[0])

this_dir = os.path.abspath(os.path.dirname(__file__))

logging.captureWarnings(False)
logging.basicConfig(level=logging.INFO, format="%(message)s")

class test_Processing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_Processing, cls).setUpClass()

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

        gdf_points, gdf_pts_crs = convert.convert_csv_to_points(in_csv, coord_columns_epsg=4326, out_epsg=28354)
        out_gdf, out_crs = clean_trim_points(gdf_points, gdf_pts_crs, 'Yld Mass(Dry)(tonne/ha)', out_csv,
                                             out_keep_shapefile=out_shp, out_removed_shapefile=out_rm_shp,
                                             boundary_polyfile=in_poly, thin_dist_m=2.5)

        self.assertTrue(os.path.exists(out_csv))
        self.assertTrue(gdf_pts_crs, out_crs)
        self.assertEqual(out_gdf.crs, {'init': 'epsg:28354', 'no_defs': True})
        self.assertEqual(len(out_gdf), 554)

    def test_createPolygonFromPointTrail(self):
        in_csv = os.path.join(this_dir + "/data/area2_yield_ISO-8859-1.csv")

        # out_ptsfile = os.path.join(TmpDir, os.path.splitext(os.path.basename(in_csv))[0] + '_points.shp')
        out_polyfile = os.path.join(TmpDir, os.path.splitext(os.path.basename(in_csv))[0] + '_poly.shp')

        gdf_points, gdf_pts_crs = convert.convert_csv_to_points(in_csv, None, coord_columns_epsg=4326, out_epsg=28354)

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
        out_shp = os.path.join(TmpDir, os.path.basename(raster_file).replace('.tif', '_randpts.shp'))
        rast_crs = pyprecag_crs.getCRSfromRasterFile(raster_file)

        with rasterio.open(os.path.normpath(raster_file)) as raster:
            rand_gdf, rand_crs = random_pixel_selection(raster, rast_crs, 50, out_shapefile=out_shp)

        self.assertEqual(len(rand_gdf), 50)
        self.assertTrue(os.path.exists(out_shp))
        self.assertEqual(rand_crs, rast_crs)

    def test_kmeansCluster(self):
        raster_files = [os.path.realpath(this_dir + '/data/area1_onebox_NDRE_250cm.tif'),
                        os.path.realpath(this_dir + '/data/area1_rgbi_jan_50cm_84sutm54.tif'),
                        os.path.realpath(this_dir + '/data/area1_onebox_NDVI_250cm.tif'),
                        os.path.realpath(this_dir + '/data/area1_onebox_PCD_250cm.tif'),
                        os.path.realpath(this_dir + '/data/area1_onebox_Yield_250cm.tif')
                        ]

        out_img = os.path.join(TmpDir, 'kmeans-cluster_3cluster_3rasters.tif')
        with self.assertRaises(TypeError) as msg:
            _ = kmeans_clustering(raster_files, out_img)
            self.assertEqual("Pixel Sizes Don't Match - [(0.5, 0.5), (2.5, 2.5)]", str(msg.exception))
        raster_files.pop(0)

        with self.assertRaises(TypeError) as msg:
            _ = kmeans_clustering(raster_files, out_img)
            self.assertEqual("1 raster(s) don't have coordinates systems assigned", str(msg.exception))
        raster_files.pop(0)

        out_df = kmeans_clustering(raster_files, out_img)

        self.assertTrue(os.path.exists(out_img))
        self.assertTrue(os.path.exists(out_img.replace('.tif', '_statistics.csv')))
        self.assertEqual(3, len(out_df['zone'].unique()))

        with rasterio.open(out_img) as src:
            self.assertEqual(1, src.count)
            self.assertEqual(src.crs.to_string(), '+init=epsg:28354')
            self.assertEqual(0, src.nodata)
            band1 = src.read(1, masked=True)
            self.assertItemsEqual(np.array([0, 1, 2, 3]), np.unique(band1.data))


class test_extractRasterStatisticsForPoints(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_extractRasterStatisticsForPoints, cls).setUpClass()

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
        out_csv = os.path.join(TmpDir, os.path.basename(raster_file).replace('.tif', '_b1grdext.csv'))
        rast_crs = pyprecag_crs.getCRSfromRasterFile(raster_file)

        with rasterio.open(os.path.normpath(raster_file)) as raster:
            pts_gdf, pts_crs = random_pixel_selection(raster, rast_crs, 50)

            out_gdf, out_crs = extract_pixel_statistics_for_points(pts_gdf, pts_crs, [raster_file],
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
        out_csv = os.path.join(TmpDir, os.path.basename(raster_file).replace('.tif', '_b1grdextwgs84.csv'))
        rast_crs = pyprecag_crs.getCRSfromRasterFile(raster_file)

        with rasterio.open(os.path.normpath(raster_file)) as raster:
            pts_gdf, pts_crs = random_pixel_selection(raster, rast_crs, 50)

        pts_crs.getFromEPSG(4326)
        pts_gdf.to_crs(epsg=4326, inplace=True)

        out_gdf, out_crs = extract_pixel_statistics_for_points(pts_gdf, pts_crs, [raster_file],
                                                               function_list=[np.nanmean, raster_ops.nancv],
                                                               size_list=[1, 3, 7], output_csvfile=out_csv)

        self.assertTrue(os.path.exists(out_csv))
        self.assertTrue(len(out_gdf), len(pts_gdf))
        self.assertTrue(len(pts_gdf.columns), len(out_gdf.columns) + 5)
        self.assertTrue(out_crs, rast_crs)
        self.assertTrue(pts_crs, out_crs)
        self.assertEqual(out_gdf['mean7x7_test_singleband_94mga54'].isnull().sum(), 0)


class test_CalculateImageIndices(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_CalculateImageIndices, cls).setUpClass()
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

        out_dir = os.path.join(TmpDir, 'test_CalculateImageIndices', 'all-opts')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/area1_rgbi_jan_50cm_84sutm54.tif')
        polygon_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        indices = CalculateIndices(red=3, infrared=4, mask=5).valid_indices()
        files = calc_indices_for_block(image_file, 2,
                                       BandMapping(red=3, green=2, infrared=4, rededge=1, mask=5),
                                       out_dir, indices=indices, image_epsg=32754, image_nodata=0,
                                       polygon_shapefile=polygon_shapefile,
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
        out_dir = os.path.join(TmpDir, 'test_CalculateImageIndices', 'no-nonvine')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/area1_rgbi_jan_50cm_84sutm54.tif')
        polygon_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        indices = CalculateIndices(red=3, infrared=4, mask=5).valid_indices()
        files = calc_indices_for_block(image_file, 2,
                                       BandMapping(red=3, green=2, infrared=4, rededge=1),
                                       out_dir, indices=indices, image_epsg=32754, image_nodata=0,
                                       polygon_shapefile=polygon_shapefile,
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
        out_dir = os.path.join(TmpDir, 'test_CalculateImageIndices', 'no-shapefile')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/area1_rgbi_jan_50cm_84sutm54.tif')
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
        out_dir = os.path.join(TmpDir, 'test_CalculateImageIndices', 'no-groupby')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        bm = BandMapping(green=2, infrared=4, rededge=1, mask=5)

        image_file = os.path.realpath(this_dir + '/data/area1_rgbi_jan_50cm_84sutm54.tif')
        polygon_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        indices = CalculateIndices(**bm).valid_indices()
        files = calc_indices_for_block(image_file, 2,
                                       BandMapping(red=3, green=2, infrared=4, rededge=1),
                                       out_dir, indices=indices, image_epsg=32754, image_nodata=0,
                                       polygon_shapefile=polygon_shapefile,
                                       out_epsg=28354)

        self.assertEqual(len(indices), 3)
        self.assertEqual(len(files), len(indices))
        with rasterio.open(files[0]) as src:
            self.assertEqual(src.nodata, -9999)
            self.assertEqual(src.crs, rasterio.crs.CRS.from_epsg(28354))
            self.assertEqual(src.dtypes, ('float32',))
            self.assertEqual(src.res, (2.0, 2.0))
            self.assertEqual(src.count, 1)


class test_ResampleToBlock(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # 'https://stackoverflow.com/a/34065561'
        super(test_ResampleToBlock, cls).setUpClass()

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

        out_dir = os.path.join(TmpDir, 'test_ResampleToBlock', 'all-opts')

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/area1_rgbi_jan_50cm_84sutm54.tif')
        polygon_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        files = resample_bands_to_block(image_file, 2,
                                        out_dir, band_nums=[6], image_epsg=32754, image_nodata=0,
                                        polygon_shapefile=polygon_shapefile,
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
        out_dir = os.path.join(TmpDir, 'test_ResampleToBlock', 'no-shapefile')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/area1_rgbi_jan_50cm_84sutm54.tif')
        files = resample_bands_to_block(image_file, 2, out_dir,
                                        band_nums=[6], image_epsg=32754, image_nodata=0, out_epsg=28354)

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
        out_dir = os.path.join(TmpDir, 'test_ResampleToBlock', 'no-groupby')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        image_file = os.path.realpath(this_dir + '/data/area1_rgbi_jan_50cm_84sutm54.tif')
        polygon_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        files = resample_bands_to_block(image_file, 2,
                                        out_dir, band_nums=[6], image_epsg=32754, image_nodata=0,
                                        polygon_shapefile=polygon_shapefile,
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

        out_dir = os.path.join(TmpDir, 'test_ResampleToBlock', 'nonstandard-nodata-notset')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # change nodata to 7777 but don't set it in the output
        new_image = os.path.join(out_dir, 'geotif_32754_no-nodata-7777.tif')
        crs_sutm54 = rasterio.crs.CRS.from_epsg(32754)

        image_file = os.path.realpath(this_dir + '/data/area1_rgbi_jan_50cm_84sutm54.tif')
        with rasterio.open(image_file) as src:
            meta = src.meta.copy()
            del meta['crs']

            with rasterio.open(new_image, 'w', crs=crs_sutm54, **meta) as dest:
                for i in range(1, src.count + 1):
                    data = src.read(i)
                    np.putmask(data, data == 0, 7777)
                    dest.write(data, i)

        polygon_shapefile = os.path.realpath(this_dir + '/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp')
        files = resample_bands_to_block(new_image, 2.5,
                                        out_dir, band_nums=[6], image_epsg=32754, image_nodata=7777,
                                        polygon_shapefile=polygon_shapefile,
                                        out_epsg=28354)

        self.assertEqual(len(files), 1)

        with rasterio.open(files[0]) as src:
            self.assertEqual(src.nodata, 7777)
            self.assertEqual(src.crs, rasterio.crs.CRS.from_epsg(28354))
            self.assertEqual(src.dtypes, ('float32',))
            self.assertEqual(src.res, (2.5, 2.5))
            self.assertEqual(src.count, 1)
