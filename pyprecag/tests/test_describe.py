from unittest import TestCase
import time
import os
from pyprecag.describe import predictCoordinateColumnNames, CsvDescribe, VectorDescribe

this_dir = os.path.abspath(os.path.dirname(__file__))

class general_describe(TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_predictCoordinateFieldnames(self):
        self.assertEquals(
            predictCoordinateColumnNames(['longitude', 'latitude', 'area', 'weight', 'tonnes/ha', 'Air Temp(degC)']),
            ['longitude', 'latitude'])
        self.assertEquals(predictCoordinateColumnNames(['northing', 'name', 'easting', ' test type']),
                          ['easting', 'northing'])
        self.assertEquals(predictCoordinateColumnNames(['name', 'lon_dms', 'lat_dms', ' test type']),
                          ['lon_dms', 'lat_dms'])
        self.assertEquals(predictCoordinateColumnNames(['row', 'column', 'name', ' test type']), [None, None])

        self.assertEquals(
            predictCoordinateColumnNames(['northing', 'easting', 'name', ' test type', 'longitude', 'latitude']),
            ['longitude', 'latitude'])

    def test_from_ISO_8859_1_csv(self):
        file = os.path.realpath(this_dir + "/data/area2_yield_file_ISO-8859-1.csv")
        descCSV = CsvDescribe(file)
        self.assertEquals(predictCoordinateColumnNames(descCSV.get_column_names()), ['Longitude', 'Latitude'])


class TestVectorDescribe_QGIS(TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_wgs84_mixedPartLine_MZ_qgisprj(self):
        vDesc = VectorDescribe(os.path.realpath(this_dir + "/data/LineMZ_wgs84_MixedPartFieldsTypes_exportedqgis.shp"))
        self.assertEqual(vDesc.crs.epsg_number, 4326)
        self.assertTrue(vDesc.is_mz_aware)
        self.assertEquals(vDesc.geometry_type, 'MultiLineString')

    def test_mga_singlePartPoly_qgisprj(self):
        vDesc = VectorDescribe(os.path.realpath(this_dir + "/data/Poly_mga54_SinglePartFieldsTypes_qgis-prj.shp"))
        self.assertEqual(vDesc.crs.epsg_number, 28354)
        self.assertFalse(vDesc.is_mz_aware)
        self.assertEquals(vDesc.geometry_type, 'Polygon')


class TestVectorDescribe_ESRI(TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_noprojection_mixedPartPoint_MZ(self):
        vDesc = VectorDescribe(os.path.realpath(this_dir + "/data/PointMZ_mga54_MixedPartFieldsTypes_noPrj.shp"))
        self.assertIsNone(vDesc.crs.srs)
        self.assertIsNone(vDesc.crs.epsg_number)
        self.assertEquals(vDesc.geometry_type, 'MultiPoint')

    def test_wgs84_mixedPartLine_MZ_esriprj(self):
        vDesc = VectorDescribe(os.path.realpath(this_dir + "/data/LineMZ_wgs84_MixedPartFieldsTypes_esri.shp"))
        self.assertEqual(vDesc.crs.epsg_number, 4326)
        self.assertTrue(vDesc.is_mz_aware)
        self.assertEquals(vDesc.geometry_type, 'MultiLineString')

    def test_wgs84_mixedPartPoly_MZ_esriprj(self):
        vDesc = VectorDescribe(os.path.realpath(this_dir + "/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp"))
        self.assertEqual(vDesc.crs.epsg_number, 4326)
        self.assertTrue(vDesc.is_mz_aware)
        self.assertEquals(vDesc.geometry_type, 'MultiPolygon')


class TestCsvDescribe(TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_csvfile_UTF8(self):
        csvDesc = CsvDescribe(os.path.realpath(this_dir + "/data/area2_yield_file_ISO-8859-1.csv"))

        self.assertEquals(csvDesc.file_encoding, 'ISO-8859-1')
        self.assertEquals(csvDesc.row_count, 10000)
        self.assertEquals(csvDesc.column_count, 18)
        self.assertEquals(predictCoordinateColumnNames(csvDesc.get_column_names()), ['Longitude', 'Latitude'])
        self.assertTrue(csvDesc.has_column_header)

    def test_csvfile_ascii(self):
        csvDesc = CsvDescribe(os.path.realpath(this_dir + "/data/area1_yield_file_ascii_wgs84.csv"))

        self.assertEquals(csvDesc.file_encoding, 'ascii')
        self.assertEquals(csvDesc.row_count, 34309)
        self.assertEquals(csvDesc.column_count, 24)

        # see: https://jira.csiro.au/browse/PA-42
        # self.assertEquals(predictCoordinateColumnNames(csvDesc.get_column_names()), ['Long', 'Lat'])
        self.assertTrue(csvDesc.has_column_header)
