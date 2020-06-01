from unittest import TestCase
import time
import os

from pkg_resources import parse_version, get_distribution

from pyprecag.describe import predictCoordinateColumnNames, CsvDescribe, VectorDescribe

this_dir = os.path.abspath(os.path.dirname(__file__))


class test_GeneralDescribe(TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_predictCoordinateFieldnames(self):
        self.assertEqual(
            predictCoordinateColumnNames(['longitude', 'latitude', 'area', 'weight', 'tonnes/ha', 'Air Temp(degC)']),
            ['longitude', 'latitude'])
        self.assertEqual(predictCoordinateColumnNames(['northing', 'name', 'easting', ' test type']),
                          ['easting', 'northing'])
        self.assertEqual(predictCoordinateColumnNames(['name', 'lon_dms', 'lat_dms', ' test type']),
                          ['lon_dms', 'lat_dms'])
        self.assertEqual(predictCoordinateColumnNames(['row', 'column', 'name', ' test type']), [None, None])

        self.assertEqual(
            predictCoordinateColumnNames(['northing', 'easting', 'name', ' test type', 'longitude', 'latitude']),
            ['longitude', 'latitude'])

    def test_from_ISO_8859_1_csv(self):
        file = os.path.realpath(this_dir + "/data/area2_yield_ISO-8859-1.csv")
        descCSV = CsvDescribe(file)
        self.assertEqual(predictCoordinateColumnNames(descCSV.get_column_names()), ['Longitude', 'Latitude'])


class test_VectorDescribe_QGIS(TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_wgs84_mixedPartLine_MZ_qgisprj(self):
        vDesc = VectorDescribe(os.path.realpath(this_dir + "/data/LineMZ_wgs84_MixedPartFieldsTypes_exportedqgis.shp"))
        self.assertEqual(vDesc.crs.epsg_number, 4326)
        self.assertTrue(vDesc.is_mz_aware)
        self.assertEqual(vDesc.geometry_type, 'MultiLineString')

    def test_mga_singlePartPoly_qgisprj(self):
        vDesc = VectorDescribe(os.path.realpath(this_dir + "/data/Poly_mga54_SinglePartFieldsTypes_qgis-prj.shp"))
        self.assertEqual(vDesc.crs.epsg_number, 28354)
        self.assertFalse(vDesc.is_mz_aware)
        self.assertEqual(vDesc.geometry_type, 'Polygon')


class test_VectorDescribe_ESRI(TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_noprojection_mixedPartPoint_MZ(self):
        vDesc = VectorDescribe(os.path.realpath(this_dir + "/data/PointMZ_mga54_MixedPartFieldsTypes_noPrj.shp"))
        self.assertIsNone(vDesc.crs.srs)
        self.assertIsNone(vDesc.crs.epsg_number)
        self.assertEqual(vDesc.geometry_type, 'MultiPoint')

    def test_wgs84_mixedPartLine_MZ_esriprj(self):
        vDesc = VectorDescribe(os.path.realpath(this_dir + "/data/LineMZ_wgs84_MixedPartFieldsTypes_esri.shp"))
        self.assertEqual(vDesc.crs.epsg_number, 4326)
        self.assertTrue(vDesc.is_mz_aware)
        self.assertEqual(vDesc.geometry_type, 'MultiLineString')
        self.assertEqual(vDesc.feature_count,8)

        from collections import OrderedDict
        self.maxDiff
        if parse_version(get_distribution('geopandas').version) > parse_version('0.6.0'):
            self.assertDictEqual(vDesc.column_properties, OrderedDict(
                [(u'Id', {'shapefile': u'Id', 'alias': 'Id', 'type': 'int', 'dtype': 'int64'}),
                 (u'float_3dLe',
                  {'shapefile': u'float_3dLe', 'alias': 'float_3dLe', 'type': 'float', 'dtype': 'float64'}),
                 (u'double_Len',
                  {'shapefile': u'double_Len', 'alias': 'double_Len', 'type': 'float', 'dtype': 'float64'}),
                 (u'short_id', {'shapefile': u'short_id', 'alias': 'short_id', 'type': 'int', 'dtype': 'int64'}),
                 (u'Long_vert', {'shapefile': u'Long_vert', 'alias': 'Long_vert', 'type': 'int', 'dtype': 'int64'}),
                 (u'Date_Creat', {'shapefile': u'Date_Creat', 'alias': 'Date_Creat', 'type': 'str', 'dtype': 'object'}),
                 (u'part_type', {'shapefile': u'part_type', 'alias': 'part_type', 'type': 'str', 'dtype': 'object'}),
                 ('geometry', {'shapefile': 'geometry', 'alias': 'geometry', 'type': 'geometry', 'dtype': 'geometry'})]))

        else:
            self.assertDictEqual(vDesc.column_properties,OrderedDict([(u'Id', {'shapefile': u'Id', 'alias': 'Id', 'type': 'int', 'dtype': 'int64'}),
                                                                  (u'float_3dLe', {'shapefile': u'float_3dLe', 'alias': 'float_3dLe', 'type': 'float', 'dtype': 'float64'}),
                                                                  (u'double_Len', {'shapefile': u'double_Len', 'alias': 'double_Len', 'type': 'float', 'dtype': 'float64'}),
                                                                  (u'short_id', {'shapefile': u'short_id', 'alias': 'short_id', 'type': 'int', 'dtype': 'int64'}),
                                                                  (u'Long_vert', {'shapefile': u'Long_vert', 'alias': 'Long_vert', 'type': 'int', 'dtype': 'int64'}),
                                                                  (u'Date_Creat', {'shapefile': u'Date_Creat', 'alias': 'Date_Creat', 'type': 'str', 'dtype': 'object'}),
                                                                  (u'part_type', {'shapefile': u'part_type', 'alias': 'part_type', 'type': 'str', 'dtype': 'object'}),
                                                                  ('geometry', {'shapefile': 'geometry', 'alias': 'geometry', 'type': 'geometry', 'dtype': 'object'})]))

    def test_wgs84_mixedPartPoly_MZ_esriprj(self):
        vDesc = VectorDescribe(os.path.realpath(this_dir + "/data/PolyMZ_wgs84_MixedPartFieldsTypes.shp"))
        self.assertEqual(vDesc.crs.epsg_number, 4326)
        self.assertTrue(vDesc.is_mz_aware)
        self.assertEqual(vDesc.geometry_type, 'MultiPolygon')


class test_CsvDescribe(TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print("%s: %.3f secs" % (self.id(), t))

    def test_csvfile_UTF8(self):
        csvDesc = CsvDescribe(os.path.realpath(this_dir + "/data/area2_yield_ISO-8859-1.csv"))

        self.assertEqual(csvDesc.file_encoding, 'ISO-8859-1')
        self.assertEqual(csvDesc.row_count, 1543)
        self.assertEqual(csvDesc.column_count, 18)
        self.assertEqual(predictCoordinateColumnNames(csvDesc.get_column_names()), ['Longitude', 'Latitude'])
        self.assertTrue(csvDesc.has_column_header)

    def test_csvfile_ascii(self):
        csvDesc = CsvDescribe(os.path.realpath(this_dir + "/data/area1_yield_ascii_wgs84.csv"))

        self.assertEqual(csvDesc.file_encoding, 'ascii')
        self.assertEqual(csvDesc.row_count, 13756)
        self.assertEqual(csvDesc.column_count, 24)
        self.assertTrue(csvDesc.has_column_header)


