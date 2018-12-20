"""  Custom Error Types """


class DriverError(ValueError):
    """Encapsulates unsupported driver and driver mode errors."""


class GeometryError(ValueError):
    """Encapsulates unexpected geometry type.
    eg. Expects Polygon Geometry, receives Line geometry"""


class SpatialReferenceError(ValueError):
    """Encapsulates spatial reference Error
        eg. expecting projected CRS but receives geographic
    """


class TransformationError(ValueError):
    """Encapsulates errors when projecting features
        eg. expecting projected CRS but receives geographic
    """


class OutOfRangeError(ValueError):
    pass
