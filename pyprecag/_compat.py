import contextlib
from packaging.version import Version
import warnings

import numpy as np
import shapely
import geopandas as gpd

SHAPELY_GE_18 = Version(shapely.__version__) >= Version("1.8")
SHAPELY_GE_182 = Version(shapely.__version__) >= Version("1.8.2")
SHAPELY_GE_20 = Version(shapely.__version__) >= Version("2.0.0.dev0")
SHAPELY_G_20a1 = Version(shapely.__version__) > Version("2.0a1")

GEOPANDAS_GE_10 = Version(gpd.__version__) > Version("0.10.0")
GEOPANDAS_GE_12 = Version(gpd.__version__) > Version("0.12.0")


# compat related to deprecation warnings introduced in Shapely 1.8
# -> creating a numpy array from a list-like of Multi-part geometries,
# although doing the correct thing (not expanding in its parts), still raises
# the warning about iteration being deprecated
# This adds a context manager to explicitly ignore this warning
try:
    from shapely.errors import ShapelyDeprecationWarning as shapely_warning
except ImportError:
    shapely_warning = None

if shapely_warning is not None and not SHAPELY_GE_20:
    @contextlib.contextmanager
    def ignore_shapely2_warnings():
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Iteration|The array interface|__len__", shapely_warning)
            yield

elif (Version(np.__version__) >= Version("1.21")) and not SHAPELY_GE_20:

    @contextlib.contextmanager
    def ignore_shapely2_warnings():
        with warnings.catch_warnings():
            # warning from numpy for existing Shapely releases (this is fixed with Shapely 1.8)
            warnings.filterwarnings("ignore", "An exception was ignored while fetching", DeprecationWarning)
            yield
else:
    @contextlib.contextmanager
    def ignore_shapely2_warnings():
        yield
