import numpy as np
import os
import rasterio


def make_dummy_tif_files(out_dir):
    """ Create fake TIFF to be used as test datasets for other components of pyprecag

    A single band tif and multi band tif will be created. The single band  and Band 1 of the multi-band contains the
    same data. Bands 2 and 3 are created by applying a formula to band 1. They are both for demonstration purposes
    only and represent no specific data and/or location in real life.

    Args:
        out_dir (str): The folder to save the output files too.

    Returns:
        List[str]: A list of TIFF's , the first being a single band tif, and the second a multi band tif.
    """

    # create and apply a transform based on epsg_number:28354
    xmin, ymax = [350110, 6060000]  # UL  coordinate required for transform
    from rasterio.transform import from_origin
    transform = from_origin(xmin, ymax, 1, 1)

    # create a dummy tif in mga54 using gausian as it will look like real data
    # source http://rasterio.readthedocs.io/en/latest/quickstart.html?highlight=new
    x = np.linspace(-4.0, 4.0, 240)
    y = np.linspace(-3.0, 3.0, 180)
    X, Y = np.meshgrid(x, y)
    Z1 = np.exp(-2 * np.log(2) * ((X - 0.5) ** 2 + (Y - 0.5) ** 2) / 1 ** 2)
    Z2 = np.exp(-3 * np.log(2) * ((X + 0.5) ** 2 + (Y + 0.5) ** 2) / 2.5 ** 2)
    Z = 10.0 * (Z2 - Z1)

    # Create a square of nodata values - https://stackoverflow.com/a/10840019
    tr1 = np.triu_indices(min(Z.shape) / 3, -min(Z.shape) / 3)

    # assign nodata to UL Square and set it no np.nan so calculations will exclude it
    Z[tr1] = np.nan

    # rescale to typical yeild values
    smin, smax = [0.01, 7]
    band1 = (Z - np.nanmin(Z)) * (smax - smin) / (np.nanmax(Z) - np.nanmin(Z)) + smin

    # otherbands is just for testing and doesn't mean anything
    band2 = Z
    band3 = Z * 1000 / 850

    # So the output file holds the nodata correctly, change nan's to -9999 nodata value
    band1[np.isnan(band1)] = -9999
    band2[np.isnan(band2)] = -9999
    band3[np.isnan(band3)] = -9999

    in_tif_single = os.path.join(out_dir, 'test_singleband_94mga54.tif')
    with rasterio.open(os.path.normpath(in_tif_single), 'w', crs='EPSG:28354', driver='GTiff',
                       height=Z.shape[0], width=Z.shape[1],
                       dtype=np.float32, count=1, nodata=-9999,
                       transform=transform) as dst:
        dst.write(band1.astype(np.float32), 1)

    in_tif_multi = os.path.join(out_dir, 'test_3band_94mga54.tif')
    with rasterio.open(os.path.normpath(in_tif_multi), 'w', crs='EPSG:28354', driver='GTiff',
                       height=Z.shape[0], width=Z.shape[1],
                       dtype=np.float32, count=3, nodata=-9999,
                       transform=transform) as dst:
        dst.write(band1.astype(np.float32), 1)
        dst.write(band2.astype(np.float32), 2)
        dst.write(band3.astype(np.float32), 3)

    return [in_tif_single, in_tif_multi]
