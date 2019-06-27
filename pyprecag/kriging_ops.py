import collections
import datetime
import glob
import inspect
import logging
import os
import platform
import re
import shutil
import subprocess
import time
import warnings

import geopandas as gpd
import pandas as pd
import rasterio
from rasterio import features

from unidecode import unidecode
from collections import OrderedDict

from . import config
from .convert import add_point_geometry_to_dataframe, numeric_pixelsize_to_string
from .describe import predictCoordinateColumnNames
from .raster_ops import raster_snap_extent, create_raster_transform

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())
# DEBUG = config.get_debug_mode()  # LOGGER.isEnabledFor(logging.DEBUG)

# Set default value for vesper_exe
vesper_exe = r"C:\Program Files (x86)\Vesper\Vesper1.6.exe"


def test_for_windows():
    if platform.system() != 'Windows' and os.path.exists(vesper_exe):
        raise IOError("Kriging currently only available on Windows with Vesper installed")


# consider enums - https://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
# The codes for the variogram. The key is the tag used in the control file_csv
VESPER_OPTIONS = {"jigraph": {"Don't Show": 0, "Show": 1},
                  "jimap": {"Don't Show": 0, "Show": 1},
                  "jlockrg": {"Local": 1, "Global": 0},
                  "jpntkrg": {"Block": 0, "Point": 1, "Punctual": 1},
                  "jsetrad": {"Calculate Radius": 0, "Set Radius": 1},
                  "jcomvar": {"Define Variogram Parameter": 0, "Compute Variogram": 1},
                  "modtyp": {"Spherical": 1, "Exponential": 2, "Gaussian": 3,
                             "Linear with sill": 4, "Stable": 5, "Generalised Cauchy": 6,
                             "Matern": 7, "Double spherical": 8, "Double exponential": 9},
                  "iwei": {"unity": 0, "No. of pairs": 1, "1/variance": 2, "no_pairs/variance": 3,
                           "no_pairs/fitted": 4}}


class VesperControl(collections.MutableMapping, dict):
    """A dictionary used to manage vesper control keys and values.

    The list of keys is confined to those in __defaults, values are checked against the default
    and must be the same type.

    Attributes:
        __defaults = values to use as defaults. Defaults are set for type as well ie 0.00 = float
        key_order = This doubles as the list of allowed keys and the order required when writing
               to file with the exception of epsg.
    """

    __defaults = dict(ivers=16121, title='kriging of control.txt configured by pyprecag epsg=',
                      datfil='vesperdata.csv', gridfile='vespergrid.txt', outdir='', epsg=0,
                      repfil='report.txt', outfil='kriged.txt', parfil='parameter.txt', numcol=4,
                      icol_x=1, icol_y=2, icol_z=3, jordkrg=1, jpntkrg=0, jlockrg=1, nest=10,
                      dstinc=1, valmis=-9999, jsetint=0, xlint=0, xhint=0, ylint=0, yhint=0,
                      jsetrad=0, radius=100, minpts=90, maxpts=100, sigsqr=0, isomod=1, modtyp=2,
                      isearch=0, igeos=0, icircs=0, phi=0, psin=0, pcos=0, jcomvar=1, nlag=20,
                      hmax=0, tolag=10, iwei=1, jigraph=0, jimap=0,
                      CO=0.0, C1=1.0, A1=10.0, C2=1.0, A2=1.0, Alfa=1.0,
                      xside=10, yside=10, lognorm=0, itrend=0, iconvex=0, igrids=1)

    # key_order is the order the tags will be written to file. epsg is excluded and written first
    key_order = ['ivers', 'title', 'datfil', 'outdir', 'repfil', 'outfil', 'parfil',
                 'igrids', 'gridfile', 'jigraph', 'jimap', 'numcol', 'icol_x', 'icol_y', 'icol_z',
                 'valmis', 'jpntkrg', 'xside', 'yside',
                 'nest', 'dstinc', 'jsetint', 'xlint', 'xhint', 'ylint', 'yhint',
                 'jsetrad', 'radius', 'minpts', 'maxpts', 'sigsqr', 'lognorm', 'itrend', 'iconvex',
                 'jlockrg', 'jcomvar', 'nlag', 'tolag', 'hmax', 'iwei',
                 'modtyp', 'CO', 'C1', 'A1', 'C2', 'A2', 'Alfa', 'jordkrg',
                 'isomod', 'isearch', 'igeos', 'icircs', 'phi', 'psin', 'pcos']

    def __init__(self, *args, **kwargs):
        self.update(**self.__defaults)
        self.update(*args, **kwargs)

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        """ check if key is from a predefined list and the value of the correct type

        epsg is an non-vesper key so it is excluded from the keyorder and written at top of file
        """

        if key not in self.key_order + ['epsg']:
            raise AttributeError('VesperControl has no attribute {}'.format(key))

        if key in VESPER_OPTIONS.keys():
            if isinstance(value, basestring):
                try:
                    value = VESPER_OPTIONS[key][value]
                except KeyError:
                    raise ValueError('{} is an invalid option for {}. Options are {}'.format(
                        value, key, VESPER_OPTIONS[key]))

            elif isinstance(value, (int, float, long)):
                if value not in VESPER_OPTIONS[key].values():
                    raise ValueError('{} is an invalid option for {}. Options are {}'.format(
                        value, key, VESPER_OPTIONS[key]))

        if isinstance(self.__defaults[key], float):
            if not isinstance(value, (int, long, float)):
                raise ValueError('value for key {} should be a float or int - Got {}'.format(
                    key, type(value.__name__)))

        elif type(self.__defaults[key]) != type(value):
            raise ValueError('value for key {} should be a {} - Got {}'.format(
                key, type(self.__defaults[key]).__name__, type(value).__name__))

        dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        dict.__setitem__(self, key, self.__defaults[key])

    def __iter__(self):
        return dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __contains__(self, x):
        return dict.__contains__(self, x)

    def updated_keys(self):
        updated_keys = dict(set(self.items()) - set(self.__defaults.items()))
        return updated_keys

    def write_to_file(self, output_file):
        if output_file == '':
            raise ValueError('Output file not specified')

        if not os.path.exists(os.path.dirname(output_file)):
            raise ValueError('Output folder {} does not exist'.format(os.path.dirname(output_file)))

        with open(output_file, "w") as w_out:
            # add epsg tag outside the standard vesper file denoted by $vsl and $end
            w_out.write('epsg={}\n'.format(self['epsg']))
            w_out.write('$vsl\n')
            for key in self.key_order:
                # need to quote the strings
                if key in ['title', 'datfil', 'gridfile', 'outdir', 'repfil', 'outfil', 'parfil']:
                    w_out.write("{}='{}'\n".format(key, self[key]))
                else:
                    w_out.write('{}={}\n'.format(key, self[key]))
            w_out.write('$end\n')


def vesper_text_to_raster(control_textfile, krig_epsg=0, nodata_value=-9999):
    """Convert an vesper kriged text file output to a prediction and a standard error (SE) tif
    raster, and create a confidence interval (CI) metadata file. If the output files already
    exists, they will be overwritten.

    Kriging using Vesper creates 3 new files:
        *_kriged_*.txt, *_report_*.txt, *_parameter_*.txt.

      There is a 100 character file path limitation set within vesper for the kriged and
      report files. By using the control filename as the input, the truncated files can be renamed
      to their correct names and is used as a base for derived tiff and metadata

    The *_kriged_*.txt file contains the data to create the desired outputs.

    Args:
        control_textfile (str): The control file for the associated vesper krige raster
                               (see note above) to convert.
        krig_epsg (int): The EPSG number to assign to the output raster
        nodata_value (int): The value to assign to nodata pixels.

    Outputs:
        *_CI_*.txt      The confidence interval (CI) metadata file.
        *_PRED_*.tif    The kriged column prediction tif
        *_SE_*.tif      The kriged column standard error (SE) tif

    """

    if not os.path.exists(control_textfile):
        raise IOError("Invalid path: {}".format(control_textfile))

    for argCheck in [('krig_epsg', krig_epsg), ('nodata_value', nodata_value)]:
        if not isinstance(argCheck[1], (int, long)):
            raise TypeError('{} must be a integer.'.format(argCheck[0]))

    start_time = time.time()

    # There is a 100 character file path limitation set for the kriged and report outputs from
    #  vesper. By using the control filename we can find these truncated files and correct names.
    if len(control_textfile) > 100:
        search_dir = os.path.dirname(control_textfile)
        search_file = os.path.basename(control_textfile[:101]).replace('control', '*')
        suffix = control_textfile[101:]
        for ea_file in glob.glob(os.path.join(search_dir, search_file)):
            # only rename if there is no extension on the file
            if os.path.splitext(ea_file)[-1] == '':
                os.rename(ea_file, ea_file + suffix)
                logging.debug('Renaming file {} to {}'.format(ea_file, ea_file + suffix))

    krige_textfile = control_textfile.replace('control', 'kriged')
    out_ci_txt = control_textfile.replace('control', 'CI')

    df_krige = pd.read_csv(krige_textfile, delim_whitespace=True)

    median_val = df_krige['SE_Pred'].median()
    with open(out_ci_txt, 'w') as ci_file:
        ci_file.writelines("Median Prediction SE    : {:.5f}\n".format(median_val))
        ci_file.writelines("95% Confidence Interval : {:.5f}\n\n".format(2 * 1.96 * median_val))
        ci_file.writelines(
            "Date/time : " + datetime.datetime.now().strftime("%d-%b-%Y %H:%M") + "\n")
        userhome = os.path.expanduser('~')
        ci_file.writelines("Username  : " + os.path.split(userhome)[-1] + "\n")

    LOGGER.debug("CI File contents : \n\tMedian Prediction SE    : {:.5f}\n"
                 "\t95% Confidence Interval : {:.5f}".format(median_val, 2 * 1.96 * median_val))

    x_field, y_field = predictCoordinateColumnNames(df_krige.columns.tolist())
    gdf_krig, gdf_crs = add_point_geometry_to_dataframe(df_krige, [x_field, y_field], krig_epsg)

    cellsize_x = float(df_krige[x_field].sort_values().drop_duplicates().diff(1).mode())
    cellsize_y = float(df_krige[y_field].sort_values().drop_duplicates().diff(1).mode())
    pixel_size = min(cellsize_x, cellsize_y)

    pixel_size_str = numeric_pixelsize_to_string(pixel_size)

    LOGGER.debug("Cellsize: {}  X: {}  Y: {}".format(pixel_size, cellsize_x, cellsize_y))
    out_se_tif = control_textfile.replace('control', 'SE_{}'.format(pixel_size_str)).replace('.txt',
                                                                                             '.tif')
    out_pred_tif = control_textfile.replace('control', 'PRED_{}'.format(pixel_size_str)).replace(
        '.txt', '.tif')

    # create an affine transformation matrix to associate the array to the coordinates.
    transform, x_cols, y_rows, new_bbox = create_raster_transform(list(gdf_krig.total_bounds),
                                                                  pixel_size=pixel_size)
    LOGGER.debug('Width (xCols):     {}   Height (yRows):     {}'.format(x_cols, y_rows))

    # create the two tifs and populate with data. This method is by far the quickest of all
    # 3 methods trialed.
    with rasterio.open(os.path.normpath(out_pred_tif), 'w', driver='GTiff',
                       width=x_cols, height=y_rows, count=1,
                       crs=rasterio.crs.CRS.from_epsg(krig_epsg),
                       dtype='float32', nodata=-9999, transform=transform) as out_pred:

        # uses points to burn values into the corresponding raster pixel.
        shapes = ((geom, value) for geom, value in zip(gdf_krig['geometry'], gdf_krig['Predicted']))
        burned = features.rasterize(shapes=shapes, out_shape=(y_rows, x_cols),
                                    transform=out_pred.transform, fill=nodata_value)

        out_pred.write(burned, indexes=1)

    with rasterio.open(os.path.normpath(out_se_tif), 'w', driver='GTiff', width=x_cols,
                       height=y_rows, count=1,
                       crs=rasterio.crs.CRS.from_epsg(krig_epsg), dtype='float32', nodata=-9999,
                       transform=transform) as out_se:

        # uses points to burn values into the corresponding raster pixel.
        shapes = ((geom, value) for geom, value in zip(gdf_krig['geometry'], gdf_krig['SE_Pred']))
        burned = features.rasterize(shapes=shapes, out_shape=(y_rows, x_cols),
                                    transform=out_se.transform, fill=nodata_value)
        out_se.write(burned, indexes=1)

    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=datetime.timedelta(
                                                   seconds=time.time() - start_time)))
    return out_pred_tif, out_se_tif, out_ci_txt


def prepare_for_vesper_krige(in_dataframe, krig_column, grid_filename, out_folder,
                             control_textfile='', coord_columns=[], epsg=0, display_graphics=False,
                             control_options=VesperControl(),
                             block_size=10,
                             vesper_exe=vesper_exe):
    """Prepare data for vesper kriging and create a windows batch file to run outside the
    python/pyprecag environment.

    The input arguments will be used to update the control_options.

    A description of the keys and their options can be found in the VESPER user manual at
    https://sydney.edu.au/agriculture/pal/documents/Vesper_1.6_User_Manual.pdf

    block_size will be replaced by the xside and yside keys in the VesperControl Object.

    Outputs:  The following files will be added to the vesper sub-folder in the output folder.
        *_control_*.txt  - The vesper control file
        *_vesperdata_*.csv - The subset of data to krige. All non-required columns are deleted.
        Do_Vesper.bat - A windows batch file to launch vesper kriging for all control files

        Any files matching the derived kriged output names will be deleted to ensure that the
        files in the folder belong to the newly create control file. If old version are to
        be kept, manually rename first.

    DEPENDANCY: Vesper must already be installed on the PC and the install directory added
    to the config.json file.

    Args:

        in_dataframe (geopandas.geodataframe.GeoDataFrame, pandas.core.frame.DataFrame):
        krig_column (str): The column containing the data to krige.
        grid_filename (str): The vesper grid file.
        control_textfile (str): The name of the control text file without the path
        out_folder (str): The folder to add outputs too. A 'Vesper' sub directory will be created
        block_size (int): The size to apply for block Kriging. Units are from the coordinate system.
        coord_columns (List): The columns representing the X and Y coordinates.
        epsg (int) : The epsg_number number for the data. If 0 the en_epsg or
                     enepsg column (if exists) will be used.
        display_graphics (bool): Option to display graphics while running vesper kriging.
        vesper_exe (str): The path for the location of the Vesper executable
        control_options (pyprecag.kriging_ops.VesperControl): Vesper control settings parameters
    Returns:
       vesper_batfile, vesper_ctrlfile: The paths to the generated batch file and control file.
    """

    # Vesper only works with Windows
    # test_for_windows()

    if not isinstance(in_dataframe, (gpd.GeoDataFrame, pd.DataFrame)):
        raise TypeError('Invalid input data :in_dataframe')

    if not os.path.exists(grid_filename):
        raise IOError("Invalid path: {}".format(grid_filename))

    if not isinstance(block_size, (int, long)):
        raise TypeError('block_size must be an integer'.format(block_size))

    warnings.warn('block_size is deprecated, use VesperControl xsize and ysize instead',
                  PendingDeprecationWarning)

    if out_folder.strip() in [None, '']:
        raise TypeError('Please specify an output folder')

    if not os.path.exists(os.path.dirname(out_folder)):
        raise IOError('Output directory {} does not exist'.format(os.path.dirname(out_folder)))

    if not isinstance(coord_columns, list):
        raise TypeError('Coordinate columns should be a list.'.format(coord_columns))

    if len(coord_columns) == 0:
        coord_columns = predictCoordinateColumnNames(in_dataframe.columns)

    for eaFld in [krig_column] + coord_columns:
        if eaFld not in in_dataframe.columns:
            raise TypeError('Column {} does not exist'.format(eaFld))

    if not isinstance(epsg, (int, long)):
        raise TypeError('EPSG {} must be a integer.'.format(epsg))

    if not isinstance(control_options, VesperControl):
        raise TypeError('control_options must of type VesperControl')

    start_time = time.time()

    # Create a filename compatible copy of the krig_column
    krig_col_file = re.sub('[^A-Za-z0-9_-]+', '', unidecode(unicode(krig_column)))

    out_sub_name = os.path.basename(grid_filename)[:20]

    if control_textfile.strip() in [None, '']:
        control_textfile = "{}_control_{}.txt".format(out_sub_name, krig_col_file)
        data_file = "{}_vesperdata_{}.csv".format(out_sub_name, krig_col_file)
        grid_file = "{}_vespergrid_{}.txt".format(out_sub_name, krig_col_file)
        report_file = "{}_report_{}.txt".format(out_sub_name, krig_col_file)
        kriged_file = "{}_kriged_{}.txt".format(out_sub_name, krig_col_file)
        param_file = "{}_parameter_{}.txt".format(out_sub_name, krig_col_file)

    elif 'control' in control_textfile:
        data_file = os.path.splitext(control_textfile)[0].replace('control', 'vesperdata') + '.csv'
        grid_file = os.path.splitext(control_textfile)[0].replace('control', 'vespergrid') + '.txt'
        report_file = os.path.splitext(control_textfile)[0].replace('control', 'report') + '.txt'
        kriged_file = os.path.splitext(control_textfile)[0].replace('control', 'kriged') + '.txt'
        param_file = os.path.splitext(control_textfile)[0].replace('control', 'parameter') + '.txt'

    else:
        data_file = os.path.splitext(control_textfile)[0] + '_vesperdata.csv'
        report_file = os.path.splitext(control_textfile)[0] + '_report.txt'
        kriged_file = os.path.splitext(control_textfile)[0] + '_kriged.txt'
        param_file = os.path.splitext(control_textfile)[0] + '_parameter.txt'
        grid_file = os.path.splitext(control_textfile)[0] + '_vespergrid.txt'

    vesper_outdir = out_folder
    if not os.path.basename(out_folder) == 'Vesper':
        vesper_outdir = os.path.join(out_folder, 'Vesper')

    vesper_ctrlfile = os.path.join(vesper_outdir, control_textfile)
    vesper_datafile = os.path.join(vesper_outdir, data_file)
    vesper_gridfile = os.path.join(vesper_outdir, grid_file)

    if not os.path.exists(vesper_outdir):
        os.mkdir(vesper_outdir)

    # check Permission Denied or WindowsError: [Error 32]) errors which may indicate that the files
    # in use either by vesper or another app.

    # Start always start with the control file. this is what gets used within QGIS
    files_list = glob.glob(
        os.path.join(vesper_outdir, "{}_*_{}.*".format(out_sub_name, krig_col_file)))
    if vesper_ctrlfile in files_list:
        files_list.insert(0, files_list.pop(files_list.index(vesper_ctrlfile)))

    for eaFile in files_list:
        try:
            with open(eaFile, "a") as f:
                pass
        except IOError as e:
            raise IOError('File(s) in use - {}'.format(eaFile))

    # if control file already exists then replace it, and delete all matching kriging
    # outputs and tiff files.
    for eaFile in glob.glob(
            os.path.join(vesper_outdir, "{}_*_{}.*".format(out_sub_name, krig_col_file))):
        os.remove(eaFile)
        LOGGER.debug('Deleted file {}'.format(eaFile))

    x_field, y_field = coord_columns

    keep_cols = [findCol for findCol in ['EN_EPSG', 'ENEPSG', 'EPSG']
                 if findCol in in_dataframe.columns]
    if epsg == 0 and len(keep_cols) > 0:
        for col in keep_cols:
            if in_dataframe.iloc[0][col] > 0:
                epsg = in_dataframe.iloc[0][col]
                break

    keep_cols += [x_field, y_field, krig_column]

    df_csv = in_dataframe[keep_cols].copy()
    df_csv.to_csv(vesper_datafile, index=False)

    shutil.copy2(grid_filename, vesper_gridfile)

    x_field, y_field = coord_columns
    i_x_column = df_csv.columns.tolist().index(x_field) + 1
    i_y_column = df_csv.columns.tolist().index(y_field) + 1
    i_k_column = df_csv.columns.tolist().index(krig_column) + 1

    # update the control_options for the new csv file
    control_options.update(epsg=epsg,
                           title="kriging of {} configured by pyprecag".format(control_textfile),
                           datfil="{}".format(data_file),
                           gridfile="{}".format(grid_file),
                           outdir="",  # blank writes to control file_csv folder
                           repfil="{}".format(report_file),
                           outfil="{}".format(kriged_file),
                           parfil="{}".format(param_file),
                           numcol=len(df_csv.columns),
                           icol_x=i_x_column,
                           icol_y=i_y_column,
                           icol_z=i_k_column,
                           jigraph=int(display_graphics),  # 1, show graph, otherwise 0
                           jimap=int(display_graphics),  # 1, show map, otherwise 0)
                           )

    # high density kriging. fix for pre VesperControl Class
    if block_size != 10:  # only update if not the default variable value
        control_options.update({"xside": block_size,
                                "yside": block_size})

    # ---------------------------------------------------------------------------------------------
    # Write a VESPER control file
    control_options.write_to_file(vesper_ctrlfile)
    if not os.path.exists(vesper_exe):
        vesper_batfile = ''
        warnings.warn('Vesper*.exe at "{exe}" does not exist. '
                      'Batch file not created'.format(exe=vesper_exe))

    else:
        vesper_batfile = os.path.join(vesper_outdir, "Do_Vesper.bat")

        bat_file_string = ("@echo off\n"
                           # "REM enable delayed variable expansion, t\n"
                           "setlocal enabledelayedexpansion\n"
                           "echo Processing control files in %CD%\n"
                           "echo.\n"
                           "set icount=0\n"
                           "set tcount=0\n"
                           "\nREM count the number of control files\n"
                           "for %%x in (*control*.txt) do set /a tcount+=1'\n"
                           "echo Found %tcount% control file_csv(s) to process....\n"
                           "echo.\n"
                           # suppress newline char and enable 'FINISHED' to be added
                           # setlocal enabledelayedexpansion and ! instead of % will
                           # evaluate a variable when they appear instead of when the
                           # batch is parsed
                           "FOR %%f IN (*control*) DO (\n"
                           "   set /a icount+=%icount%+1\n"
                           "    <nul set /p mystr=Kriging !icount! of %tcount% - %%~nxf \n"
                           '   start "" /HIGH /WAIT /SHARED \"{exe}\" %%~nxf\n'
                           "   echo - Finished \n"
                           ")\n"  # close bracket for for /f do loop
                           # "echo.\n"
                           # "pause\n"
                           .format(exe=vesper_exe))

        with open(vesper_batfile, 'w') as wBatFile:
            wBatFile.write(bat_file_string)

    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=datetime.timedelta(
                                                   seconds=time.time() - start_time)))
    return vesper_batfile, vesper_ctrlfile


def run_vesper(ctrl_file, bMinimiseWindow=True, vesper_exe=vesper_exe):
    """
    Run VESPER for selected control file.

    By default the current working directory will be set to the folder containing the control file.

    By default processing will wait for the current VESPER control file to complete.
    Args:
        ctrl_file (str):  the control file to run as VESPER argument.
        bMinimiseWindow (bool):  Option to automatically minimise the VESPER window on launch.
        vesper_exe (str): The path for the location of the Vesper executable
    """

    # Vesper only works with Windows
    test_for_windows()

    if not os.path.exists(vesper_exe):
        raise IOError('Vesper*.exe at "{}"'.format(vesper_exe)
                      + ' does not exist. Please install and configure for kriging to occur')

    task_time = time.time()

    info = subprocess.STARTUPINFO()
    if bMinimiseWindow:
        # run vesper minimised.
        info.dwFlags = subprocess.STARTF_USESHOWWINDOW
        info.wShowWindow = 6  # SW_MINIMIZE

    process = subprocess.Popen([vesper_exe, ctrl_file], cwd=os.path.dirname(ctrl_file),
                               startupinfo=info)
    process.wait()
    logging.info('{:<30}\t{dur:<15}'.format('Vesper Kriging', dur=datetime.timedelta(
        seconds=time.time() - task_time)))
