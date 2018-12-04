import datetime
import glob
import inspect
import logging
import os
import re
import shutil
import subprocess
import time

import geopandas as gpd
import pandas as pd
import rasterio
from rasterio import features
from unidecode import unidecode

from . import config
from .convert import addPointGeometryToDataframe, numeric_pixelsize_to_string
from .describe import predictCoordinateColumnNames
from .raster_ops import RasterSnapExtent

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())  # Handle logging, no logging has been configured
DEBUG = config.get_config_key('debug_mode')  # LOGGER.isEnabledFor(logging.DEBUG)


def vesperTextToRaster(control_textfile, krig_epsg=0, nodata_value=-9999):
    """Convert an vesper kriged text file output to a prediction and a standard error (SE) tif raster, and create a
    confidence interval (CI) metadata file. If the output files already exists, they will be overwritten.

    Kriging using Vesper creates 3 new files:
        *_kriged_*.txt, *_report_*.txt, *_parameter_*.txt.

      There is a 100 character file path limitation set within vesper for the kriged and report files. By using the
      control filename as the input, the truncated files can be renamed to their correct names and is used as a base for
      derived tiff and metadata

    The *_kriged_*.txt file contains the data to create the desired outputs.

    Args:
        control_textfile (str): The control file for the associated vesper krige raster (see note above) to convert.
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
    # There is a 100 character file path limitation set for the kriged and report outputs from vesper. By using the
    # control filename we can find these truncated files and correct their names.
    if len(control_textfile) > 100:
        search_dir = os.path.dirname(control_textfile)
        search_file = os.path.basename(control_textfile[:101]).replace('control', '*')
        suffix = control_textfile[101:]
        for ea_file in glob.glob(os.path.join(search_dir, search_file)):
            os.rename(ea_file, ea_file + suffix)
            logging.debug('Renaming file {} to {}'.format(ea_file, ea_file + suffix))

    krige_textfile = control_textfile.replace('control', 'kriged')
    out_CITxt = control_textfile.replace('control', 'CI')

    # using pd.read_table takes care of the multiple white space delimiter
    dfKrige = pd.read_table(krige_textfile, delim_whitespace=True)

    median_val = dfKrige['SE_Pred'].median()
    with open(out_CITxt, 'w') as ci_file:
        ci_file.writelines("Median Prediction SE    : {:.5f}\n".format(median_val))
        ci_file.writelines("95% Confidence Interval : {:.5f}\n\n".format(2 * 1.96 * median_val))
        ci_file.writelines("Date/time : " + datetime.datetime.now().strftime("%d-%b-%Y %H:%M") + "\n")
        userhome = os.path.expanduser('~')
        ci_file.writelines("Username  : " + os.path.split(userhome)[-1] + "\n")

    LOGGER.debug("CI File contents : \n\tMedian Prediction SE    : {:.5f}\n"
                 "\t95% Confidence Interval : {:.5f}".format(median_val, 2 * 1.96 * median_val))

    x_field, y_field = predictCoordinateColumnNames(dfKrige.columns.tolist())
    gdfKrig, gdfCRS = addPointGeometryToDataframe(dfKrige, [x_field, y_field], krig_epsg)

    cellsize_x = float(dfKrige[x_field].sort_values().drop_duplicates().diff(1).mode())
    cellsize_y = float(dfKrige[y_field].sort_values().drop_duplicates().diff(1).mode())
    pixel_size = min(cellsize_x, cellsize_y)
    
    pixel_size_str = numeric_pixelsize_to_string(pixel_size)
    
    LOGGER.debug("Cellsize: {}  X: {}  Y: {}".format(pixel_size, cellsize_x, cellsize_y))
    out_SETif = control_textfile.replace('control', 'SE_{}'.format(pixel_size_str)).replace('.txt', '.tif')
    out_PredTif = control_textfile.replace('control', 'PRED_{}'.format(pixel_size_str)).replace('.txt', '.tif')
    

    x_min, y_min, x_max, y_max = RasterSnapExtent(*gdfKrig.total_bounds, pixel_size=pixel_size)

    # get the number of rows/cols
    x_cols = int((x_max - x_min) / pixel_size) + 1
    y_rows = int((y_max - y_min) / pixel_size) + 1
    LOGGER.debug('Width (xCols):     {}   Height (yRows):     {}'.format(x_cols, y_rows))

    # create an affine transformation matrix to associate the array to the coordinates.
    from rasterio.transform import from_origin
    transform = from_origin(x_min, y_max, pixel_size, pixel_size)

    # create the two tifs and populate with data. This method is by far the quickest of all 3 methods trialed.
    with rasterio.open(os.path.normpath(out_PredTif), 'w', driver='GTiff',
                       width=x_cols, height=y_rows, count=1, crs=rasterio.crs.CRS.from_epsg(krig_epsg),
                       dtype='float32', nodata=-9999, transform=transform) as outPred:

        # uses points to burn values into the corresponding raster pixel.
        shapes = ((geom, value) for geom, value in zip(gdfKrig['geometry'], gdfKrig['Predicted']))
        burned = features.rasterize(shapes=shapes, out_shape=(y_rows,x_cols), transform=outPred.transform, fill=nodata_value)

        outPred.write(burned, indexes=1)

    with rasterio.open(os.path.normpath(out_SETif), 'w', driver='GTiff', width=x_cols, height=y_rows, count=1,
                       crs= rasterio.crs.CRS.from_epsg(krig_epsg), dtype='float32', nodata=-9999, transform=transform) as outSE:

        # uses points to burn values into the corresponding raster pixel.
        shapes = ((geom, value) for geom, value in zip(gdfKrig['geometry'], gdfKrig['SE_Pred']))
        burned = features.rasterize(shapes=shapes, out_shape=(y_rows,x_cols), transform=outSE.transform, fill=nodata_value)
        outSE.write(burned, indexes=1)

    # noinspection PyStringFormat
    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=datetime.timedelta(seconds=time.time() - start_time)))
    return out_PredTif, out_SETif, out_CITxt


def prepareForVesperKrig(in_dataframe, krig_column, grid_filename, out_folder, control_textfile='',
                         block_size=10, coord_columns=[], epsg=0, display_graphics=False):
    """Prepare data for vesper kriging and create a windows batch file to run outside the python/pyprecag environment.

    Outputs:  The following files will be added to the vesper sub-folder in the output folder.
        *_control_*.txt  - The vesper control file
        *_vesperdata_*.csv - The subset of data to krige. All non-required columns are deleted.
        Do_Vesper.bat - A windows batch file to launch vesper kriging for all control files in the vesper sub-folder
                        outside the python/pyprecag environment.

        Any files matching the derived kriged output names will be deleted to ensure that the files in the folder
        belong to the newly create control file. If old version are to be kept, manually rename first.

    DEPENDANCY: Vesper must already be installed on the PC and the install directory added to the config.json file.

    Args:
        in_dataframe (geopandas.geodataframe.GeoDataFrame, pandas.core.frame.DataFrame):
        krig_column (str): The column containing the data to krige.
        grid_filename (str): The vesper grid file.
        control_textfile (str): The name of the control text file without the path
        out_folder (str): The folder to add outputs too. A 'Vesper' sub directory will be created
        block_size (int): The size to apply for block Kriging. This is the same units as the coordinate system.
        coord_columns (List): The columns representing the X and Y coordinates.
        epsg (int) : The epsg_number number for the data. If 0 the en_epsg or enepsg column (if exists) will be used.
        display_graphics (bool): Option to display graphics while running vesper kriging.

    Returns:
       vesper_batfile, vesper_ctrlfile: The paths to the generated batch file and control file.
    """

    if not isinstance(in_dataframe, (gpd.GeoDataFrame, pd.DataFrame)):
        raise TypeError('Invalid input data :in_dataframe')

    if not os.path.exists(grid_filename):
        raise IOError("Invalid path: {}".format(grid_filename))

    if not isinstance(block_size, (int, long)):
        raise TypeError('block_size must be an integer'.format(block_size))

    if out_folder.strip() in [None, '']:
        raise TypeError('Please specify an output folder')

    if not os.path.exists(os.path.dirname(out_folder)):
        raise IOError('Output directory {} does not exist'.format(os.path.dirname(out_folder)))

    if not isinstance(coord_columns, list):
        raise TypeError('Coordinate columns should be a list.'.format(coord_columns))

    vesper_exe = config.get_config_key('vesperEXE')

    if len(coord_columns) == 0:
        coord_columns = predictCoordinateColumnNames(in_dataframe.columns)

    for eaFld in [krig_column] + coord_columns:
        if eaFld not in in_dataframe.columns:
            raise TypeError('Column {} does not exist'.format(eaFld))

    if not isinstance(epsg, (int, long)):
        raise TypeError('EPSG {} must be a integer.'.format(epsg))

    start_time = time.time()

    # Create a filename compatible copy of the krig_column
    krig_col_file = re.sub('[^A-Za-z0-9_-]+', '', unidecode(unicode(krig_column)))

    outSubName = os.path.basename(grid_filename)[:20]

    if control_textfile.strip() in [None, '']:
        control_textfile = "{}_control_{}.txt".format(outSubName, krig_col_file)
        data_file = "{}_vesperdata_{}.csv".format(outSubName, krig_col_file)
        grid_file = "{}_vespergrid_{}.txt".format(outSubName, krig_col_file)
        report_file = "{}_report_{}.txt".format(outSubName, krig_col_file)
        kriged_file = "{}_kriged_{}.txt".format(outSubName, krig_col_file)
        param_file = "{}_parameter_{}.txt".format(outSubName, krig_col_file)

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
    vesper_batfile = os.path.join(vesper_outdir, "Do_Vesper.bat")

    if not os.path.exists(vesper_outdir):
        os.mkdir(vesper_outdir)

    # check Permission Denied or WindowsError: [Error 32]) errors which may indicate that the files
    # in use either by vesper or another app.

    # Start always start with the control file. this is what gets used within QGIS
    filesList = glob.glob(os.path.join(vesper_outdir, "{}_*_{}.*".format(outSubName, krig_col_file)))
    if vesper_ctrlfile in filesList:
        filesList.insert(0, filesList.pop(filesList.index(vesper_ctrlfile)))

    for eaFile in filesList:
        try:
            with open(eaFile, "a") as f:
                pass
        except IOError as e:
            raise IOError('File(s) in use - {}'.format(eaFile))

    # if control file already exists then replace it, and delete all matching kriging outputs and tiff files.
    for eaFile in glob.glob(os.path.join(vesper_outdir, "{}_*_{}.*".format(outSubName, krig_col_file))):
        os.remove(eaFile)
        LOGGER.debug('Deleted file {}'.format(eaFile))

    x_field, y_field = coord_columns

    keep_cols = [findCol for findCol in ['EN_EPSG', 'ENEPSG'] if findCol in in_dataframe.columns]
    if epsg == 0 and len(keep_cols) > 0:
        for col in keep_cols:
            if in_dataframe.iloc[0][col] > 0:
                epsg = in_dataframe.iloc[0][col]
                break

    keep_cols += [x_field, y_field, krig_column]

    dfCSV = in_dataframe[keep_cols].copy()
    dfCSV.to_csv(vesper_datafile, index=False)

    shutil.copy2(grid_filename, vesper_gridfile)

    x_field, y_field = coord_columns
    i_x_column = dfCSV.columns.tolist().index(x_field) + 1
    i_y_column = dfCSV.columns.tolist().index(y_field) + 1
    i_k_column = dfCSV.columns.tolist().index(krig_column) + 1

    # ------------------------------------------------------------------------------------------------------------------
    # Write a VESPER control file for LOCAL KRIGING ONLY
    with open(vesper_ctrlfile, 'w') as wCntrlFile:
        wCntrlFile.write("$vsl\n")
        wCntrlFile.write("ivers=16121\n")
        wCntrlFile.write(
            "title='kriging of {} configured by pyprecag epsg={}'\n".format(control_textfile, epsg))

        wCntrlFile.write("datfil='{}'\n".format(data_file))
        wCntrlFile.write("gridfile='{}'\n".format(grid_file))
        wCntrlFile.write("outdir=''\n".format(vesper_outdir))
        wCntrlFile.write("repfil='{}'\n".format(report_file))
        wCntrlFile.write("outfil='{}'\n".format(kriged_file))
        wCntrlFile.write("parfil='{}'\n".format(param_file))

        wCntrlFile.write("numcol={}\n".format(len(dfCSV.columns)))
        wCntrlFile.write("icol_x={}\n".format(i_x_column))
        wCntrlFile.write("icol_y={}\n".format(i_y_column))
        wCntrlFile.write("icol_z={}\n".format(i_k_column))
        wCntrlFile.write("jordkrg=1\n")
        wCntrlFile.write("jpntkrg=0\n")             # 1 = point kriging , 0 = block kriging
        wCntrlFile.write("jlockrg=1\n")             # 1 = local variogram kriging , 0 = global variogram
        wCntrlFile.write("nest=10\n")
        wCntrlFile.write("dstinc=1\n")
        wCntrlFile.write("valmis=-9999\n")          # missing value
        wCntrlFile.write("jsetint=0\n")
        wCntrlFile.write("xlint=0\n")
        wCntrlFile.write("xhint=0\n")
        wCntrlFile.write("ylint=0\n")
        wCntrlFile.write("yhint=0\n")
        wCntrlFile.write("jsetrad=0\n")             # 1 = set radius, 0 = calculate radius
        wCntrlFile.write("radius=100\n")            # search radius (when jsetrad=1)
        wCntrlFile.write("minpts=90\n")             # min. no. of points for interpolation
        wCntrlFile.write("maxpts=100\n")            # max. no. of points for interpolation
        wCntrlFile.write("sigsqr=0\n")
        wCntrlFile.write("isomod=1\n")
        wCntrlFile.write("modtyp=2\n")
        wCntrlFile.write("isearch=0\n")
        wCntrlFile.write("igeos=0\n")
        wCntrlFile.write("icircs=0\n")
        wCntrlFile.write("phi=0\n")
        wCntrlFile.write("psin=0\n")
        wCntrlFile.write("pcos=0\n")
        wCntrlFile.write("jcomvar=1\n")
        wCntrlFile.write("nlag=20\n")
        wCntrlFile.write("hmax=0\n")
        wCntrlFile.write("tolag=10\n")
        wCntrlFile.write("iwei=1\n")
        wCntrlFile.write("jigraph={}\n".format(int(display_graphics)))  # 1=show graph of variogram, otherwise 0
        wCntrlFile.write("jimap={}\n".format(int(display_graphics)))    # 1=show map of interpolation, otherwise 0
        wCntrlFile.write("CO=0\n")
        wCntrlFile.write("C1=1\n")
        wCntrlFile.write("A1=10\n")
        wCntrlFile.write("C2=1\n")
        wCntrlFile.write("A2=1\n")
        wCntrlFile.write("Alfa=1\n")
        wCntrlFile.write("xside={}\n".format(block_size))           # Block size (in x direction) for block kriging
        wCntrlFile.write("yside={}\n".format(block_size))           # Block size (in y direction) for block kriging
        wCntrlFile.write("lognorm=0\n")
        wCntrlFile.write("itrend=0\n")
        wCntrlFile.write("iconvex=0\n")
        wCntrlFile.write("igrids=1\n")
        wCntrlFile.write("$end\n")

    with open(vesper_batfile, 'w') as wBatFile:
        wBatFile.write("@echo off\n")
        # enable delayed variable expansion, t
        wBatFile.write('setlocal enabledelayedexpansion\n')
        wBatFile.write("echo Processing control files in %CD%\n")
        wBatFile.write("echo.\n")
        wBatFile.write('set icount=0\n')
        wBatFile.write('set tcount=0\n')
        wBatFile.write(r'for %%x in (*control*.txt) do set /a tcount+=1' + '\n')  # count the number of control files
        wBatFile.write('echo Found %tcount% control file(s) to process....\n')
        wBatFile.write("echo.\n")
        wBatFile.write("FOR %%f IN (*control*) DO (" + "\n")
        wBatFile.write(r'set /a icount+=%icount%+1' + '\n')

        # this will suppress newline char and enable 'FINISHED' to be added to the same line.
        # setlocal enabledelayedexpansion and ! instead of % will evaluate a variable when they appear
        # instead of when the batch is parsed
        wBatFile.write("<nul set /p mystr=Kriging !icount! of %tcount% - %%~nxf " + "\n")
        wBatFile.write(r'   start "" /HIGH /WAIT /SHARED "{}" %%~nxf'.format(vesper_exe) + "\n")
        wBatFile.write("    echo - Finished" + "\n")
        wBatFile.write(")\n")  # close bracket for for /f do loop
        # wBatFile.write("echo.\n")
        # wBatFile.write("pause\n")

    LOGGER.info('{:<30}\t{dur:<15}\t{}'.format(inspect.currentframe().f_code.co_name, '',
                                               dur=datetime.timedelta(seconds=time.time() - start_time)))
    return vesper_batfile, vesper_ctrlfile


def run_vesper(ctrl_file, bMinimiseWindow=True):
    """
    Run VESPER for selected control file.

    By default the current working directory will be set to the folder containing the control file.

    By default processing will wait for the current VESPER control file to complete.
    Args:
        ctrl_file (str):  the control file to run as VESPER argument.
        bMinimiseWindow (bool):  Option to automatically minimise the VESPER window on launch.
    """

    task_time = time.time()
    vesper_exe = config.get_config_key('vesperEXE')

    info = subprocess.STARTUPINFO()
    if bMinimiseWindow:
        # run vesper minimised.

        info.dwFlags = subprocess.STARTF_USESHOWWINDOW
        info.wShowWindow = 6  # SW_MINIMIZE


    # need to change to the path of the control file so use cwd=
    process = subprocess.Popen([vesper_exe, ctrl_file], cwd=os.path.dirname(ctrl_file), startupinfo=info)
    process.wait()
    logging.info('{:<30}\t{dur:<15}'.format('Vesper Kriging', dur=datetime.timedelta(seconds=time.time() - task_time)))
