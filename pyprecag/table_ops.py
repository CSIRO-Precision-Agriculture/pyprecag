import logging

import pandas as pd
from geopandas import GeoDataFrame
from scipy.stats import stats

from . import config, TEMPDIR

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())


def t_test(arr):
    """ Moving window t-test function.
        to get two columns in to the function, one must be set as the index.

        Example:
            input_table = input_table.set_index(treatment_column, drop=False)

            input_table['p_value'] = input_table['controls_mean'].rolling(
                window=size, center=True).apply(t_test, raw=False)
    """

    tstat, pvalue = stats.ttest_rel(arr.index, arr)
    return pvalue


def response_index(arr):
    """ Moving window response index function
        to get two columns in to the function, one must be set as the index.
        Example:
            input_table = input_table.set_index(treatment_column, drop=False)
            input_table['RI'] = input_table['controls_mean'].rolling(
                window=size, center=True).apply(response_index, raw=False)
    """

    idx_s = arr.index.to_series()
    return idx_s.mean() / arr.mean()


def calculate_strip_stats(input_table, treatment_column, control_columns=[], size=5):
    """Calculate statistics for a strip

    A moving window is used for some of the statistics. This window is centred so for a window
    size of 5, 2 NAN or blanks will be added to start and end of the output column.

    Statistics include (output column names):
        controls_mean  - row by row mean of the control columns
        treat_diff -   row by row difference between the treatment and controls_mean columns
        av_treat_diff - calculate mean of values using a moving window using the treat_diff column
        p_value - calculate p_value using a moving window using treatment and controls_mean columns
        RI  - Response Index using the treatment and controls_mean columns

    Args:
        input_table (pandas.core.frame.DataFrame): the table to calculate statistics for
        treatment_column (str): The column containing the treatment values
        control_columns (List[str]): The column containing the control values.
                                     This can be one or two columns
        size (int):The size of the moving window.

    Returns:
        pandas.core.frame.DataFrame: The output table containing new statistics columns
        control_mean (str): The column used as the control mean.
    """

    if isinstance(input_table, GeoDataFrame):
        # drop geometry etc. and create flat table.
        input_table = pd.DataFrame(input_table.drop(columns='geometry', inplace=False))

    if not isinstance(control_columns, list):
        raise TypeError("control_columns should be a list.")

    if not isinstance(treatment_column, basestring):
        raise TypeError("treatment_column should be a string.")

    if treatment_column is None or treatment_column == '':
        raise ValueError('Invalid treatment column')

    missing = [ea for ea in [treatment_column] + control_columns
               if ea and ea not in input_table.columns]
    if len(missing) > 0:
        raise ValueError('columns not found - {}'.format(len(missing), ','.join(missing)))

    input_table = input_table.copy()

    ''' Statistics ----------------------------------------------------------------- '''
    if len(control_columns) > 1:
        # calculate the mean for the column(s)
        control_mean = '-'.join(control_columns)
        control_mean = control_mean.replace(' Strip Value', '')
        control_mean = control_mean.replace(' Strip Control', '')
        control_mean = '{}_mean'.format(control_mean.strip())
        input_table[control_mean] = input_table[control_columns].mean(axis=1)

    else:
        control_mean = control_columns[0]

        # calculate the difference
    input_table['treat_diff'] = input_table[treatment_column] - input_table[control_mean]

    # Moving Mean for values diff
    input_table['av_treat_dif'] = input_table['treat_diff'].rolling(size, center=True).mean()

    ''' Rolling window using two-tailed paired student t-test
        https://pythonfordatascience.org/paired-samples-t-test-python/ 
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_rel.html 
        https://stackoverflow.com/a/52029516
    '''

    input_table = input_table.set_index(treatment_column, drop=False)

    input_table['p_value'] = input_table[control_mean].rolling(
        window=size, center=True).apply(t_test, raw=False)

    input_table['RI'] = input_table[control_mean].rolling(
        window=size, center=True).apply(response_index, raw=False)

    # reset index
    input_table.set_index('TrialPtID', drop=False, inplace=True)

    return input_table, control_mean
