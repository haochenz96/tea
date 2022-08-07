import re
import pandas as pd
import numpy as np

CONDENSED_FORMAT = 'chr.*:(\d+):.*/.*$'

def check_matrix_format(input_array, string_format) -> bool:
    """
    Check if EVERY INDEX of the input is in the correct format for the given format type.
    """
    try:
        # every cell must adhere to the format for it to return True
        return input_array.index.map(lambda x: re.search(string_format, x)).isna().sum() == 0 
    except TypeError:
        print('[ERROR] --- Input array is not a pandas.Series or a numpy.ndarray.')
        return False

def isNaN(num):
    # check for float('nan') values
    return num != num
    