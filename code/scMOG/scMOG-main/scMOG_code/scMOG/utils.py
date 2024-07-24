# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 09:14:43 2023

@author: LCY
"""
import numpy as np
import pandas as pd
import scipy


def ensure_arr(x) -> np.ndarray:
    """Return x as a np.array"""
    if isinstance(x, np.matrix):
        return np.squeeze(np.asarray(x))
    elif isinstance(x, np.ndarray):
        return x
    elif isinstance(x, (scipy.sparse.csr_matrix, scipy.sparse.csc_matrix)):
        return x.toarray()
    elif isinstance(x, (pd.Series, pd.DataFrame)):
        return x.values
    else:
        raise TypeError(f"Unrecognized type: {type(x)}")














