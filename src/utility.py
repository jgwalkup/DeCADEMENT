"""
This module, utility.py, contains three functions facilitating the calculation.

    1) LHS():    perform 1-D latin hypercube sampling
    2) expand(): put a community on the spatial grid
    3) export(): export the output object to the local disk
"""

import pickle

import numpy as np
import pandas as pd
from scipy.stats import distributions


def LHS(n, loc, upc, dist):
    """
    Latin hypercube sampling.

    Parameters:
        n:    integer; size of desired sampling
        loc:  scalar; lower bound of desired distribution
        upc:  scalar; upper bound of desired distribution
        dist: string; either 'uniform' or 'normal'
    Returns:
        lhs: 1D array
    """

    lower_limits = np.arange(0, n) / n
    higher_limits = np.arange(1, n + 1) / n

    points = np.random.uniform(low=lower_limits, high=higher_limits, size=n)
    np.random.shuffle(points)

    scale = upc - loc

    # BUGFIX (2026-07-22): when loc == upc (a zero-width min/max range -- e.g.
    # Osmo_Consti_Prod_min == Osmo_Consti_Prod_max in a real parameter file),
    # scale is exactly 0. scipy's rv_continuous.ppf() internally gates every
    # output behind `cond0 = self._argcheck(*args) & (scale > 0)` (strict
    # inequality), so with scale == 0, cond0 is False for every input and
    # ppf() returns its NaN fill value unconditionally -- regardless of the
    # sampled `points`. That silently made LHS() return an all-NaN array
    # whenever a trait's min/max were set equal, which is a well-defined
    # degenerate case (every draw should just be the constant `loc`), not
    # an error. This NaN then propagated into taxon-level production-rate
    # tables (e.g. Consti_Osmo_C/Induci_Osmo_C) at initialization -- present
    # from day 0, in every taxon and every grid cell, regardless of any
    # simulation logic downstream. `points`/np.random are still consumed
    # identically above so the RNG stream is unaffected for every other
    # (non-degenerate) call to LHS() elsewhere in the model; only the final
    # mapping to a distribution value changes here.
    if scale == 0:
        return np.full(n, loc, dtype=np.asarray(points).dtype)

    if dist == "uniform":
        rv = distributions.uniform(loc=loc, scale=scale)
    elif dist == "normal":
        rv = distributions.norm(loc=loc, scale=scale)

    lhs = rv.ppf(points)

    return lhs


def expand(df, gridsize):
    """
    Put data (df/series) on a spatial grid.

    Parameters:
         df:       dataframe/series
         gridsize: integer
    Return:
         df_expanded: dataframe/series;
    """
    df_expanded = pd.concat([df] * gridsize)

    return df_expanded


def export(output, name):
    """
    Save the output object as a .pickle file.

    Parameters:
        output: object of the class, Output
        name:   naming output file
    """

    with open(str(name) + ".pickle", "wb") as f:
        pickle.dump(output, f, pickle.HIGHEST_PROTOCOL)
