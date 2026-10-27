#!/usr/bin/env python
# coding: utf-8

"""
hi_kinematics_functions
=======================

Shared helpers extracted from the notebook exports in this directory.

Every function here was previously copy-pasted into many notebooks. The five
marked ZERO-DRIFT below were byte-identical (modulo whitespace and comments)
in every copy, so consolidating them is behaviour-preserving:

    get_subcube      9 copies  - 148, 152, 154, 156, 157, 158, 159, 164, 165
    get_theta        8 copies  - 148, 154, 156, 157, 158, 159, 164, 165
    pearson          8 copies  - 148, 152, 154, 155, 156, 158, 162, 165
    second_moment    8 copies  - 148, 154, 156, 157, 158, 159, 164, 165
    separate_clouds  7 copies  - 044, 045, 046, 047, 048, 049, 050

`first_moment` is included because `second_moment` calls it, but it is NOT
zero-drift -- see the note on that function before relying on it.

Bodies are reproduced verbatim from the notebooks; nothing has been
"fixed" in the move. Known sharp edges are documented in the docstrings
rather than silently changed, so results stay reproducible against the
existing figures. See the module notes at the bottom.
"""

import numpy as np
from astropy import units as u


# =====================================================================
#  cube slicing
# =====================================================================

def get_subcube(bigcube, coords, dx=-0.125, dy=0.125, x0=360, y0=-50):
    """Cut a (lon, lat) box out of a cube, leaving the velocity axis whole.

    ZERO-DRIFT: identical in 148, 152, 154, 156, 157, 158, 159, 164, 165.

    Parameters
    ----------
    bigcube : ndarray, shape (nv, nb, nl)
        Velocity-first cube.
    coords : sequence of 4 floats
        (lon_outer, lon_inner, lat_outer, lat_inner) in degrees.
    dx, dy, x0, y0 : float
        Grid parameters, inverting the pixel->world convention
        lon = px*dx + x0 and lat = px*dy + y0. The defaults are the
        HI4PI/DHT grid used throughout these notebooks.

        NOTE: other parameterisations appear elsewhere in the corpus
        (lon offsets of 720 and 180.125, lat offset of -30.125). If you
        are working off a differently-cropped cube, pass them explicitly
        rather than relying on the defaults.

    Returns
    -------
    ndarray
        View into `bigcube` covering the requested box.
    """
    lo, li, bo, bi = coords  # assume in degrees
    xo, xi = ((np.array([lo, li]) - x0) / dx).astype('int')
    yo, yi = ((np.array([bo, bi]) - y0) / dy).astype('int')

    return bigcube[:, yo:yi, xo:xi]


# =====================================================================
#  moment maps
# =====================================================================

def first_moment(cube, arr, oned=False):
    """Intensity-weighted mean along axis 0, plus a log column-density map.

    *** NOT ZERO-DRIFT -- read this before using. ***

    Three variants exist in the notebooks:

      * 148, 154, 158, 159, 165 - this body (has the `oned` flag)
      * 156, 157, 164           - identical maths, no `oned` flag
      * 152                     - uses 1.97e21 instead of 1.97e20

    This is the first of those. It is a strict superset of the
    156/157/164 body (`oned` defaults to False, so old calls are
    unaffected). File 152's 1.97e21 is a 1-dex difference in the second
    return value and is treated here as a typo -- if 152's numbers are
    ones you published, reconcile before using this version.

    Returns
    -------
    tuple of (ndarray, ndarray)
        (first moment, log10 of zeroth * 1.3 * 1.97e20)

    The second element is a log column density, not a moment; it rides
    along because every caller in the notebooks expects the 2-tuple.
    The 1.3 is the channel width in km/s and 1.97e20 is X_CO.
    """
    zeroth = np.nansum(cube, axis=0)

    if oned:
        num = np.nansum(cube * arr)
    else:
        num = np.nansum(cube * arr[:, np.newaxis, np.newaxis], axis=0)

    return num / zeroth, np.log10(zeroth * 1.3 * 1.97e20)


def second_moment(cube, arr):
    """Intensity-weighted variance along axis 0.

    ZERO-DRIFT: identical in 148, 154, 156, 157, 158, 159, 164, 165.

    Note this returns the *variance*, not the velocity dispersion --
    take the square root if you want sigma.

    Returns
    -------
    tuple of (ndarray, ndarray)
        (second moment, first moment) -- the first moment is returned
        too so callers don't have to recompute it.
    """
    zeroth = np.nansum(cube, axis=0)
    first = first_moment(cube, arr)[0]

    num = np.nansum(cube * (arr[:, np.newaxis, np.newaxis]
                            - first[np.newaxis, :]) ** 2, axis=0)
    return num / zeroth, first


# =====================================================================
#  geometry
# =====================================================================

def get_theta(xy):
    """Position angle of a point in degrees, wrapped to (-180, 180].

    ZERO-DRIFT: identical in 148, 154, 156, 157, 158, 159, 164, 165.

    Parameters
    ----------
    xy : sequence of 2 floats
        (x, y). Scalars, not arrays -- the `if` statements below are
        scalar tests and will raise on array input.

    NOTE: x == 0 divides by zero. With numpy scalars that yields inf
    (and a RuntimeWarning) so arctan returns 90 degrees, which is
    correct; with plain Python floats it raises ZeroDivisionError.
    Preserved as-is from the notebooks.
    """
    x, y = xy
    thet = np.arctan(y / x) * u.rad.to(u.deg)

    if x < 0:
        thet += 180

    if thet > 180:
        thet -= 360

    return thet


# =====================================================================
#  statistics
# =====================================================================

def pearson(xs, ys):
    """Pearson correlation coefficient, NaN-tolerant in the sums.

    ZERO-DRIFT: identical in 148, 152, 154, 155, 156, 158, 162, 165.

    NOTE: `N = len(xs)` counts NaNs even though the sums skip them, so
    the result is only exact when xs and ys have no NaNs. This is the
    behaviour the existing figures were made with, so it is preserved.
    For NaN-heavy inputs, mask first or use scipy.stats.pearsonr.
    """
    N = len(xs)
    xmn = np.nanmean(xs)
    ymn = np.nanmean(ys)

    num = np.nansum(xs * ys) - N * xmn * ymn
    den = (np.sqrt(np.nansum(xs * xs) - N * xmn ** 2)
           * np.sqrt(np.nansum(ys * ys) - N * ymn ** 2))

    return num / den


# =====================================================================
#  segmentation
# =====================================================================

def separate_clouds(labels_map, data, N):
    """Split a labelled cube into one masked cube per label.

    ZERO-DRIFT: identical in 044, 045, 046, 047, 048, 049, 050.

    Parameters
    ----------
    labels_map : ndarray
        Integer label cube, same shape as `data`. Labels run 1..N;
        0 is treated as background and dropped.
    data : ndarray, shape (x, y, z)
    N : int
        Number of labels to extract.

    Returns
    -------
    ndarray, shape (N, x, y, z)
        out[n-1] holds `data` masked to label n, zeros elsewhere.

    NOTE: allocates N full-size copies of `data`. For a large cube with
    many labels this will be big -- N * data.nbytes.
    """
    x, y, z = data.shape
    out = np.zeros((N, x, y, z))
    for n in range(1, N + 1):
        m = (labels_map == n)
        mskd = m * data
        out[n - 1] = mskd

    return out


# =====================================================================
#  module notes
# =====================================================================
#
# Provenance
# ----------
# Bodies were taken from file 148 (get_subcube, get_theta, pearson,
# second_moment, first_moment) and file 044 (separate_clouds), which
# carry the fullest comments. All other copies are textually identical
# after stripping whitespace and comments.
#
# Not included
# ------------
# Functions that are duplicated but have drifted, and so need a
# reconciliation decision before they can be shared:
#
#   fit_gradient    16 files, 11 distinct bodies  (the worst)
#   plane_rvld      16 files,  2 distinct bodies  (nearly stable)
#   read_DHT        12 files,  4 distinct bodies
#   moment_maps      8 files,  6 distinct bodies
#   make_zxz         7 files,  5 distinct bodies
#   calc_mass        8 files,  4 distinct bodies
#
# Also note `surface_density`, `transform_1D`, `transform_2D`, `dist`,
# `rms_noise` and the `l_deg`/`b_deg`/`v_kms` coordinate helpers already
# exist in ../../pys/ -- prefer those over re-adding them here, but check
# the constants first: the pys/reid_analysis.py coordinate helpers use
# different offsets than the inlined notebook versions.
