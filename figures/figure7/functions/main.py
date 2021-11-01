# -*- coding: utf-8 -*-
"""
Helper functions for plotting gridded phi measurements

Author: Conor Bacon
Date: 29 February 2020
"""

import pathlib

import numpy as np
import pandas as pd
import pyproj


def grid_data(summary, grid_on, lbound, ubound, outfile, cproj, gproj,
              min_cell_size, n_min, n_max):
    """
    Recursive function to perform the sub-gridding.

    Parameters
    ----------
    summary : pandas DataFrame object
        A summary file in the MFAST format.

    lbound : array-like [float, float]
        Lower-left bound for the cell.

    ubound : array-like [float, float]
        Upper-right bound for the cell.

    outfile : str
        Name of file to write data to.

    """

    # Calculate current cell size
    cell_size = ubound - lbound

    # Calculate new cell size
    new_cell_size = cell_size / 2

    # Test if the cell size has reached the minimum scale
    divide = True
    if np.all(np.less(new_cell_size / 2, min_cell_size)):
        print("Reached minimum cell size.")
        divide = False

    # Divide the cell into 4 cells
    print("Divide into 4 cells...")
    for i in range(2):
        for j in range(2):
            ll = lbound + np.array([new_cell_size[0]*i, new_cell_size[1]*j])
            ur = ll + new_cell_size

            # Filter summary file
            df = summary[(summary["midx"] > ll[0]) & (summary["midx"] <= ur[0]) &
                         (summary["midy"] > ll[1]) & (summary["midy"] <= ur[1])]

            # If the number of measurements exceeds the cell max count, subdivide
            if len(df) > n_max and divide:
                print("Subdividing...")
                grid_data(df, grid_on, ll, ur, outfile, cproj, gproj,
                          min_cell_size, n_min, n_max)
            elif len(df) > n_min:
                # Calculate the average value of phi for the cell
                R, phi, Rh = find_R(df[grid_on].values)
                #phi_unc = np.sqrt(1 - Rh)

                # Write out the line
                with open(outfile, "a+") as f:
                    # Make the midpoint the mean location of the midpoints?
                    midpoint = (ll + ur) / 2
                    midpoint = pyproj.transform(gproj, cproj, midpoint[0], midpoint[1])
                    llpx, llpy = pyproj.transform(gproj, cproj, ll[0], ll[1])
                    urpx, urpy = pyproj.transform(gproj, cproj, ur[0], ur[1])

                    f.write((f"{llpx} {llpy} {urpx} {urpy} "
                             f"{midpoint[0]} {midpoint[1]} "
                             f"{len(df)} {phi:3.1f} {Rh}\n"))
            else:
                # Just write out the cell summary file
                p = pathlib.Path.cwd()
                outd = p / "gridded_results"
                # df.to_csv(str(outd / ))


def find_R(angles):
    """
    Calculate the resultant vector.

    Parameters
    ----------
    angles : array-like
        Orientation data from which to calculate the resultant vector.

    """

    factor = 2 * np.pi / 180.
    angles = angles * factor

    xr = np.sum(np.cos(angles))
    yr = np.sum(np.sin(angles))

    n = len(angles)

    R = np.sqrt(xr ** 2 + yr ** 2)
    phi = np.arctan2(yr, xr) * 180 / (2*np.pi)
    Rh = R / n

    return R, phi, Rh
