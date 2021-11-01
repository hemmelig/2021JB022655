# -*- coding: utf-8 -*-
"""
Helper functions for calculating the shear wave anisotropy (SWA)

Author: Conor Bacon
Date: 29 February 2020
"""

import pathlib

import numpy as np
import pandas as pd
import pyproj


def vbar(vmodel, event_depth, station_elevation):
    """
    Calculates the average velocity between source and receiver for a given
    velocity model.

    Only need the difference in the vertical axis as sines of angles will
    cancel.

    Parameters
    ----------
    vmodel : pandas DataFrame
        Contains the velocity model, ordered by deepest layer first.
    event_depth : float
        Depth of event. Units should be consistent with station_elevation.
    station_elevation : float
        Elevation of event. Units should be consistent with depth of event.

    Returns
    -------
    vbar : float
        Average velocity within the model.

    """

    average = 0.0

    for i, layer in vmodel.iterrows():
        layer_top = vmodel.iloc[i+1][0]
        layer_bottom = layer[0]

        if station_elevation < layer_bottom and station_elevation > layer_top:
            layer_top = station_elevation
        if layer_top == -100.0:
            break
        if event_depth <= layer_top:
            continue
        elif event_depth > layer_top and event_depth <= layer_bottom:
            # Handle interpolated distance
            dist_in_layer = abs(event_depth - layer_top)
        else:
            # Handle full layer
            dist_in_layer = abs(layer_bottom - layer_top)

        average += dist_in_layer * layer[2]

    return average / (event_depth - station_elevation)


def swa(dt, vbar, d):
    """
    Calculates the shear wave anisotropy, as defined in Thomas & Kendall.

    Parameters
    ----------
    dt : float
        Measured delay time.
    vbar : float
        Average velocity between source and receiver.
    d : float
        Straightline distance between source and receiver.

    Returns
    -------
    swa : float
        % anisotropy as measured along ray.

    """

    val = 2*d / (dt*vbar)

    a1 = -1*val + np.sqrt(4 + (val)**2)

    return a1 * 100


def pc_swa(summary, vmodel, stations, fixed_layer=None, deep_layer=None):
    """
    Calculates the percentage anisotropy for all events in a
    summary file.

    Parameters
    ----------
    summary : pandas DataFrame
        Contains shear-wave splitting observations.
    vmodel : pandas DataFrame
        Contains the velocity model, ordered by deepest layer first.
    stations : pandas DataFrame
        Reference file containing station information.
    fixed_layer : float, optional
        Impose a finite thickness of anisotropic layer.

    Returns
    -------
    summary : pandas DataFrame
        Original summary file with % anisotropy column added.

    """

    pc_swas = np.zeros(len(summary))
    for i, event in summary.iterrows():
        # Look up the station elevation and calculate d
        stat = stations[stations["Name"] == event["stat"]]
        selv = -1 * stat["Elevation"].values[0]
        if fixed_layer is None:
            dep = event["depthkm"]
        else:
            dep = fixed_layer if event["depthkm"] > fixed_layer else event["depthkm"]

        if deep_layer is None:
            pass
        else:
            selv = deep_layer
        h = dep - selv
        s = event["10dist(ev-stat)"]
        d = np.sqrt(h**2 + s**2)

        # Calculate average velocity
        vb = vbar(vmodel, dep, selv)

        # Calculate the shear wave anisotropy
        if deep_layer is None:
            pc_swas[i] = swa(event["tlag"], vb, d)
        else:
            diff = event["tlag"] - 0.09
            tlag = diff if diff > 0 else 0
            pc_swas[i] = swa(tlag, vb, d)

    return pc_swas


def grid_data(summary, value, x, y, bin_radius, min_count, run_name):
    """
    Performs a poor man's tomography, as described in Nowacki
    et al., 2018. See paper for details.

    Parameters
    ----------
    summary : pandas DataFrame
        Contains shear-wave splitting observations.
    value : str
        Variable to average.
    x : array
        x coordinates of grid nodes, units m
    y : array
        y coordinates of grid nodes, units m
    bin_radius : float
        Distance from node to include results, units m.
    min_count : int
        Minimum number of events required in a bin.
    run_name : str
        Unique name for the run.

    Returns
    -------
    out : array, shape (x, y)
        Average values of SWA at each grid node.

    """

    out = np.zeros((len(x), len(y)))

    biggest_df = pd.DataFrame()

    outfile = f"xyz_files/{run_name}.xyz"
    p = pathlib.Path.cwd() / outfile
    if p.is_file():
        p.unlink()
    with open(outfile, "a+") as f:
        for i, xval in enumerate(x):
            for j, yval in enumerate(y):
                # Calculate distance between all midpoints and the cell
                d = np.sqrt((summary["midx"] - xval)**2 + (summary["midy"] - yval)**2)
                df = summary[d < bin_radius]
                if len(df) > min_count:
                    out[i, j] = df[value].mean()
                f.write(f"{xval} {yval} {out[i, j]}\n")
                if len(df) > len(biggest_df):
                    biggest_df = df.copy()

    return out, biggest_df


def convert_xyz_file(file_to_convert, cproj, gproj):
    """
    Converts xyz coords to lon/lat coords for a given .xyz file.

    Parameters
    ----------
    file_to_convert : str
        .xyz file to convert.
    cproj : pyproj Proj object
        Coordinate projection.
    gproj : pyproj Proj object
        Grid projection.

    """

    df = pd.read_csv(file_to_convert, header=None, sep=" ")

    out = pd.DataFrame()
    crds = pyproj.transform(gproj, cproj, df[0].values, df[1].values)
    out[0] = crds[0]
    out[1] = crds[1]
    out[2] = df[2]

    out = out.replace(0.0, "NaN")
    out.to_csv(file_to_convert, header=None, index=False, sep=" ")


def box_filter(df, ll_corner, ur_corner):
    """
    Filter a summary file for event within a geographical box.

    Parameters
    ----------
    df : pandas DataFrame
        Summary of SWS measurements
    ll_corner : float, array
        Coordinates of lower left corner of box.
    ur_corner : float, array
        Coordinates of upper right corner of box.

    Returns
    -------
    out : pandas DataFrame
        Filtered summary file.

    """

    # Filter by longitude
    out = df[(df["evlo"] >= ll_corner[0]) & (df["evlo"] <= ur_corner[0])]

    # Filter by latitude
    out = out[(out["evla"] >= ll_corner[1]) & (out["evlo"] <= ur_corner[1])]

    # Filter by depth
    out = out[(out["depthkm"] >= ll_corner[2]) & (out["depthkm"] <= ur_corner[2])]

    return out
