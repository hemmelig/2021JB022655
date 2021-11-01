# -*- coding: utf-8 -*-
"""
Shear wave anisotropy (SWA)
The degree of splitting (i.e. the strength of anisotropy) is in part controlled
by the density and aspect ratio of cracks and aligned features within the
crust. The delay time is a path-integrated value, which means that in some
cases it can be a somewhat misleading metric against which to assess anisotropy
in the crust. The shear wave anisotropy (SWA) is a measure of the amount of
splitting normalised by the distance traversed by a ray. It has been used in
some works on a local scale in a volcanic setting (Nowacki et al., 2018).

This script will create a file containing gridded SWA values for this study in
a format that can be plotted using GMT. Specifically, this script can be used
to generate some of the input data required for Figure 6 of:

    On the origin of seismic anisotropy in the shallow crust of the Northern
    Volcanic Zone, Iceland
    Bacon, C.A., Johnson, J., White, R.S., and Rawlinson, N.

which has been submitted to the Journal of Geophysical Research - Solid Earth.

"""

# --- Import libraries ---
import pathlib

import numpy as np
import pandas as pd
import pyproj

import functions as swa


data_dir = pathlib.Path.cwd() / "../data"
(pathlib.Path.cwd() / "xyz_files").mkdir(exist_ok=True)

# --- Read in station file ---
stations = pd.read_csv(data_dir / "stations/used.stations")

# --- Read in velocity model ---
vmodel = pd.read_csv(data_dir / "velocity_model/askja_vmodel.txt", header=None,
                     delim_whitespace=True)[::-1].reset_index(drop=True)

# --- Read in filtered summary file ---
summary = pd.read_csv(data_dir / "filtered_shallow_results_sww.summ")

# Drop all stations not in the stations file
for i, row in summary.iterrows():
    stat = row["stat"]
    if stat not in stations["Name"].values:
        summary.drop(summary.loc[summary["stat"] == stat].index, inplace=True)
summary.reset_index(drop=True, inplace=True)

# --- Calculate SWA for all events shallower than 10 km ---
layer_thickness = 4
summary["pc_swa"] = swa.pc_swa(summary, vmodel, stations,
                               fixed_layer=layer_thickness)

# Filter out some of the nonsense results (station < event)
summary = summary[summary["pc_swa"] < 10000]

# Regional average SWA
print(f"Average SWA: {summary.pc_swa.mean()}")

# --- Recursively grid the SWA data ---
# Overlay a grid with nodes every 500m onto the region and calculate the
# average SWA within a 2.5 km radius of the node, excluding any node with fewer
# than 10 observations.
run_name = f"fixed_{layer_thickness}km_layer"

# Set up grid
increment = 500
x = np.arange(-18000, 18000+increment, increment)
y = np.arange(-17000, 17000+increment, increment)

swa_out, big_df = swa.grid_data(summary, "pc_swa", x, y, 2500, 10, run_name)

# Convert .xyz to lon/lat coords.
cproj = pyproj.Proj("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
gproj = pyproj.Proj("+proj=lcc +lon_0=-16.5 +lat_0=65.1 +lat_1=64.7 +lat_2=65.5 +datum=WGS84 +units=m +no_defs") 
swa.convert_xyz_file(f"xyz_files/{run_name}.xyz", cproj, gproj)
