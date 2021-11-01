# -*- coding: utf-8 -*-
"""
Gridding measurements of fast direction
We grid the measurements of fast direction by sub-gridding where there is an
excess of data. The measurements are assigned to the midpoint between the
station and event.

This script will create a file containing gridded phi values for this study in
a format that can be plotted using GMT. Specifically, this script can be used
to generate some of the input data required for Figure 7 of:

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

import functions as gridphi


data_dir = pathlib.Path.cwd() / "../data"
(pathlib.Path.cwd() / "gridded_results").mkdir(exist_ok=True)

# --- Read in station file ---
stations = pd.read_csv(data_dir / "stations/used.stations")

# --- Read in filtered summary file ---
summary = pd.read_csv(data_dir / "filtered_shallow_results.summ")

# Drop all stations not in the stations file
for i, row in summary.iterrows():
    stat = row["stat"]
    if stat not in stations["Name"].values:
        summary.drop(summary.loc[summary["stat"] == stat].index, inplace=True)
summary.reset_index(drop=True, inplace=True)

# --- Grid setup ---
# Specify the grid size, minimum cell size for the quadtree subdivision and
# the conditions under which division occurs

min_x, min_y = -35000, -35000
cell_size = np.array([70000, 70000])

# Set the minimum cell size
min_cell_size = [2e3, 2e3]

# Set min/max data count for a cell
n_min = 50
n_max = 200

# Get cell bounds
lbound = np.array([min_x, min_y])
ubound = lbound + cell_size

run_name = "figure7"
outfile = f"gridded_results/{run_name}_grid_data.bins"

p = pathlib.Path.cwd()
try:
    (p / outfile).unlink()
except FileNotFoundError:
    pass

# Grid the data
cproj = pyproj.Proj("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
gproj = pyproj.Proj("+proj=lcc +lon_0=-16.5 +lat_0=65.1 +lat_1=64.7 +lat_2=65.5 +datum=WGS84 +units=m +no_defs") 
gridphi.grid_data(summary, "fast", lbound, ubound, outfile, cproj, gproj, min_cell_size, n_min, n_max)

grid = pd.read_csv(outfile, delimiter=" ", header=None)
xy_file = f"gridded_results/{run_name}_grid.xy"
with open(xy_file, "w") as f:
    for i, row in grid.iterrows():
        f.write((">\n"
                 f"{row[0]} {row[1]}\n"
                 f"{row[0]} {row[3]}\n"
                 f"{row[2]} {row[3]}\n"
                 f"{row[2]} {row[1]}\n"
                 f"{row[0]} {row[1]}\n"))
