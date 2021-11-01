# coding: utf-8
"""
Creating raypaths file
This script will create a file containing raypaths for this study in a format
that can be plotted using GMT. Specifically, this script can be used to
generate some of the input data required for Figure 5 of:

    On the origin of seismic anisotropy in the shallow crust of the Northern
    Volcanic Zone, Iceland
    Bacon, C.A., Johnson, J., White, R.S., and Rawlinson, N.

which has been submitted to the Journal of Geophysical Research - Solid Earth.

"""

# --- Import libraries ---
import pathlib

import pandas as pd


data_dir = pathlib.Path.cwd() / "../data"
(pathlib.Path.cwd() / "raypaths").mkdir(exist_ok=True)

# --- Read in station file ---
stations = pd.read_csv(data_dir / "stations/used.stations")

# --- Read in filtered summary file ---
summary = pd.read_csv(data_dir / "filtered_shallow_results_sww.summ")

# Drop all stations not in the stations file
for i, row in summary.iterrows():
    stat = row["stat"]
    if stat not in stations["Name"].values:
        summary.drop(summary.loc[summary["stat"] == stat].index, inplace=True)
summary.reset_index(drop=True, inplace=True)

# --- Write out raypaths file ---
with open("raypaths/raypaths.xy", "w") as f:
    for i, event in summary.iterrows():
        print(">", file=f)
        print(event["evlo"], event["evla"], file=f)
        print(event["slon"], event["slat"], file=f)
    print(">", file=f)
