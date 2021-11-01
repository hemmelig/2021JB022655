# -*- coding: utf-8 -*-
"""
This script will create a summary figure for an example splitting measurement
Specifically, this script can be used to generate Figure 2 of:

    On the origin of seismic anisotropy in the shallow crust of the Northern
    Volcanic Zone, Iceland
    Bacon, C.A., Johnson, J., White, R.S., and Rawlinson, N.

which has been submitted to the Journal of Geophysical Research - Solid Earth.

"""

import pathlib

import pandas as pd

import functions as mfast_tools


results_dir = pathlib.Path.cwd()

event = "2013.101.11.50.21"
station = "FLAT"
uid = f"{event}.{station}"

# Read in results 
results = pd.read_csv(mfast_tools.resolve(results_dir / uid,
                                         f"{uid}.*.fb1.*.res")).loc[0]

mfast_tools.plot_summary(results_dir, event, station, results, results_dir)
