# coding: utf-8
"""
Modelling changes in delay as a function of depth
This script will create files containing anisotropic moels for this study in a
format that can be plotted using GMT. Specifically, this script can be used to
generate some of the input data required for Figure 3 of:

    On the origin of seismic anisotropy in the shallow crust of the Northern
    Volcanic Zone, Iceland
    Bacon, C.A., Johnson, J., White, R.S., and Rawlinson, N.

which has been submitted to the Journal of Geophysical Research - Solid Earth.

"""

# --- Import libraries ---
import pathlib

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline


data_dir = pathlib.Path.cwd() / "../data"
(pathlib.Path.cwd() / "result").mkdir(exist_ok=True)
(pathlib.Path.cwd() / "models").mkdir(exist_ok=True)

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

inspected_splits = pd.read_csv(data_dir / "inspected_splits.csv")
labels = []
for i, event in summary.iterrows():
    evid = ".".join(event["1event"].split(".")[:6])
    label = inspected_splits[inspected_splits["EventID"] == evid]
    if len(label) == 0:
        labels.append("u")
    else:
        labels.append(label["Label"].values[0])
summary["label"] = labels

labelled_results = summary[summary["label"] != "u"]
good_splits = labelled_results[labelled_results["label"] == 1]
bad_splits = labelled_results[labelled_results["label"] == 0]

good_splits_nm = good_splits[good_splits["stat"] != "MYVO"]
good_splits_nm[["tlag", "depthkm"]].to_csv("results/data_scatter.txt", index=False, header=None)

# --- Make station average elevation GMT file ---
with open("results/mean_elevation.txt", "w") as f:
    for i in [0, 0.2]:
        print(i, -stations["Elevation"].mean(), file=f)

# --- Useful functions ---
def fixed_window_rolling_average(df, thickness, dx=0.1, bounds=[0., 28.0]):
    """
    Calculates a rolling average using a window of fixed width.

    Parameters
    ----------
    df : pandas.DataFrame object
        DataFrame containing splitting results.
    thickness : float
        Size of window, in km.
    dx : float, optional
        Size of increment for rolling average, in km. Default 0.1 km.
    bounds : list of float, optional
        Min/max depth to perform rolling average over, in km. Default 0-28 km.

    Returns
    -------
    depth : list of float
        Depths for each of the assigned values.
    mean : list of float
        Mean values for the rolling window, assigned to the midpoint.
    median : list of float
        Median values for the rolling window, assigned to the midpoint.
    std : list of float
        Standard deviation values for the rolling window, assigned to the
        midpoint.

    """

    depth, amean, median, std = [], [], [], []

    top, max_depth = bounds
    bottom = top + thickness
    while bottom <= max_depth:
        depth.append(top + thickness / 2)
        vals = df[(df["depthkm"] >= top) & (df["depthkm"] <= bottom)]["tlag"]
        amean.append(vals.mean())
        median.append(vals.median())
        std.append(vals.std())

        top += dx
        bottom += dx
    
    return depth, amean, median, std

# Calculate values to plot
thickness = 2.0
bounds = [0, 11]
dx = 0.75
        
x, mean, _, std = fixed_window_rolling_average(
    good_splits,
    thickness,
    bounds=bounds,
    dx=dx)

# Create inputs for GMT
with open("results/1d_means.txt", "w") as f:
    for x_, y, err in zip(x, mean, std):
        print(y, x_, err, err, file=f)

x_nm, mean_nm, _, std_nm = fixed_window_rolling_average(
    good_splits_nm,
    thickness,
    bounds=bounds,
    dx=dx)

# Create inputs for GMT
with open("results/1d_means_nomyvo.txt", "w") as f:
    for x_, y, err in zip(x_nm, mean_nm, std_nm):
        print(y, x_, err, err, file=f)


# --- Create profiles for different anisotropic models ---
# --- Fitting exponential function to data ---
def func(x, a, b, c):
    return a - np.exp(x * b) / c

popt, pcov = curve_fit(func, x_nm, mean_nm, sigma=std_nm, p0=(0.09, 0.9, 9e5))

xs = np.linspace(-1, 10, 1000)
y_exp_nm = func(xs, *popt)
realign = np.argmin(abs(y_exp_nm[::-1]))
with open("models/exp_model.txt", "w") as f:
    for x, y in zip(y_exp_nm, xs[::-1] - xs[realign]):
        print(x, -1*y - stations["Elevation"].mean(), file=f)

# --- Anisotropic layer of finite thickness (from surface) ---
finite_lower_x = [0, popt[0], popt[0]]
finite_y = [-stations["Elevation"].mean(), 3, 10]

d = 3 + stations["Elevation"].mean()
vbar = (1.9 + 2.35 + 2.85 + 3.25 * 0.8) / 3.8

val = 2*d / (popt[0]*vbar)
a_pers = (-1*val + np.sqrt(4 + (val)**2)) * 100
print(f"Average SWA for finite thickness layer down to 3 km ={a_pers:4.1f} %")

with open("models/finite_layer_model.txt", "w") as f:
    for x, y in zip(finite_lower_x, finite_y):
        print(x, y, file=f)

# --- Spline fit to data ---
spl = UnivariateSpline(x_nm, mean_nm)
realign = np.argmin(abs(spl(xs)[::-1]))
with open("models/spline_model.txt", "w") as f:
    for x, y in zip(spl(xs), xs[::-1] - xs[realign]):
        print(x, -1*y - stations["Elevation"].mean(), file=f)
